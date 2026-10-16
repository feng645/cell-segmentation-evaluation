suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(scales)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

morphology_theme_academic <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "#4d4d4d", size = base_size - 1),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#f3f4f6", color = "#c7cbd1", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.35)
    )
}

read_json_any_mainv2 <- function(path) {
  if (grepl("\\.gz$", path)) {
    jsonlite::fromJSON(paste(readLines(gzfile(path), warn = FALSE), collapse = ""), simplifyVector = FALSE)
  } else {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
}

ring_matrix_mainv2 <- function(coords) {
  mat <- tryCatch(do.call(rbind, coords), error = function(...) NULL)
  if (is.null(mat) || ncol(mat) < 2) return(NULL)
  mat <- matrix(as.numeric(mat[, 1:2, drop = FALSE]), ncol = 2)
  if (nrow(mat) < 3 || any(is.na(mat))) return(NULL)
  mat
}

close_ring_mainv2 <- function(mat) {
  if (nrow(mat) < 2) return(mat)
  if (!isTRUE(all.equal(mat[1, ], mat[nrow(mat), ]))) {
    mat <- rbind(mat, mat[1, ])
  }
  mat
}

polygon_area_mainv2 <- function(mat) {
  mat <- close_ring_mainv2(mat)
  sum(mat[-1, 1] * mat[-nrow(mat), 2] - mat[-1, 2] * mat[-nrow(mat), 1]) / 2
}

polygon_perimeter_mainv2 <- function(mat) {
  mat <- close_ring_mainv2(mat)
  sum(sqrt(rowSums((mat[-1, , drop = FALSE] - mat[-nrow(mat), , drop = FALSE])^2)))
}

polygon_centroid_mainv2 <- function(mat) {
  mat <- close_ring_mainv2(mat)
  cross <- mat[-nrow(mat), 1] * mat[-1, 2] - mat[-1, 1] * mat[-nrow(mat), 2]
  area2 <- sum(cross)
  if (abs(area2) < 1e-12) {
    return(c(mean(mat[, 1]), mean(mat[, 2])))
  }
  cx <- sum((mat[-nrow(mat), 1] + mat[-1, 1]) * cross) / (3 * area2)
  cy <- sum((mat[-nrow(mat), 2] + mat[-1, 2]) * cross) / (3 * area2)
  c(cx, cy)
}

convex_hull_area_mainv2 <- function(points) {
  if (nrow(points) < 3) return(NA_real_)
  idx <- chull(points)
  hull <- rbind(points[idx, , drop = FALSE], points[idx[1], , drop = FALSE])
  abs(sum(hull[-1, 1] * hull[-nrow(hull), 2] - hull[-1, 2] * hull[-nrow(hull), 1]) / 2)
}

resolve_feature_cell_id_mainv2 <- function(feature) {
  cell_id <- feature$id %||% feature$properties$cell_id %||% feature$properties$cell
  if (is.null(cell_id)) return(NULL)
  as.character(cell_id)
}

normalize_cell_types_mainv2 <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  unknown_idx <- tolower(trimws(x)) %in% c("unknown", "unassigned", "none", "na")
  x[unknown_idx] <- "Unknown"
  x
}

normalize_morphology_coord_columns_mainv2 <- function(dt) {
  dt <- as.data.table(dt)

  choose_coord <- function(primary, secondary, tertiary = NULL) {
    if (!is.null(primary) && primary %in% names(dt)) return(as.numeric(dt[[primary]]))
    if (!is.null(secondary) && secondary %in% names(dt)) return(as.numeric(dt[[secondary]]))
    if (!is.null(tertiary) && tertiary %in% names(dt)) return(as.numeric(dt[[tertiary]]))
    rep(NA_real_, nrow(dt))
  }

  dt[, x := choose_coord("x.y", "x", "x.x")]
  dt[, y := choose_coord("y.y", "y", "y.x")]

  drop_cols <- intersect(c("x.x", "y.x", "x.y", "y.y"), names(dt))
  if (length(drop_cols)) {
    dt[, (drop_cols) := NULL]
  }
  dt
}

extract_geometry_metrics_mainv2 <- function(cell_id, geometries) {
  total_area <- 0
  total_perimeter <- 0
  weighted_centroid_num <- c(0, 0)
  all_points <- NULL

  for (geom in geometries) {
    if (is.null(geom$type) || is.null(geom$coordinates)) next
    polygon_list <- if (toupper(geom$type) == "POLYGON") {
      list(geom$coordinates)
    } else if (toupper(geom$type) == "MULTIPOLYGON") {
      geom$coordinates
    } else {
      next
    }

    for (poly in polygon_list) {
      if (!length(poly)) next
      ext_ring <- ring_matrix_mainv2(poly[[1]])
      if (is.null(ext_ring)) next

      ext_area <- abs(polygon_area_mainv2(ext_ring))
      hole_rings <- if (length(poly) > 1L) Filter(Negate(is.null), lapply(poly[-1], ring_matrix_mainv2)) else list()
      hole_areas <- vapply(hole_rings, function(ring) abs(polygon_area_mainv2(ring)), numeric(1))
      poly_area <- ext_area - sum(hole_areas, na.rm = TRUE)
      if (!is.finite(poly_area) || poly_area <= 0) next

      poly_perimeter <- polygon_perimeter_mainv2(ext_ring) +
        sum(vapply(hole_rings, polygon_perimeter_mainv2, numeric(1)), na.rm = TRUE)
      centroid_num <- polygon_centroid_mainv2(ext_ring) * ext_area
      if (length(hole_rings)) {
        for (hole_ring in hole_rings) {
          hole_area <- abs(polygon_area_mainv2(hole_ring))
          centroid_num <- centroid_num - polygon_centroid_mainv2(hole_ring) * hole_area
          all_points <- rbind(all_points, hole_ring)
        }
      }

      total_area <- total_area + poly_area
      total_perimeter <- total_perimeter + poly_perimeter
      weighted_centroid_num <- weighted_centroid_num + centroid_num
      all_points <- rbind(all_points, ext_ring)
    }
  }

  if (!is.finite(total_area) || total_area <= 0) return(NULL)
  centroid_final <- weighted_centroid_num / total_area
  hull_area <- convex_hull_area_mainv2(all_points)
  solidity <- ifelse(is.finite(hull_area) && hull_area > 0, total_area / hull_area, NA_real_)

  eccentricity <- NA_real_
  if (!is.null(all_points) && nrow(all_points) >= 3) {
    cov_mat <- stats::cov(all_points)
    if (!any(is.na(cov_mat))) {
      eig <- sort(eigen(cov_mat)$values, decreasing = TRUE)
      if (length(eig) >= 2 && eig[1] > 0) {
        a <- sqrt(eig[1])
        b <- sqrt(max(eig[2], 0))
        eccentricity <- sqrt(max(0, 1 - (b^2 / a^2)))
      }
    }
  }

  circularity <- ifelse(total_perimeter > 0, (4 * pi * total_area) / (total_perimeter^2), NA_real_)
  data.table(
    cell_id = cell_id,
    x = centroid_final[1],
    y = centroid_final[2],
    area = total_area,
    perimeter = total_perimeter,
    convex_area = hull_area,
    solidity = solidity,
    eccentricity = eccentricity,
    circularity = circularity
  )
}

extract_feature_metrics <- function(feature) {
  cell_id <- resolve_feature_cell_id_mainv2(feature)
  geom <- feature$geometry
  if (is.null(cell_id) || is.null(geom)) return(NULL)
  extract_geometry_metrics_mainv2(cell_id, list(geom))
}

empty_morphology_metrics_dt_mainv2 <- function() {
  data.table(
    cell_id = character(),
    x = numeric(),
    y = numeric(),
    area = numeric(),
    perimeter = numeric(),
    convex_area = numeric(),
    solidity = numeric(),
    eccentricity = numeric(),
    circularity = numeric()
  )
}

compute_morphology_metrics_raw_mainv2 <- function(geojson_path) {
  gj <- read_json_any_mainv2(geojson_path)
  features <- gj$features %||% list()
  if (!length(features)) return(empty_morphology_metrics_dt_mainv2())

  feature_ids <- vapply(features, resolve_feature_cell_id_mainv2, character(1), USE.NAMES = FALSE)
  keep <- !is.na(feature_ids) & nzchar(feature_ids)
  features <- features[keep]
  feature_ids <- feature_ids[keep]
  if (!length(features)) return(empty_morphology_metrics_dt_mainv2())

  grouped_features <- split(features, feature_ids)
  metrics <- rbindlist(Filter(Negate(is.null), lapply(names(grouped_features), function(cell_id) {
    extract_geometry_metrics_mainv2(
      cell_id = cell_id,
      geometries = lapply(grouped_features[[cell_id]], `[[`, "geometry")
    )
  })), fill = TRUE)
  if (!nrow(metrics)) return(empty_morphology_metrics_dt_mainv2())
  metrics
}

compute_morphology_metrics_from_geojson <- function(geojson_path) {
  message(
    "[morphology] raw grouped geometry backend for ",
    basename(geojson_path)
  )
  compute_morphology_metrics_raw_mainv2(geojson_path)
}

plot_morphology_method_outputs <- function(merged, method_name, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_dt <- copy(merged)[cell_type_custom_final != "Unknown"]
  if (!nrow(plot_dt)) return(invisible(NULL))

  yield_dt <- plot_dt[, .(n = .N), by = cell_type_custom_final][order(-n)]
  p_yield <- ggplot(yield_dt, aes(x = reorder(cell_type_custom_final, n), y = n)) +
    geom_col(fill = "#2C7FB8", width = 0.74) +
    coord_flip() +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste(method_name, "cell yield by type"),
      subtitle = "Annotated cells retained in morphology evaluation",
      x = "Annotated cell type",
      y = "Cell count"
    ) +
    morphology_theme_academic()

  p_area <- ggplot(plot_dt[area > 0], aes(x = cell_type_custom_final, y = area)) +
    geom_violin(fill = "#7FC97F", color = NA, alpha = 0.8, trim = FALSE) +
    geom_boxplot(width = 0.12, outlier.size = 0.35, alpha = 0.6) +
    scale_y_log10(labels = label_number()) +
    coord_flip() +
    labs(
      title = paste(method_name, "cell area distribution"),
      subtitle = "Polygon-derived cell area by annotated cell type",
      x = "Annotated cell type",
      y = "Cell area"
    ) +
    morphology_theme_academic()

  shape_long <- melt(
    plot_dt,
    measure.vars = c("eccentricity", "solidity", "circularity"),
    variable.name = "metric",
    value.name = "value"
  )
  p_shape <- ggplot(shape_long[!is.na(value)], aes(x = cell_type_custom_final, y = value)) +
    geom_violin(fill = "#FDB462", color = NA, alpha = 0.8, scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.12, outlier.size = 0.35, alpha = 0.6) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    coord_flip() +
    labs(
      title = paste(method_name, "shape descriptors"),
      subtitle = "Polygon geometry descriptors by annotated cell type",
      x = "Annotated cell type",
      y = "Metric value"
    ) +
    morphology_theme_academic()

  ggsave(file.path(outdir, "cell_yield_by_type.png"), p_yield, width = 6.6, height = 4.6, dpi = 320)
  ggsave(file.path(outdir, "cell_area_by_type.png"), p_area, width = 7.0, height = 5.2, dpi = 320)
  ggsave(file.path(outdir, "shape_descriptors_by_type.png"), p_shape, width = 7.2, height = 8.8, dpi = 320)
}

plot_morphology_cross_method_summary <- function(summary_dt, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_dt <- copy(summary_dt)[cell_type_custom_final != "Unknown"]
  if (!nrow(plot_dt)) return(invisible(NULL))

  p_area <- ggplot(plot_dt, aes(x = cell_type_custom_final, y = method, fill = area_median)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = comma(round(area_median, 0))), size = 2.5) +
    scale_fill_gradient(low = "#f7fcf5", high = "#00441b", na.value = "#f1f1f1") +
    labs(
      title = "Morphology Area Comparison",
      subtitle = "Median cell area by method and annotated cell type",
      x = "Annotated cell type",
      y = "Segmentation method",
      fill = "Median area"
    ) +
    morphology_theme_academic() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

  p_yield <- ggplot(plot_dt, aes(x = cell_type_custom_final, y = method, fill = proportion)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", proportion)), size = 2.7) +
    scale_fill_gradient(low = "#fff5f0", high = "#99000d", na.value = "#f1f1f1", limits = c(0, 1)) +
    labs(
      title = "Morphology Cell-type Composition",
      subtitle = "Fraction of evaluated cells assigned to each cell type",
      x = "Annotated cell type",
      y = "Segmentation method",
      fill = "Fraction"
    ) +
    morphology_theme_academic() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

  ggsave(file.path(outdir, "morphology_area_comparison.png"), p_area, width = 8.8, height = 4.8, dpi = 320)
  ggsave(file.path(outdir, "morphology_yield_comparison.png"), p_yield, width = 8.8, height = 4.8, dpi = 320)
}

run_morphology_eval <- function(config) {
  output_root <- config$output_root
  eval_dir <- file.path(output_root, "morphology_eval")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)

  manifest <- jsonlite::read_json(file.path(output_root, "manifest.json"), simplifyVector = FALSE)
  methods_dt <- rbindlist(lapply(manifest$methods, as.data.table), fill = TRUE)
  annotation_manifest <- fread(file.path(output_root, "annotation", "annotation_manifest.csv"))

  results <- list()
  summary_rows <- list()
  idx <- 1L

  for (i in seq_len(nrow(methods_dt))) {
    row <- methods_dt[i]
    method_name <- row$name[[1]]
    if (!identical(row$status[[1]], "ready")) next
    if (is.na(row$polygon_path[[1]]) || !file.exists(row$polygon_path[[1]])) next

    annotation_row <- annotation_manifest[method == method_name]
    if (!nrow(annotation_row)) next

    message("[morphology] running method=", method_name)
    morph_dt <- compute_morphology_metrics_from_geojson(row$polygon_path[[1]])
    annotation_dt <- fread(annotation_row$annotation_path[[1]])
    annotation_dt[, cell_id := as.character(cell_id)]
    annotation_dt[, cell_type_custom_final := normalize_cell_types_mainv2(cell_type_custom_final)]
    merged <- merge(annotation_dt, morph_dt, by = "cell_id", all.x = FALSE, all.y = FALSE)
    merged <- normalize_morphology_coord_columns_mainv2(merged)

    summary_dt <- merged[, .(
      n = .N,
      proportion = .N / nrow(merged),
      area_median = median(area, na.rm = TRUE),
      area_p10 = quantile(area, 0.1, na.rm = TRUE),
      area_p90 = quantile(area, 0.9, na.rm = TRUE),
      eccentricity_median = median(eccentricity, na.rm = TRUE),
      solidity_median = median(solidity, na.rm = TRUE),
      circularity_median = median(circularity, na.rm = TRUE)
    ), by = cell_type_custom_final]

    per_cell_path <- file.path(eval_dir, paste0(method_name, "_morphology_per_cell.csv"))
    summary_path <- file.path(eval_dir, paste0(method_name, "_morphology_by_celltype.csv"))
    fwrite(merged, per_cell_path)
    fwrite(summary_dt, summary_path)
    plot_morphology_method_outputs(merged, method_name, file.path(eval_dir, paste0(method_name, "_morphology_plots")))

    summary_rows[[idx]] <- copy(summary_dt)[, method := method_name]
    results[[idx]] <- data.table(
      method = method_name,
      n_cells = nrow(merged),
      n_cell_types = uniqueN(merged$cell_type_custom_final),
      median_area = median(merged$area, na.rm = TRUE)
    )
    idx <- idx + 1L
  }

  summary_all <- if (length(summary_rows)) rbindlist(summary_rows, fill = TRUE) else data.table()
  if (nrow(summary_all)) {
    fwrite(summary_all, file.path(eval_dir, "morphology_summary_all_methods.csv"))
    plot_morphology_cross_method_summary(summary_all, eval_dir)
  }

  manifest_path <- file.path(eval_dir, "morphology_manifest.csv")
  manifest_dt <- if (length(results)) rbindlist(results, fill = TRUE) else data.table(
    method = character(),
    n_cells = integer(),
    n_cell_types = integer(),
    median_area = numeric()
  )
  fwrite(manifest_dt, manifest_path)
  manifest_path
}
