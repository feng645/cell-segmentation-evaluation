suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
})

viz_theme_academic <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, color = "#4d4d4d"),
      axis.title = element_text(face = "plain"),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      strip.background = element_rect(fill = "#f3f4f6", color = "#c7cbd1", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.35)
    )
}

read_json_any <- function(path) {
  if (grepl("\\.gz$", path)) {
    jsonlite::fromJSON(paste(readLines(gzfile(path), warn = FALSE), collapse = ""), simplifyVector = FALSE)
  } else {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
}

build_celltype_palette <- function(cell_types) {
  cell_types <- unique(as.character(cell_types))
  cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types) & cell_types != "Unknown"]
  cell_types <- sort(cell_types)
  # Colorblind-safe qualitative palette, drawn from Paul Tol "bright" + "vibrant"
  # + "muted" (https://personal.sron.nl/~pault/). Ordered for maximum pairwise
  # distinguishability under CVD simulation; Unknown reserved for neutral grey.
  base_colors <- c(
    "#4477AA",  # blue
    "#EE6677",  # red
    "#228833",  # green
    "#AA3377",  # purple
    "#EE7733",  # orange
    "#009988",  # teal
    "#66CCEE",  # cyan
    "#CC6677",  # rose
    "#CCBB44",  # yellow
    "#332288",  # indigo
    "#999933",  # olive
    "#882255"   # wine
  )
  palette <- if (length(cell_types)) setNames(rep_len(base_colors, length(cell_types)), cell_types) else character(0)
  c(palette, Unknown = "#c7cdd3")
}

extract_polygon_paths <- function(geojson_path, method_name) {
  gj <- read_json_any(geojson_path)
  rows <- list()
  row_idx <- 1L

  for (feature_idx in seq_along(gj$features)) {
    feature <- gj$features[[feature_idx]]
    cell_id <- feature$id %||% feature$properties$cell_id %||% feature$properties$cell %||% paste0("cell_", feature_idx)
    geom <- feature$geometry
    if (is.null(geom$type) || is.null(geom$coordinates)) next

    polygon_list <- list()
    if (toupper(geom$type) == "POLYGON") {
      polygon_list <- list(geom$coordinates)
    } else if (toupper(geom$type) == "MULTIPOLYGON") {
      polygon_list <- geom$coordinates
    } else {
      next
    }

    for (poly_idx in seq_along(polygon_list)) {
      rings <- polygon_list[[poly_idx]]
      if (length(rings) == 0L) next
      ring <- rings[[1]]
      mat <- tryCatch(do.call(rbind, ring), error = function(...) NULL)
      if (is.null(mat) || ncol(mat) < 2 || nrow(mat) < 3) next
      mat <- matrix(as.numeric(mat[, 1:2, drop = FALSE]), ncol = 2)
      rows[[row_idx]] <- data.table(
        method = method_name,
        cell_id = as.character(cell_id),
        feature_id = feature_idx,
        poly_id = poly_idx,
        order = seq_len(nrow(mat)),
        x = mat[, 1],
        y = mat[, 2]
      )
      row_idx <- row_idx + 1L
    }
  }

  if (length(rows) == 0L) {
    return(data.table(method = character(), cell_id = character(), feature_id = integer(), poly_id = integer(), order = integer(), x = numeric(), y = numeric()))
  }

  rbindlist(rows, fill = TRUE)
}

load_manifest <- function(output_root) {
  jsonlite::read_json(file.path(output_root, "manifest.json"), simplifyVector = FALSE)
}

manifest_methods_dt <- function(manifest) {
  if (is.null(manifest$methods) || length(manifest$methods) == 0L) {
    return(data.table())
  }
  rbindlist(lapply(manifest$methods, as.data.table), fill = TRUE)
}

sample_transcripts_for_plot <- function(assignments_path, n_max = 40000L) {
  dt <- data.table::fread(assignments_path)
  if (nrow(dt) <= n_max) return(dt)
  set.seed(42)
  dt[sample.int(nrow(dt), n_max)]
}

compute_polygon_centroids <- function(geojson_path, method_name) {
  dt <- extract_polygon_paths(geojson_path, method_name)
  if (!nrow(dt)) {
    return(data.table(method = character(), cell_id = character(), x = numeric(), y = numeric()))
  }
  dt[, .(
    x = mean(range(x, na.rm = TRUE)),
    y = mean(range(y, na.rm = TRUE))
  ), by = .(method, cell_id)]
}

build_assignment_agreement_plot <- function(assignments_path, out_path) {
  dt <- fread(assignments_path)
  method_cols <- grep("^cell_id_", names(dt), value = TRUE)
  methods <- sub("^cell_id_", "", method_cols)
  pairs <- CJ(method_a = methods, method_b = methods)
  pairs[, agreement := 0]

  for (i in seq_len(nrow(pairs))) {
    col_a <- paste0("cell_id_", pairs$method_a[[i]])
    col_b <- paste0("cell_id_", pairs$method_b[[i]])
    a_assigned <- dt[[col_a]] != "-1"
    b_assigned <- dt[[col_b]] != "-1"
    union_n <- sum(a_assigned | b_assigned)
    intersect_n <- sum(a_assigned & b_assigned)
    pairs$agreement[[i]] <- if (union_n > 0) intersect_n / union_n else NA_real_
  }

  p <- ggplot(pairs, aes(x = method_b, y = method_a, fill = agreement)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.2f", agreement)), size = 3) +
    scale_fill_gradient(low = "#edf2f7", high = "#1f4e79", na.value = "#f4f5f7", limits = c(0, 1)) +
    labs(
      title = "Assignment Overlap Between Methods",
      subtitle = "Transcript-level Jaccard overlap based on assigned-versus-background calls",
      x = NULL,
      y = NULL,
      fill = "Jaccard"
    ) +
    viz_theme_academic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(out_path, p, width = 6.8, height = 5.8, dpi = 320)
}

build_assignment_summary_plot <- function(summary_dt, out_path) {
  summary_long <- melt(
    copy(summary_dt)[, .(method, assigned_transcripts, background_transcripts)],
    id.vars = "method",
    variable.name = "category",
    value.name = "n_transcripts"
  )
  summary_long[, category := factor(category,
                                    levels = c("assigned_transcripts", "background_transcripts"),
                                    labels = c("Assigned", "Background"))]

  p <- ggplot(summary_long, aes(x = reorder(method, -n_transcripts), y = n_transcripts, fill = category)) +
    geom_col(width = 0.72, position = "stack", color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("Assigned" = "#315c78", "Background" = "#d6d9dd")) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Transcript Assignment Summary",
      subtitle = "Assigned versus background transcripts across segmentation methods",
      x = NULL,
      y = "Transcripts",
      fill = NULL
    ) +
    viz_theme_academic()

  ggsave(out_path, p, width = 7.8, height = 4.8, dpi = 320)
}

build_cell_count_plot <- function(cell_manifest_dt, out_path) {
  plot_dt <- melt(
    copy(cell_manifest_dt)[status == "ready", .(method, n_cells_raw, n_cells_eval)],
    id.vars = "method",
    variable.name = "cell_set",
    value.name = "n_cells"
  )
  plot_dt[, cell_set := factor(cell_set, levels = c("n_cells_raw", "n_cells_eval"), labels = c("Raw", "Eval"))]

  p <- ggplot(plot_dt, aes(x = reorder(method, -n_cells), y = n_cells, fill = cell_set)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62, color = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("Raw" = "#8b9fb1", "Eval" = "#1f4e79")) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Cell Yield by Method",
      subtitle = "Raw cells and cells retained for evaluation",
      x = NULL,
      y = "Cells",
      fill = NULL
    ) +
    viz_theme_academic()

  ggsave(out_path, p, width = 7.8, height = 4.8, dpi = 320)
}

build_qc_distribution_plot <- function(cell_table_manifest_dt, out_path) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(cell_table_manifest_dt))) {
    row <- cell_table_manifest_dt[i]
    if (!identical(row$status, "ready")) next
    dt <- data.table::fread(row$cell_table_path)
    if (!nrow(dt)) next
    dt[, method := row$method]
    rows[[idx]] <- dt[, .(method, n_transcripts, n_genes, pct_mt)]
    idx <- idx + 1L
  }
  qc_dt <- rbindlist(rows, fill = TRUE)

  long_dt <- melt(
    qc_dt,
    id.vars = "method",
    measure.vars = c("n_transcripts", "n_genes", "pct_mt"),
    variable.name = "metric",
    value.name = "value"
  )
  long_dt[, value := as.numeric(value)]
  long_dt[metric %in% c("n_transcripts", "n_genes"), value := log10(value + 1)]
  long_dt[, metric := factor(metric,
                             levels = c("n_transcripts", "n_genes", "pct_mt"),
                             labels = c("log10(UMIs + 1)", "log10(Genes + 1)", "Percent MT"))]

  p <- ggplot(long_dt[is.finite(value)], aes(x = method, y = value, fill = method)) +
    geom_violin(scale = "width", color = NA, alpha = 0.65) +
    geom_boxplot(width = 0.12, outlier.size = 0.15, alpha = 0.85, color = "#2f2f2f") +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = "Per-cell Quality Distributions",
      subtitle = "Distribution of count and QC summaries across methods",
      x = NULL,
      y = NULL
    ) +
    viz_theme_academic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(out_path, p, width = 8.2, height = 9.6, dpi = 320)
}

build_cell_area_plot <- function(cell_table_manifest_dt, out_path) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(cell_table_manifest_dt))) {
    row <- cell_table_manifest_dt[i]
    if (!identical(row$status, "ready")) next
    dt <- fread(row$cell_table_path)
    if (!nrow(dt) || !"cell_area" %in% names(dt)) next
    dt <- dt[is.finite(cell_area) & cell_area > 0]
    if (!nrow(dt)) next
    dt[, method := row$method]
    rows[[idx]] <- dt[, .(method, cell_area)]
    idx <- idx + 1L
  }
  area_dt <- rbindlist(rows, fill = TRUE)
  if (!nrow(area_dt)) return(invisible(NULL))
  area_dt[, log_area := log10(cell_area)]

  p <- ggplot(area_dt, aes(x = method, y = log_area, fill = method)) +
    geom_violin(scale = "width", color = NA, alpha = 0.68) +
    geom_boxplot(width = 0.12, outlier.size = 0.15, alpha = 0.85, color = "#2f2f2f") +
    scale_fill_brewer(palette = "Pastel2", guide = "none") +
    labs(
      title = "Cell Area Distribution by Method",
      subtitle = "Polygon-derived cell area shown on a log10 scale",
      x = NULL,
      y = expression(log[10]("Area"))
    ) +
    viz_theme_academic() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(out_path, p, width = 7.8, height = 4.9, dpi = 320)
}

build_gene_composition_plot <- function(cell_table_manifest_dt, output_root, out_path, top_n = 20L) {
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(cell_table_manifest_dt))) {
    row <- cell_table_manifest_dt[i]
    if (!identical(row$status, "ready")) next
    counts_path <- file.path(output_root, "expression", paste0(row$method, "_counts_long.csv.gz"))
    if (!file.exists(counts_path)) next
    dt <- fread(counts_path)
    if (!nrow(dt)) next
    gene_dt <- dt[, .(count = sum(count)), by = gene]
    gene_dt[, method := row$method]
    rows[[idx]] <- gene_dt
    idx <- idx + 1L
  }
  gene_dt <- rbindlist(rows, fill = TRUE)
  if (!nrow(gene_dt)) return(invisible(NULL))

  top_genes <- gene_dt[, .(count = sum(count)), by = gene][order(-count)][seq_len(min(top_n, .N)), gene]
  plot_dt <- gene_dt[gene %in% top_genes]
  plot_dt[, fraction := count / sum(count), by = method]

  p <- ggplot(plot_dt, aes(x = gene, y = method, fill = fraction)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient(low = "#f7fbff", high = "#084594") +
    labs(
      title = "Top-gene Composition by Method",
      subtitle = "Relative abundance of the most prevalent genes across method-specific assigned counts",
      x = NULL,
      y = NULL,
      fill = "Fraction"
    ) +
    viz_theme_academic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(out_path, p, width = 9.0, height = 4.8, dpi = 320)
}

build_centroid_overview_plot <- function(manifest, out_path, sample_n = 2500L) {
  methods_df <- manifest_methods_dt(manifest)
  methods_df <- methods_df[status == "ready" & mode == "polygon" & !is.na(polygon_path) & nzchar(polygon_path)]
  rows <- list()
  idx <- 1L
  for (i in seq_len(nrow(methods_df))) {
    row <- methods_df[i]
    dt <- compute_polygon_centroids(row$polygon_path, row$name)
    if (!nrow(dt)) next
    if (nrow(dt) > sample_n) {
      set.seed(900 + i)
      dt <- dt[sample.int(nrow(dt), sample_n)]
    }
    rows[[idx]] <- dt
    idx <- idx + 1L
  }
  cent_dt <- rbindlist(rows, fill = TRUE)
  if (!nrow(cent_dt)) return(invisible(NULL))

  p <- ggplot(cent_dt, aes(x = x, y = y)) +
    geom_point(size = 0.22, alpha = 0.55, color = "#1f4e79") +
    coord_equal() +
    facet_wrap(~method, ncol = 2) +
    labs(
      title = "Cell Centroid Overview",
      subtitle = "Sampled polygon centroids by segmentation method",
      x = "x (um)",
      y = "y (um)"
    ) +
    viz_theme_academic()

  ggsave(out_path, p, width = 10.4, height = 8.2, dpi = 320)
}

build_polygon_overview_plot <- function(manifest, assignments_path, out_path) {
  methods_df <- manifest_methods_dt(manifest)
  methods_df <- methods_df[status == "ready" & mode == "polygon" & !is.na(polygon_path) & nzchar(polygon_path)]
  if (!nrow(methods_df)) {
    warning("No polygon methods available for overview plot.")
    return(invisible(NULL))
  }

  transcripts_dt <- sample_transcripts_for_plot(assignments_path, n_max = 35000L)
  panel_list <- list()
  idx <- 1L
  for (i in seq_len(nrow(methods_df))) {
    row <- methods_df[i]
    sub_dt <- extract_polygon_paths(row$polygon_path, row$name)
    if (!nrow(sub_dt)) next
    sampled_cells <- unique(sub_dt$cell_id)
    if (length(sampled_cells) > 250L) {
      set.seed(100 + i)
      sampled_cells <- sample(sampled_cells, 250L)
      sub_dt <- sub_dt[cell_id %in% sampled_cells]
    }
    panel_list[[idx]] <- 
      ggplot() +
      geom_point(
        data = transcripts_dt,
        aes(x = x, y = y),
        color = "#b8bec5",
        size = 0.10,
        alpha = 0.30
      ) +
      geom_path(
        data = sub_dt,
        aes(x = x, y = y, group = interaction(feature_id, poly_id)),
        color = "#1f4e79",
        linewidth = 0.18,
        alpha = 0.75
      ) +
      coord_equal() +
      labs(
        title = row$name,
        x = "x (um)",
        y = "y (um)"
      ) +
      viz_theme_academic()
    idx <- idx + 1L
  }

  if (length(panel_list) == 0L) {
    warning("No polygon coordinates available for overview plot.")
    return(invisible(NULL))
  }

  combined <- wrap_plots(panel_list, ncol = 2) +
    plot_annotation(
      title = "Boundary Overview Across Methods",
      subtitle = "Sampled polygon outlines overlaid on the canonical transcript universe"
    )

  ggsave(out_path, combined, width = 10.6, height = 8.4, dpi = 320)
}

build_cellpose_comparison_plot <- function(assignments_path, manifest, out_path) {
  methods_df <- manifest_methods_dt(manifest)
  cellpose_row <- methods_df[name == "cellpose" & status == "ready" & !is.na(polygon_path) & nzchar(polygon_path)]
  if (!nrow(cellpose_row)) {
    warning("Cellpose polygon artifact not available for comparison plot.")
    return(invisible(NULL))
  }

  transcripts_dt <- sample_transcripts_for_plot(assignments_path, n_max = 45000L)
  transcripts_dt[, cellpose_status := ifelse(cell_id_cellpose == "-1", "Background", "Assigned")]
  poly_dt <- extract_polygon_paths(cellpose_row$polygon_path[[1]], "cellpose")
  sampled_cells <- unique(poly_dt$cell_id)
  if (length(sampled_cells) > 260L) {
    set.seed(84)
    sampled_cells <- sample(sampled_cells, 260L)
    poly_dt <- poly_dt[cell_id %in% sampled_cells]
  }

  p <- ggplot() +
    geom_point(
      data = transcripts_dt,
      aes(x = x, y = y, color = cellpose_status),
      size = 0.12,
      alpha = 0.38
    ) +
    geom_path(
      data = poly_dt,
      aes(x = x, y = y, group = interaction(feature_id, poly_id)),
      color = "#111827",
      linewidth = 0.20,
      alpha = 0.82
    ) +
    scale_color_manual(values = c("Assigned" = "#006d77", "Background" = "#c7cdd3")) +
    coord_equal() +
    labs(
      title = "Cellpose Polygon Overlay",
      subtitle = "User-supplied Cellpose polygons over canonical transcripts",
      x = "x (um)",
      y = "y (um)",
      color = NULL
    ) +
    viz_theme_academic()

  ggsave(out_path, p, width = 7.2, height = 6.6, dpi = 320)
}

load_annotation_per_method <- function(annotation_manifest_path) {
  if (is.null(annotation_manifest_path) || !file.exists(annotation_manifest_path)) return(list())
  am <- fread(annotation_manifest_path)
  out <- list()
  for (i in seq_len(nrow(am))) {
    p <- am$annotation_path[[i]]
    if (is.null(p) || !file.exists(p)) next
    dt <- fread(p, select = c("cell_id", "cell_type_custom_final"))
    dt[, cell_id := as.character(cell_id)]
    dt[, cell_type_custom_final := ifelse(is.na(cell_type_custom_final) | !nzchar(cell_type_custom_final), "Unknown", cell_type_custom_final)]
    out[[am$method[[i]]]] <- dt
  }
  out
}

load_marker_gene_map <- function(markers_path) {
  if (is.null(markers_path) || !file.exists(markers_path)) return(data.table(gene = character(), marker_cell_type = character()))
  m <- fread(markers_path)
  if (!all(c("cell_type", "marker") %in% names(m))) return(data.table(gene = character(), marker_cell_type = character()))
  m <- unique(m[, .(gene = as.character(marker), marker_cell_type = as.character(cell_type))])
  m <- m[!is.na(gene) & nzchar(gene) & !is.na(marker_cell_type) & nzchar(marker_cell_type)]
  m
}

build_zoom_region_plots <- function(manifest, assignments_path, out_dir, zoom_regions,
                                    transcript_sample_n = 50000L,
                                    annotation_manifest_path = NULL,
                                    markers_path = NULL) {
  methods_df <- manifest_methods_dt(manifest)
  methods_df <- methods_df[status == "ready" & mode == "polygon" & !is.na(polygon_path) & nzchar(polygon_path)]
  if (!length(zoom_regions) || !nrow(methods_df)) return(invisible(NULL))

  transcripts_dt <- sample_transcripts_for_plot(assignments_path, n_max = transcript_sample_n)
  annot_by_method <- load_annotation_per_method(annotation_manifest_path)
  marker_map <- load_marker_gene_map(markers_path)

  all_celltypes <- unique(c(
    unlist(lapply(annot_by_method, function(d) unique(d$cell_type_custom_final))),
    marker_map$marker_cell_type
  ))
  palette <- build_celltype_palette(all_celltypes)

  if (nrow(marker_map)) {
    tx <- merge(transcripts_dt, marker_map, by = "gene", all.x = TRUE)
    tx[is.na(marker_cell_type), marker_cell_type := NA_character_]
  } else {
    tx <- copy(transcripts_dt)
    tx[, marker_cell_type := NA_character_]
  }

  for (i in seq_along(zoom_regions)) {
    zr <- zoom_regions[[i]]
    bbox <- zr$bbox %||% NULL
    name <- zr$name %||% sprintf("region_%02d", i)
    if (is.null(bbox) || length(bbox) != 4L) next

    tx_zoom <- tx[x >= bbox[1] & x <= bbox[3] & y >= bbox[2] & y <= bbox[4]]
    tx_bg <- tx_zoom[is.na(marker_cell_type)]
    tx_marker <- tx_zoom[!is.na(marker_cell_type)]

    panels_markers <- list()
    panels_annot <- list()
    idx <- 1L
    for (j in seq_len(nrow(methods_df))) {
      row <- methods_df[j]
      poly_dt <- extract_polygon_paths(row$polygon_path, row$name)
      if (!nrow(poly_dt)) next
      poly_dt <- poly_dt[x >= bbox[1] & x <= bbox[3] & y >= bbox[2] & y <= bbox[4]]

      # Variant A: transcripts colored by their marker cell type
      p_m <- ggplot() +
        geom_point(data = tx_bg, aes(x = x, y = y), color = "#d0d5da", size = 0.18, alpha = 0.35) +
        geom_point(data = tx_marker, aes(x = x, y = y, color = marker_cell_type), size = 0.35, alpha = 0.85) +
        geom_path(
          data = poly_dt,
          aes(x = x, y = y, group = interaction(feature_id, poly_id)),
          color = "#2b2b2b",
          linewidth = 0.22,
          alpha = 0.75
        ) +
        scale_color_manual(values = palette, na.value = "#d0d5da", name = "Marker cell type") +
        coord_equal(xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)]) +
        labs(title = row$name, x = "x (um)", y = "y (um)") +
        viz_theme_academic() +
        guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))
      panels_markers[[idx]] <- p_m

      # Variant B: polygons colored by cell's annotated cell type (same palette)
      annot_dt <- annot_by_method[[row$name]]
      if (!is.null(annot_dt)) {
        poly_annot <- merge(poly_dt, annot_dt, by = "cell_id", all.x = TRUE)
        poly_annot[is.na(cell_type_custom_final), cell_type_custom_final := "Unknown"]
      } else {
        poly_annot <- copy(poly_dt)
        poly_annot[, cell_type_custom_final := "Unknown"]
      }
      p_a <- ggplot() +
        geom_polygon(
          data = poly_annot,
          aes(x = x, y = y,
              group = interaction(feature_id, poly_id),
              fill = cell_type_custom_final),
          color = "#2b2b2b",
          linewidth = 0.15,
          alpha = 0.9
        ) +
        scale_fill_manual(values = palette, na.value = "#c7cdd3", name = "Cell type") +
        coord_equal(xlim = bbox[c(1, 3)], ylim = bbox[c(2, 4)]) +
        labs(title = row$name, x = "x (um)", y = "y (um)") +
        viz_theme_academic() +
        guides(fill = guide_legend(override.aes = list(alpha = 1)))
      panels_annot[[idx]] <- p_a

      idx <- idx + 1L
    }
    if (length(panels_markers) == 0L) next

    sub_bbox <- sprintf("Bounding box [%.1f, %.1f] to [%.1f, %.1f]", bbox[1], bbox[2], bbox[3], bbox[4])
    combined_m <- wrap_plots(panels_markers, ncol = 2, guides = "collect") +
      plot_annotation(title = paste("Zoom-in Boundary + Marker Transcripts:", name), subtitle = sub_bbox)
    combined_a <- wrap_plots(panels_annot, ncol = 2, guides = "collect") +
      plot_annotation(title = paste("Zoom-in Annotated Boundaries:", name), subtitle = sub_bbox)

    ggsave(file.path(out_dir, paste0("zoom_", name, "_by_markers.png")), combined_m, width = 11.2, height = 8.8, dpi = 320)
    ggsave(file.path(out_dir, paste0("zoom_", name, "_by_annotation.png")), combined_a, width = 11.2, height = 8.8, dpi = 320)
  }
}

run_visualization_eval <- function(output_root, visual_cfg = list(), markers_path = NULL) {
  manifest <- load_manifest(output_root)
  viz_dir <- file.path(output_root, "visualization")
  dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)

  assignment_summary_dt <- fread(file.path(output_root, "assignments", "assignment_summary_by_method.csv"))
  cell_manifest_dt <- fread(file.path(output_root, "cell_tables", "method_cell_output_manifest.csv"))
  assignments_path <- file.path(output_root, "assignments", "transcripts_with_assignments.csv.gz")
  annotation_manifest_path <- file.path(output_root, "annotation", "annotation_manifest.csv")
  transcript_sample_n <- as.integer(visual_cfg$transcript_sample_n %||% 35000L)
  centroid_sample_n <- as.integer(visual_cfg$centroid_sample_n %||% 2500L)

  build_assignment_summary_plot(assignment_summary_dt, file.path(viz_dir, "assignment_summary.png"))
  build_assignment_agreement_plot(assignments_path, file.path(viz_dir, "assignment_agreement_heatmap.png"))
  build_cell_count_plot(cell_manifest_dt, file.path(viz_dir, "cell_yield_summary.png"))
  build_qc_distribution_plot(cell_manifest_dt, file.path(viz_dir, "qc_distributions.png"))
  build_cell_area_plot(cell_manifest_dt, file.path(viz_dir, "cell_area_distributions.png"))
  build_gene_composition_plot(cell_manifest_dt, output_root, file.path(viz_dir, "top_gene_composition_heatmap.png"))
  build_centroid_overview_plot(manifest, file.path(viz_dir, "centroid_overview.png"), sample_n = centroid_sample_n)
  build_polygon_overview_plot(manifest, assignments_path, file.path(viz_dir, "boundary_overview.png"))
  build_cellpose_comparison_plot(assignments_path, manifest, file.path(viz_dir, "cellpose_overlay.png"))
  build_zoom_region_plots(
    manifest = manifest,
    assignments_path = assignments_path,
    out_dir = viz_dir,
    zoom_regions = visual_cfg$zoom_regions %||% list(),
    transcript_sample_n = transcript_sample_n,
    annotation_manifest_path = annotation_manifest_path,
    markers_path = markers_path
  )

  invisible(viz_dir)
}
