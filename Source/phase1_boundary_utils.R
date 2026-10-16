suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(tiff)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

read_json_any <- function(path) {
  if (grepl("\\.gz$", path)) {
    jsonlite::fromJSON(paste(readLines(gzfile(path), warn = FALSE), collapse = ""), simplifyVector = FALSE)
  } else {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
}

write_json_any <- function(payload, path) {
  if (grepl("\\.gz$", path)) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"), con = con)
  } else {
    jsonlite::write_json(payload, path, auto_unbox = TRUE, pretty = FALSE, null = "null")
  }
}

resolve_feature_cell_id <- function(feature) {
  cid <- feature$id %||% feature$properties$cell_id %||% feature$properties$cell %||% feature$properties$id
  if (is.null(cid)) {
    stop("Could not resolve cell id from feature.")
  }
  as.character(cid)
}

ring_matrix <- function(coords) {
  mat <- tryCatch(do.call(rbind, coords), error = function(...) NULL)
  if (is.null(mat) || ncol(mat) < 2) return(NULL)
  mat <- matrix(as.numeric(mat[, 1:2, drop = FALSE]), ncol = 2)
  if (nrow(mat) < 3 || any(is.na(mat))) return(NULL)
  mat
}

close_ring_mat <- function(mat) {
  if (nrow(mat) < 2) return(mat)
  if (!identical(unname(mat[1, ]), unname(mat[nrow(mat), ]))) {
    mat <- rbind(mat, mat[1, ])
  }
  mat
}

polygon_bbox <- function(mat) {
  c(min(mat[, 1]), min(mat[, 2]), max(mat[, 1]), max(mat[, 2]))
}

polygon_area_mat <- function(mat) {
  mat <- close_ring_mat(mat)
  abs(sum(mat[-1, 1] * mat[-nrow(mat), 2] - mat[-1, 2] * mat[-nrow(mat), 1]) / 2)
}

polygon_signed_area_mat <- function(mat) {
  mat <- close_ring_mat(mat)
  sum(mat[-1, 1] * mat[-nrow(mat), 2] - mat[-1, 2] * mat[-nrow(mat), 1]) / 2
}

polygon_centroid_mat <- function(mat) {
  mat <- close_ring_mat(mat)
  area2 <- sum(mat[-1, 1] * mat[-nrow(mat), 2] - mat[-1, 2] * mat[-nrow(mat), 1])
  if (!is.finite(area2) || abs(area2) < 1e-12) {
    return(c(mean(mat[, 1], na.rm = TRUE), mean(mat[, 2], na.rm = TRUE)))
  }
  cx <- sum((mat[-nrow(mat), 1] + mat[-1, 1]) * (mat[-nrow(mat), 1] * mat[-1, 2] - mat[-1, 1] * mat[-nrow(mat), 2])) / (3 * area2)
  cy <- sum((mat[-nrow(mat), 2] + mat[-1, 2]) * (mat[-nrow(mat), 1] * mat[-1, 2] - mat[-1, 1] * mat[-nrow(mat), 2])) / (3 * area2)
  c(cx, cy)
}

coords_from_mat <- function(mat) {
  lapply(seq_len(nrow(mat)), function(i) list(unname(mat[i, 1]), unname(mat[i, 2])))
}

polygon_record_from_mats <- function(cell_id, ring_mats) {
  ring_mats <- Filter(Negate(is.null), ring_mats)
  if (!length(ring_mats)) return(NULL)

  ext_ring <- close_ring_mat(ring_mats[[1]])
  hole_rings <- if (length(ring_mats) > 1L) lapply(ring_mats[-1], close_ring_mat) else list()
  ext_area <- polygon_area_mat(ext_ring)
  ext_centroid <- polygon_centroid_mat(ext_ring)
  hole_areas <- vapply(hole_rings, polygon_area_mat, numeric(1))
  hole_centroids <- lapply(hole_rings, polygon_centroid_mat)
  area <- ext_area - sum(hole_areas, na.rm = TRUE)
  if (!is.finite(area) || area <= 0) return(NULL)

  centroid_num <- ext_centroid * ext_area
  if (length(hole_centroids)) {
    for (idx in seq_along(hole_centroids)) {
      centroid_num <- centroid_num - hole_centroids[[idx]] * hole_areas[[idx]]
    }
  }
  centroid <- centroid_num / area

  list(
    cell_id = as.character(cell_id),
    vx = ext_ring[, 1],
    vy = ext_ring[, 2],
    hole_vx = lapply(hole_rings, function(mat) mat[, 1]),
    hole_vy = lapply(hole_rings, function(mat) mat[, 2]),
    bbox = polygon_bbox(ext_ring),
    area = area,
    centroid = centroid,
    geometry = lapply(c(list(ext_ring), hole_rings), coords_from_mat)
  )
}

polygon_record_from_rings <- function(cell_id, rings) {
  polygon_record_from_mats(cell_id, lapply(rings, ring_matrix))
}

flatten_geojson_polygons <- function(gj) {
  res <- list()
  for (feature in gj$features) {
    cid <- resolve_feature_cell_id(feature)
    geom <- feature$geometry
    if (is.null(geom$type) || is.null(geom$coordinates)) next
    polygons <- list()
    if (toupper(geom$type) == "POLYGON") {
      polygons <- list(geom$coordinates)
    } else if (toupper(geom$type) == "MULTIPOLYGON") {
      polygons <- geom$coordinates
    } else {
      next
    }

    for (poly_rings in polygons) {
      poly <- polygon_record_from_rings(cid, poly_rings)
      if (is.null(poly)) next
      res[[length(res) + 1L]] <- poly
    }
  }
  res
}

read_bidcell_affine_summary <- function(path) {
  dt <- data.table::fread(path, header = FALSE, sep = "\t", col.names = c("key", "value"))
  vals <- as.list(dt$value)
  names(vals) <- dt$key
  numeric_keys <- setdiff(names(vals), "origin")
  for (k in numeric_keys) vals[[k]] <- as.numeric(vals[[k]])
  vals$origin <- isTRUE(as.logical(vals$origin))
  vals
}

transform_bidcell_polygon <- function(poly, transform_cfg, method_cfg) {
  transform_type <- transform_cfg$type %||% "identity"
  if (identical(transform_type, "identity")) {
    return(poly)
  }
  if (!identical(transform_type, "bidcell_crop_tiff")) {
    stop("Unsupported BIDCell transform type: ", transform_type)
  }

  affine <- read_bidcell_affine_summary(method_cfg$input$affine_summary)
  seg_img <- tiff::readTIFF(method_cfg$input$seg_path, as.is = TRUE)
  img_width <- ncol(seg_img)
  img_height <- nrow(seg_img)
  x_scale <- affine$size_x / img_width
  y_scale <- affine$size_y / img_height
  min_x <- affine$min_x + (affine$global_shift_x %||% 0)
  min_y <- affine$min_y + (affine$global_shift_y %||% 0)
  ts_x <- affine$scale_ts_x %||% 1
  ts_y <- affine$scale_ts_y %||% 1

  out <- poly
  ring_mats <- c(
    list(cbind(out$vx, out$vy)),
    Map(cbind, out$hole_vx %||% list(), out$hole_vy %||% list())
  )
  transformed_rings <- lapply(ring_mats, function(mat) {
    mat <- close_ring_mat(mat)
    out_mat <- mat
    out_mat[, 1] <- min_x + out_mat[, 1] * x_scale * ts_x
    out_mat[, 2] <- min_y + out_mat[, 2] * y_scale * ts_y
    out_mat
  })
  polygon_record_from_mats(out$cell_id, transformed_rings)
}

summarize_polygons_by_cell <- function(polygons) {
  area_dt <- rbindlist(lapply(polygons, function(poly) data.table(cell_id = poly$cell_id, area = poly$area)))
  area_dt <- area_dt[, .(area = sum(area, na.rm = TRUE)), by = cell_id]

  centroid_dt <- rbindlist(lapply(polygons, function(poly) data.table(
    cell_id = poly$cell_id,
    area = poly$area,
    x = poly$centroid[[1]],
    y = poly$centroid[[2]]
  )))
  centroid_dt <- centroid_dt[, .(
    x = sum(x * area, na.rm = TRUE) / sum(area, na.rm = TRUE),
    y = sum(y * area, na.rm = TRUE) / sum(area, na.rm = TRUE)
  ), by = cell_id]

  list(area_dt = area_dt, centroid_dt = centroid_dt)
}

geojson_from_polygons <- function(polygons, method_name = NULL, dataset_id = NULL) {
  by_cell <- split(polygons, vapply(polygons, `[[`, character(1), "cell_id"))
  features <- lapply(names(by_cell), function(cell_id) {
    cell_polys <- by_cell[[cell_id]]
    coords <- lapply(cell_polys, `[[`, "geometry")
    geom <- if (length(coords) == 1L) {
      list(type = "Polygon", coordinates = coords[[1]])
    } else {
      list(type = "MultiPolygon", coordinates = coords)
    }
    list(
      type = "Feature",
      id = cell_id,
      properties = Filter(Negate(is.null), list(
        cell_id = cell_id,
        method = method_name,
        dataset_id = dataset_id
      )),
      geometry = geom
    )
  })
  list(type = "FeatureCollection", features = features)
}

load_geojson_method <- function(method_cfg, polygons_dir) {
  gj <- read_json_any(method_cfg$input$path)
  polygons <- flatten_geojson_polygons(gj)
  if (length(polygons) == 0L) {
    return(list(name = method_cfg$name, status = "empty", reason = "No usable polygons were loaded."))
  }
  polygon_path <- file.path(polygons_dir, paste0(method_cfg$name, ".geojson.gz"))
  write_json_any(geojson_from_polygons(polygons, method_name = method_cfg$name, dataset_id = method_cfg$dataset_id %||% NULL), polygon_path)

  summary <- summarize_polygons_by_cell(polygons)
  list(
    name = method_cfg$name,
    status = "ready",
    mode = "polygon",
    loader = method_cfg$input$loader,
    polygons = polygons,
    polygon_path = polygon_path,
    n_cells_raw = uniqueN(summary$area_dt$cell_id),
    cell_areas = setNames(summary$area_dt$area, summary$area_dt$cell_id),
    cell_centroids = summary$centroid_dt
  )
}

load_xenium_vertex_csv_method <- function(method_cfg, polygons_dir) {
  dt <- data.table::fread(method_cfg$input$path)
  required <- c("cell_id", "vertex_x", "vertex_y")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing Xenium boundary columns: ", paste(missing, collapse = ", "))
  }

  polygons <- dt[, {
    mat <- as.matrix(.SD[, .(vertex_x, vertex_y)])
    if (nrow(mat) < 3) return(NULL)
    poly <- polygon_record_from_mats(as.character(cell_id[1]), list(mat))
    if (is.null(poly)) return(NULL)
    list(list(poly))
  }, by = cell_id]$V1

  if (length(polygons) == 0L) {
    return(list(name = method_cfg$name, status = "empty", reason = "No usable polygons were loaded."))
  }

  polygon_path <- file.path(polygons_dir, paste0(method_cfg$name, ".geojson.gz"))
  write_json_any(geojson_from_polygons(polygons, method_name = method_cfg$name, dataset_id = method_cfg$dataset_id %||% NULL), polygon_path)

  summary <- summarize_polygons_by_cell(polygons)
  list(
    name = method_cfg$name,
    status = "ready",
    mode = "polygon",
    loader = method_cfg$input$loader,
    polygons = polygons,
    polygon_path = polygon_path,
    n_cells_raw = uniqueN(summary$area_dt$cell_id),
    cell_areas = setNames(summary$area_dt$area, summary$area_dt$cell_id),
    cell_centroids = summary$centroid_dt
  )
}

load_bidcell_geojson_method <- function(method_cfg, polygons_dir) {
  gj <- read_json_any(method_cfg$input$path)
  raw_polygons <- flatten_geojson_polygons(gj)
  transform_cfg <- method_cfg$input$transform %||% list(type = "identity")
  polygons <- lapply(raw_polygons, transform_bidcell_polygon, transform_cfg = transform_cfg, method_cfg = method_cfg)
  if (length(polygons) == 0L) {
    return(list(name = method_cfg$name, status = "empty", reason = "No usable polygons were loaded."))
  }
  polygon_path <- file.path(polygons_dir, paste0(method_cfg$name, ".geojson.gz"))
  write_json_any(geojson_from_polygons(polygons, method_name = method_cfg$name, dataset_id = method_cfg$dataset_id %||% NULL), polygon_path)

  summary <- summarize_polygons_by_cell(polygons)
  list(
    name = method_cfg$name,
    status = "ready",
    mode = "polygon",
    loader = method_cfg$input$loader,
    polygons = polygons,
    polygon_path = polygon_path,
    n_cells_raw = uniqueN(summary$area_dt$cell_id),
    cell_areas = setNames(summary$area_dt$area, summary$area_dt$cell_id),
    cell_centroids = summary$centroid_dt
  )
}

read_affine_matrix_any <- function(path) {
  A <- as.matrix(data.table::fread(path, header = FALSE, sep = " ", fill = TRUE))
  if (!all(dim(A) == c(3, 3))) {
    A <- as.matrix(data.table::fread(path, header = FALSE))
  }
  if (!all(dim(A) == c(3, 3))) {
    stop("Affine transform must be a 3x3 matrix: ", path)
  }
  storage.mode(A) <- "double"
  A
}

load_mask_method <- function(method_cfg) {
  mask <- tiff::readTIFF(method_cfg$input$path, as.is = TRUE)
  if (length(dim(mask)) != 2L) {
    stop("Mask must be 2D: ", method_cfg$input$path)
  }
  labels <- sort(unique(as.integer(round(as.vector(mask)))))
  labels <- labels[labels > 0]
  list(
    name = method_cfg$name,
    status = "ready",
    mode = "mask",
    loader = method_cfg$input$loader,
    mask = mask,
    affine = read_affine_matrix_any(method_cfg$input$affine_matrix),
    n_cells_raw = length(labels),
    cell_areas = NULL,
    cell_centroids = NULL,
    polygon_path = NULL
  )
}

prepare_method_state <- function(method_cfg, polygons_dir) {
  if (!isTRUE(method_cfg$enabled)) {
    return(list(
      name = method_cfg$name,
      status = "disabled",
      reason = method_cfg$reason %||% "Disabled by run config.",
      mode = NULL,
      loader = NULL,
      n_cells_raw = 0L,
      polygon_path = NULL,
      cell_areas = NULL,
      cell_centroids = NULL
    ))
  }

  if (is.null(method_cfg$input$path) || !file.exists(method_cfg$input$path)) {
    return(list(
      name = method_cfg$name,
      status = "missing_input",
      reason = paste0("Missing input path: ", method_cfg$input$path),
      mode = NULL,
      loader = method_cfg$input$loader %||% NULL,
      n_cells_raw = 0L,
      polygon_path = NULL,
      cell_areas = NULL,
      cell_centroids = NULL
    ))
  }

  loader <- method_cfg$input$loader
  if (identical(loader, "geojson")) return(load_geojson_method(method_cfg, polygons_dir))
  if (identical(loader, "xenium_vertex_csv")) return(load_xenium_vertex_csv_method(method_cfg, polygons_dir))
  if (identical(loader, "bidcell_geojson")) return(load_bidcell_geojson_method(method_cfg, polygons_dir))
  if (identical(loader, "mask_tiff_affine")) return(load_mask_method(method_cfg))
  stop("Unsupported method loader: ", loader)
}
