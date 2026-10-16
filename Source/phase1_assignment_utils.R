suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

point_in_ring <- function(px, py, vx, vy) {
  n <- length(vx)
  inside <- rep(FALSE, length(px))
  j <- n
  for (i in seq_len(n)) {
    cond <- ((vy[i] > py) != (vy[j] > py)) &
      (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i] + 1e-12) + vx[i])
    inside <- xor(inside, cond)
    j <- i
  }
  inside
}

assign_points_to_polygons <- function(transcripts_dt, polygons, x_col = "x", y_col = "y", grid_n = 200L) {
  px <- transcripts_dt[[x_col]]
  py <- transcripts_dt[[y_col]]
  N <- length(px)
  assign_vec <- rep("-1", N)

  x_min <- min(px, na.rm = TRUE); x_max <- max(px, na.rm = TRUE)
  y_min <- min(py, na.rm = TRUE); y_max <- max(py, na.rm = TRUE)
  x_span <- x_max - x_min + 1e-9
  y_span <- y_max - y_min + 1e-9

  # Build grid index: assign each transcript to a grid cell (integer key)
  gx <- pmin(as.integer(floor((px - x_min) / x_span * grid_n)), grid_n - 1L)
  gy <- pmin(as.integer(floor((py - y_min) / y_span * grid_n)), grid_n - 1L)
  grid_key <- gx * grid_n + gy
  grid_idx <- split(seq_len(N), grid_key)

  for (poly in polygons) {
    # Grid cells overlapping this polygon's bbox
    gx_lo <- max(0L, as.integer(floor((poly$bbox[1] - x_min) / x_span * grid_n)))
    gx_hi <- min(grid_n - 1L, as.integer(floor((poly$bbox[3] - x_min) / x_span * grid_n)))
    gy_lo <- max(0L, as.integer(floor((poly$bbox[2] - y_min) / y_span * grid_n)))
    gy_hi <- min(grid_n - 1L, as.integer(floor((poly$bbox[4] - y_min) / y_span * grid_n)))

    gx_seq <- seq.int(gx_lo, gx_hi)
    gy_seq <- seq.int(gy_lo, gy_hi)
    keys <- as.character(outer(gx_seq, gy_seq, function(a, b) a * grid_n + b))
    candidates <- unlist(grid_idx[keys], use.names = FALSE)
    if (length(candidates) == 0L) next

    candidates <- candidates[assign_vec[candidates] == "-1"]
    if (length(candidates) == 0L) next

    inside <- point_in_ring(px[candidates], py[candidates], poly$vx, poly$vy)
    if (any(inside) && length(poly$hole_vx %||% list())) {
      for (hole_idx in seq_along(poly$hole_vx)) {
        hole_inside <- point_in_ring(
          px[candidates],
          py[candidates],
          poly$hole_vx[[hole_idx]],
          poly$hole_vy[[hole_idx]]
        )
        inside <- inside & !hole_inside
        if (!any(inside)) break
      }
    }
    if (any(inside)) {
      assign_vec[candidates[inside]] <- poly$cell_id
    }
  }

  assign_vec
}

assign_points_to_mask <- function(transcripts_dt, mask, affine, x_col = "x", y_col = "y") {
  H <- nrow(mask)
  W <- ncol(mask)
  pts <- cbind(transcripts_dt[[x_col]], transcripts_dt[[y_col]], 1.0)
  px_py <- pts %*% t(affine)
  px <- floor(px_py[, 1])
  py <- floor(px_py[, 2])

  in_bounds <- !is.na(px) & !is.na(py) & px >= 0 & px < W & py >= 0 & py < H
  out <- rep("-1", nrow(transcripts_dt))
  if (any(in_bounds)) {
    label_vals <- mask[cbind(py[in_bounds] + 1L, px[in_bounds] + 1L)]
    label_vals <- as.character(as.integer(round(label_vals)))
    label_vals[is.na(label_vals) | label_vals == "0"] <- "-1"
    out[in_bounds] <- label_vals
  }
  out
}

lightweight_method_record <- function(method_state) {
  list(
    name = method_state$name,
    status = method_state$status %||% NULL,
    reason = method_state$reason %||% NULL,
    mode = method_state$mode %||% NULL,
    loader = method_state$loader %||% NULL,
    polygon_path = method_state$polygon_path %||% NULL,
    n_cells_raw = method_state$n_cells_raw %||% 0L,
    cell_areas = method_state$cell_areas %||% NULL,
    cell_centroids = method_state$cell_centroids %||% NULL
  )
}

append_method_assignments_dt <- function(transcripts_dt, method_state) {
  if (!identical(method_state$status, "ready")) {
    assignments <- rep("-1", nrow(transcripts_dt))
  } else if (identical(method_state$mode, "polygon")) {
    assignments <- assign_points_to_polygons(transcripts_dt, method_state$polygons, x_col = "x", y_col = "y")
  } else if (identical(method_state$mode, "mask")) {
    assignments <- assign_points_to_mask(transcripts_dt, method_state$mask, method_state$affine, x_col = "x", y_col = "y")
  } else {
    stop("Unsupported method mode: ", method_state$mode)
  }

  summary <- data.table(
    method = method_state$name,
    status = method_state$status %||% "unknown",
    n_cells_raw = as.integer(method_state$n_cells_raw %||% 0L),
    assigned_transcripts = sum(assignments != "-1"),
    background_transcripts = sum(assignments == "-1"),
    assignment_rate = mean(assignments != "-1"),
    polygon_path = method_state$polygon_path %||% NA_character_,
    mode = method_state$mode %||% NA_character_,
    reason = method_state$reason %||% NA_character_
  )
  summary[, low_assignment_warning := assignment_rate < 0.05]

  list(assignments = assignments, summary = summary)
}

build_assignments_for_phase1 <- function(config, transcript_meta, polygons_dir, assignments_dir) {
  dir.create(polygons_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(assignments_dir, recursive = TRUE, showWarnings = FALSE)

  transcripts_dt <- copy(transcript_meta$transcripts_eval_dt)
  n_methods <- length(config$methods)
  requested_cores <- suppressWarnings(as.integer(config$runtime$phase1_cores %||% NA_integer_))
  if (!is.finite(requested_cores) || requested_cores < 1L) {
    requested_cores <- parallel::detectCores()
  }
  n_cores <- min(n_methods, max(1L, requested_cores))

  # Each method reads the same transcripts and writes to its own polygon file — fully independent.
  # Workers get 1 DT thread each to avoid oversubscription across forks.
  results <- parallel::mclapply(config$methods, function(method_cfg) {
    data.table::setDTthreads(1L)
    message("[phase1] preparing method=", method_cfg$name)
    method_state <- prepare_method_state(method_cfg, polygons_dir)
    res <- append_method_assignments_dt(transcripts_dt, method_state)
    message(
      "[phase1] completed method=", method_cfg$name,
      " status=", res$summary$status[[1]],
      " assigned=", res$summary$assigned_transcripts[[1]]
    )
    list(
      name = method_cfg$name,
      assignments = res$assignments,
      summary = res$summary,
      record = lightweight_method_record(method_state)
    )
  }, mc.cores = n_cores)

  method_summaries <- list()
  method_records <- list()
  for (r in results) {
    col_name <- paste0("cell_id_", r$name)
    transcripts_dt[, (col_name) := r$assignments]
    method_summaries[[length(method_summaries) + 1L]] <- r$summary
    method_records[[length(method_records) + 1L]] <- r$record
  }

  assignments_path <- file.path(assignments_dir, "transcripts_with_assignments.csv.gz")
  assignment_summary_path <- file.path(assignments_dir, "assignment_summary_by_method.csv")
  data.table::fwrite(transcripts_dt, assignments_path)
  summary_dt <- rbindlist(method_summaries, fill = TRUE)
  data.table::fwrite(summary_dt, assignment_summary_path)

  list(
    assignments_path = assignments_path,
    summary_path = assignment_summary_path,
    assignments_dt = transcripts_dt,
    method_records = method_records
  )
}
