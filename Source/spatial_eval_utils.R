suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(scales)
})

source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/generic_banksy.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

spatial_theme_academic <- function(base_size = 10) {
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

normalize_cell_types_spatial <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  idx <- tolower(trimws(x)) %in% c("unknown", "unassigned", "none", "na")
  x[idx] <- "Unknown"
  x
}

read_markers_spatial <- function(path) {
  markers <- fread(path)
  required <- c("cell_type", "marker")
  missing <- setdiff(required, names(markers))
  if (length(missing) > 0) {
    stop("Markers CSV is missing required columns: ", paste(missing, collapse = ", "))
  }
  markers[, cell_type := as.character(cell_type)]
  markers[, marker := as.character(marker)]
  unique(markers[!is.na(cell_type) & !is.na(marker) & cell_type != "" & marker != ""])
}

safe_quantile_mainv2 <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

safe_mean_mainv2 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x, na.rm = TRUE)
}

choose2_mainv2 <- function(x) {
  x <- as.numeric(x)
  x * (x - 1) / 2
}

compute_nmi_mainv2 <- function(labels_a, labels_b) {
  keep <- !is.na(labels_a) & !is.na(labels_b)
  labels_a <- as.character(labels_a[keep])
  labels_b <- as.character(labels_b[keep])
  if (!length(labels_a)) return(NA_real_)

  tab <- table(labels_a, labels_b)
  n <- sum(tab)
  if (n <= 0) return(NA_real_)

  p_ij <- tab / n
  p_i <- rowSums(p_ij)
  p_j <- colSums(p_ij)

  mi <- 0
  nz <- which(p_ij > 0, arr.ind = TRUE)
  if (nrow(nz) > 0) {
    mi <- sum(vapply(seq_len(nrow(nz)), function(idx) {
      i <- nz[idx, 1]
      j <- nz[idx, 2]
      p_ij[i, j] * log(p_ij[i, j] / (p_i[i] * p_j[j]))
    }, numeric(1)))
  }

  h_i <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  h_j <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))
  denom <- sqrt(h_i * h_j)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  mi / denom
}

compute_ari_mainv2 <- function(labels_a, labels_b) {
  keep <- !is.na(labels_a) & !is.na(labels_b)
  labels_a <- as.character(labels_a[keep])
  labels_b <- as.character(labels_b[keep])
  if (!length(labels_a)) return(NA_real_)

  tab <- table(labels_a, labels_b)
  n <- sum(tab)
  if (n < 2) return(NA_real_)

  sum_comb <- sum(choose2_mainv2(tab))
  row_comb <- sum(choose2_mainv2(rowSums(tab)))
  col_comb <- sum(choose2_mainv2(colSums(tab)))
  total_comb <- choose2_mainv2(n)
  if (!is.finite(total_comb) || total_comb <= 0) return(NA_real_)

  expected <- (row_comb * col_comb) / total_comb
  max_index <- 0.5 * (row_comb + col_comb)
  denom <- max_index - expected
  if (!is.finite(denom) || abs(denom) < 1e-12) return(NA_real_)
  (sum_comb - expected) / denom
}

compute_knn_metrics_mainv2 <- function(cell_dt, coord_cols = c("x", "y"), type_col = "cell_type", k = 10L) {
  cell_dt <- as.data.table(cell_dt)
  if (!all(coord_cols %in% names(cell_dt))) {
    stop("Missing coordinate columns: ", paste(setdiff(coord_cols, names(cell_dt)), collapse = ", "))
  }
  if (!type_col %in% names(cell_dt)) {
    stop("Missing type column: ", type_col)
  }

  cell_dt <- unique(cell_dt[, .(
    cell_id = as.character(cell_id),
    x = as.numeric(get(coord_cols[1])),
    y = as.numeric(get(coord_cols[2])),
    cell_type = normalize_cell_types_spatial(get(type_col))
  )], by = "cell_id")
  cell_dt <- cell_dt[is.finite(x) & is.finite(y)]
  if (nrow(cell_dt) < 3) {
    stop("At least three cells with coordinates are required for spatial evaluation.")
  }

  coords <- as.matrix(cell_dt[, .(x, y)])
  k_use <- min(as.integer(k), nrow(cell_dt) - 1L)
  knn <- FNN::get.knn(coords, k = k_use)

  types_vec <- cell_dt$cell_type
  neighbor_types <- matrix(types_vec[knn$nn.index], nrow = nrow(knn$nn.index), ncol = ncol(knn$nn.index))
  s_same <- rowMeans(neighbor_types == types_vec, na.rm = TRUE)

  d_same <- rep(NA_real_, nrow(cell_dt))
  for (tp in unique(types_vec)) {
    idx <- which(types_vec == tp)
    if (length(idx) < 2) next
    sub_coords <- coords[idx, , drop = FALSE]
    sub_knn <- FNN::get.knn(sub_coords, k = 1L)
    d_same[idx] <- sub_knn$nn.dist[, 1]
  }

  per_cell <- data.table(
    cell_id = cell_dt$cell_id,
    cell_type = types_vec,
    x = coords[, 1],
    y = coords[, 2],
    s_same = s_same,
    d_same = d_same
  )

  all_types <- sort(unique(types_vec))
  comp_list <- vector("list", nrow(cell_dt))
  for (i in seq_len(nrow(cell_dt))) {
    nt <- neighbor_types[i, ]
    tab <- as.data.table(table(neighbor_type = factor(nt, levels = all_types)))
    tab[, cell_type := types_vec[i]]
    tab[, frac := N / k_use]
    comp_list[[i]] <- tab[, .(cell_type, neighbor_type, frac)]
  }
  comp_dt <- rbindlist(comp_list, fill = TRUE)
  composition <- comp_dt[, .(pi = mean(frac)), by = .(cell_type, neighbor_type)]

  neighbor_ids <- matrix(cell_dt$cell_id[knn$nn.index], nrow = nrow(knn$nn.index), ncol = ncol(knn$nn.index))
  list(
    per_cell = per_cell,
    composition = composition,
    neighbor_ids = neighbor_ids,
    neighbor_dists = knn$nn.dist,
    k = k_use
  )
}

build_domain_backend_inputs_mainv2 <- function(method_name, eval_dir, cell_meta, counts_dt, domain_cfg) {
  input_dir <- file.path(eval_dir, "domain_backend_inputs", method_name)
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  cells_dt <- unique(as.data.table(cell_meta)[, .(
    cell_id = as.character(cell_id),
    x = as.numeric(x),
    y = as.numeric(y),
    cell_type = normalize_cell_types_spatial(cell_type)
  )], by = "cell_id")

  keep_ids <- cells_dt$cell_id
  counts_sub <- as.data.table(counts_dt)
  counts_sub[, cell_id := as.character(cell_id)]
  counts_sub <- counts_sub[cell_id %in% keep_ids]

  cells_path <- file.path(input_dir, "cells_eval.csv")
  counts_path <- file.path(input_dir, "counts_eval_long.csv.gz")
  params_path <- file.path(input_dir, "backend_params.json")
  manifest_path <- file.path(input_dir, "input_manifest.json")

  fwrite(cells_dt, cells_path)
  fwrite(counts_sub, counts_path)

  params <- list(
    backend = as.character(domain_cfg$backend %||% "banksy"),
    strict_labels = isTRUE(domain_cfg$strict_labels %||% TRUE),
    allow_unlabeled = isTRUE(domain_cfg$allow_unlabeled %||% FALSE),
    banksy = domain_cfg$banksy %||% list()
  )
  jsonlite::write_json(params, params_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

  manifest <- list(
    method = method_name,
    cells_path = cells_path,
    counts_path = counts_path,
    params_path = params_path
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")

  list(
    input_dir = input_dir,
    manifest_path = manifest_path,
    cells_path = cells_path,
    counts_path = counts_path,
    params_path = params_path
  )
}

validate_domain_labels_mainv2 <- function(labels_dt, cells_dt, allow_unlabeled = FALSE) {
  labels_dt <- as.data.table(labels_dt)
  cells_dt <- as.data.table(cells_dt)
  required <- c("cell_id", "domain_label")
  missing <- setdiff(required, names(labels_dt))
  if (length(missing) > 0) {
    stop("Domain label output is missing required columns: ", paste(missing, collapse = ", "))
  }

  labels_dt[, cell_id := as.character(cell_id)]
  cells_dt[, cell_id := as.character(cell_id)]
  labels_dt[, domain_label := as.character(domain_label)]
  labels_dt <- labels_dt[!is.na(cell_id) & nzchar(cell_id)]
  labels_dt <- labels_dt[!is.na(domain_label) & nzchar(trimws(domain_label))]
  if (anyDuplicated(labels_dt$cell_id) > 0) {
    dup <- unique(labels_dt$cell_id[duplicated(labels_dt$cell_id)])
    stop("Domain labels contain duplicate cell_id values: ", paste(head(dup, 5), collapse = ", "))
  }

  extra_ids <- setdiff(labels_dt$cell_id, cells_dt$cell_id)
  if (length(extra_ids) > 0) {
    stop("Domain labels contain cell_id values outside the current eval universe: ", paste(head(extra_ids, 5), collapse = ", "))
  }

  missing_ids <- setdiff(cells_dt$cell_id, labels_dt$cell_id)
  if (length(missing_ids) > 0 && !isTRUE(allow_unlabeled)) {
    stop("Domain labels are missing eval cells. Example missing cell_id values: ", paste(head(missing_ids, 5), collapse = ", "))
  }

  merged <- merge(cells_dt, labels_dt, by = "cell_id", all.x = TRUE, all.y = FALSE)
  merged[, domain_label := fifelse(is.na(domain_label), NA_character_, as.character(domain_label))]
  merged
}

evaluate_domain_labels_mainv2 <- function(labels_merged, backend_name, backend_cfg = list(), cluster_col = "domain_label") {
  score_dt <- labels_merged[cell_type != "Unknown" & !is.na(domain_label)]
  composition <- labels_merged[!is.na(domain_label), .N, by = .(domain_label, cell_type)][
    , frac := N / sum(N), by = domain_label
  ][order(domain_label, -frac, cell_type)]
  domain_summary <- composition[, .SD[1], by = domain_label][
    labels_merged[!is.na(domain_label), .(n_cells = .N), by = domain_label],
    on = "domain_label"
  ][, .(
    domain_label,
    n_cells,
    dominant_cell_type = cell_type,
    dominant_frac = frac
  )]

  data.table(
    domain_method = backend_name,
    cluster_col = cluster_col,
    n_cells = nrow(labels_merged),
    n_cells_labeled = sum(!is.na(labels_merged$domain_label)),
    label_coverage = mean(!is.na(labels_merged$domain_label)),
    n_cells_scored = nrow(score_dt),
    n_domains = uniqueN(labels_merged$domain_label[!is.na(labels_merged$domain_label)]),
    nmi = compute_nmi_mainv2(score_dt$cell_type, score_dt$domain_label),
    ari = compute_ari_mainv2(score_dt$cell_type, score_dt$domain_label),
    backend_resolution = paste(as.character(backend_cfg$resolutions %||% NA), collapse = ";"),
    backend_k_geom = as.integer(backend_cfg$k_geom %||% NA_integer_),
    backend_lambda = as.numeric(backend_cfg$lambda %||% NA_real_)
  ) -> metrics

  list(
    per_cell = labels_merged,
    composition = composition,
    domain_summary = domain_summary,
    metrics = metrics
  )
}

run_domain_backend_banksy_mainv2 <- function(input_info, output_dir, domain_cfg) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  labels_path <- file.path(output_dir, "domain_labels.csv")
  plot_dir <- file.path(output_dir, "plots")
  banksy_cfg <- domain_cfg$banksy %||% list()
  res <- run_banksy_from_standardized_inputs_mainv2(
    cells_path = input_info$cells_path,
    counts_path = input_info$counts_path,
    labels_path = labels_path,
    plot_dir = plot_dir,
    k_geom = as.integer(banksy_cfg$k_geom %||% 20L),
    lambda = as.numeric(banksy_cfg$lambda %||% 0.6),
    resolutions = as.numeric(banksy_cfg$resolutions %||% c(0.2)),
    compute_agf = isTRUE(banksy_cfg$compute_agf %||% TRUE),
    use_agf = isTRUE(banksy_cfg$use_agf %||% TRUE)
  )
  list(
    labels_path = labels_path,
    backend_name = "banksy",
    cluster_col = res$cluster_col,
    metadata = res$metadata
  )
}

run_domain_backend_script_mainv2 <- function(input_info, output_dir, script_cfg) {
  command <- script_cfg$command %||% character()
  if (!length(command)) {
    stop("Domain backend `script` requires `spatial$domain_clustering$script$command`.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  labels_path <- file.path(output_dir, "domain_labels.csv")
  log_path <- file.path(output_dir, "script_backend.log")
  args <- c(
    command[-1],
    as.character(script_cfg$extra_args %||% character()),
    input_info$manifest_path,
    labels_path
  )

  workdir <- script_cfg$workdir %||% NULL
  if (!is.null(workdir) && nzchar(as.character(workdir[[1]]))) {
    workdir <- as.character(workdir[[1]])
    if (!dir.exists(workdir)) {
      stop("Domain backend script workdir does not exist: ", workdir)
    }
    old_wd <- getwd()
    setwd(workdir)
    on.exit(setwd(old_wd), add = TRUE)
  }

  out <- system2(
    command = command[[1]],
    args = args,
    stdout = TRUE,
    stderr = TRUE,
    wait = TRUE
  )
  writeLines(out, log_path)
  status <- attr(out, "status") %||% 0L
  if (!identical(as.integer(status), 0L)) {
    stop("Domain backend script failed with exit status ", status, ". See ", log_path)
  }
  if (!file.exists(labels_path)) {
    stop("Domain backend script did not create expected labels file: ", labels_path)
  }
  list(
    labels_path = labels_path,
    backend_name = "script",
    cluster_col = "domain_label",
    metadata = data.table(
      script_command = paste(command, collapse = " "),
      script_log = log_path
    )
  )
}

run_domain_backend_mainv2 <- function(method_name, eval_dir, cell_meta, counts_dt, domain_cfg) {
  if (!isTRUE(domain_cfg$enabled %||% TRUE)) {
    return(list(status = "skipped", reason = "domain backend disabled"))
  }
  backend <- tolower(as.character(domain_cfg$backend %||% "banksy"))
  if (identical(backend, "none")) {
    return(list(status = "skipped", reason = "domain backend disabled"))
  }

  input_info <- build_domain_backend_inputs_mainv2(
    method_name = method_name,
    eval_dir = eval_dir,
    cell_meta = cell_meta,
    counts_dt = counts_dt,
    domain_cfg = domain_cfg
  )
  backend_dir <- file.path(eval_dir, "domain_backend_outputs", method_name)

  backend_run <- tryCatch({
    if (identical(backend, "banksy")) {
      run_domain_backend_banksy_mainv2(input_info, backend_dir, domain_cfg)
    } else if (identical(backend, "script")) {
      run_domain_backend_script_mainv2(input_info, backend_dir, domain_cfg$script %||% list())
    } else {
      stop("Unsupported domain backend: ", backend)
    }
  }, error = function(e) {
    list(status = "failed", reason = conditionMessage(e), input_info = input_info)
  })

  if (!is.null(backend_run$status) && backend_run$status %in% c("failed", "skipped")) {
    return(backend_run)
  }

  cells_dt <- fread(input_info$cells_path)
  labels_dt <- fread(backend_run$labels_path)
  merged <- validate_domain_labels_mainv2(
    labels_dt = labels_dt,
    cells_dt = cells_dt,
    allow_unlabeled = isTRUE(domain_cfg$allow_unlabeled %||% FALSE)
  )
  eval_res <- evaluate_domain_labels_mainv2(
    labels_merged = merged,
    backend_name = backend_run$backend_name,
    backend_cfg = if (identical(backend, "banksy")) domain_cfg$banksy %||% list() else list(),
    cluster_col = backend_run$cluster_col %||% "domain_label"
  )

  list(
    status = "ok",
    input_info = input_info,
    backend_output_dir = backend_dir,
    labels_path = backend_run$labels_path,
    backend_name = backend_run$backend_name,
    metadata = backend_run$metadata %||% NULL,
    per_cell = eval_res$per_cell,
    composition = eval_res$composition,
    domain_summary = eval_res$domain_summary,
    metrics = eval_res$metrics
  )
}

resolve_leakage_weights_mainv2 <- function(distances, weighting_scheme = "uniform", distance_scale = NULL) {
  weighting_scheme <- tolower(as.character(weighting_scheme %||% "uniform"))
  distances <- as.numeric(distances)
  if (!weighting_scheme %in% c("uniform", "distance_decay")) {
    weighting_scheme <- "uniform"
  }

  if (identical(weighting_scheme, "uniform")) {
    return(list(
      weights = rep(1, length(distances)),
      weighting_scheme = weighting_scheme,
      distance_scale = NA_real_
    ))
  }

  inferred_scale <- suppressWarnings(as.numeric(distance_scale[[1]] %||% NA_real_))
  if (!length(inferred_scale) || !is.finite(inferred_scale) || inferred_scale <= 0) {
    positive_dist <- distances[is.finite(distances) & distances > 0]
    inferred_scale <- if (length(positive_dist)) stats::median(positive_dist) else 1
  }
  if (!length(inferred_scale) || !is.finite(inferred_scale) || inferred_scale <= 0) {
    inferred_scale <- 1
  }

  weights <- exp(-distances / inferred_scale)
  weights[!is.finite(weights)] <- 0

  list(
    weights = weights,
    weighting_scheme = weighting_scheme,
    distance_scale = inferred_scale
  )
}

compute_marker_leakage_mainv2 <- function(
  counts_dt,
  knn_neighbor_ids,
  knn_neighbor_dists,
  cell_meta,
  markers_path,
  type_col = "cell_type",
  leakage_cfg = list()
) {
  markers <- read_markers_spatial(markers_path)
  counts_dt <- as.data.table(counts_dt)
  cell_meta <- as.data.table(cell_meta)

  counts_dt[, cell_id := as.character(cell_id)]
  counts_dt[, gene := as.character(gene)]
  counts_dt[, count := as.numeric(count)]
  cell_meta[, cell_id := as.character(cell_id)]
  cell_meta[, (type_col) := normalize_cell_types_spatial(get(type_col))]

  marker_scaling <- tolower(as.character(leakage_cfg$scaling %||% "raw"))
  if (!marker_scaling %in% c("raw", "gene_rate")) {
    stop("Unknown spatial marker_leakage scaling '", marker_scaling, "'. Expected one of: raw, gene_rate.")
  }

  marker_counts <- counts_dt[gene %in% markers$marker]
  if (!nrow(marker_counts)) {
    return(NULL)
  }
  if (marker_scaling == "gene_rate") {
    gene_totals <- counts_dt[, .(gene_total = sum(count, na.rm = TRUE)), by = gene]
    marker_counts <- merge(marker_counts, gene_totals, by = "gene", all.x = TRUE)
    marker_counts[is.na(gene_total) | gene_total <= 0, gene_total := 1]
    marker_counts[, count := count / gene_total]
    marker_counts[, gene_total := NULL]
  }
  marker_counts <- merge(marker_counts, markers, by.x = "gene", by.y = "marker", allow.cartesian = TRUE)
  marker_counts <- marker_counts[, .(marker_sum = sum(count, na.rm = TRUE)), by = .(cell_id, marker_type = cell_type)]

  cell_ids <- cell_meta$cell_id
  self_dt <- data.table(
    focal_id = cell_ids,
    neighbor_id = cell_ids,
    distance = 0
  )
  knn_dt <- data.table(
    focal_id = rep(cell_ids, each = ncol(knn_neighbor_ids)),
    neighbor_id = as.character(as.vector(t(knn_neighbor_ids))),
    distance = as.numeric(as.vector(t(knn_neighbor_dists)))
  )
  neighborhood_dt <- rbindlist(list(self_dt, knn_dt), use.names = TRUE, fill = TRUE)

  weight_res <- resolve_leakage_weights_mainv2(
    distances = neighborhood_dt$distance,
    weighting_scheme = leakage_cfg$weighting %||% "uniform",
    distance_scale = leakage_cfg$distance_scale %||% NULL
  )
  neighborhood_dt[, weight := weight_res$weights]

  neigh_counts <- merge(
    neighborhood_dt,
    marker_counts,
    by.x = "neighbor_id",
    by.y = "cell_id",
    all.x = FALSE,
    all.y = FALSE,
    allow.cartesian = TRUE
  )
  neigh_counts[, weighted_marker_sum := weight * marker_sum]

  totals <- neigh_counts[, .(
    marker_sum = sum(weighted_marker_sum, na.rm = TRUE)
  ), by = .(focal_id, marker_type)]
  totals_wide <- dcast(totals, focal_id ~ marker_type, value.var = "marker_sum", fill = 0)
  setnames(totals_wide, "focal_id", "cell_id")

  merged <- merge(
    cell_meta[, .(cell_id, focal_type = get(type_col))],
    totals_wide,
    by = "cell_id",
    all.x = TRUE
  )
  marker_cols <- setdiff(names(merged), c("cell_id", "focal_type"))
  if (!length(marker_cols)) {
    return(NULL)
  }
  merged[, (marker_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, as.numeric(x))), .SDcols = marker_cols]
  merged[, total_markers := rowSums(.SD, na.rm = TRUE), .SDcols = marker_cols]

  merged[, target_sum := 0]
  for (mc in marker_cols) {
    merged[focal_type == mc, target_sum := get(mc)]
  }
  merged[, purity_nbr := fifelse(total_markers > 0, target_sum / total_markers, NA_real_)]

  contam_long <- melt(
    merged,
    id.vars = c("cell_id", "focal_type", "total_markers"),
    measure.vars = marker_cols,
    variable.name = "neighbor_marker_type",
    value.name = "marker_sum"
  )
  contam_long <- contam_long[neighbor_marker_type != focal_type]
  contam_long[, contam_frac := fifelse(total_markers > 0, marker_sum / total_markers, NA_real_)]

  purity_summary <- merged[focal_type != "Unknown", .(
    purity_p25 = safe_quantile_mainv2(purity_nbr, 0.25),
    purity_median = safe_quantile_mainv2(purity_nbr, 0.5),
    purity_p75 = safe_quantile_mainv2(purity_nbr, 0.75),
    purity_mean = safe_mean_mainv2(purity_nbr)
  ), by = focal_type]

  contam_summary <- contam_long[!is.na(contam_frac), .(
    contam_p25 = safe_quantile_mainv2(contam_frac, 0.25),
    contam_median = safe_quantile_mainv2(contam_frac, 0.5),
    contam_p75 = safe_quantile_mainv2(contam_frac, 0.75),
    contam_mean = safe_mean_mainv2(contam_frac)
  ), by = .(focal_type, neighbor_marker_type)]

  metadata <- data.table(
    weighting_scheme = weight_res$weighting_scheme,
    distance_scale = weight_res$distance_scale,
    marker_scaling = marker_scaling
  )

  list(
    per_cell = merged[, .(cell_id, focal_type, purity_nbr, total_markers)],
    purity_summary = purity_summary,
    contam_summary = contam_summary,
    metadata = metadata,
    weighting_scheme = weight_res$weighting_scheme,
    distance_scale = weight_res$distance_scale
  )
}

plot_knn_metrics_mainv2 <- function(knn_res, method_name, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  same_dt <- copy(knn_res$per_cell)[cell_type != "Unknown" & !is.na(s_same)]
  same_summary <- same_dt[, .(
    s_same_p25 = safe_quantile_mainv2(s_same, 0.25),
    s_same_median = safe_quantile_mainv2(s_same, 0.5),
    s_same_p75 = safe_quantile_mainv2(s_same, 0.75)
  ), by = cell_type]

  p_same <- ggplot(same_dt, aes(x = s_same, y = cell_type)) +
    geom_count(shape = 21, fill = "#9ecae1", color = "white", alpha = 0.85, stroke = 0.25) +
    geom_segment(
      data = same_summary,
      aes(x = s_same_p25, xend = s_same_p75, y = cell_type, yend = cell_type),
      inherit.aes = FALSE,
      linewidth = 1.1,
      color = "#2b4c7e"
    ) +
    geom_point(
      data = same_summary,
      aes(x = s_same_median, y = cell_type),
      inherit.aes = FALSE,
      shape = 23,
      size = 2.7,
      fill = "#08306b",
      color = "white",
      stroke = 0.25
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.1),
      limits = c(0, 1)
    ) +
    scale_size_area(max_size = 7, name = "Cells") +
    labs(
      title = paste(method_name, "homotypic neighborhood fraction"),
      subtitle = sprintf("Discrete same-type neighbor fractions for k = %d; diamonds mark medians and bars show IQR", knn_res$k),
      x = "Fraction of nearest neighbors sharing the focal cell type",
      y = "Focal cell type"
    ) +
    spatial_theme_academic()

  p_dist <- ggplot(knn_res$per_cell[cell_type != "Unknown" & !is.na(d_same) & d_same > 0], aes(x = cell_type, y = d_same)) +
    geom_violin(fill = "#fd8d3c", alpha = 0.65, color = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.7) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(sigma = 1),
      labels = label_number()
    ) +
    coord_flip() +
    labs(
      title = paste(method_name, "nearest same-type distance"),
      subtitle = "Distance to the nearest cell with the same annotation; pseudo-log scale suppresses long-tail distortion",
      x = "Focal cell type",
      y = "Nearest same-type distance"
    ) +
    spatial_theme_academic()

  p_comp <- ggplot(knn_res$composition, aes(x = neighbor_type, y = cell_type, fill = pi)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", pi)), size = 2.5) +
    scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "#f1f1f1", limits = c(0, 1)) +
    labs(
      title = paste(method_name, "neighborhood composition"),
      subtitle = "Average neighbor-type composition by focal cell type",
      x = "Neighbor cell type",
      y = "Focal cell type",
      fill = expression(pi)
    ) +
    spatial_theme_academic() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

  ggsave(file.path(plot_dir, "knn_violin.png"), p_same / p_dist, width = 8.4, height = 9.6, dpi = 320)
  ggsave(file.path(plot_dir, "composition_heatmap.png"), p_comp, width = 7.2, height = 5.2, dpi = 320)
}

plot_domain_clustering_mainv2 <- function(domain_res, method_name, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  map_dt <- copy(domain_res$per_cell)[!is.na(domain_label)]
  comp_dt <- copy(domain_res$composition)
  metrics <- domain_res$metrics[1]
  metric_subtitle <- sprintf(
    "%s; NMI = %.3f, ARI = %.3f",
    metrics$domain_method[[1]],
    metrics$nmi[[1]],
    metrics$ari[[1]]
  )

  p_map <- ggplot(map_dt, aes(x = x, y = y, color = domain_label)) +
    geom_point(size = 0.7, alpha = 0.85) +
    coord_equal() +
    labs(
      title = paste(method_name, "spatial domains"),
      subtitle = metric_subtitle,
      x = "x",
      y = "y",
      color = "Domain"
    ) +
    spatial_theme_academic()

  p_comp <- ggplot(comp_dt, aes(x = cell_type, y = domain_label, fill = frac)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", frac)), size = 2.5) +
    scale_fill_gradient(low = "#fff7ec", high = "#7f2704", na.value = "#f1f1f1", limits = c(0, 1)) +
    labs(
      title = paste(method_name, "domain composition"),
      subtitle = "Within-domain cell-type fractions",
      x = "Annotated cell type",
      y = "Domain",
      fill = "Fraction"
    ) +
    spatial_theme_academic() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))

  ggsave(file.path(plot_dir, "domain_map.png"), p_map, width = 7.0, height = 5.8, dpi = 320)
  ggsave(file.path(plot_dir, "domain_composition_heatmap.png"), p_comp, width = 8.0, height = 5.2, dpi = 320)
}

plot_leakage_metrics_mainv2 <- function(leakage, method_name, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  purity_dt <- leakage$per_cell[focal_type != "Unknown" & !is.na(purity_nbr)]
  contam_dt <- leakage$contam_summary
  subtitle_suffix <- if (identical(leakage$weighting_scheme, "distance_decay")) {
    sprintf("Distance-decay weighting (scale = %.2f)", leakage$distance_scale)
  } else {
    "Uniform closed-neighborhood weights"
  }

  if (nrow(purity_dt)) {
    p_purity <- ggplot(purity_dt, aes(x = focal_type, y = purity_nbr)) +
      geom_violin(fill = "#74c476", alpha = 0.7, color = NA) +
      geom_boxplot(width = 0.12, outlier.size = 0.35, alpha = 0.6) +
      coord_flip() +
      labs(
        title = paste(method_name, "neighborhood marker purity"),
        subtitle = subtitle_suffix,
        x = "Focal cell type",
        y = "Purity_nbr"
      ) +
      spatial_theme_academic()
    ggsave(file.path(plot_dir, "marker_purity_violin.png"), p_purity, width = 7.0, height = 4.8, dpi = 320)
  }

  if (nrow(contam_dt)) {
    p_contam <- ggplot(contam_dt, aes(x = neighbor_marker_type, y = focal_type, fill = contam_median)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", contam_median)), size = 2.5) +
      scale_fill_gradient(low = "#fff5f0", high = "#99000d", na.value = "#f1f1f1", limits = c(0, 1)) +
      labs(
        title = paste(method_name, "marker contamination"),
        subtitle = subtitle_suffix,
        x = "Neighbor marker cell type",
        y = "Focal cell type",
        fill = "Median contam"
      ) +
      spatial_theme_academic() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggsave(file.path(plot_dir, "marker_contamination_heatmap.png"), p_contam, width = 7.8, height = 5.2, dpi = 320)
  }
}

plot_spatial_cross_method_summary <- function(knn_dt, purity_dt, domain_dt, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(knn_dt)) {
    plot_dt <- copy(knn_dt)[cell_type != "Unknown"]
    p_knn <- ggplot(plot_dt, aes(x = cell_type, y = method, fill = s_same_median)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", s_same_median)), size = 2.7) +
      scale_fill_gradient(low = "#f7fbff", high = "#08519c", na.value = "#f1f1f1", limits = c(0, 1)) +
      labs(
        title = "Spatial Homotypic Mixing Comparison",
        subtitle = "Median same-type neighbor fraction by method and annotated cell type",
        x = "Annotated cell type",
        y = "Segmentation method",
        fill = "Median s_same"
      ) +
      spatial_theme_academic() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggsave(file.path(outdir, "spatial_knn_comparison.png"), p_knn, width = 8.8, height = 4.8, dpi = 320)
  }

  if (nrow(purity_dt)) {
    plot_dt <- copy(purity_dt)[focal_type != "Unknown"]
    p_purity <- ggplot(plot_dt, aes(x = focal_type, y = method, fill = purity_median)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", purity_median)), size = 2.7) +
      scale_fill_gradient(low = "#f7fcf5", high = "#006d2c", na.value = "#f1f1f1", limits = c(0, 1)) +
      labs(
        title = "Spatial Marker Purity Comparison",
        subtitle = "Median neighborhood marker purity by method and annotated cell type",
        x = "Annotated cell type",
        y = "Segmentation method",
        fill = "Median purity"
      ) +
      spatial_theme_academic() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggsave(file.path(outdir, "spatial_marker_purity_comparison.png"), p_purity, width = 8.8, height = 4.8, dpi = 320)
  }

  if (nrow(domain_dt)) {
    domain_plot_dt <- melt(
      copy(domain_dt[, .(method, nmi, ari)]),
      id.vars = "method",
      variable.name = "metric",
      value.name = "value"
    )
    p_domain <- ggplot(domain_plot_dt, aes(x = metric, y = method, fill = value)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", value)), size = 2.8) +
      scale_fill_gradient(low = "#f7fbff", high = "#54278f", na.value = "#f1f1f1", limits = c(0, 1)) +
      labs(
        title = "Spatial Domain Consistency Comparison",
        subtitle = "Agreement between inferred domains and annotated cell types",
        x = NULL,
        y = "Segmentation method",
        fill = "Score"
      ) +
      spatial_theme_academic()
    ggsave(file.path(outdir, "spatial_domain_consistency_comparison.png"), p_domain, width = 6.6, height = 4.2, dpi = 320)
  }
}

run_spatial_eval <- function(config) {
  output_root <- config$output_root
  eval_dir <- file.path(output_root, "spatial_eval")
  morph_dir <- file.path(output_root, "morphology_eval")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)

  morph_manifest <- fread(file.path(morph_dir, "morphology_manifest.csv"))
  methods_to_run <- morph_manifest$method

  spatial_cfg <- config$spatial %||% list()
  domain_cfg <- spatial_cfg$domain_clustering %||% list(backend = "banksy")
  leakage_cfg <- spatial_cfg$marker_leakage %||% list()
  k_neighbors <- as.integer(spatial_cfg$k_neighbors %||% 10L)

  results <- list()
  knn_rows <- list()
  purity_rows <- list()
  domain_rows <- list()
  idx <- 1L

  for (method_name in methods_to_run) {
    meta_path <- file.path(morph_dir, paste0(method_name, "_morphology_per_cell.csv"))
    counts_path <- file.path(output_root, "expression", paste0(method_name, "_counts_long.csv.gz"))
    if (!file.exists(meta_path) || !file.exists(counts_path)) next

    message("[spatial] running method=", method_name)
    meta_dt <- fread(meta_path)
    counts_dt <- fread(counts_path)
    meta_dt <- meta_dt[is.finite(x) & is.finite(y)]
    if (!nrow(meta_dt)) next
    meta_dt[, cell_type := normalize_cell_types_spatial(cell_type_custom_final)]

    knn_res <- compute_knn_metrics_mainv2(
      cell_dt = meta_dt[, .(cell_id, cell_type, x, y)],
      coord_cols = c("x", "y"),
      type_col = "cell_type",
      k = k_neighbors
    )

    fwrite(knn_res$per_cell, file.path(eval_dir, paste0(method_name, "_spatial_knn_per_cell.csv")))
    fwrite(knn_res$composition, file.path(eval_dir, paste0(method_name, "_spatial_composition_matrix.csv")))

    knn_summary <- knn_res$per_cell[, .(
      s_same_p25 = safe_quantile_mainv2(s_same, 0.25),
      s_same_median = safe_quantile_mainv2(s_same, 0.5),
      s_same_p75 = safe_quantile_mainv2(s_same, 0.75),
      d_same_p25 = safe_quantile_mainv2(d_same, 0.25),
      d_same_median = safe_quantile_mainv2(d_same, 0.5),
      d_same_p75 = safe_quantile_mainv2(d_same, 0.75)
    ), by = cell_type]
    fwrite(knn_summary, file.path(eval_dir, paste0(method_name, "_spatial_knn_by_celltype.csv")))
    plot_knn_metrics_mainv2(knn_res, method_name, file.path(eval_dir, paste0(method_name, "_spatial_plots")))

    domain_status <- "skipped"
    domain_reason <- NA_character_
    domain_method <- tolower(as.character(domain_cfg$backend %||% "banksy"))
    domain_nmi <- NA_real_
    domain_ari <- NA_real_
    domain_n_clusters <- NA_integer_

    domain_res <- tryCatch(
      run_domain_backend_mainv2(
        method_name = method_name,
        eval_dir = eval_dir,
        cell_meta = meta_dt[, .(cell_id, cell_type, x, y)],
        counts_dt = counts_dt,
        domain_cfg = domain_cfg
      ),
      error = function(e) list(status = "failed", reason = conditionMessage(e))
    )

    if (identical(domain_res$status, "ok")) {
      fwrite(domain_res$per_cell, file.path(eval_dir, paste0(method_name, "_spatial_domain_per_cell.csv")))
      fwrite(domain_res$composition, file.path(eval_dir, paste0(method_name, "_spatial_domain_composition.csv")))
      fwrite(domain_res$domain_summary, file.path(eval_dir, paste0(method_name, "_spatial_domain_by_domain.csv")))
      fwrite(domain_res$metrics, file.path(eval_dir, paste0(method_name, "_spatial_domain_metrics.csv")))
      if (!is.null(domain_res$metadata) && inherits(domain_res$metadata, "data.table")) {
        fwrite(domain_res$metadata, file.path(eval_dir, paste0(method_name, "_spatial_domain_backend_metadata.csv")))
      }
      plot_domain_clustering_mainv2(domain_res, method_name, file.path(eval_dir, paste0(method_name, "_spatial_plots")))
      domain_rows[[idx]] <- copy(domain_res$metrics)[, method := method_name]
      domain_status <- "ok"
      domain_method <- domain_res$backend_name
      domain_nmi <- domain_res$metrics$nmi[[1]]
      domain_ari <- domain_res$metrics$ari[[1]]
      domain_n_clusters <- as.integer(domain_res$metrics$n_domains[[1]])
    } else {
      domain_status <- domain_res$status %||% "failed"
      domain_reason <- domain_res$reason %||% "unknown domain backend failure"
      backend_meta_path <- file.path(eval_dir, paste0(method_name, "_spatial_domain_backend_status.csv"))
      fwrite(data.table(status = domain_status, reason = domain_reason, backend = domain_method), backend_meta_path)
    }

    leakage <- compute_marker_leakage_mainv2(
      counts_dt = counts_dt,
      knn_neighbor_ids = knn_res$neighbor_ids,
      knn_neighbor_dists = knn_res$neighbor_dists,
      cell_meta = meta_dt[, .(cell_id, cell_type)],
      markers_path = config$markers_path,
      type_col = "cell_type",
      leakage_cfg = leakage_cfg
    )
    leakage_weighting <- leakage_cfg$weighting %||% "uniform"
    leakage_scaling <- tolower(as.character(leakage_cfg$scaling %||% "raw"))
    if (!is.null(leakage)) {
      fwrite(leakage$per_cell, file.path(eval_dir, paste0(method_name, "_spatial_marker_purity_per_cell.csv")))
      fwrite(leakage$purity_summary, file.path(eval_dir, paste0(method_name, "_spatial_marker_purity_by_celltype.csv")))
      fwrite(leakage$contam_summary, file.path(eval_dir, paste0(method_name, "_spatial_marker_contamination_by_pair.csv")))
      fwrite(leakage$metadata, file.path(eval_dir, paste0(method_name, "_spatial_marker_leakage_metadata.csv")))
      plot_leakage_metrics_mainv2(leakage, method_name, file.path(eval_dir, paste0(method_name, "_spatial_plots")))
      purity_rows[[idx]] <- copy(leakage$purity_summary)[, method := method_name]
      leakage_weighting <- leakage$weighting_scheme
      if (!is.null(leakage$metadata$marker_scaling)) {
        leakage_scaling <- leakage$metadata$marker_scaling[[1]]
      }
    }

    knn_rows[[idx]] <- copy(knn_summary)[, method := method_name]
    results[[idx]] <- data.table(
      method = method_name,
      n_cells = nrow(knn_res$per_cell),
      k_neighbors = knn_res$k,
      domain_backend = domain_method,
      domain_status = domain_status,
      domain_reason = domain_reason,
      domain_n_clusters = domain_n_clusters,
      domain_nmi = domain_nmi,
      domain_ari = domain_ari,
      leakage_weighting = leakage_weighting,
      leakage_scaling = leakage_scaling
    )
    idx <- idx + 1L
  }

  knn_dt <- if (length(knn_rows)) rbindlist(knn_rows, fill = TRUE) else data.table()
  purity_dt <- if (length(purity_rows)) rbindlist(purity_rows, fill = TRUE) else data.table()
  domain_dt <- if (length(domain_rows)) rbindlist(domain_rows, fill = TRUE) else data.table()

  if (nrow(knn_dt)) fwrite(knn_dt, file.path(eval_dir, "spatial_knn_summary_all_methods.csv"))
  if (nrow(purity_dt)) fwrite(purity_dt, file.path(eval_dir, "spatial_marker_purity_summary_all_methods.csv"))
  if (nrow(domain_dt)) fwrite(domain_dt, file.path(eval_dir, "spatial_domain_consistency_summary_all_methods.csv"))
  plot_spatial_cross_method_summary(knn_dt, purity_dt, domain_dt, eval_dir)

  manifest_path <- file.path(eval_dir, "spatial_manifest.csv")
  manifest_dt <- if (length(results)) rbindlist(results, fill = TRUE) else data.table(
    method = character(),
    n_cells = integer(),
    k_neighbors = integer(),
    domain_backend = character(),
    domain_status = character(),
    domain_reason = character(),
    domain_n_clusters = integer(),
    domain_nmi = numeric(),
    domain_ari = numeric(),
    leakage_weighting = character(),
    leakage_scaling = character()
  )
  fwrite(manifest_dt, manifest_path)
  manifest_path
}
