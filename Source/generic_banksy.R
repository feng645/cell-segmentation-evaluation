suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(jsonlite)
})

ensure_banksy_runtime_mainv2 <- function() {
  required <- c("Banksy", "SpatialExperiment", "SummarizedExperiment", "scater", "cowplot", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(
      "BANKSY backend requires the following packages: ",
      paste(missing, collapse = ", "),
      ". Install them or switch to the `script` backend."
    )
  }
  invisible(TRUE)
}

build_sparse_counts_from_long_mainv2 <- function(counts_dt, keep_cells) {
  counts_dt <- as.data.table(counts_dt)
  counts_dt[, cell_id := as.character(cell_id)]
  counts_dt[, gene := as.character(gene)]
  counts_dt[, count := as.numeric(count)]
  counts_dt <- counts_dt[cell_id %in% keep_cells]
  if (!nrow(counts_dt)) {
    stop("No assignment-derived counts remain after restricting to eval cells.")
  }

  genes <- sort(unique(counts_dt$gene))
  cells <- as.character(keep_cells)
  counts_dt[, gene_idx := match(gene, genes)]
  counts_dt[, cell_idx := match(cell_id, cells)]
  counts_dt <- counts_dt[!is.na(gene_idx) & !is.na(cell_idx)]

  Matrix::sparseMatrix(
    i = counts_dt$gene_idx,
    j = counts_dt$cell_idx,
    x = counts_dt$count,
    dims = c(length(genes), length(cells)),
    dimnames = list(genes, cells)
  )
}

build_spatial_experiment_from_standardized_inputs_mainv2 <- function(cells_path, counts_path) {
  ensure_banksy_runtime_mainv2()

  cells_dt <- fread(cells_path)
  required_cell_cols <- c("cell_id", "x", "y")
  missing_cell_cols <- setdiff(required_cell_cols, names(cells_dt))
  if (length(missing_cell_cols) > 0) {
    stop("cells_eval.csv is missing required columns: ", paste(missing_cell_cols, collapse = ", "))
  }

  cells_dt[, cell_id := as.character(cell_id)]
  cells_dt <- unique(cells_dt, by = "cell_id")
  cells_dt <- cells_dt[is.finite(x) & is.finite(y)]
  if (nrow(cells_dt) < 3) {
    stop("BANKSY backend requires at least three cells with finite spatial coordinates.")
  }

  counts_dt <- fread(counts_path)
  expr_mat <- build_sparse_counts_from_long_mainv2(counts_dt, keep_cells = cells_dt$cell_id)
  coords <- as.matrix(cells_dt[, .(x, y)])
  rownames(coords) <- cells_dt$cell_id
  colnames(coords) <- c("x", "y")

  se <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = expr_mat),
    colData = as.data.frame(cells_dt),
    spatialCoords = coords
  )
  SummarizedExperiment::colData(se)$sdimx <- SpatialExperiment::spatialCoords(se)[, 1]
  SummarizedExperiment::colData(se)$sdimy <- SpatialExperiment::spatialCoords(se)[, 2]
  se
}

open_plot_device_mainv2 <- function(path, width = 10, height = 5, units = "in", res = 300) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height, units = units, res = res)
  } else {
    png(filename = path, width = width, height = height, units = units, res = res, type = "cairo")
  }
  if (identical(names(dev.cur()), "null device")) stop("Failed to open graphics device for: ", path)
}

run_banksy_on_standardized_se_mainv2 <- function(
  se,
  assay_name = "counts",
  k_geom = 20L,
  lambda = 0.6,
  resolutions = c(0.2),
  compute_agf = TRUE,
  use_agf = TRUE
) {
  se <- Banksy::computeBanksy(se, assay_name = assay_name, compute_agf = compute_agf, k_geom = as.integer(k_geom))
  se <- Banksy::runBanksyPCA(se, use_agf = use_agf, lambda = lambda)

  for (res in as.numeric(resolutions)) {
    se <- Banksy::clusterBanksy(
      se,
      use_agf = use_agf,
      lambda = lambda,
      k_neighbors = as.integer(k_geom),
      resolution = res
    )
  }

  cluster_cols <- grep("^clust_", colnames(SummarizedExperiment::colData(se)), value = TRUE)
  if (!length(cluster_cols)) {
    stop("BANKSY did not produce any cluster columns.")
  }
  cluster_col <- cluster_cols[[1]]
  SummarizedExperiment::colData(se)$banksy_cluster <- SummarizedExperiment::colData(se)[[cluster_col]]
  list(se = se, cluster_col = cluster_col, cluster_cols = cluster_cols)
}

write_banksy_labels_mainv2 <- function(se, labels_path, cluster_col, cluster_cols) {
  labels_dt <- data.table(
    cell_id = colnames(se),
    domain_label = as.character(SummarizedExperiment::colData(se)[[cluster_col]])
  )
  labels_dt[, backend_name := "banksy"]
  labels_dt[, cluster_col := cluster_col]
  for (cnm in cluster_cols) {
    labels_dt[[cnm]] <- SummarizedExperiment::colData(se)[[cnm]]
  }
  fwrite(labels_dt, labels_path)
  labels_dt
}

plot_banksy_clusters_mainv2 <- function(se, cluster_cols, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  plot_path <- file.path(plot_dir, "banksy_clusters.png")
  plots <- lapply(cluster_cols, function(cnm) {
    scater::plotColData(se, x = "sdimx", y = "sdimy", colour_by = cnm, point_size = 0.25) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = cnm) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.title = ggplot2::element_blank())
  })
  open_plot_device_mainv2(plot_path, width = 10, height = 5, units = "in", res = 300)
  on.exit(dev.off(), add = TRUE)
  print(cowplot::plot_grid(plotlist = plots, ncol = 2))
  invisible(plot_path)
}

run_banksy_from_standardized_inputs_mainv2 <- function(
  cells_path,
  counts_path,
  labels_path,
  plot_dir = NULL,
  k_geom = 20L,
  lambda = 0.6,
  resolutions = c(0.2),
  compute_agf = TRUE,
  use_agf = TRUE
) {
  se <- build_spatial_experiment_from_standardized_inputs_mainv2(cells_path, counts_path)
  res <- run_banksy_on_standardized_se_mainv2(
    se = se,
    assay_name = "counts",
    k_geom = k_geom,
    lambda = lambda,
    resolutions = resolutions,
    compute_agf = compute_agf,
    use_agf = use_agf
  )

  dir.create(dirname(labels_path), recursive = TRUE, showWarnings = FALSE)
  labels_dt <- write_banksy_labels_mainv2(
    se = res$se,
    labels_path = labels_path,
    cluster_col = res$cluster_col,
    cluster_cols = res$cluster_cols
  )

  if (!is.null(plot_dir)) {
    plot_banksy_clusters_mainv2(res$se, res$cluster_cols, plot_dir)
  }

  list(
    labels_path = labels_path,
    cluster_col = res$cluster_col,
    cluster_cols = res$cluster_cols,
    metadata = data.table(
      backend_name = "banksy",
      cluster_col = res$cluster_col,
      n_cells = nrow(labels_dt),
      k_geom = as.integer(k_geom),
      lambda = as.numeric(lambda),
      resolutions = paste(as.character(resolutions), collapse = ";"),
      compute_agf = isTRUE(compute_agf),
      use_agf = isTRUE(use_agf)
    )
  )
}
