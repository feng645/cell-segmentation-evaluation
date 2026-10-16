suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(cellmarkeraccordion)
})

annotation_theme_academic <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "#4d4d4d", size = base_size - 1),
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.35)
    )
}

read_counts_long <- function(path) {
  dt <- data.table::fread(path)
  required <- c("cell_id", "gene", "count")
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Counts file is missing required columns: ", paste(missing, collapse = ", "))
  }
  dt[, cell_id := as.character(cell_id)]
  dt[, gene := as.character(gene)]
  dt[, count := as.numeric(count)]
  dt
}

read_marker_table <- function(path) {
  markers <- data.table::fread(path)
  required <- c("cell_type", "marker")
  missing <- setdiff(required, names(markers))
  if (length(missing) > 0) {
    stop("Marker file is missing required columns: ", paste(missing, collapse = ", "))
  }
  markers[, cell_type := as.character(cell_type)]
  markers[, marker := as.character(marker)]
  markers <- unique(markers[!is.na(cell_type) & !is.na(marker) & cell_type != "" & marker != ""])
  if (!nrow(markers)) {
    stop("Marker file does not contain any valid marker rows after cleanup.")
  }
  as.data.frame(markers)
}

normalize_annotation_labels <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  unknown_idx <- tolower(trimws(x)) %in% c("unknown", "unassigned", "na", "none")
  x[unknown_idx] <- "Unknown"
  x
}

build_sparse_counts_matrix <- function(counts_dt, keep_cells = NULL) {
  dt <- copy(counts_dt)
  if (!is.null(keep_cells)) {
    keep_cells <- as.character(keep_cells)
    dt <- dt[cell_id %in% keep_cells]
  }
  if (!nrow(dt)) {
    stop("No counts remain after filtering cells for annotation.")
  }

  genes <- sort(unique(dt$gene))
  cells <- sort(unique(dt$cell_id))
  dt[, gene_idx := match(gene, genes)]
  dt[, cell_idx := match(cell_id, cells)]

  mat <- sparseMatrix(
    i = dt$gene_idx,
    j = dt$cell_idx,
    x = dt$count,
    dims = c(length(genes), length(cells)),
    dimnames = list(genes, cells)
  )
  mat
}

qc_filter_obj <- function(obj, assay = DefaultAssay(obj), min_feat = 30, max_feat = 8000, max_mt = 15, mt_pattern = "^MT") {
  DefaultAssay(obj) <- assay
  feat_col <- paste0("nFeature_", assay)
  mt_col <- "percent.mt"
  obj[[mt_col]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  md <- obj@meta.data
  keep <- md[[feat_col]] >= min_feat & md[[feat_col]] <= max_feat & md[[mt_col]] <= max_mt
  keep[is.na(keep)] <- FALSE
  subset(obj, cells = rownames(md)[keep])
}

normalize_obj_for_annotation <- function(obj, assay = DefaultAssay(obj), hvgs = 600L) {
  DefaultAssay(obj) <- assay
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = hvgs, verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj
}

cluster_obj_for_annotation <- function(obj, assay = DefaultAssay(obj), pc_use = 10L, resolution = 0.2) {
  DefaultAssay(obj) <- assay
  obj <- RunPCA(obj, npcs = pc_use, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = seq_len(pc_use), verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(
    obj,
    dims = seq_len(pc_use),
    umap.method = "uwot",
    metric = "cosine",
    n.neighbors = 35,
    min.dist = 0.15,
    spread = 1.0,
    verbose = FALSE
  )
  obj
}

annotate_seurat_obj <- function(obj, markers_path, assay = DefaultAssay(obj), min_n_marker = 2L) {
  DefaultAssay(obj) <- assay
  markers <- read_marker_table(markers_path)
  seurat_version <<- as.character(packageVersion("Seurat"))
  annotated <- accordion_custom(
    obj,
    assay = assay,
    annotation_resolution = "cell",
    marker_table = markers,
    category_column = "cell_type",
    marker_column = "marker",
    min_n_marker = min_n_marker,
    plot = TRUE,
    annotation_name = "cell_type_custom"
  )
  per_cell <- normalize_annotation_labels(annotated$cell_type_custom_per_cell)
  annotated$cell_type_custom_final <- per_cell
  annotated$cell_type_custom_per_cell <- per_cell
  annotated
}

build_annotation_table <- function(obj, method_name) {
  umap <- Embeddings(obj, "umap")
  meta <- as.data.table(obj@meta.data, keep.rownames = "cell_id")
  cluster_col <- grep("_snn_res\\.", names(meta), value = TRUE)
  cluster_col <- if (length(cluster_col)) cluster_col[[1]] else NA_character_

  out <- data.table(
    cell_id = meta$cell_id,
    method = method_name,
    nCount = meta[[paste0("nCount_", DefaultAssay(obj))]],
    nFeature = meta[[paste0("nFeature_", DefaultAssay(obj))]],
    percent_mt = meta[["percent.mt"]],
    cluster = if (!is.na(cluster_col)) as.character(meta[[cluster_col]]) else NA_character_,
    cell_type_custom_per_cell = as.character(meta[["cell_type_custom_per_cell"]]),
    cell_type_custom_final = as.character(meta[["cell_type_custom_final"]]),
    umap_1 = umap[, 1],
    umap_2 = umap[, 2]
  )
  out[, cell_type_custom_per_cell := normalize_annotation_labels(cell_type_custom_per_cell)]
  out[, cell_type_custom_final := normalize_annotation_labels(cell_type_custom_final)]
  out
}

plot_annotation_outputs <- function(annotation_dt, outdir, method_name) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  plot_dt <- copy(annotation_dt)
  plot_dt[, cell_type_custom_final := factor(cell_type_custom_final)]

  p_umap <- ggplot(plot_dt, aes(x = umap_1, y = umap_2, color = cell_type_custom_final)) +
    geom_point(size = 0.45, alpha = 0.85) +
    labs(
      title = paste(method_name, "annotation UMAP"),
      subtitle = "Marker-driven cell-type calls on unified per-method counts",
      x = "UMAP1",
      y = "UMAP2",
      color = "Cell type"
    ) +
    annotation_theme_academic()

  yield_dt <- plot_dt[, .N, by = cell_type_custom_final][order(-N)]
  p_yield <- ggplot(yield_dt, aes(x = reorder(cell_type_custom_final, N), y = N)) +
    geom_col(fill = "#315c78", width = 0.72) +
    coord_flip() +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = paste(method_name, "annotated cell yield"),
      subtitle = "Per-cell annotation counts after unified QC and clustering",
      x = NULL,
      y = "Cells"
    ) +
    annotation_theme_academic()

  ggsave(file.path(outdir, paste0(method_name, "_annotation_umap.png")), p_umap, width = 7.4, height = 5.8, dpi = 320)
  ggsave(file.path(outdir, paste0(method_name, "_annotation_yield.png")), p_yield, width = 6.8, height = 4.8, dpi = 320)
}

run_annotation_for_method <- function(method_name, counts_path, cell_table_path, markers_path, output_dir, annotation_cfg = list()) {
  counts_dt <- read_counts_long(counts_path)
  cell_table_dt <- data.table::fread(cell_table_path)
  cell_table_dt[, cell_id := as.character(cell_id)]
  keep_cells <- as.character(cell_table_dt[keep_for_eval == TRUE, cell_id])

  assay_name <- "originalexp"
  mat <- build_sparse_counts_matrix(counts_dt, keep_cells = keep_cells)
  obj <- CreateSeuratObject(counts = mat, assay = assay_name, min.cells = 0, min.features = 0)
  obj <- qc_filter_obj(
    obj,
    assay = assay_name,
    min_feat = as.integer(annotation_cfg$min_genes %||% 30L),
    max_feat = as.integer(annotation_cfg$max_genes %||% 8000L),
    max_mt = as.numeric(annotation_cfg$max_pct_mt %||% 15),
    mt_pattern = annotation_cfg$mt_pattern %||% "^MT"
  )
  obj <- normalize_obj_for_annotation(
    obj,
    assay = assay_name,
    hvgs = as.integer(annotation_cfg$hvgs %||% 600L)
  )
  obj <- cluster_obj_for_annotation(
    obj,
    assay = assay_name,
    pc_use = as.integer(annotation_cfg$pc_use %||% 10L),
    resolution = as.numeric(annotation_cfg$resolution %||% 0.2)
  )
  obj <- annotate_seurat_obj(
    obj,
    markers_path = markers_path,
    assay = assay_name,
    min_n_marker = as.integer(annotation_cfg$min_n_marker %||% 2L)
  )

  annotation_dt <- build_annotation_table(obj, method_name = method_name)
  extra_cols <- intersect(c("cell_id", "x", "y", "z", "cell_area", "keep_for_eval", "drop_reason"), names(cell_table_dt))
  if (length(extra_cols) > 1L) {
    annotation_dt <- merge(
      annotation_dt,
      unique(cell_table_dt[, ..extra_cols], by = "cell_id"),
      by = "cell_id",
      all.x = TRUE
    )
  }
  out_csv <- file.path(output_dir, paste0(method_name, "_annotation_per_cell.csv"))
  data.table::fwrite(annotation_dt, out_csv)
  plot_annotation_outputs(annotation_dt, output_dir, method_name)

  list(
    method = method_name,
    annotation_path = out_csv,
    cell_table_path = cell_table_path,
    n_cells_annotated = nrow(annotation_dt),
    n_cell_types = uniqueN(annotation_dt$cell_type_custom_final)
  )
}

run_annotation_eval <- function(config) {
  output_root <- config$output_root
  annotation_dir <- file.path(output_root, "annotation")
  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

  cell_manifest <- data.table::fread(file.path(output_root, "cell_tables", "method_cell_output_manifest.csv"))
  results <- list()
  idx <- 1L

  annotation_cfg <- config$annotation %||% list()
  if (is.null(annotation_cfg$mt_pattern)) {
    annotation_cfg$mt_pattern <- config$platform_protocol$mt_gene_regex %||% "^MT"
  }
  if (is.null(annotation_cfg$min_genes)) {
    annotation_cfg$min_genes <- config$platform_protocol$cell_filters$min_genes %||% 30L
  }
  if (is.null(annotation_cfg$max_genes)) {
    annotation_cfg$max_genes <- config$platform_protocol$cell_filters$max_genes %||% 8000L
  }
  if (is.null(annotation_cfg$max_pct_mt)) {
    annotation_cfg$max_pct_mt <- config$platform_protocol$cell_filters$max_pct_mt %||% 15
  }

  for (i in seq_len(nrow(cell_manifest))) {
    row <- cell_manifest[i]
    if (!identical(row$status, "ready")) next
    method_name <- row$method[[1]]
    counts_path <- file.path(output_root, "expression", paste0(method_name, "_counts_long.csv.gz"))
    cell_table_path <- row$eval_cell_table_path[[1]] %||% row$cell_table_path[[1]]
    if (!file.exists(counts_path) || !file.exists(cell_table_path)) next
    message("[annotation] running method=", method_name)
    res <- run_annotation_for_method(
      method_name = method_name,
      counts_path = counts_path,
      cell_table_path = cell_table_path,
      markers_path = config$markers_path,
      output_dir = annotation_dir,
      annotation_cfg = annotation_cfg
    )
    results[[idx]] <- as.data.table(res)
    idx <- idx + 1L
  }

  manifest_dt <- if (length(results)) rbindlist(results, fill = TRUE) else data.table(
    method = character(),
    annotation_path = character(),
    cell_table_path = character(),
    n_cells_annotated = integer(),
    n_cell_types = integer()
  )
  manifest_path <- file.path(annotation_dir, "annotation_manifest.csv")
  fwrite(manifest_dt, manifest_path)
  manifest_path
}
