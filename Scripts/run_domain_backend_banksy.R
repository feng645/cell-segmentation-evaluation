suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
  library(cowplot)
  library(scater)
})

source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/generic_banksy.R")

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_domain_backend_banksy.R <input_manifest.json> <output_labels.csv>")
}

manifest_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_labels_path <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)

manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
params_path <- manifest$params_path
params <- if (!is.null(params_path) && file.exists(params_path)) {
  jsonlite::read_json(params_path, simplifyVector = TRUE)
} else {
  list()
}
banksy_cfg <- params$banksy %||% list()

res <- run_banksy_from_standardized_inputs_mainv2(
  cells_path = manifest$cells_path,
  counts_path = manifest$counts_path,
  labels_path = output_labels_path,
  plot_dir = file.path(dirname(output_labels_path), "plots"),
  k_geom = as.integer(banksy_cfg$k_geom %||% 20L),
  lambda = as.numeric(banksy_cfg$lambda %||% 0.6),
  resolutions = as.numeric(banksy_cfg$resolutions %||% c(0.2)),
  compute_agf = isTRUE(banksy_cfg$compute_agf %||% TRUE),
  use_agf = isTRUE(banksy_cfg$use_agf %||% TRUE)
)

meta_path <- file.path(dirname(output_labels_path), "banksy_backend_metadata.csv")
data.table::fwrite(res$metadata, meta_path)
cat(output_labels_path, "\n")
