sys.source(
  "/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/merscope_lung_roi.R",
  envir = environment()
)

config$output_root <- "/users/9/feng0421/Code/R/Evaluation/MainV2/Outputs/merscope_lung_roi_gene_rate"

config$expression <- list(marker_scaling = "gene_rate")
config$spatial <- list(
  domain_clustering = list(
    enabled = TRUE,
    backend = "banksy"
  ),
  marker_leakage = list(
    weighting = "uniform",
    scaling = "gene_rate"
  )
)
