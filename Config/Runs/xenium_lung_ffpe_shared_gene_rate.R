sys.source(
  "/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/xenium_lung_ffpe_shared.R",
  envir = environment()
)

config$output_root <- "/projects/standard/zhangl4/shared/Data_feng0421/Evaluation_MainV2/xenium_lung_ffpe_gene_rate"

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
