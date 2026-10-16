sys.source(
  "/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/merscope_lung_roi.R",
  envir = environment()
)

config$output_root <- "/users/9/feng0421/Code/R/Evaluation/MainV2/Outputs/merscope_lung_roi_banksy_backend"
config$spatial <- list(
  domain_clustering = list(
    enabled = TRUE,
    backend = "banksy"
  )
)
