suppressPackageStartupMessages({
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/config_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/spatial_eval_utils.R")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_spatial_eval.R <run_config.R>")
}

compiled <- compile_run_config(args[[1]])
manifest_path <- run_spatial_eval(compiled)
message("Spatial evaluation outputs written to: ", manifest_path)
