suppressPackageStartupMessages({
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/config_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/morphology_eval_utils.R")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_morphology_eval.R <run_config.R>")
}

compiled <- compile_run_config(args[[1]])
manifest_path <- run_morphology_eval(compiled)
message("Morphology evaluation outputs written to: ", manifest_path)
