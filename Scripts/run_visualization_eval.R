suppressPackageStartupMessages({
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/config_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/visualization_utils.R")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_visualization_eval.R <run_config.R | output_root>")
}

input_arg <- args[[1]]
markers_path <- NULL
if (file.exists(input_arg) && grepl("\\.R$", input_arg)) {
  compiled <- compile_run_config(input_arg)
  output_root <- compiled$output_root
  visual_cfg <- compiled$visualization %||% list()
  markers_path <- compiled$markers_path
} else {
  output_root <- normalizePath(input_arg, winslash = "/", mustWork = TRUE)
  visual_cfg <- list()
}

viz_dir <- run_visualization_eval(output_root, visual_cfg = visual_cfg, markers_path = markers_path)
message("Visualization outputs written to: ", viz_dir)
