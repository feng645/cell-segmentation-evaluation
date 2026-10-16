suppressPackageStartupMessages({
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/config_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_transcript_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_boundary_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_assignment_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_cell_output_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_runner_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/annotation_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/expression_eval_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/morphology_eval_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/spatial_eval_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/visualization_utils.R")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_all_eval.R <run_config.R>")
}

compiled <- compile_run_config(args[[1]])
message("[run_all] phase1")
run_phase1_eval_r(compiled)
message("[run_all] annotation")
run_annotation_eval(compiled)
message("[run_all] expression")
run_expression_eval(compiled)
message("[run_all] morphology")
run_morphology_eval(compiled)
message("[run_all] spatial")
run_spatial_eval(compiled)
message("[run_all] visualization")
run_visualization_eval(
  compiled$output_root,
  visual_cfg = compiled$visualization %||% list(),
  markers_path = compiled$markers_path
)
message("[run_all] completed: ", compiled$output_root)
