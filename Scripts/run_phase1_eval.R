suppressPackageStartupMessages({
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/config_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_transcript_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_boundary_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_assignment_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_cell_output_utils.R")
  source("/users/9/feng0421/Code/R/Evaluation/MainV2/Source/phase1_runner_utils.R")
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_phase1_eval.R <run_config.R>")
}

run_config_path <- args[[1]]
compiled <- compile_run_config(run_config_path)

dir.create(compiled$output_root, recursive = TRUE, showWarnings = FALSE)
compiled_json_path <- file.path(compiled$output_root, "compiled_run_config.json")
write_compiled_config(compiled, compiled_json_path)

message("Compiled config: ", compiled_json_path)
manifest_path <- run_phase1_eval_r(compiled)

message("Phase1 evaluation completed for dataset: ", compiled$dataset_id)
message("Manifest written to: ", manifest_path)
