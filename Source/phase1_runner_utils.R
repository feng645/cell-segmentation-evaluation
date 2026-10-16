suppressPackageStartupMessages({
  library(jsonlite)
})

build_manifest <- function(config, transcript_meta, assignments_meta, cell_manifest_path) {
  list(
    dataset_id = config$dataset_id,
    platform = config$platform,
    run_config_path = config$run_config_path,
    profile_path = config$profile_path,
    canonical_transcripts = list(
      input_path = config$canonical_transcripts$path,
      transcripts_source_path = transcript_meta$source_path,
      transcripts_eval_path = transcript_meta$eval_path,
      n_source_transcripts = transcript_meta$n_source_transcripts,
      n_eval_transcripts = transcript_meta$n_eval_transcripts
    ),
    assignments = list(
      transcripts_with_assignments = assignments_meta$assignments_path,
      assignment_summary = assignments_meta$summary_path
    ),
    cell_output_manifest = cell_manifest_path,
    methods = lapply(assignments_meta$method_records, function(method) {
      list(
        name = method$name,
        status = method$status %||% NULL,
        reason = method$reason %||% NULL,
        mode = method$mode %||% NULL,
        loader = method$loader %||% NULL,
        polygon_path = method$polygon_path %||% NULL,
        n_cells_raw = method$n_cells_raw %||% 0L
      )
    })
  )
}

run_phase1_eval_r <- function(config) {
  output_root <- config$output_root
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

  transcripts_dir <- file.path(output_root, "transcripts")
  polygons_dir <- file.path(output_root, "polygons")
  assignments_dir <- file.path(output_root, "assignments")
  expression_dir <- file.path(output_root, "expression")
  cell_tables_dir <- file.path(output_root, "cell_tables")
  dir.create(transcripts_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(polygons_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(assignments_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(expression_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cell_tables_dir, recursive = TRUE, showWarnings = FALSE)

  transcript_meta <- standardize_transcripts_for_phase1(config, transcripts_dir)
  assignments_meta <- build_assignments_for_phase1(config, transcript_meta, polygons_dir, assignments_dir)
  cell_manifest_path <- build_cell_outputs_for_phase1(config, assignments_meta, expression_dir, cell_tables_dir)

  manifest <- build_manifest(config, transcript_meta, assignments_meta, cell_manifest_path)
  manifest_path <- file.path(output_root, "manifest.json")
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  manifest_path
}
