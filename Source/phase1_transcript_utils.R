suppressPackageStartupMessages({
  library(data.table)
})

normalize_assignment_id <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "0" | toupper(x) == "UNASSIGNED" | tolower(x) %in% c("none", "null", "nan")] <- "-1"
  x
}

normalize_boolish <- function(x) {
  out <- rep(FALSE, length(x))
  if (is.logical(x)) {
    out[!is.na(x)] <- x[!is.na(x)]
    return(out)
  }
  txt <- trimws(tolower(as.character(x)))
  out[txt %in% c("1", "true", "t", "yes", "y")] <- TRUE
  out
}

optional_column <- function(dt, col_name, default = NA) {
  if (is.null(col_name) || !nzchar(col_name) || !col_name %in% names(dt)) {
    if (length(default) == nrow(dt)) {
      return(default)
    }
    return(rep(default, nrow(dt)))
  }
  dt[[col_name]]
}

write_table_auto <- function(dt, path) {
  data.table::fwrite(dt, path)
}

validate_transcript_ids <- function(dt) {
  if (!"transcript_id" %in% names(dt)) {
    stop("Standardized transcript table is missing `transcript_id`.")
  }
  if (anyNA(dt$transcript_id) || any(trimws(dt$transcript_id) == "")) {
    stop("Transcript IDs must be non-missing and non-empty after standardization.")
  }
  dup_ids <- unique(dt$transcript_id[duplicated(dt$transcript_id)])
  if (length(dup_ids) > 0) {
    stop("Transcript IDs must be unique after standardization. Example duplicates: ", paste(head(dup_ids, 5), collapse = ", "))
  }
  invisible(TRUE)
}

standardize_transcript_table <- function(raw_dt, column_cfg) {
  required_cols <- c(
    gene = column_cfg$gene,
    x = column_cfg$x,
    y = column_cfg$y
  )
  missing <- required_cols[is.na(required_cols) | !required_cols %in% names(raw_dt)]
  if (length(missing) > 0) {
    stop("Missing required transcript columns: ", paste(unique(missing), collapse = ", "))
  }

  out <- data.table(
    transcript_id = as.character(optional_column(raw_dt, column_cfg$transcript_id, default = seq_len(nrow(raw_dt)))),
    gene = as.character(raw_dt[[column_cfg$gene]]),
    x = suppressWarnings(as.numeric(raw_dt[[column_cfg$x]])),
    y = suppressWarnings(as.numeric(raw_dt[[column_cfg$y]])),
    z = suppressWarnings(as.numeric(optional_column(raw_dt, column_cfg$z, default = NA_real_))),
    cell_id_prior = normalize_assignment_id(optional_column(raw_dt, column_cfg$prior_cell_id, default = NA_character_)),
    fov = as.character(optional_column(raw_dt, column_cfg$fov, default = "")),
    qv = suppressWarnings(as.numeric(optional_column(raw_dt, column_cfg$qv, default = NA_real_))),
    overlaps_nucleus = normalize_boolish(optional_column(raw_dt, column_cfg$overlaps_nucleus, default = FALSE))
  )

  if (is.null(column_cfg$transcript_id) || !column_cfg$transcript_id %in% names(raw_dt)) {
    out[, transcript_id := paste0("generated_", seq_len(.N))]
  }

  out
}

apply_transcript_filters <- function(dt, filter_cfg) {
  keep <- rep(TRUE, nrow(dt))
  reasons <- vector("list", nrow(dt))

  add_reason <- function(mask, reason) {
    hit <- which(mask)
    if (!length(hit)) return(invisible(NULL))
    for (idx in hit) {
      reasons[[idx]] <<- c(reasons[[idx]], reason)
    }
    invisible(NULL)
  }

  if (isTRUE(filter_cfg$drop_blank_genes)) {
    bad <- is.na(dt$gene) | !nzchar(trimws(dt$gene))
    add_reason(bad, "blank_gene")
    keep <- keep & !bad
  }

  if (!is.null(filter_cfg$min_qv)) {
    bad <- is.na(dt$qv) | dt$qv < as.numeric(filter_cfg$min_qv)
    add_reason(bad, "low_qv")
    keep <- keep & !bad
  }

  if (!is.null(filter_cfg$exclude_gene_regex) && nzchar(filter_cfg$exclude_gene_regex)) {
    bad <- grepl(filter_cfg$exclude_gene_regex, dt$gene, perl = TRUE)
    bad[is.na(bad)] <- FALSE
    add_reason(bad, "excluded_gene_regex")
    keep <- keep & !bad
  }

  bad_coords <- !is.finite(dt$x) | !is.finite(dt$y)
  add_reason(bad_coords, "nonfinite_xy")
  keep <- keep & !bad_coords

  dt[, keep_for_eval := keep]
  dt[, drop_reason := vapply(reasons, function(x) {
    if (!length(x)) return(NA_character_)
    paste(unique(x), collapse = ";")
  }, character(1))]
  dt
}

standardize_transcripts_for_phase1 <- function(config, transcripts_dir) {
  dir.create(transcripts_dir, recursive = TRUE, showWarnings = FALSE)

  input_path <- config$canonical_transcripts$path
  column_cfg <- config$platform_protocol$transcript_columns
  filter_cfg <- config$platform_protocol$transcript_filters %||% list()

  message("[phase1] standardizing transcripts from ", input_path)
  raw_dt <- data.table::fread(input_path)
  source_dt <- standardize_transcript_table(raw_dt, column_cfg)
  validate_transcript_ids(source_dt)
  source_dt <- apply_transcript_filters(copy(source_dt), filter_cfg)
  eval_dt <- source_dt[keep_for_eval == TRUE]

  source_path <- file.path(transcripts_dir, "transcripts_source.csv.gz")
  eval_path <- file.path(transcripts_dir, "transcripts_eval.csv.gz")
  write_table_auto(source_dt, source_path)
  write_table_auto(eval_dt, eval_path)

  message(
    "[phase1] transcripts standardized: source_rows=", nrow(source_dt),
    " eval_rows=", nrow(eval_dt)
  )

  list(
    source_path = source_path,
    eval_path = eval_path,
    n_source_transcripts = nrow(source_dt),
    n_eval_transcripts = nrow(eval_dt),
    transcripts_eval_dt = eval_dt
  )
}
