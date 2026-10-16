suppressPackageStartupMessages({
  library(data.table)
})

build_method_counts <- function(assignments_dt, method_name) {
  cell_col <- paste0("cell_id_", method_name)
  if (!cell_col %in% names(assignments_dt)) {
    stop("Assignment column not found: ", cell_col)
  }
  dt <- assignments_dt[get(cell_col) != "-1", .(
    cell_id = as.character(get(cell_col)),
    gene = as.character(gene)
  )]
  if (nrow(dt) == 0L) {
    return(data.table(cell_id = character(), gene = character(), count = integer()))
  }
  dt[, .(count = .N), by = .(cell_id, gene)]
}

compute_cell_drop_reason <- function(cell_dt, min_genes, max_genes, max_pct_mt) {
  reasons <- vector("list", nrow(cell_dt))

  add_reason <- function(mask, reason) {
    hit <- which(mask)
    if (!length(hit)) return(invisible(NULL))
    for (idx in hit) {
      reasons[[idx]] <<- c(reasons[[idx]], reason)
    }
    invisible(NULL)
  }

  add_reason(cell_dt$n_genes < min_genes, "low_genes")
  add_reason(cell_dt$n_genes > max_genes, "high_genes")
  add_reason(fifelse(is.na(cell_dt$pct_mt), FALSE, cell_dt$pct_mt > max_pct_mt), "high_pct_mt")

  vapply(reasons, function(x) {
    if (!length(x)) return(NA_character_)
    paste(unique(x), collapse = ";")
  }, character(1))
}

build_cell_outputs_for_phase1 <- function(config, assignments_meta, expression_dir, cell_tables_dir) {
  dir.create(expression_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cell_tables_dir, recursive = TRUE, showWarnings = FALSE)

  assignments_dt <- assignments_meta$assignments_dt
  method_records <- assignments_meta$method_records
  filter_cfg <- config$platform_protocol$cell_filters %||% list()
  mt_regex <- config$platform_protocol$mt_gene_regex %||% "(?i)^mt[-_]"
  min_genes <- as.integer(filter_cfg$min_genes %||% 30L)
  max_genes <- as.integer(filter_cfg$max_genes %||% 8000L)
  max_pct_mt <- as.numeric(filter_cfg$max_pct_mt %||% 15)

  manifest_rows <- list()

  for (method in method_records) {
    method_name <- method$name
    if (!identical(method$status, "ready")) {
      manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
        method = method_name,
        status = method$status %||% "unknown",
        reason = method$reason %||% NA_character_,
        n_cells_raw = 0L,
        n_cells_eval = 0L,
        counts_long_path = NA_character_,
        cell_table_path = NA_character_,
        raw_cell_table_path = NA_character_,
        eval_cell_table_path = NA_character_
      )
      next
    }

    message("[phase1] aggregating cell outputs for method=", method_name)
    counts_dt <- build_method_counts(assignments_dt, method_name)
    counts_path <- file.path(expression_dir, paste0(method_name, "_counts_long.csv.gz"))
    data.table::fwrite(counts_dt, counts_path)

    if (nrow(counts_dt) == 0L) {
      cell_dt <- data.table(
        cell_id = character(),
        n_transcripts = integer(),
        n_genes = integer(),
        n_mt_transcripts = integer(),
        pct_mt = numeric(),
        cell_area = numeric(),
        x = numeric(),
        y = numeric(),
        z = numeric(),
        keep_for_eval = logical(),
        drop_reason = character()
      )
    } else {
      counts_dt[, is_mt := grepl(mt_regex, gene, perl = TRUE)]
      cell_dt <- counts_dt[, .(
        n_transcripts = sum(count),
        n_genes = .N,
        n_mt_transcripts = sum(count[is_mt], na.rm = TRUE)
      ), by = cell_id]
      cell_dt[, pct_mt := fifelse(n_transcripts > 0, 100 * n_mt_transcripts / n_transcripts, NA_real_)]
      if (!is.null(method$cell_areas)) {
        cell_dt[, cell_area := as.numeric(method$cell_areas[cell_id])]
      } else {
        cell_dt[, cell_area := NA_real_]
      }
      if (!is.null(method$cell_centroids) && inherits(method$cell_centroids, "data.table")) {
        cell_dt <- merge(cell_dt, method$cell_centroids, by = "cell_id", all.x = TRUE)
      } else {
        cell_dt[, `:=`(x = NA_real_, y = NA_real_)]
      }
      cell_dt[, z := NA_real_]
      cell_dt[, keep_for_eval := n_genes >= min_genes & n_genes <= max_genes & fifelse(is.na(pct_mt), 0, pct_mt) <= max_pct_mt]
      cell_dt[, drop_reason := compute_cell_drop_reason(.SD, min_genes, max_genes, max_pct_mt)]
      cell_dt[keep_for_eval == TRUE, drop_reason := NA_character_]
    }

    raw_path <- file.path(cell_tables_dir, paste0(method_name, "_cells_raw.csv.gz"))
    eval_path <- file.path(cell_tables_dir, paste0(method_name, "_cells_eval.csv.gz"))
    data.table::fwrite(cell_dt, raw_path)
    data.table::fwrite(cell_dt[keep_for_eval == TRUE], eval_path)

    manifest_rows[[length(manifest_rows) + 1L]] <- data.table(
      method = method_name,
      status = "ready",
      reason = NA_character_,
      n_cells_raw = nrow(cell_dt),
      n_cells_eval = sum(cell_dt$keep_for_eval, na.rm = TRUE),
      counts_long_path = counts_path,
      cell_table_path = raw_path,
      raw_cell_table_path = raw_path,
      eval_cell_table_path = eval_path
    )
  }

  manifest_dt <- rbindlist(manifest_rows, fill = TRUE)
  manifest_path <- file.path(cell_tables_dir, "method_cell_output_manifest.csv")
  data.table::fwrite(manifest_dt, manifest_path)
  manifest_path
}
