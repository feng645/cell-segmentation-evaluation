profile <- list(
  platform = "cosmx",
  platform_protocol = list(
    transcript_format = "csv",
    transcript_columns = list(
      transcript_id = NULL,
      gene = NULL,
      x = NULL,
      y = NULL,
      z = NULL,
      prior_cell_id = NULL,
      fov = NULL,
      qv = NULL,
      overlaps_nucleus = NULL
    ),
    transcript_filters = list(
      min_qv = NULL,
      drop_blank_genes = TRUE,
      exclude_gene_regex = NULL
    ),
    cell_filters = list(
      min_genes = 30L,
      max_genes = 8000L,
      max_pct_mt = 15.0
    ),
    mt_gene_regex = "(?i)^mt[-_]"
  ),
  runtime = list(
    python = "/users/9/feng0421/.conda/envs/segger/bin/python",
    chunk_size = 250000L
  ),
  notes = "Placeholder profile reserved for future CosMx support."
)
