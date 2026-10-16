profile <- list(
  platform = "xenium",
  platform_protocol = list(
    transcript_format = "csv",
    transcript_columns = list(
      transcript_id = "transcript_id",
      gene = "feature_name",
      x = "x_location",
      y = "y_location",
      z = "z_location",
      prior_cell_id = "cell_id",
      fov = "fov_name",
      qv = "qv",
      overlaps_nucleus = "overlaps_nucleus"
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
  spatial = list(
    k_neighbors = 10L,
    domain_clustering = list(
      enabled = TRUE,
      backend = "banksy",
      strict_labels = TRUE,
      allow_unlabeled = FALSE,
      banksy = list(
        k_geom = 20L,
        lambda = 0.6,
        resolutions = c(0.2),
        compute_agf = TRUE,
        use_agf = TRUE
      ),
      script = list(
        command = character(),
        extra_args = character(),
        workdir = NULL
      )
    ),
    marker_leakage = list(
      weighting = "uniform",
      distance_scale = NULL,
      scaling = "gene_rate"
    )
  ),
  expression = list(
    marker_scaling = "gene_rate"
  )
)
