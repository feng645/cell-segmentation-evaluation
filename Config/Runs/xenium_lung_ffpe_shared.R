config <- list(
  dataset_id = "xenium_lung_ffpe",
  platform = "xenium",
  markers_path = "/users/9/feng0421/Code/R/Evaluation/Main/Extra/cell_type_markers_Xenium.csv",
  output_root = "/projects/standard/zhangl4/shared/Data_feng0421/Evaluation_MainV2/xenium_lung_ffpe",
  runtime = list(
    chunk_size = 100000L
  ),
  visualization = list(
    transcript_sample_n = 500000L,
    centroid_sample_n = 20000L,
    zoom_regions = list(
      list(name = "roi_a", bbox = c(3000, 1500, 3500, 2000)),
      list(name = "roi_b", bbox = c(3000, 1500, 4000, 2500))
    )
  ),
  canonical_transcripts = list(
    path = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_V1_humanLung_Cancer_FFPE_outs/transcripts_cleaned.csv",
    description = "Recommended cleaned canonical transcript universe for Xenium."
  ),
  methods = list(
    list(
      name = "baysor",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_Baysor/segmentation_polygons_2d.json"
      )
    ),
    list(
      name = "bidcell",
      enabled = TRUE,
      input = list(
        loader = "bidcell_geojson",
        path = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_BIDCell/model_outputs/2026_04_04_00_11_19/test_output/cells_output.json",
        transform = list(type = "identity")
      )
    ),
    list(
      name = "proseg",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_Proseg/cell-polygons.geojson.gz"
      )
    ),
    list(
    name = "segger",
    enabled = TRUE,
    input = list(
      loader = "geojson",
      path = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_segger/predictions/lung_cancer_v1_pred_0.5_True_4_6.0_5_5.0_20260406/segger_cell_boundaries.geojson.gz"
    )
),
    list(
      name = "cellpose",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/projects/standard/zhangl4/shared/Data_feng0421/Xenium_Cellpose/cellpose_boundaries.geojson.gz"
      )
    )
  )
)
