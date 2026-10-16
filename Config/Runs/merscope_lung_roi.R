config <- list(
  dataset_id = "merscope_lung_roi",
  platform = "merscope",
  markers_path = "/users/9/feng0421/Code/R/Evaluation/Main/Extra/cell_type_markers_Merfish.csv",
  output_root = "/users/9/feng0421/Code/R/Evaluation/MainV2/Outputs/merscope_lung_roi",
  visualization = list(
    transcript_sample_n = 35000L,
    centroid_sample_n = 2500L,
    zoom_regions = list(
      list(name = "roi_a", bbox = c(620, 585, 690, 655)),
      list(name = "roi_b", bbox = c(550, 550, 800, 800))
    )
  ),
  canonical_transcripts = list(
    path = "/users/9/feng0421/BIDCell/data/MerFish/roi_detected_transcripts_filtered.csv",
    description = "Recommended cleaned canonical transcript universe for the ROI."
  ),
  methods = list(
    list(
      name = "baysor",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/users/9/feng0421/Code/R/Evaluation/Baysor_out/MERFISH_Baysor_polygons_2d.json"
      )
    ),
    list(
      name = "bidcell",
      enabled = TRUE,
      input = list(
        loader = "bidcell_geojson",
        path = "/users/9/feng0421/Code/R/Evaluation/BIDCell_out/model_outputs/2025_12_22_04_16_49/test_output/cells_output.json",
        affine_summary = "/users/9/feng0421/Code/R/Evaluation/BIDCell_out/affine.csv",
        seg_path = "/users/9/feng0421/Code/R/Evaluation/BIDCell_out/model_outputs/2025_12_22_04_16_49/test_output/epoch_1_step_4000_connected.tif",
        transform = list(type = "bidcell_crop_tiff")
      )
    ),
    list(
      name = "proseg",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/users/9/feng0421/Code/R/Evaluation/Proseg_out/merfish_roi/cell-polygons.geojson.gz"
      )
    ),
    list(
      name = "segger",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/users/9/feng0421/Code/R/Evaluation/Segger_out/segger_cell_boundaries.geojson.gz"
      )
    ),
    list(
      name = "cellpose",
      enabled = TRUE,
      input = list(
        loader = "geojson",
        path = "/users/9/feng0421/Code/R/Evaluation/Cellpose_out/cellpose_boundaries.geojson.gz"
      )
    )
  )
)
