# Same as merscope_lung_roi_gene_rate.R, but with Baysor rerun using the correct
# unassigned prior label ("-1"). Output goes to the shared data directory.
sys.source(
  "/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/merscope_lung_roi_gene_rate.R",
  envir = environment()
)

config$output_root <- "/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix/Evaluation/merscope_lung_roi_gene_rate"

for (i in seq_along(config$methods)) {
  if (identical(config$methods[[i]]$name, "baysor")) {
    config$methods[[i]]$input$path <- "/projects/standard/zhangl4/shared/Data_feng0421/Baysor_rerun_unassigned_fix/MERSCOPE_ROI/segmentation_polygons_2d.json"
  }
}
