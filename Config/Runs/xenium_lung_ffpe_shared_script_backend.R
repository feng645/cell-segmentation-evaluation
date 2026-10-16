sys.source(
  "/users/9/feng0421/Code/R/Evaluation/MainV2/Config/Runs/xenium_lung_ffpe_shared.R",
  envir = environment()
)

config$output_root <- "/projects/standard/zhangl4/shared/Data_feng0421/Evaluation_MainV2/xenium_lung_ffpe_script_backend"
config$spatial <- list(
  domain_clustering = list(
    enabled = TRUE,
    backend = "script",
    script = list(
      command = c(
        "bash",
        "/users/9/feng0421/Code/R/Evaluation/MainV2/Scripts/example_domain_backend_banksy.sh"
      ),
      extra_args = character(),
      workdir = "/users/9/feng0421/Code/R/Evaluation/MainV2"
    )
  )
)
