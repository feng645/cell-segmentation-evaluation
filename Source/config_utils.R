suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

load_config_object <- function(path, object_name = "config") {
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Expected object `", object_name, "` in ", path)
  }
  get(object_name, envir = env, inherits = FALSE)
}

deep_merge_lists <- function(base, override) {
  if (is.null(base)) return(override)
  if (is.null(override)) return(base)
  if (!is.list(base) || !is.list(override)) return(override)

  base_names <- names(base)
  override_names <- names(override)

  if (is.null(base_names) || is.null(override_names)) {
    return(override)
  }

  merged <- base
  for (nm in override_names) {
    merged[[nm]] <- deep_merge_lists(base[[nm]], override[[nm]])
  }
  merged
}

path_if_present <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) {
    return(NULL)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

normalize_method_config <- function(method_cfg) {
  method_cfg$name <- tolower(method_cfg$name)
  if (is.null(method_cfg$enabled)) {
    method_cfg$enabled <- TRUE
  }
  if (!isTRUE(method_cfg$enabled)) {
    return(method_cfg)
  }
  if (is.null(method_cfg$input)) {
    stop("Method `", method_cfg$name, "` is missing `input` configuration.")
  }

  input_fields <- c("path", "seg_path", "affine_summary", "affine_matrix")
  for (field in input_fields) {
    if (!is.null(method_cfg$input[[field]])) {
      method_cfg$input[[field]] <- path_if_present(method_cfg$input[[field]])
    }
  }

  method_cfg
}

validate_markers_csv <- function(path) {
  marker_cols <- names(data.table::fread(path, nrows = 0))
  required <- c("cell_type", "marker")
  missing <- setdiff(required, marker_cols)
  if (length(missing) > 0) {
    stop("Markers CSV is missing required columns: ", paste(missing, collapse = ", "))
  }
  marker_dt <- data.table::fread(path)
  dup_n <- marker_dt[, .N, by = .(cell_type, marker)][N > 1, sum(N - 1)]
  if (is.finite(dup_n) && dup_n > 0) {
    warning("Markers CSV contains duplicated (cell_type, marker) rows: ", dup_n, " duplicates will be ignored downstream.")
  }
  invisible(TRUE)
}

compile_run_config <- function(run_config_path) {
  run_config_path <- normalizePath(run_config_path, winslash = "/", mustWork = TRUE)
  mainv2_root <- normalizePath(file.path(dirname(run_config_path), "..", ".."), winslash = "/", mustWork = TRUE)

  run_cfg <- load_config_object(run_config_path, "config")
  if (is.null(run_cfg$platform) || !nzchar(run_cfg$platform)) {
    stop("Run config must define `platform`.")
  }

  profile_path <- file.path(mainv2_root, "Config", "Profiles", paste0(tolower(run_cfg$platform), "_profile.R"))
  if (!file.exists(profile_path)) {
    stop("Platform profile not found: ", profile_path)
  }
  profile_cfg <- load_config_object(profile_path, "profile")

  compiled <- deep_merge_lists(profile_cfg, run_cfg)
  compiled$mainv2_root <- mainv2_root
  compiled$run_config_path <- run_config_path
  compiled$profile_path <- normalizePath(profile_path, winslash = "/", mustWork = TRUE)

  required_top_level <- c("dataset_id", "platform", "markers_path", "canonical_transcripts", "methods", "output_root")
  missing_top_level <- required_top_level[!vapply(required_top_level, function(nm) !is.null(compiled[[nm]]), logical(1))]
  if (length(missing_top_level) > 0) {
    stop("Run config is missing required fields: ", paste(missing_top_level, collapse = ", "))
  }

  compiled$output_root <- path_if_present(compiled$output_root)
  compiled$markers_path <- path_if_present(compiled$markers_path)
  compiled$canonical_transcripts$path <- path_if_present(compiled$canonical_transcripts$path)

  if (!file.exists(compiled$markers_path)) {
    stop("Markers CSV not found: ", compiled$markers_path)
  }
  validate_markers_csv(compiled$markers_path)
  if (!file.exists(compiled$canonical_transcripts$path)) {
    stop("Canonical transcripts file not found: ", compiled$canonical_transcripts$path)
  }

  if (is.null(compiled$runtime)) compiled$runtime <- list()
  if (is.null(compiled$runtime$python)) {
    compiled$runtime$python <- "/users/9/feng0421/.conda/envs/segger/bin/python"
  }
  compiled$runtime$python <- path_if_present(compiled$runtime$python)
  if (is.null(compiled$runtime$chunk_size)) {
    compiled$runtime$chunk_size <- 250000L
  }

  compiled$methods <- lapply(compiled$methods, normalize_method_config)
  method_names <- vapply(compiled$methods, `[[`, character(1), "name")
  if (anyDuplicated(method_names) > 0) {
    stop("Method names must be unique. Duplicates: ", paste(unique(method_names[duplicated(method_names)]), collapse = ", "))
  }

  compiled
}

write_compiled_config <- function(config, path) {
  jsonlite::write_json(
    x = config,
    path = path,
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null",
    na = "null"
  )
}
