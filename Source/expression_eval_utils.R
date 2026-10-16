suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

expression_theme_academic <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "#4d4d4d", size = base_size - 1),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "#f3f4f6", color = "#c7cbd1", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.35)
    )
}

read_markers_table <- function(path) {
  markers <- data.table::fread(path)
  required <- c("cell_type", "marker")
  missing <- setdiff(required, names(markers))
  if (length(missing) > 0) {
    stop("Markers CSV is missing required columns: ", paste(missing, collapse = ", "))
  }
  markers[, cell_type := as.character(cell_type)]
  markers[, marker := as.character(marker)]
  unique(markers[!is.na(cell_type) & !is.na(marker) & cell_type != "" & marker != ""])
}

normalize_cell_types <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- "Unknown"
  unknown_idx <- tolower(trimws(x)) %in% c("unknown", "unassigned", "none", "na")
  x[unknown_idx] <- "Unknown"
  x
}

resolve_marker_scaling_mainv2 <- function(scaling) {
  scaling <- tolower(as.character(scaling %||% "raw"))
  if (!scaling %in% c("raw", "gene_rate")) {
    stop("Unknown marker_scaling '", scaling, "'. Expected one of: raw, gene_rate.")
  }
  scaling
}

apply_marker_gene_rate_mainv2 <- function(marker_counts_dt, gene_universe_counts_dt,
                                          count_col = "count", gene_col = "gene") {
  gene_totals <- gene_universe_counts_dt[, .(..gene_total = sum(get(count_col), na.rm = TRUE)), by = c(gene_col)]
  out <- merge(marker_counts_dt, gene_totals, by = gene_col, all.x = TRUE)
  out[is.na(..gene_total) | ..gene_total <= 0, ..gene_total := 1]
  out[[count_col]] <- out[[count_col]] / out[["..gene_total"]]
  out[, ..gene_total := NULL]
  out
}

expression_metrics_mainv2 <- function(transcripts, cell_meta, markers_path,
                                      annotation_col = "cell_type_custom_final",
                                      marker_scaling = "raw") {
  marker_scaling <- resolve_marker_scaling_mainv2(marker_scaling)
  cell_meta <- as.data.table(cell_meta)
  transcripts <- as.data.table(transcripts)

  if (!all(c("cell_id", "gene", "count") %in% names(transcripts))) {
    stop("Transcripts table must contain columns: cell_id, gene, count")
  }
  if (!"cell_id" %in% names(cell_meta)) {
    stop("Cell metadata must contain column: cell_id")
  }
  if (!annotation_col %in% names(cell_meta)) {
    stop("Cell metadata is missing annotation column: ", annotation_col)
  }

  markers <- read_markers_table(markers_path)
  marker_union <- unique(markers$marker)

  transcripts[, cell_id := as.character(cell_id)]
  transcripts[, gene := as.character(gene)]
  transcripts[, count := as.numeric(count)]
  transcripts[cell_id %in% c("-1", ""), cell_id := NA_character_]
  transcripts[, assigned := !is.na(cell_id) & nzchar(cell_id)]

  meta_dt <- unique(cell_meta[, .(cell_id = as.character(cell_id), annot = normalize_cell_types(get(annotation_col)))])
  eval_cell_ids <- unique(meta_dt$cell_id)

  total_transcripts <- transcripts[, sum(count, na.rm = TRUE)]
  assigned_transcripts <- transcripts[assigned == TRUE & cell_id %in% eval_cell_ids, sum(count, na.rm = TRUE)]
  assignment <- data.table(
    total_transcripts = total_transcripts,
    assigned_transcripts = assigned_transcripts,
    background_transcripts = total_transcripts - assigned_transcripts,
    assignment_rate = fifelse(total_transcripts > 0, assigned_transcripts / total_transcripts, NA_real_)
  )

  assigned_dt <- transcripts[assigned == TRUE]
  assigned_dt <- assigned_dt[cell_id %in% eval_cell_ids]
  per_cell_counts <- assigned_dt[, .(n_total = sum(count, na.rm = TRUE)), by = cell_id]

  marker_dt <- assigned_dt[gene %in% marker_union, .(count = sum(count, na.rm = TRUE)), by = .(cell_id, gene)]

  marker_by_type_raw <- merge(marker_dt, markers, by.x = "gene", by.y = "marker", allow.cartesian = TRUE)
  marker_by_type_raw <- marker_by_type_raw[, .(marker_sum = sum(count, na.rm = TRUE)), by = .(cell_id, marker_type = cell_type)]
  marker_totals_raw <- marker_by_type_raw[, .(marker_total = sum(marker_sum, na.rm = TRUE)), by = cell_id]

  marker_dt_scaled <- if (marker_scaling == "gene_rate") {
    apply_marker_gene_rate_mainv2(marker_dt, assigned_dt[, .(gene, count)])
  } else {
    marker_dt
  }
  marker_by_type <- merge(marker_dt_scaled, markers, by.x = "gene", by.y = "marker", allow.cartesian = TRUE)
  marker_by_type <- marker_by_type[, .(marker_sum = sum(count, na.rm = TRUE)), by = .(cell_id, marker_type = cell_type)]
  marker_totals <- marker_by_type[, .(marker_total = sum(marker_sum, na.rm = TRUE)), by = cell_id]

  per_cell <- copy(meta_dt)
  per_cell <- merge(per_cell, per_cell_counts, by = "cell_id", all.x = TRUE)
  per_cell <- merge(per_cell, marker_totals_raw, by = "cell_id", all.x = TRUE)
  per_cell[is.na(n_total), n_total := 0]
  per_cell[is.na(marker_total), marker_total := 0]
  per_cell[, marker_zero := marker_total <= 0]
  per_cell[, n_marker := marker_total]
  per_cell[, marker_total := NULL]
  per_cell <- merge(per_cell, marker_totals, by = "cell_id", all.x = TRUE)
  per_cell[is.na(marker_total), marker_total := 0]
  per_cell[, annot := normalize_cell_types(annot)]
  setnames(per_cell, "annot", annotation_col)

  purity_base <- per_cell[get(annotation_col) != "Unknown"]
  purity_dt <- merge(
    purity_base[, .(cell_id, marker_total, annot = get(annotation_col))],
    marker_by_type,
    by = "cell_id",
    all.x = TRUE
  )
  purity_dt[is.na(marker_sum), marker_sum := 0]
  purity_dt[, purity_on_target := fifelse(marker_total > 0 & marker_type == annot, marker_sum / marker_total, NA_real_)]
  purity_dt <- purity_dt[marker_type == annot, .(cell_id, purity_on_target)]

  per_cell <- merge(per_cell, purity_dt, by = "cell_id", all.x = TRUE)
  per_cell[marker_zero == TRUE, purity_on_target := NA_real_]

  contam_base <- per_cell[get(annotation_col) != "Unknown"]
  contam_long <- merge(
    contam_base[, .(cell_id, annot = get(annotation_col), marker_total)],
    marker_by_type,
    by = "cell_id",
    all.x = TRUE
  )
  contam_long[is.na(marker_sum), marker_sum := 0]
  contam_long[, contam_fraction := fifelse(marker_total > 0, marker_sum / marker_total, NA_real_)]
  contam_long <- contam_long[marker_type != annot]
  setnames(contam_long, c("annot", "marker_type"), c(annotation_col, "contam_type"))

  counts_summary <- per_cell[, .(
    n_total_p25 = quantile(n_total, 0.25, na.rm = TRUE),
    n_total_median = quantile(n_total, 0.5, na.rm = TRUE),
    n_total_p75 = quantile(n_total, 0.75, na.rm = TRUE),
    n_marker_p25 = quantile(n_marker, 0.25, na.rm = TRUE),
    n_marker_median = quantile(n_marker, 0.5, na.rm = TRUE),
    n_marker_p75 = quantile(n_marker, 0.75, na.rm = TRUE),
    marker_zero_frac = mean(marker_zero, na.rm = TRUE)
  ), by = annotation_col]

  purity_summary <- per_cell[get(annotation_col) != "Unknown", .(
    purity_p25 = quantile(purity_on_target, 0.25, na.rm = TRUE),
    purity_median = quantile(purity_on_target, 0.5, na.rm = TRUE),
    purity_p75 = quantile(purity_on_target, 0.75, na.rm = TRUE),
    purity_mean = mean(purity_on_target, na.rm = TRUE),
    marker_zero_frac = mean(marker_zero, na.rm = TRUE)
  ), by = annotation_col]

  contam_summary <- contam_long[!is.na(contam_fraction), .(
    contam_p25 = quantile(contam_fraction, 0.25, na.rm = TRUE),
    contam_median = quantile(contam_fraction, 0.5, na.rm = TRUE),
    contam_p75 = quantile(contam_fraction, 0.75, na.rm = TRUE),
    contam_mean = mean(contam_fraction, na.rm = TRUE)
  ), by = c(annotation_col, "contam_type")]

  list(
    assignment = assignment,
    per_cell = per_cell,
    counts_summary = counts_summary,
    purity_summary = purity_summary,
    contam_long = contam_long,
    contam_summary = contam_summary,
    metadata = data.table(marker_scaling = marker_scaling)
  )
}

plot_expression_method_outputs <- function(metrics, method_name, annotation_col, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  per_cell <- copy(metrics$per_cell)
  contam_long <- copy(metrics$contam_long)

  per_cell[[annotation_col]] <- factor(per_cell[[annotation_col]])
  contam_long[[annotation_col]] <- factor(contam_long[[annotation_col]])

  p_counts <- ggplot(per_cell, aes(x = log1p(n_total))) +
    geom_histogram(bins = 40, fill = "#2c7fb8", alpha = 0.75, color = "white") +
    facet_wrap(as.formula(paste("~", annotation_col)), scales = "free_y") +
    labs(
      title = paste(method_name, "per-cell transcript counts"),
      subtitle = "Assigned transcript burden after unified transcript-to-polygon reassignment",
      x = "log1p transcripts per cell",
      y = "Cells"
    ) +
    expression_theme_academic()

  p_marker <- ggplot(per_cell, aes(x = log1p(n_marker))) +
    geom_histogram(bins = 40, fill = "#7fcdbb", alpha = 0.75, color = "white") +
    facet_wrap(as.formula(paste("~", annotation_col)), scales = "free_y") +
    labs(
      title = paste(method_name, "per-cell marker counts"),
      subtitle = "Marker-only transcript burden using user-defined marker panel",
      x = "log1p marker transcripts per cell",
      y = "Cells"
    ) +
    expression_theme_academic()

  purity_df <- per_cell[!is.na(purity_on_target) & get(annotation_col) != "Unknown"]
  if (nrow(purity_df)) {
    p_purity <- ggplot(purity_df, aes(x = .data[[annotation_col]], y = purity_on_target)) +
      geom_violin(fill = "#feb24c", alpha = 0.85, color = NA) +
      geom_boxplot(width = 0.12, outlier.size = 0.4, alpha = 0.6) +
      coord_flip() +
      labs(
        title = paste(method_name, "marker purity"),
        subtitle = "On-target marker fraction per annotated cell",
        x = NULL,
        y = "On-target purity"
      ) +
      expression_theme_academic()
    ggsave(file.path(outdir, "marker_purity.png"), p_purity, width = 7.2, height = 4.8, dpi = 320)
  }

  contam_df <- contam_long[!is.na(contam_fraction)]
  if (nrow(contam_df)) {
    p_contam <- ggplot(contam_df, aes(x = .data[[annotation_col]], y = contam_fraction)) +
      geom_violin(fill = "#f1695b", alpha = 0.65, color = NA) +
      geom_boxplot(width = 0.12, outlier.size = 0.35, alpha = 0.6) +
      facet_wrap(~contam_type) +
      coord_flip() +
      labs(
        title = paste(method_name, "marker contamination"),
        subtitle = "Off-target marker fraction stratified by contaminating cell-type marker set",
        x = NULL,
        y = "Off-target contamination"
      ) +
      expression_theme_academic()
    ggsave(file.path(outdir, "marker_contamination.png"), p_contam, width = 10.4, height = 6.4, dpi = 320)
  }

  ggsave(file.path(outdir, "per_cell_counts_log1p.png"), p_counts, width = 8.4, height = 5.2, dpi = 320)
  ggsave(file.path(outdir, "per_cell_marker_counts_log1p.png"), p_marker, width = 8.4, height = 5.2, dpi = 320)
}

plot_expression_cross_method_summary <- function(purity_dt, counts_dt, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(purity_dt)) {
    plot_dt <- purity_dt[cell_type_custom_final != "Unknown"]
    p_purity <- ggplot(plot_dt, aes(x = cell_type_custom_final, y = method, fill = purity_median)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.2f", purity_median)), size = 2.7) +
      scale_fill_gradient(low = "#f7fbff", high = "#084594", na.value = "#f1f1f1", limits = c(0, 1)) +
      labs(
        title = "Expression Purity Comparison",
        subtitle = "Median on-target marker purity by method and annotated cell type",
        x = NULL,
        y = NULL,
        fill = "Median purity"
      ) +
      expression_theme_academic() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggsave(file.path(outdir, "expression_purity_comparison.png"), p_purity, width = 8.6, height = 4.8, dpi = 320)
  }

  if (nrow(counts_dt)) {
    plot_dt <- counts_dt[cell_type_custom_final != "Unknown"]
    p_marker <- ggplot(plot_dt, aes(x = cell_type_custom_final, y = method, fill = n_marker_median)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = comma(round(n_marker_median, 0))), size = 2.5) +
      scale_fill_gradient(low = "#fff5eb", high = "#b30000", na.value = "#f1f1f1") +
      labs(
        title = "Marker Count Comparison",
        subtitle = "Median marker transcripts per cell by method and annotated cell type",
        x = NULL,
        y = NULL,
        fill = "Median markers"
      ) +
      expression_theme_academic() +
      theme(axis.text.x = element_text(angle = 40, hjust = 1))
    ggsave(file.path(outdir, "expression_marker_count_comparison.png"), p_marker, width = 8.8, height = 4.8, dpi = 320)
  }
}

run_expression_eval <- function(config) {
  output_root <- config$output_root
  eval_dir <- file.path(output_root, "expression_eval")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)

  cell_manifest <- fread(file.path(output_root, "cell_tables", "method_cell_output_manifest.csv"))
  annotation_manifest <- fread(file.path(output_root, "annotation", "annotation_manifest.csv"))

  results <- list()
  purity_rows <- list()
  count_rows <- list()
  idx <- 1L

  for (i in seq_len(nrow(cell_manifest))) {
    row <- cell_manifest[i]
    method_name <- row$method[[1]]
    if (!identical(row$status[[1]], "ready")) next

    annotation_row <- annotation_manifest[method == method_name]
    if (!nrow(annotation_row)) next

    counts_path <- file.path(output_root, "expression", paste0(method_name, "_counts_long.csv.gz"))
    annotation_path <- annotation_row$annotation_path[[1]]
    if (!file.exists(counts_path) || !file.exists(annotation_path)) next

    message("[expression] running method=", method_name)
    counts_dt <- fread(counts_path)
    annotation_dt <- fread(annotation_path)
    expression_cfg <- config$expression %||% list()
    metrics <- expression_metrics_mainv2(
      transcripts = counts_dt,
      cell_meta = annotation_dt,
      markers_path = config$markers_path,
      annotation_col = "cell_type_custom_final",
      marker_scaling = expression_cfg$marker_scaling %||% "raw"
    )
    assignment_out <- metrics$assignment

    prefix <- file.path(eval_dir, paste0(method_name, "_expression"))
    fwrite(assignment_out, paste0(prefix, "_assignment_rates.csv"))
    fwrite(metrics$per_cell, paste0(prefix, "_per_cell.csv"))
    fwrite(metrics$counts_summary, paste0(prefix, "_counts_by_celltype.csv"))
    fwrite(metrics$purity_summary, paste0(prefix, "_purity_by_celltype.csv"))
    fwrite(metrics$contam_summary, paste0(prefix, "_contamination_by_pair.csv"))
    fwrite(metrics$metadata, paste0(prefix, "_metadata.csv"))

    plot_expression_method_outputs(
      metrics = metrics,
      method_name = method_name,
      annotation_col = "cell_type_custom_final",
      outdir = paste0(prefix, "_plots")
    )

    purity_rows[[idx]] <- copy(metrics$purity_summary)[, method := method_name]
    count_rows[[idx]] <- copy(metrics$counts_summary)[, method := method_name]
    results[[idx]] <- data.table(
      method = method_name,
      assignment_rate = assignment_out$assignment_rate[[1]],
      n_cells = nrow(metrics$per_cell),
      n_cell_types = uniqueN(metrics$per_cell$cell_type_custom_final),
      marker_scaling = metrics$metadata$marker_scaling[[1]]
    )
    idx <- idx + 1L
  }

  purity_dt <- if (length(purity_rows)) rbindlist(purity_rows, fill = TRUE) else data.table()
  counts_dt <- if (length(count_rows)) rbindlist(count_rows, fill = TRUE) else data.table()
  if (nrow(purity_dt)) fwrite(purity_dt, file.path(eval_dir, "expression_purity_summary_all_methods.csv"))
  if (nrow(counts_dt)) fwrite(counts_dt, file.path(eval_dir, "expression_counts_summary_all_methods.csv"))
  plot_expression_cross_method_summary(purity_dt, counts_dt, eval_dir)

  manifest_dt <- if (length(results)) rbindlist(results, fill = TRUE) else data.table(
    method = character(),
    assignment_rate = numeric(),
    n_cells = integer(),
    n_cell_types = integer()
  )
  manifest_path <- file.path(eval_dir, "expression_manifest.csv")
  fwrite(manifest_dt, manifest_path)
  manifest_path
}
