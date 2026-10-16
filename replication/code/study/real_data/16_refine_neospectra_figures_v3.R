# Q1-oriented, reporting-only refinements of the NeoSpectra v3 figures.
#
# This script deliberately reads the frozen analysis object and its exported
# CSV tables. It neither refits a model nor modifies any analysis result. The
# original v3 figures are retained; refined outputs use the `_v3_q1` suffix.

script_dir <- local({
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argument)) {
    dirname(normalizePath(sub("^--file=", "", argument[1L])))
  } else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required.", call. = FALSE)
}

overwrite <- "--overwrite" %in% commandArgs(trailingOnly = TRUE)
results_dir <- file.path(
  project_root, "real_data", "results", "v3", "neospectra"
)
figures_dir <- file.path(
  project_root, "real_data", "figures", "v3", "neospectra"
)

input_paths <- c(
  analysis = file.path(results_dir, "analysis_v3.rds"),
  prediction = file.path(results_dir, "prediction_metrics_by_task_v3.csv"),
  coefficients = file.path(
    results_dir, "full_data_coefficient_curves_v3.csv"
  ),
  stability = file.path(results_dir, "stability_summary_v3.csv"),
  simulation = file.path(
    results_dir, "real_design_simulation_paired_summary_v3.csv"
  )
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop(
    "Required frozen NeoSpectra v3 inputs are missing: ",
    paste(basename(missing_inputs), collapse = ", "),
    call. = FALSE
  )
}

analysis <- readRDS(input_paths[["analysis"]])
if (!identical(analysis$version, "neospectra_v3") ||
    !identical(analysis$mode, "full")) {
  stop("analysis_v3.rds is not the full NeoSpectra v3 analysis.",
       call. = FALSE)
}

read_export <- function(path) {
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
prediction_task <- read_export(input_paths[["prediction"]])
full_curves <- read_export(input_paths[["coefficients"]])
stability_summary <- read_export(input_paths[["stability"]])
simulation_paired_summary <- read_export(input_paths[["simulation"]])

assert_export_matches <- function(label, exported, frozen) {
  if (!identical(names(exported), names(frozen)) ||
      nrow(exported) != nrow(frozen) ||
      !isTRUE(all.equal(
        exported, frozen,
        tolerance = 1e-12, check.attributes = FALSE
      ))) {
    stop(label, " CSV does not match analysis_v3.rds.", call. = FALSE)
  }
  invisible(TRUE)
}

assert_export_matches(
  "Prediction", prediction_task, analysis$prediction$task
)
assert_export_matches(
  "Coefficient", full_curves, analysis$full$curves
)
assert_export_matches(
  "Stability", stability_summary, analysis$stability$summary
)
assert_export_matches(
  "Semi-simulation", simulation_paired_summary,
  analysis$simulation_paired_summary
)

required_columns <- list(
  prediction = c("evaluation", "task_id", "method", "RMSE"),
  coefficients = c("method", "wavelength", "beta"),
  stability = c(
    "source", "method", "wavelength", "successful_fits", "beta_q10",
    "beta_q90", "zero_selection_frequency"
  ),
  simulation = c(
    "scenario", "n_train", "snr_label", "ise_advantage_mean",
    "ise_advantage_se"
  )
)
input_tables <- list(
  prediction = prediction_task,
  coefficients = full_curves,
  stability = stability_summary,
  simulation = simulation_paired_summary
)
for (table_name in names(required_columns)) {
  missing_columns <- setdiff(
    required_columns[[table_name]], names(input_tables[[table_name]])
  )
  if (length(missing_columns)) {
    stop(
      table_name, " input is missing columns: ",
      paste(missing_columns, collapse = ", "), call. = FALSE
    )
  }
}

bootstrap_summary <- stability_summary[
  stability_summary$source == "bootstrap", , drop = FALSE
]
if (!nrow(bootstrap_summary)) {
  stop("No bootstrap stability rows are available.", call. = FALSE)
}
coefficient_plot_data <- merge(
  full_curves,
  bootstrap_summary[c(
    "method", "wavelength", "beta_q10", "beta_q90",
    "zero_selection_frequency"
  )],
  by = c("method", "wavelength"), all.x = TRUE, sort = TRUE
)
if (anyNA(coefficient_plot_data[c("beta_q10", "beta_q90")])) {
  stop("Coefficient bands are incomplete after the reporting merge.",
       call. = FALSE)
}

method_colours <- c(FMCP = "#0072B2", fSCAD = "#D55E00")
coefficient_plot <- ggplot2::ggplot(
  coefficient_plot_data,
  ggplot2::aes(x = wavelength, colour = method, fill = method)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = beta_q10, ymax = beta_q90),
    alpha = 0.18, colour = NA
  ) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  ggplot2::geom_line(ggplot2::aes(y = beta), linewidth = 0.8) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::scale_fill_manual(values = method_colours) +
  ggplot2::labs(
    x = "NIR wavelength (nm)",
    y = "Coefficient per nm (log1p(SOC) / reflectance)",
    colour = NULL, fill = NULL,
    title = "Full-data coefficient functions with 80% bootstrap bands",
    subtitle = "Bands are conditional on the fixed nested-CV tuning choice"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    plot.title.position = "plot",
    plot.subtitle = ggplot2::element_text(colour = "grey35"),
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

zero_values <- bootstrap_summary$zero_selection_frequency
all_zero_stability <- length(zero_values) > 0L &&
  all(is.finite(zero_values)) &&
  all(abs(zero_values) <= sqrt(.Machine$double.eps))
zero_plot <- ggplot2::ggplot(
  bootstrap_summary,
  ggplot2::aes(
    x = wavelength, y = zero_selection_frequency, colour = method
  )
) +
  ggplot2::geom_hline(
    yintercept = c(0.5, 0.8), linetype = 3, colour = "grey65"
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::coord_cartesian(ylim = c(0, 1), clip = "off") +
  ggplot2::labs(
    x = "NIR wavelength (nm)",
    y = "Conditional bootstrap zero-selection frequency",
    colour = NULL,
    title = "Local-zero stability conditional on fixed CV tuning"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "top",
    plot.title.position = "plot",
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )
if (all_zero_stability) {
  zero_plot <- zero_plot + ggplot2::annotate(
    "label",
    x = mean(range(bootstrap_summary$wavelength, finite = TRUE)),
    y = 0.88,
    label = paste0(
      "All reported zero-selection frequencies are 0\n",
      "for FMCP and fSCAD under fixed CV tuning"
    ),
    size = 3.6, colour = "grey20", fill = "white",
    linewidth = 0.25
  )
}

cv_plot_data <- prediction_task[
  prediction_task$evaluation == "outer_cv", , drop = FALSE
]
if (!nrow(cv_plot_data)) {
  stop("No outer-CV prediction rows are available.", call. = FALSE)
}
cv_plot <- ggplot2::ggplot(
  cv_plot_data,
  ggplot2::aes(x = method, y = RMSE, colour = method)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62) +
  ggplot2::geom_point(
    ggplot2::aes(group = task_id),
    position = ggplot2::position_jitter(
      width = 0.08, height = 0,
      seed = analysis$configuration$master_seed
    ),
    alpha = 0.65, size = 1.5
  ) +
  ggplot2::labs(
    x = NULL, y = "Held-out RMSE (log1p(SOC))", colour = NULL,
    title = "Grouped nested-CV prediction error by outer fold"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "none")

simulation_plot <- ggplot2::ggplot(
  simulation_paired_summary,
  ggplot2::aes(
    x = n_train, y = ise_advantage_mean, group = snr_label,
    colour = snr_label
  )
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = ise_advantage_mean - 1.96 * ise_advantage_se,
      ymax = ise_advantage_mean + 1.96 * ise_advantage_se
    ),
    width = 0
  ) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ scenario, scales = "free_y") +
  ggplot2::labs(
    x = "Training sample size",
    y = "ISE(fSCAD) - ISE(FMCP)", colour = "SNR",
    title = "Real-design semi-simulation: positive values favour FMCP"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "top")

output_spec <- data.frame(
  stem = c(
    "coefficient_bootstrap", "zero_selection_stability",
    "nested_cv_rmse", "simulation_ise_advantage"
  ),
  width = c(7.4, 7.2, 6.6, 8.2),
  height = c(4.9, 4.4, 4.4, 5.4),
  stringsAsFactors = FALSE
)
plots <- list(
  coefficient_bootstrap = coefficient_plot,
  zero_selection_stability = zero_plot,
  nested_cv_rmse = cv_plot,
  simulation_ise_advantage = simulation_plot
)
output_paths <- unlist(lapply(output_spec$stem, function(stem) {
  file.path(figures_dir, paste0(stem, "_v3_q1.", c("pdf", "png")))
}), use.names = FALSE)
if (!overwrite && any(file.exists(output_paths))) {
  stop(
    "Refined Q1 output already exists; use --overwrite to replace only ",
    "the `_v3_q1` derivatives.", call. = FALSE
  )
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
for (index in seq_len(nrow(output_spec))) {
  stem <- output_spec$stem[index]
  plot <- plots[[stem]]
  pdf_path <- file.path(figures_dir, paste0(stem, "_v3_q1.pdf"))
  png_path <- file.path(figures_dir, paste0(stem, "_v3_q1.png"))
  ggplot2::ggsave(
    pdf_path, plot,
    width = output_spec$width[index], height = output_spec$height[index],
    bg = "white"
  )
  ggplot2::ggsave(
    png_path, plot,
    width = output_spec$width[index], height = output_spec$height[index],
    dpi = 320, bg = "white"
  )
}

if (any(!file.exists(output_paths)) || any(file.info(output_paths)$size <= 0)) {
  stop("One or more refined figure outputs are missing or empty.",
       call. = FALSE)
}

cat("NeoSpectra v3 Q1 reporting figures written without refitting.\n")
cat("Frozen CSV/RDS consistency checks: 4/4 passed.\n")
cat(
  "All zero-selection frequencies equal zero: ",
  if (all_zero_stability) "yes (annotation added)" else "no", "\n",
  sep = ""
)
cat(paste0(normalizePath(output_paths), "\n"))
