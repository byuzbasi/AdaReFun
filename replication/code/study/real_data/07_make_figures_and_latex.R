script_dir <- local({
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) dirname(normalizePath(sub("^--file=", "", file_argument[1L]))) else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
for (file in c("utils.R", "config.R", "data.R")) {
  source(file.path(script_dir, "R", file))
}
activate_fmcp_libraries()
assert_packages("ggplot2")
font_cache <- file.path(tempdir(), "fmcp-real-fontconfig")
dir.create(font_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(XDG_CACHE_HOME = font_cache)

mode <- real_mode()
config <- apply_resolved_grid(real_data_config(mode, project_root))
datasets <- load_processed_datasets(config)
analysis_path <- file.path(
  config$results_dir, versioned_name("analysis", mode, "rds")
)
full_path <- file.path(
  config$results_dir, versioned_name("full_models", mode, "rds")
)
manifest_path <- file.path(config$results_dir, "data_manifest_v1.csv")
if (!all(file.exists(c(analysis_path, full_path, manifest_path)))) {
  stop("Analysis, full-model, or manifest input is missing.", call. = FALSE)
}
analysis <- readRDS(analysis_path)
full <- readRDS(full_path)
manifest <- utils::read.csv(manifest_path, check.names = FALSE)

dir.create(config$figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$outputs_dir, recursive = TRUE, showWarnings = FALSE)

figure_stems <- c(
  "real_data_profiles", "cv_performance", "primary_coefficients",
  "stability_null_probability", "paired_fmcp_fscad",
  "permutation_evidence", "gamma_sensitivity"
)
figure_paths <- unlist(lapply(figure_stems, function(stem) {
  c(
    file.path(config$figures_dir, versioned_name(stem, mode, "pdf")),
    file.path(config$figures_dir, versioned_name(stem, mode, "png"))
  )
}))
table_stems <- c(
  "table_data", "table_cv_performance", "table_paired_fscad",
  "table_full_sparse", "table_stable_segments", "table_permutation",
  "table_gamma", "table_validation"
)
table_paths <- setNames(vapply(
  table_stems,
  function(stem) file.path(config$outputs_dir,
                           versioned_name(stem, mode, "tex")),
  FUN.VALUE = character(1L)
), table_stems)
tex_path <- file.path(
  config$outputs_dir,
  versioned_name("real_data_application", mode, "tex")
)
pdf_path <- file.path(
  config$outputs_dir,
  versioned_name("real_data_application", mode, "pdf")
)
report_readme <- file.path(
  config$outputs_dir,
  versioned_name("README", mode, "md")
)
ensure_new_output(
  c(figure_paths, table_paths, tex_path, pdf_path, report_readme),
  overwrite = real_overwrite()
)

method_levels <- config$methods
method_colors <- c(
  FMCP = "#0072B2", fSCAD = "#D55E00", Smooth = "#009E73",
  PCR = "#CC79A7", PenS = "#E69F00", NGR = "#56B4E9",
  FLiRTI = "#8C564B", OLS = "#7F7F7F"
)
theme_real <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(color = "grey35"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank()
    )
}

save_figure <- function(plot, stem, width = 8.2, height = 5.2) {
  pdf <- file.path(config$figures_dir, versioned_name(stem, mode, "pdf"))
  png <- file.path(config$figures_dir, versioned_name(stem, mode, "png"))
  ggplot2::ggsave(pdf, plot, width = width, height = height,
                  units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(png, plot, width = width, height = height,
                  units = "in", dpi = 320, bg = "white")
}

profile_rows <- list()
profile_index <- 1L
for (dataset_index in seq_along(datasets)) {
  dataset <- datasets[[dataset_index]]
  set.seed(task_seed(config, "analysis", dataset_index, 900L))
  chosen <- sort(sample(
    seq_along(dataset$y), min(24L, length(dataset$y)), replace = FALSE
  ))
  profile_rows[[profile_index]] <- data.frame(
    dataset = dataset$name,
    subject_id = rep(dataset$subject_id[chosen], each = ncol(dataset$x)),
    argval_physical = rep(dataset$argvals_physical, times = length(chosen)),
    value = as.vector(t(dataset$x[chosen, , drop = FALSE])),
    stringsAsFactors = FALSE
  )
  profile_index <- profile_index + 1L
}
profiles <- bind_rows_base(profile_rows)
profile_plot <- ggplot2::ggplot(
  profiles,
  ggplot2::aes(argval_physical, value, group = subject_id)
) +
  ggplot2::geom_line(color = "#336B87", alpha = 0.25, linewidth = 0.45) +
  ggplot2::stat_summary(
    ggplot2::aes(group = 1), fun = mean, geom = "line",
    color = "#D55E00", linewidth = 1.05
  ) +
  ggplot2::facet_wrap(~dataset, scales = "free", ncol = 2) +
  ggplot2::labs(
    title = "Functional predictors in the real-data applications",
    subtitle = "Twenty-four reproducibly sampled curves; orange line is the sample mean",
    x = "Physical argument", y = "Functional predictor"
  ) + theme_real()
save_figure(profile_plot, "real_data_profiles")

performance <- analysis$performance
performance$method <- factor(performance$method, levels = rev(method_levels))
cv_plot <- ggplot2::ggplot(
  performance,
  ggplot2::aes(RMSE, method, color = method)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = RMSE_lower, xmax = RMSE_upper),
    orientation = "y", width = 0.15, linewidth = 0.65
  ) +
  ggplot2::geom_point(size = 2.6) +
  ggplot2::facet_wrap(~dataset, scales = "free_x") +
  ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
  ggplot2::labs(
    title = "Repeated cross-validation prediction error",
    subtitle = paste0(
      config$cv_repeats, " repeated ", config$cv_folds,
      "-fold CV; bars are 95% subject-bootstrap intervals"
    ),
    x = "RMSE (lower is better)", y = NULL
  ) + theme_real() + ggplot2::theme(legend.position = "none")
save_figure(cv_plot, "cv_performance")

primary_curves <- full$coefficients[
  full$coefficients$context == "full" &
    full$coefficients$method %in% c("FMCP", "fSCAD", "Smooth"),
  , drop = FALSE
]
primary_curves$method <- factor(
  primary_curves$method, levels = c("FMCP", "fSCAD", "Smooth")
)
coefficient_plot <- ggplot2::ggplot(
  primary_curves,
  ggplot2::aes(argval_physical, beta_physical, color = method)
) +
  ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~dataset, scales = "free", ncol = 2) +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::labs(
    title = "Full-data coefficient functions",
    subtitle = "Coefficients are back-transformed to response per physical argument unit",
    x = "Physical argument", y = expression(hat(beta)(t))
  ) + theme_real()
save_figure(coefficient_plot, "primary_coefficients")

stability <- analysis$stability_probabilities
stability$method <- factor(stability$method, levels = c("FMCP", "fSCAD"))
stability_plot <- ggplot2::ggplot(
  stability,
  ggplot2::aes(argval_physical, null_selection_probability, color = method)
) +
  ggplot2::geom_hline(
    yintercept = config$primary_stability_threshold,
    linetype = 2, color = "grey45"
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~dataset, scales = "free_x") +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(
    title = "Bootstrap stability of estimated null regions",
    subtitle = paste0(
      config$stability_replicates,
      " subject bootstraps; dashed line marks the primary 0.80 threshold"
    ),
    x = "Physical argument", y = "Null-selection probability"
  ) + theme_real()
save_figure(stability_plot, "stability_null_probability")

paired_primary <- analysis$paired[
  analysis$paired$comparator == "fSCAD", , drop = FALSE
]
paired_primary$metric <- factor(
  paired_primary$metric, levels = c("RMSE", "MAE", "R2")
)
paired_plot <- ggplot2::ggplot(
  paired_primary,
  ggplot2::aes(estimate, dataset, color = dataset)
) +
  ggplot2::geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(xmin = lower, xmax = upper),
    orientation = "y", width = 0.15, linewidth = 0.7
  ) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~metric, scales = "free_x", nrow = 1) +
  ggplot2::scale_x_continuous(
    expand = ggplot2::expansion(mult = c(0.10, 0.10))
  ) +
  ggplot2::scale_color_manual(values = c(Truck = "#336B87", Tecator = "#D55E00")) +
  ggplot2::labs(
    title = "Paired FMCP minus fSCAD contrasts",
    subtitle = paste0(
      "Percentile intervals from ", config$subject_bootstrap_replicates,
      " subject-level paired bootstraps"
    ),
    x = "Difference (negative favors FMCP for RMSE/MAE; positive for R2)",
    y = NULL
  ) + theme_real() + ggplot2::theme(
    legend.position = "none",
    panel.spacing.x = grid::unit(2.5, "lines"),
    axis.text.x = ggplot2::element_text(size = 8)
  )
save_figure(paired_plot, "paired_fmcp_fscad", height = 4.6)

permutation_primary <- analysis$permutation_summary[
  analysis$permutation_summary$primary_statistic, , drop = FALSE
]
permutation_long <- bind_rows_base(list(
  data.frame(
    permutation_primary[c("dataset", "method")],
    quantity = "Observed BIC gain", value = permutation_primary$observed
  ),
  data.frame(
    permutation_primary[c("dataset", "method")],
    quantity = "Permutation 95th percentile", value = permutation_primary$null_q95
  )
))
permutation_long$quantity <- factor(
  permutation_long$quantity,
  levels = c("Permutation 95th percentile", "Observed BIC gain")
)
permutation_plot <- ggplot2::ggplot(
  permutation_long,
  ggplot2::aes(quantity, value, color = method, group = method)
) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~dataset, scales = "free_y") +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::labs(
    title = "Permutation evidence for functional association",
    subtitle = paste0(
      "Observed intercept-to-model BIC gain versus the 95th percentile from ",
      config$permutation_replicates, " response permutations"
    ),
    x = NULL, y = "BIC gain over intercept-only model"
  ) + theme_real() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 15, hjust = 1))
save_figure(permutation_plot, "permutation_evidence")

primary_gamma_curve <- full$coefficients[
  full$coefficients$context == "full" &
    full$coefficients$method == "FMCP", , drop = FALSE
]
primary_gamma_curve$shape <- 3
primary_gamma_curve$context <- "gamma_3_primary"
gamma_curves <- bind_rows_base(list(full$gamma_curves, primary_gamma_curve))
gamma_curves$gamma <- ifelse(
  gamma_curves$shape == 3, "gamma = 3 (primary)",
  paste0("gamma = ", gamma_curves$shape)
)
gamma_plot <- ggplot2::ggplot(
  gamma_curves,
  ggplot2::aes(argval_physical, beta_physical, color = gamma)
) +
  ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::facet_wrap(~dataset, scales = "free") +
  ggplot2::scale_color_manual(values = c("gamma = 2.7" = "#0072B2",
                                        "gamma = 3 (primary)" = "#009E73",
                                        "gamma = 3.7" = "#D55E00")) +
  ggplot2::labs(
    title = "FMCP concavity sensitivity",
    subtitle = "Pre-specified primary and alternative MCP concavity values",
    x = "Physical argument", y = expression(hat(beta)(t))
  ) + theme_real()
save_figure(gamma_plot, "gamma_sensitivity")

format_number <- function(x, digits = 3L) {
  ifelse(is.na(x), "--", trimws(formatC(
    x, digits = digits, format = "g"
  )))
}

format_interval <- function(estimate, lower, upper, digits = 3L) {
  paste0(
    format_number(estimate, digits), " [",
    format_number(lower, digits), ", ",
    format_number(upper, digits), "]"
  )
}

write_latex_tabular <- function(data, path, align = NULL) {
  data[] <- lapply(data, function(column) latex_escape(as.character(column)))
  if (is.null(align)) align <- paste0("l", paste(rep("r", ncol(data) - 1L), collapse = ""))
  header <- paste(latex_escape(names(data)), collapse = " & ")
  body <- apply(data, 1L, function(row) paste0(paste(row, collapse = " & "), " \\\\"))
  lines <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule", paste0(header, " \\\\"), "\\midrule",
    body, "\\bottomrule", "\\end{tabular}"
  )
  atomic_write_lines(lines, path)
}

data_table <- data.frame(
  Dataset = manifest$dataset, Subjects = manifest$n,
  Grid = manifest$grid_points, Response = manifest$response,
  `Physical domain` = paste0(
    format_number(manifest$physical_grid_min, 5), "--",
    format_number(manifest$physical_grid_max, 5)
  ), check.names = FALSE, stringsAsFactors = FALSE
)
write_latex_tabular(data_table, table_paths[["table_data"]], "lrrll")

performance_table <- data.frame(
  Dataset = performance$dataset,
  Method = as.character(performance$method),
  `RMSE [95% CI]` = format_interval(
    performance$RMSE, performance$RMSE_lower, performance$RMSE_upper
  ),
  `MAE [95% CI]` = format_interval(
    performance$MAE, performance$MAE_lower, performance$MAE_upper
  ),
  `R2 [95% CI]` = format_interval(
    performance$R2, performance$R2_lower, performance$R2_upper
  ), check.names = FALSE, stringsAsFactors = FALSE
)
write_latex_tabular(performance_table, table_paths[["table_cv_performance"]], "llrrr")

paired_table <- data.frame(
  Dataset = paired_primary$dataset, Metric = paired_primary$metric,
  `FMCP - fSCAD [95% CI]` = format_interval(
    paired_primary$estimate, paired_primary$lower, paired_primary$upper
  ),
  `Pr(FMCP better)` = format_number(
    paired_primary$probability_fmcp_better, 3
  ),
  `Pr(numerical tie)` = format_number(
    paired_primary$probability_numerical_tie, 3
  ),
  `Bootstrap p` = format_number(paired_primary$bootstrap_two_sided_p, 3),
  check.names = FALSE, stringsAsFactors = FALSE
)
write_latex_tabular(paired_table, table_paths[["table_paired_fscad"]], "llrrrr")

full_sparse <- analysis$full_summary[
  analysis$full_summary$method %in% c("FMCP", "fSCAD"), , drop = FALSE
]
full_sparse_table <- data.frame(
  Dataset = full_sparse$dataset, Method = full_sparse$method,
  M = full_sparse$M, Lambda = format_number(full_sparse$lambda, 3),
  Kappa = format_number(full_sparse$kappa, 3),
  BIC = format_number(full_sparse$bic, 4),
  `Null grid fraction` = format_number(
    full_sparse$estimated_null_grid_fraction, 3
  ), check.names = FALSE, stringsAsFactors = FALSE
)
write_latex_tabular(full_sparse_table, table_paths[["table_full_sparse"]], "llrrrrr")

stable_primary <- analysis$stable_null_segments[
  abs(analysis$stable_null_segments$threshold -
        config$primary_stability_threshold) < 1e-12, , drop = FALSE
]
if (nrow(stable_primary)) {
  segment_table <- data.frame(
    Dataset = stable_primary$dataset, Method = stable_primary$method,
    Segment = stable_primary$segment,
    Start = format_number(stable_primary$physical_start, 5),
    End = format_number(stable_primary$physical_end, 5),
    `Grid points` = stable_primary$grid_points,
    check.names = FALSE, stringsAsFactors = FALSE
  )
} else {
  segment_table <- data.frame(
    Dataset = "All", Method = "FMCP/fSCAD", Segment = "--",
    Start = "--", End = "--", `Grid points` =
      "No segment reached 0.80 stability",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}
write_latex_tabular(segment_table, table_paths[["table_stable_segments"]], "llrrrr")

permutation_table <- data.frame(
  Dataset = permutation_primary$dataset, Method = permutation_primary$method,
  `Observed BIC gain` = format_number(permutation_primary$observed, 4),
  `Null 95%` = format_number(permutation_primary$null_q95, 4),
  `Permutation p` = format_number(permutation_primary$p_value_upper, 4),
  `Successful / planned` = paste0(
    permutation_primary$successful_permutations, "/",
    permutation_primary$planned_permutations
  ), check.names = FALSE, stringsAsFactors = FALSE
)
write_latex_tabular(permutation_table, table_paths[["table_permutation"]], "llrrrr")

gamma_table <- data.frame(
  Dataset = analysis$gamma_summary$dataset,
  Gamma = format_number(analysis$gamma_summary$gamma, 3),
  Role = ifelse(analysis$gamma_summary$gamma == 3, "Primary", "Sensitivity"),
  M = analysis$gamma_summary$M,
  Lambda = format_number(analysis$gamma_summary$lambda, 3),
  Kappa = format_number(analysis$gamma_summary$kappa, 3),
  BIC = format_number(analysis$gamma_summary$bic, 4),
  `Null grid fraction` = format_number(
    analysis$gamma_summary$estimated_null_grid_fraction, 3
  ), check.names = FALSE, stringsAsFactors = FALSE
)
write_latex_tabular(gamma_table, table_paths[["table_gamma"]], "lrlrrrrr")

validation_table <- data.frame(
  Check = analysis$validation$check,
  Expected = analysis$validation$expected,
  Observed = analysis$validation$observed,
  Passed = ifelse(analysis$validation$passed, "Yes", "No"),
  stringsAsFactors = FALSE
)
write_latex_tabular(validation_table, table_paths[["table_validation"]], "lrrl")

best_rows <- do.call(rbind, lapply(split(performance, performance$dataset), function(data) {
  data[which.min(data$RMSE), , drop = FALSE]
}))
best_text <- paste(vapply(seq_len(nrow(best_rows)), function(index) {
  paste0(
    best_rows$dataset[index], ": ", as.character(best_rows$method[index]),
    " (RMSE ", format_number(best_rows$RMSE[index], 4), ")"
  )
}, character(1L)), collapse = "; ")

primary_text <- paste(vapply(seq_len(nrow(paired_primary)), function(index) {
  paste0(
    paired_primary$dataset[index], " ", paired_primary$metric[index], " ",
    format_interval(
      paired_primary$estimate[index], paired_primary$lower[index],
      paired_primary$upper[index], 3
    )
  )
}, character(1L)), collapse = "; ")

figure_file <- function(stem) versioned_name(stem, mode, "pdf")
table_file <- function(stem) versioned_name(stem, mode, "tex")
report_title <- if (identical(mode, "smoke")) {
  "Real-data application of functional MCP regression (smoke validation)"
} else {
  "Real-data application of functional MCP regression: prediction and stable local sparsity"
}

tex_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[utf8]{inputenc}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage{lmodern}",
  "\\usepackage{booktabs}",
  "\\usepackage{graphicx}",
  "\\usepackage{geometry}",
  "\\usepackage{float}",
  "\\usepackage{xcolor}",
  "\\usepackage[hidelinks]{hyperref}",
  "\\geometry{margin=1in}",
  "\\graphicspath{{../../../real_data/figures/v1/}}",
  paste0("\\title{", latex_escape(report_title), "}"),
  "\\author{FMCP real-data analysis pipeline}",
  "\\date{August 2026}",
  "\\begin{document}",
  "\\maketitle",
  "\\begin{abstract}",
  paste(
    "We evaluate functional minimax concave penalization (FMCP) on two",
    "scalar-on-function regression applications: particulate-matter emissions",
    "from the Truck data and fat content from Tecator near-infrared spectra.",
    "The design separates predictive validation from local-null interpretation,",
    "uses identical subject splits for every method, and compares FMCP directly",
    "with functional SCAD (fSCAD) and six additional benchmarks."
  ),
  "\\end{abstract}",
  "\\section{Data}",
  paste(
    "Truck contains 108 observational rows and 90 second-indexed acceleration",
    "measurements. The distributed data do not provide truck or test-cycle IDs;",
    "rows are therefore the resampling units, which is a limitation if hidden",
    "within-truck dependence exists. Tecator contains 215 meat samples, 100",
    "near-infrared absorbance channels over the documented 850--1050 nm range,",
    "and fat percentage as the scalar response. In accordance with the Tecator",
    "permission note, the instrument and company name (Tecator) are stated here."
  ),
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Real-data analysis samples. No missing-value imputation, outlier removal, or observation exclusion was applied.}",
  paste0("\\input{", table_file("table_data"), "}"),
  "\\end{table}",
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=\\textwidth]{", figure_file("real_data_profiles"), "}"),
  "\\caption{Representative functional predictors and their sample means.}",
  "\\end{figure}",
  "\\section{Methods}",
  paste(
    "All observed grids were linearly mapped to [0,1] for estimation while",
    "physical labels were retained for interpretation. Within each training",
    "set only, the response was divided by its standard deviation and the",
    "functional predictor by one global integrated RMS scale computed after",
    "curve centering. Test observations used the corresponding training scales.",
    "All functional integrals used trapezoid weights. Coefficient functions were",
    "back-transformed by the response scale, predictor scale, and physical-domain",
    "length. No pointwise standardization was used."
  ),
  paste(
    "FMCP and controlled fSCAD share the same cubic B-spline design, second-derivative",
    "roughness penalty, numerical integration, BIC rule, convergence controls,",
    "and refitting rule; only the nonconvex penalty shape differs. The primary MCP",
    "concavity is 3 and the SCAD shape is 3.7. Smooth, PCR, PenS, NGR, FLiRTI,",
    "and OLS form additional benchmarks. PFLR-style comparators were implemented",
    "as isolated trapezoid adaptations because both real-data grids have an even",
    "number of channels and the distributed PFLR routines assume Simpson grids."
  ),
  paste(
    "The initial basis grid was M in {10,20,30,40,50}. A single pre-specified",
    "boundary audit could expand a boundary when its selection rate exceeded 10%;",
    "no second expansion was permitted. Predictive performance used",
    paste0(config$cv_repeats, " repeated ", config$cv_folds, "-fold CV."),
    paste0(config$subject_bootstrap_replicates, " paired subject bootstraps"),
    "formed percentile intervals. Local-null stability used",
    paste0(config$stability_replicates, " subject bootstraps with 0.80 as the"),
    "primary probability threshold (0.70 and 0.90 sensitivity). Functional",
    "association was tested with",
    paste0(config$permutation_replicates, " response permutations; each permutation"),
    "repeated tuning, and the primary statistic was the BIC gain over the",
    "intercept-only model. Differences no larger than the square root of machine",
    "precision on the metric scale were reported as numerical ties."
  ),
  "\\section{Prediction results}",
  paste0("The lowest pooled CV RMSE by data set was: ", latex_escape(best_text), "."),
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=\\textwidth]{", figure_file("cv_performance"), "}"),
  "\\caption{Repeated-CV RMSE with paired subject-bootstrap intervals.}",
  "\\end{figure}",
  "\\begin{table}[H]\\centering\\scriptsize",
  "\\caption{Repeated-CV performance. Intervals resample subjects and retain all repeat-specific predictions for a sampled subject.}",
  "\\resizebox{\\textwidth}{!}{%",
  paste0("\\input{", table_file("table_cv_performance"), "}}"),
  "\\end{table}",
  "\\clearpage",
  paste0(
    "The pre-specified paired FMCP-minus-fSCAD contrasts were ",
    latex_escape(primary_text),
    ". For RMSE and MAE, negative differences favor FMCP; for R2, positive differences favor FMCP."
  ),
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=\\textwidth]{", figure_file("paired_fmcp_fscad"), "}"),
  "\\caption{Paired FMCP--fSCAD contrasts.}",
  "\\end{figure}",
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Primary paired FMCP--fSCAD comparisons.}",
  paste0("\\input{", table_file("table_paired_fscad"), "}"),
  "\\end{table}",
  "\\section{Coefficient functions and local-null stability}",
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=\\textwidth]{", figure_file("primary_coefficients"), "}"),
  "\\caption{Back-transformed full-data coefficient functions.}",
  "\\end{figure}",
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Full-data sparse-model selections. A grid point is numerically null at tolerance $\\sqrt{\\epsilon_{machine}}$.}",
  paste0("\\input{", table_file("table_full_sparse"), "}"),
  "\\end{table}",
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=\\textwidth]{", figure_file("stability_null_probability"), "}"),
  "\\caption{Bootstrap probability that each coefficient-grid location is estimated null.}",
  "\\end{figure}",
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Contiguous stable-null segments at the primary 0.80 threshold. Endpoints are observed-grid endpoints, not interpolated change points.}",
  paste0("\\input{", table_file("table_stable_segments"), "}"),
  "\\end{table}",
  "\\section{Association and sensitivity analyses}",
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=0.92\\textwidth]{", figure_file("permutation_evidence"), "}"),
  "\\caption{Observed and permutation-reference BIC gains.}",
  "\\end{figure}",
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Permutation tests for functional association. Monte Carlo p-values use $(1+r)/(1+B)$.}",
  paste0("\\input{", table_file("table_permutation"), "}"),
  "\\end{table}",
  "\\begin{figure}[H]\\centering",
  paste0("\\includegraphics[width=\\textwidth]{", figure_file("gamma_sensitivity"), "}"),
  "\\caption{FMCP coefficient sensitivity to the primary and alternative MCP concavity values.}",
  "\\end{figure}",
  "\\begin{table}[H]\\centering\\small",
  "\\caption{FMCP concavity sensitivity.}",
  paste0("\\input{", table_file("table_gamma"), "}"),
  "\\end{table}",
  "\\section{Interpretation and limitations}",
  paste(
    "Prediction and local sparsity answer different questions. A smooth mean",
    "coefficient can predict well without identifying reproducible null regions;",
    "FMCP and fSCAD explicitly target local zeros. Accordingly, superiority is",
    "not inferred from a single full-data curve or training BIC. The repeated-CV",
    "contrasts quantify prediction, while bootstrap null-selection probabilities",
    "quantify stability of the local interpretation."
  ),
  paste(
    "The Truck analysis may understate uncertainty if rows share unrecorded trucks",
    "or cycles. Tecator wavelengths are nominally mapped across the documented",
    "range because channel-specific labels are absent. The real-data truth is",
    "unknown, so stable-null regions are reproducible estimates rather than",
    "verified biological or mechanical zeros. Boundary selections remaining after",
    "the single allowed expansion are reported rather than triggering data-driven",
    "additional searches."
  ),
  "\\section{Reproducibility checks}",
  "\\begin{table}[H]\\centering\\small",
  "\\caption{Pipeline validation checks.}",
  paste0("\\input{", table_file("table_validation"), "}"),
  "\\end{table}",
  "\\begin{thebibliography}{9}",
  "\\bibitem{lin2017} Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally sparse estimator for functional linear regression models. \\emph{Journal of Computational and Graphical Statistics}, 26(2), 306--318.",
  "\\bibitem{zhang2010} Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave penalty. \\emph{The Annals of Statistics}, 38(2), 894--942.",
  "\\bibitem{guan2020} Guan, T., Lin, Z., and Cao, J. (2020). Estimating truncated functional linear models with a nested group bridge approach. \\emph{Journal of Computational and Graphical Statistics}, 29(3), 620--628. doi:10.1080/10618600.2020.1713797.",
  "\\bibitem{james2009} James, G., Wang, J., and Zhu, J. (2009). Functional linear regression that's interpretable. \\emph{The Annals of Statistics}, 37(5A), 2083--2108.",
  "\\bibitem{clark2007} Clark, N. et al. (2007). Heavy-duty vehicle chassis dynamometer testing for emissions inventory, air quality modeling, source apportionment and air toxics emissions inventory: E55/59 all phases. Coordinating Research Council.",
  "\\end{thebibliography}",
  "\\end{document}"
)
atomic_write_lines(tex_lines, tex_path)

readme_lines <- c(
  paste0("# ", report_title), "",
  paste0("- LaTeX source: `", basename(tex_path), "`"),
  paste0("- Compiled PDF: `", basename(pdf_path), "`"),
  paste0("- Mode: `", mode, "`"),
  paste0("- Figure source directory: `", config$figures_dir, "`"),
  "",
  "All table cells are generated from versioned RDS/CSV analysis outputs."
)
atomic_write_lines(readme_lines, report_readme)

old_directory <- getwd()
on.exit(setwd(old_directory), add = TRUE)
setwd(config$outputs_dir)
latex_command <- Sys.which("pdflatex")
if (!nzchar(latex_command)) {
  stop("pdflatex is required to compile the manuscript.", call. = FALSE)
}
for (pass in 1:2) {
  output <- system2(
    latex_command,
    c("-interaction=nonstopmode", "-halt-on-error", basename(tex_path)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      "pdflatex failed on pass ", pass, ":\n",
      paste(tail(output, 30L), collapse = "\n"), call. = FALSE
    )
  }
}
if (!file.exists(pdf_path)) {
  stop("Expected compiled PDF was not created.", call. = FALSE)
}
cat("Figures, LaTeX tables, and compiled manuscript created at:\n", pdf_path, "\n")
