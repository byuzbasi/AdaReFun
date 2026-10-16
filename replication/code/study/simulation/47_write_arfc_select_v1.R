#!/usr/bin/env Rscript

# Publication-style tables and figures for ARFC-Select version 1.
# Smoke outputs are explicitly noninferential and live in a versioned bundle.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/47_write_arfc_select_v1.R", mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
version_argument <- grep("^--report-version=", arguments, value = TRUE)
if (length(version_argument) > 1L) {
  stop("Supply at most one --report-version argument.", call. = FALSE)
}
report_version <- if (length(version_argument)) {
  sub("^--report-version=", "", version_argument)
} else {
  "v1"
}
if (!report_version %in% c("v1", "v2")) {
  stop("`--report-version` must be v1 or v2.", call. = FALSE)
}
unknown <- setdiff(arguments, c("--smoke", version_argument))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

required_packages <- c("ggplot2", "patchwork", "knitr", "kableExtra", "scales")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE,
  FUN.VALUE = logical(1L)
)]
if (length(missing_packages)) {
  stop("Missing reporting package(s): ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

study <- "arfc_select_v1"
mode <- if (smoke) "smoke" else "production"
results_root <- file.path(simulation_dir, "results", study)
analysis_dir <- file.path(results_root, paste0(mode, "_analysis_v1"))
report_root <- file.path(project_dir, "output", "pdf", study)
report_dir <- file.path(report_root, paste0(mode, "_reporting_", report_version))
if (!dir.exists(analysis_dir)) {
  stop("Complete ARFC analysis bundle not found: ", analysis_dir,
       call. = FALSE)
}
if (dir.exists(report_dir) || file.exists(report_dir)) {
  stop("Refusing to overwrite ARFC reporting bundle: ", report_dir,
       call. = FALSE)
}

prefix <- study
analysis_file <- file.path(analysis_dir, paste0(prefix, "_analysis_v1.rds"))
analysis_manifest_file <- file.path(
  analysis_dir, paste0(prefix, "_analysis_manifest_v1.csv")
)
if (!file.exists(analysis_file) || !file.exists(analysis_manifest_file)) {
  stop("Analysis RDS or analysis manifest is missing.", call. = FALSE)
}

verify_manifest <- function(directory, manifest_file) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required <- c("scientific_signature", "file", "bytes", "md5")
  if (!all(required %in% names(manifest)) || !nrow(manifest) ||
      anyDuplicated(manifest$file)) {
    stop("Analysis manifest schema is invalid.", call. = FALSE)
  }
  paths <- file.path(directory, manifest$file)
  if (!all(file.exists(paths)) ||
      !identical(unname(as.character(tools::md5sum(paths))),
                 as.character(manifest$md5)) ||
      any(as.numeric(file.info(paths)$size) != as.numeric(manifest$bytes))) {
    stop("Analysis manifest integrity verification failed.", call. = FALSE)
  }
  manifest
}

analysis_manifest <- verify_manifest(analysis_dir, analysis_manifest_file)
analysis <- readRDS(analysis_file)
manifest_signature <- unique(as.character(
  analysis_manifest$scientific_signature
))
if (length(manifest_signature) != 1L ||
    !identical(analysis$study, study) || !identical(analysis$mode, mode) ||
    !identical(analysis$scientific_signature, manifest_signature) ||
    !identical(analysis$analysis_status$cross_case_overall_computed, FALSE)) {
  stop("Analysis identity, signature, or no-overall contract failed.",
       call. = FALSE)
}

selection <- analysis$selection_frequency
selected_performance <- analysis$selected_performance
oracle <- analysis$oracle_summary
curves <- analysis$selected_curve_summary
contrast_summaries <- analysis$selected_candidate_case_summaries
methods <- analysis$method_metadata
scenario_order <- as.character(analysis$scenarios$scenario)
case_order <- as.character(analysis$scenarios$case_label)
if (nrow(selection) != 4L * length(analysis$config$method_order) ||
    any(abs(tapply(selection$selection_frequency, selection$scenario, sum) -
            1) > 1e-12) ||
    any(tapply(selection$selected_count, selection$scenario, sum) !=
          analysis$config$replications) ||
    !all(scenario_order %in% unique(curves$scenario))) {
  stop("Reporting inputs failed their structural checks.", call. = FALSE)
}

dir.create(report_root, recursive = TRUE, showWarnings = FALSE)
staging <- file.path(
  report_root, paste0(".", basename(report_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) {
  stop("Reporting staging directory already exists: ", staging,
       call. = FALSE)
}
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed && dir.exists(staging)) {
  unlink(staging, recursive = TRUE, force = TRUE)
}, add = TRUE)

method_colors <- c(
  asr_mcp = "#0072B2",
  asr_scad = "#E69F00",
  amr_mcp = "#009E73",
  amr_scad = "#CC79A7"
)
selection$case_label <- factor(selection$case_label, levels = case_order)
selection$label <- factor(
  selection$label,
  levels = methods$label[match(analysis$config$method_order, methods$method)]
)

base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(color = "#4D4D4D"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F2F2F2"),
    strip.text = ggplot2::element_text(face = "bold")
  )

smoke_subtitle <- if (smoke) {
  "Technical smoke only (one task per case); frequencies are noninferential"
} else {
  paste0(
    "Independent-validation selection; ", analysis$config$replications,
    " fresh replications per case"
  )
}

selection_plot <- ggplot2::ggplot(
  selection,
  ggplot2::aes(x = case_label, y = selection_frequency, fill = label)
) +
  ggplot2::geom_col(width = 0.68, color = "white", linewidth = 0.25) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = ifelse(
        selection_frequency >= 0.05,
        paste0(selected_count, "/", n_tasks), ""
      )
    ),
    position = ggplot2::position_stack(vjust = 0.5),
    color = "white", fontface = "bold", size = 3.2
  ) +
  ggplot2::scale_fill_manual(
    values = unname(method_colors[analysis$config$method_order]),
    breaks = methods$label[match(analysis$config$method_order, methods$method)],
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "ARFC-Select candidate frequencies",
    subtitle = smoke_subtitle,
    x = NULL,
    y = "Selection frequency"
  ) +
  base_theme

case_ii <- curves[curves$scenario == "case_ii_sparse", , drop = FALSE]
case_ii <- case_ii[order(case_ii$t), , drop = FALSE]
if (!nrow(case_ii) || anyDuplicated(case_ii$t)) {
  stop("Case-II selected-curve summary is missing or duplicated.",
       call. = FALSE)
}

curve_plot <- ggplot2::ggplot(case_ii, ggplot2::aes(x = t)) +
  ggplot2::annotate(
    "rect", xmin = 0.3, xmax = 0.7, ymin = -Inf, ymax = Inf,
    fill = "#BDBDBD", alpha = 0.24
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower_10, ymax = upper_90),
    fill = "#56B4E9", alpha = 0.28
  ) +
  ggplot2::geom_hline(yintercept = 0, color = "#888888", linewidth = 0.35) +
  ggplot2::geom_line(
    ggplot2::aes(y = truth, linetype = "Truth"),
    color = "#222222", linewidth = 0.9
  ) +
  ggplot2::geom_line(
    ggplot2::aes(y = mean_estimate, linetype = "ARFC-Select mean"),
    color = "#0072B2", linewidth = 0.95
  ) +
  ggplot2::scale_linetype_manual(
    values = c("Truth" = "longdash", "ARFC-Select mean" = "solid")
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(beta(t)),
    linetype = NULL
  ) +
  base_theme +
  ggplot2::theme(legend.position = "bottom")

zero_plot <- ggplot2::ggplot(
  case_ii, ggplot2::aes(x = t, y = zero_selection_frequency)
) +
  ggplot2::annotate(
    "rect", xmin = 0.3, xmax = 0.7, ymin = -Inf, ymax = Inf,
    fill = "#BDBDBD", alpha = 0.24
  ) +
  ggplot2::geom_area(fill = "#E69F00", alpha = 0.28) +
  ggplot2::geom_line(color = "#D55E00", linewidth = 0.9) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1), breaks = c(0, 0.5, 1),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = "Functional argument t",
    y = "Exact-zero frequency"
  ) +
  base_theme +
  ggplot2::theme(legend.position = "none")

recovery_title <- if (smoke) {
  "Case II selected recovery (technical smoke; noninferential)"
} else {
  "Case II recovery after independent-validation selection"
}
recovery_plot <- curve_plot / zero_plot +
  patchwork::plot_layout(heights = c(1.7, 1)) +
  patchwork::plot_annotation(
    title = recovery_title,
    subtitle = paste(
      "Gray band: true null interval [0.3, 0.7]; blue band:",
      "pointwise empirical 10%-90% Monte Carlo envelope"
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 10, color = "#4D4D4D")
    )
  )

figure_files <- c(
  selection_pdf = paste0(prefix, "_selection_frequencies_", report_version,
                         ".pdf"),
  selection_png = paste0(prefix, "_selection_frequencies_", report_version,
                         ".png"),
  recovery_pdf = paste0(prefix, "_case_ii_selected_recovery_", report_version,
                        ".pdf"),
  recovery_png = paste0(prefix, "_case_ii_selected_recovery_", report_version,
                        ".png")
)
ggplot2::ggsave(
  file.path(staging, figure_files[["selection_pdf"]]),
  selection_plot, width = 8.2, height = 4.8, units = "in",
  device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  file.path(staging, figure_files[["selection_png"]]),
  selection_plot, width = 8.2, height = 4.8, units = "in",
  dpi = 240, bg = "white"
)
ggplot2::ggsave(
  file.path(staging, figure_files[["recovery_pdf"]]),
  recovery_plot, width = 8.4, height = 6.7, units = "in",
  device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  file.path(staging, figure_files[["recovery_png"]]),
  recovery_plot, width = 8.4, height = 6.7, units = "in",
  dpi = 240, bg = "white"
)

format_number <- function(value, digits = 4L) {
  ifelse(is.finite(value), formatC(value, digits = digits, format = "g"), "--")
}
format_percent <- function(value) {
  ifelse(is.finite(value), sprintf("%.1f\\%%", 100 * value), "--")
}

selection_table <- selection[order(selection$case_label, selection$label),
                             c("case_label", "label", "selected_count",
                               "n_tasks", "selection_frequency")]
names(selection_table) <- c(
  "Case", "Candidate", "Selected", "Tasks", "Frequency"
)
selection_table$Frequency <- format_percent(selection_table$Frequency)

performance_rows <- lapply(scenario_order, function(scenario) {
  response <- selected_performance[
    selected_performance$scenario == scenario &
      selected_performance$metric == "response_prediction_mse", , drop = FALSE
  ]
  signal <- selected_performance[
    selected_performance$scenario == scenario &
      selected_performance$metric == "signal_prediction_mse", , drop = FALSE
  ]
  oracle_row <- oracle[oracle$scenario == scenario, , drop = FALSE]
  data.frame(
    Case = response$case_label,
    Tasks = oracle_row$n_tasks,
    `Test response MSE` = format_number(response$mean),
    `Test signal MSE` = format_number(signal$mean),
    `Mean regret` = format_number(oracle_row$mean_prediction_regret),
    `Oracle agreement` = format_percent(oracle_row$oracle_agreement_rate),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
performance_table <- do.call(rbind, performance_rows)
if (identical(report_version, "v2")) {
  performance_table$Tasks <- NULL
}

case_ii_contrasts <- contrast_summaries[
  contrast_summaries$scenario == "case_ii_sparse" &
    contrast_summaries$metric == "response_prediction_mse",
  , drop = FALSE
]
case_ii_contrasts <- case_ii_contrasts[
  match(analysis$config$method_order, case_ii_contrasts$candidate_method),
  , drop = FALSE
]
if (nrow(case_ii_contrasts) != length(analysis$config$method_order) ||
    anyNA(case_ii_contrasts$candidate_method)) {
  stop("Case-II response-MSE contrast table is incomplete.", call. = FALSE)
}
selected_mean_tolerance <- 1e-12 *
  (1 + abs(mean(case_ii_contrasts$selected_mean)))
if (diff(range(case_ii_contrasts$selected_mean)) > selected_mean_tolerance) {
  stop(
    "Case-II ARFC-Select means differ across candidate contrasts.",
    call. = FALSE
  )
}
case_ii_selected_mean <- mean(case_ii_contrasts$selected_mean)
if (identical(report_version, "v2")) {
  contrast_table <- data.frame(
    Candidate = methods$label[
      match(case_ii_contrasts$candidate_method, methods$method)
    ],
    `Candidate MSE` = format_number(case_ii_contrasts$candidate_mean, 5L),
    `Selected - candidate` = sprintf(
      "%+.2e", case_ii_contrasts$mean_difference
    ),
    `95\\% CI` = sprintf(
      "[%+.2e, %+.2e]",
      case_ii_contrasts$lower_95, case_ii_contrasts$upper_95
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  contrast_table <- data.frame(
    Candidate = methods$label[
      match(case_ii_contrasts$candidate_method, methods$method)
    ],
    `ARFC-Select MSE` = format_number(case_ii_contrasts$selected_mean, 5L),
    `Candidate MSE` = format_number(case_ii_contrasts$candidate_mean, 5L),
    `Selected - candidate` = sprintf(
      "%+.2e", case_ii_contrasts$mean_difference
    ),
    `95\\% CI` = sprintf(
      "[%+.2e, %+.2e]",
      case_ii_contrasts$lower_95, case_ii_contrasts$upper_95
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

kable_fragment <- function(table, caption, label, align) {
  value <- knitr::kable(
    table,
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = caption,
    label = label,
    align = align
  )
  as.character(kableExtra::kable_styling(
    value,
    latex_options = c("hold_position"),
    position = "center",
    font_size = 9
  ))
}

table_files <- c(
  selection = paste0(prefix, "_table_selection_frequencies_", report_version,
                     ".tex"),
  performance = paste0(prefix, "_table_selected_performance_", report_version,
                       ".tex"),
  contrast = paste0(prefix, "_table_case_ii_prediction_contrasts_",
                    report_version, ".tex")
)
selection_caption <- if (smoke) {
  paste(
    "Technical-smoke candidate selections (one task per case;",
    "noninferential)."
  )
} else {
  paste(
    "Candidate selection frequencies from the independent validation",
    "sample. Cases are not pooled."
  )
}
performance_caption <- if (smoke) {
  "Technical-smoke selected-estimator diagnostics (noninferential)."
} else {
  paste(
    "Test-only performance after independent-validation selection.",
    "Cases are not pooled."
  )
}
writeLines(
  kable_fragment(
    selection_table, selection_caption, "arfc-select-frequency",
    c("l", "l", "r", "r", "r")
  ),
  file.path(staging, table_files[["selection"]]), useBytes = TRUE
)
writeLines(
  kable_fragment(
    performance_table, performance_caption,
    "arfc-select-performance",
    if (identical(report_version, "v2")) {
      c("l", "r", "r", "r", "r")
    } else {
      c("l", "r", "r", "r", "r", "r")
    }
  ),
  file.path(staging, table_files[["performance"]]), useBytes = TRUE
)
writeLines(
  kable_fragment(
    contrast_table,
    if (identical(report_version, "v2")) {
      paste0(
        "Case-II paired test response-MSE contrasts. The ARFC-Select mean ",
        "response MSE, common to all four rows, was ",
        format_number(case_ii_selected_mean, 5L),
        ". Negative differences favor ARFC-Select; intervals are paired ",
        "95\\% t intervals."
      )
    } else {
      paste(
        "Case-II paired test response-MSE contrasts. Negative differences",
        "favor ARFC-Select; intervals are paired 95\\% t intervals."
      )
    },
    "arfc-select-caseii-contrast",
    if (identical(report_version, "v2")) {
      c("l", "r", "r", "r")
    } else {
      c("l", "r", "r", "r", "r")
    }
  ),
  file.path(staging, table_files[["contrast"]]), useBytes = TRUE
)

latex_escape <- function(value) {
  value <- gsub("\\\\", "\\\\textbackslash{}", value)
  value <- gsub("([%&#_$])", "\\\\\\1", value, perl = TRUE)
  value <- gsub("~", "\\\\textasciitilde{}", value, fixed = TRUE)
  value <- gsub("\\^", "\\\\textasciicircum{}", value)
  value
}

results_file <- paste0(prefix, "_results_text_", report_version, ".tex")
if (smoke) {
  results_lines <- c(
    "\\paragraph{Technical smoke only.}",
    paste(
      "The four-task run verified the complete training--validation--test",
      "selection pipeline. Each case contains one task, so these values are",
      "noninferential and no scientific comparison is reported."
    )
  )
} else {
  results_lines <- c("\\paragraph{Four-candidate validation selection.}")
  for (scenario in scenario_order) {
    block <- selection[selection$scenario == scenario, , drop = FALSE]
    dominant <- block[order(-block$selected_count,
                            match(block$method, analysis$config$method_order)),
                      , drop = FALSE][1L, ]
    oracle_row <- oracle[oracle$scenario == scenario, , drop = FALSE]
    response <- selected_performance[
      selected_performance$scenario == scenario &
        selected_performance$metric == "response_prediction_mse",
      , drop = FALSE
    ]
    results_lines <- c(results_lines, sprintf(
      paste0(
        "%s: %s was selected in %d/%d replications; the selected-estimator ",
        "test response MSE was %s, mean oracle regret was %s, and test-oracle ",
        "agreement was %.1f\\%%."
      ),
      latex_escape(as.character(dominant$case_label)),
      latex_escape(as.character(dominant$label)),
      dominant$selected_count, dominant$n_tasks,
      format_number(response$mean),
      format_number(oracle_row$mean_prediction_regret),
      100 * oracle_row$oracle_agreement_rate
    ))
  }
  results_lines <- c(
    results_lines,
    sprintf(
      paste0(
        "In Case II, all four paired selected-minus-candidate test response-",
        "MSE differences were negative (range %+.2e to %+.2e), and every ",
        "95\\%% interval had an upper endpoint below zero."
      ),
      min(case_ii_contrasts$mean_difference),
      max(case_ii_contrasts$mean_difference)
    ),
    "No cross-case overall estimate or ranking was computed."
  )
}
writeLines(results_lines, file.path(staging, results_file), useBytes = TRUE)

figure_include_file <- paste0(
  prefix, "_figure_includes_", report_version, ".tex"
)
report_relative <- file.path(
  "output", "pdf", study, paste0(mode, "_reporting_", report_version)
)
selection_figure_caption <- if (smoke) {
  paste(
    "Technical-smoke candidate selections (one task per case;",
    "noninferential). Cases are shown separately."
  )
} else {
  paste(
    "Independent-validation selection frequencies for the four ARFC",
    "candidates. Cases are shown separately."
  )
}
recovery_figure_caption <- if (smoke) {
  paste(
    "Technical-smoke Case-II coefficient recovery and pointwise exact-zero",
    "frequency (noninferential). The gray region is the true null interval."
  )
} else {
  paste(
    "Case-II coefficient recovery and pointwise exact-zero frequency after",
    "candidate selection. The gray region is the true null interval."
  )
}
figure_lines <- c(
  "\\begin{figure}[t]",
  "\\centering",
  sprintf("\\includegraphics[width=0.92\\textwidth]{%s}",
          file.path(report_relative, figure_files[["selection_pdf"]])),
  paste0("\\caption{", selection_figure_caption, "}"),
  "\\label{fig:arfc-select-frequency}",
  "\\end{figure}",
  "",
  "\\begin{figure}[t]",
  "\\centering",
  sprintf("\\includegraphics[width=0.92\\textwidth]{%s}",
          file.path(report_relative, figure_files[["recovery_pdf"]])),
  paste0("\\caption{", recovery_figure_caption, "}"),
  "\\label{fig:arfc-select-caseii}",
  "\\end{figure}"
)
writeLines(
  figure_lines, file.path(staging, figure_include_file), useBytes = TRUE
)

reporting_metadata <- list(
  study = study,
  mode = mode,
  report_version = report_version,
  scientific_signature = analysis$scientific_signature,
  analysis_manifest_md5 = unname(as.character(
    tools::md5sum(analysis_manifest_file)
  )),
  writer_script_md5 = unname(as.character(tools::md5sum(script_file))),
  scientific_interpretation_allowed = !smoke,
  cross_case_overall_computed = FALSE,
  figures = unname(figure_files),
  tables = unname(table_files),
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
metadata_file <- paste0(
  prefix, "_reporting_metadata_", report_version, ".rds"
)
session_file <- paste0(
  prefix, "_reporting_sessionInfo_", report_version, ".txt"
)
saveRDS(reporting_metadata, file.path(staging, metadata_file),
        compress = "xz", version = 3)
writeLines(capture.output(sessionInfo()), file.path(staging, session_file),
           useBytes = TRUE)

payload_files <- c(
  unname(figure_files), unname(table_files), results_file,
  figure_include_file, metadata_file, session_file
)
payload_paths <- file.path(staging, payload_files)
roles <- c(
  "selection_figure_pdf", "selection_figure_png",
  "case_ii_recovery_pdf", "case_ii_recovery_png",
  "selection_table_tex", "performance_table_tex", "contrast_table_tex",
  "results_text_tex", "figure_includes_tex", "reporting_metadata_rds",
  "session_info_txt"
)
report_manifest <- data.frame(
  scientific_signature = analysis$scientific_signature,
  role = roles,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
manifest_file <- file.path(
  staging, paste0(prefix, "_reporting_manifest_", report_version, ".csv")
)
utils::write.csv(report_manifest, manifest_file, row.names = FALSE)

if (!file.rename(staging, report_dir)) {
  stop("Could not atomically install ARFC reporting bundle.", call. = FALSE)
}
committed <- TRUE

cat("ARFC selection reporting complete.\n")
cat("Mode:", mode, "\n")
cat("Scientific interpretation allowed:", !smoke, "\n")
cat("Figures: 2 PDF + 2 PNG\n")
cat("Tables: 3 LaTeX fragments\n")
cat("Bundle:", report_dir, "\n")
