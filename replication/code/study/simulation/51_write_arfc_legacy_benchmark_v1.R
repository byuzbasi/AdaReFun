#!/usr/bin/env Rscript

# Versioned LaTeX tables and publication figures for the ARFC benchmark.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/51_write_arfc_legacy_benchmark_v1.R",
                mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
report_argument <- grep("^--report-version=", arguments, value = TRUE)
if (length(report_argument) > 1L) {
  stop("Specify --report-version at most once.", call. = FALSE)
}
report_version <- if (length(report_argument)) {
  sub("^--report-version=", "", report_argument)
} else {
  "v1"
}
if (!report_version %in% c("v1", "v2")) {
  stop("--report-version must be v1 or v2.", call. = FALSE)
}
unknown <- setdiff(arguments, c("--smoke", report_argument))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}
mode <- if (smoke) "smoke" else "production"
prefix <- "arfc_legacy_benchmark_v1"
results_root <- file.path(simulation_dir, "results", prefix)
analysis_dir <- file.path(results_root, paste0(mode, "_analysis_v1"))
analysis_file <- file.path(analysis_dir, paste0(prefix, "_analysis_v1.rds"))
analysis_manifest_file <- file.path(
  analysis_dir, paste0(prefix, "_analysis_manifest_v1.csv")
)
generated_parent <- file.path(
  project_dir, "config", "generated", prefix
)
reporting_dir <- file.path(
  generated_parent, paste0(mode, "_reporting_", report_version)
)
pdf_parent <- file.path(project_dir, "output", "pdf", prefix)
pdf_dir <- file.path(pdf_parent, paste0(mode, "_reporting_", report_version))

if (!file.exists(analysis_file) || !file.exists(analysis_manifest_file)) {
  stop("Validated analysis bundle is missing.", call. = FALSE)
}
if (dir.exists(reporting_dir) || dir.exists(pdf_dir)) {
  stop("Reporting output already exists; refusing overwrite.", call. = FALSE)
}

verify_manifest <- function(directory, manifest_file) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required <- c("scientific_signature", "file", "bytes", "md5")
  if (!all(required %in% names(manifest)) || anyDuplicated(manifest$file)) {
    stop("Analysis manifest schema is invalid.", call. = FALSE)
  }
  paths <- file.path(directory, manifest$file)
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == manifest$bytes) ||
      !all(unname(as.character(tools::md5sum(paths))) == manifest$md5)) {
    stop("Analysis manifest integrity failure.", call. = FALSE)
  }
  manifest
}

analysis_manifest <- verify_manifest(analysis_dir, analysis_manifest_file)
analysis <- readRDS(analysis_file)
if (!identical(analysis$study, prefix) || !identical(analysis$mode, mode) ||
    length(unique(analysis_manifest$scientific_signature)) != 1L ||
    !identical(unique(analysis_manifest$scientific_signature),
               analysis$scientific_signature) ||
    !all(analysis$validation$passed)) {
  stop("Analysis identity or validation gate failed.", call. = FALSE)
}

required_packages <- c("ggplot2", "knitr", "kableExtra", "png")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE,
  FUN.VALUE = logical(1L)
)]
if (length(missing_packages)) {
  stop("Missing reporting package(s): ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

config <- analysis$config
reporting_order <- config$reporting_order
mode_suffix <- if (smoke) "_smoke" else ""
pmse_digits <- if (identical(report_version, "v2")) 6L else 4L
mode_caption <- if (smoke) {
  "Technical smoke (noninferential)"
} else {
  "Production benchmark"
}

dir.create(generated_parent, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_parent, recursive = TRUE, showWarnings = FALSE)
reporting_staging <- file.path(
  generated_parent,
  paste0(".", basename(reporting_dir), "_staging_", Sys.getpid())
)
pdf_staging <- file.path(
  pdf_parent, paste0(".", basename(pdf_dir), "_staging_", Sys.getpid())
)
if (dir.exists(reporting_staging) || dir.exists(pdf_staging)) {
  stop("Reporting staging directory already exists.", call. = FALSE)
}
dir.create(file.path(reporting_staging, "figures"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_staging, recursive = TRUE, showWarnings = FALSE)
on.exit({
  if (dir.exists(reporting_staging)) unlink(reporting_staging, recursive = TRUE)
  if (dir.exists(pdf_staging)) unlink(pdf_staging, recursive = TRUE)
}, add = TRUE)

format_number <- function(value, digits = 4L) {
  ifelse(
    is.finite(value),
    ifelse(abs(value) > 0 & abs(value) < 10^(-digits),
           formatC(value, format = "e", digits = 2L),
           formatC(value, format = "f", digits = digits)),
    "--"
  )
}
format_percent <- function(value) {
  ifelse(is.finite(value), paste0(formatC(100 * value, format = "f",
                                          digits = 1L), "\\%"), "--")
}

write_kable <- function(data, path, caption, label, align = NULL,
                        font_size = 8L) {
  table <- knitr::kable(
    data, format = "latex", booktabs = TRUE, escape = FALSE,
    caption = caption, label = label, align = align, linesep = ""
  )
  table <- kableExtra::kable_styling(
    table, latex_options = c("hold_position"),
    font_size = font_size, position = "center"
  )
  writeLines(as.character(table), path, useBytes = TRUE)
}

case_ii <- analysis$case_ii_performance
case_ii$method <- factor(case_ii$method, levels = reporting_order,
                         ordered = TRUE)
case_ii <- case_ii[order(case_ii$method), , drop = FALSE]
case_ii_table <- data.frame(
  Method = case_ii$label,
  PMSE = format_number(case_ii$response_prediction_mse, pmse_digits),
  `$\\mathrm{ISE}_0$` = format_number(case_ii$ise_null),
  `$\\mathrm{ISE}_1$` = format_number(case_ii$ise_nonnull),
  `True zero` = format_percent(case_ii$true_zero_rate),
  `False zero` = format_percent(case_ii$false_zero_rate),
  check.names = FALSE, stringsAsFactors = FALSE
)
write_kable(
  case_ii_table,
  file.path(reporting_staging, paste0("table_case_ii_benchmark_",
                                      report_version,
                                      mode_suffix, ".tex")),
  paste0(
    mode_caption,
    ". Case-II performance on the common test samples. PMSE is response ",
    "prediction MSE; exact-zero entries are not applicable to nonsparse ",
    "estimators."
  ),
  "arfc-legacy-caseii", align = c("l", rep("r", 5L))
)

pmse_table <- reshape(
  analysis$safety_summary[c(
    "method", "label", "case_label", "response_prediction_mse"
  )],
  idvar = c("method", "label"), timevar = "case_label",
  direction = "wide"
)
pmse_table$method <- factor(pmse_table$method, levels = reporting_order,
                            ordered = TRUE)
pmse_table <- pmse_table[order(pmse_table$method), , drop = FALSE]
case_columns <- paste0("response_prediction_mse.",
                       c("Case I", "Case II", "Case III", "Case IV"))
pmse_display <- data.frame(
  Method = pmse_table$label,
  `Case I` = format_number(pmse_table[[case_columns[1L]]], pmse_digits),
  `Case II` = format_number(pmse_table[[case_columns[2L]]], pmse_digits),
  `Case III` = format_number(pmse_table[[case_columns[3L]]], pmse_digits),
  `Case IV` = format_number(pmse_table[[case_columns[4L]]], pmse_digits),
  check.names = FALSE, stringsAsFactors = FALSE
)
write_kable(
  pmse_display,
  file.path(reporting_staging, paste0("table_four_case_pmse_",
                                      report_version,
                                      mode_suffix, ".tex")),
  paste0(mode_caption,
         ". Response prediction MSE by case; no cross-case pooling is used."),
  "arfc-legacy-pmse", align = c("l", rep("r", 4L))
)

structural_methods <- reporting_order[reporting_order %in% c(
  "arfc_select", "asr_mcp", "amr_mcp", "asr_scad", "amr_scad",
  "fixed_fmcp", "fixed_fscad", "flirti", "oracle"
)]
safety <- analysis$safety_summary[
  analysis$safety_summary$method %in% structural_methods &
    analysis$safety_summary$scenario %in%
      c("case_i_null", "case_iii_crossing", "case_iv_away"), , drop = FALSE
]
safety_wide <- reshape(
  safety[c("method", "label", "case_label", "safety_value")],
  idvar = c("method", "label"), timevar = "case_label", direction = "wide"
)
safety_wide$method <- factor(safety_wide$method, levels = structural_methods,
                             ordered = TRUE)
safety_wide <- safety_wide[order(safety_wide$method), , drop = FALSE]
safety_columns <- paste0(
  "safety_value.", c("Case I", "Case III", "Case IV")
)
safety_display <- data.frame(
  Method = safety_wide$label,
  `Case I false active` = format_percent(safety_wide[[safety_columns[1L]]]),
  `Case III false zero` = format_percent(safety_wide[[safety_columns[2L]]]),
  `Case IV false zero` = format_percent(safety_wide[[safety_columns[3L]]]),
  check.names = FALSE, stringsAsFactors = FALSE
)
write_kable(
  safety_display,
  file.path(reporting_staging, paste0("table_structural_safety_",
                                      report_version,
                                      mode_suffix, ".tex")),
  paste0(
    mode_caption,
    ". Structural safety metrics. Oracle is unavailable in Cases III--IV."
  ),
  "arfc-legacy-safety", align = c("l", rep("r", 3L))
)

curve_data <- analysis$curve_summary
coefficient_methods <- c(
  "arfc_select", "fixed_fmcp", "fixed_fscad", "smooth", "flirti", "oracle"
)
coefficient_data <- curve_data[
  curve_data$method %in% coefficient_methods, , drop = FALSE
]
coefficient_data$label <- factor(
  coefficient_data$label,
  levels = unname(config$method_labels[coefficient_methods]), ordered = TRUE
)
coefficient_colors <- c(
  "ARFC--Select" = "#0072B2", "fMCP" = "#D55E00",
  "fSCAD" = "#009E73", "Smooth" = "#777777",
  "FLiRTI" = "#CC79A7", "Oracle" = "#E69F00"
)
coefficient_types <- c(
  "ARFC--Select" = "solid", "fMCP" = "longdash",
  "fSCAD" = "dotdash", "Smooth" = "dotted",
  "FLiRTI" = "twodash", "Oracle" = "dashed"
)
arfc_band <- coefficient_data[coefficient_data$method == "arfc_select", ]
coefficient_plot <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect", xmin = 0.3, xmax = 0.7, ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.20
  ) +
  ggplot2::geom_ribbon(
    data = arfc_band,
    ggplot2::aes(x = t, ymin = q05, ymax = q95),
    fill = coefficient_colors[["ARFC--Select"]], alpha = 0.12
  ) +
  ggplot2::geom_line(
    data = coefficient_data,
    ggplot2::aes(x = t, y = mean, colour = label, linetype = label),
    linewidth = 0.8, na.rm = TRUE
  ) +
  ggplot2::geom_line(
    data = unique(coefficient_data[c("t", "truth")]),
    ggplot2::aes(x = t, y = truth), colour = "black",
    linewidth = 0.9, linetype = "dashed"
  ) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
  ggplot2::scale_colour_manual(values = coefficient_colors, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = coefficient_types, drop = FALSE) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::labs(
    title = "Case II coefficient recovery",
    subtitle = paste0(
      mode_caption,
      "; grey = true null; blue ribbon = ARFC-Select 90% envelope"
    ),
    x = "t", y = expression(hat(beta)(t)), colour = NULL, linetype = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom", legend.box = "horizontal",
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 9)
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  )

structural_plot_methods <- structural_methods
zero_data <- curve_data[
  curve_data$method %in% structural_plot_methods &
    is.finite(curve_data$zero_frequency), , drop = FALSE
]
zero_data$panel <- ifelse(
  zero_data$method %in%
    c("arfc_select", "asr_mcp", "amr_mcp", "asr_scad", "amr_scad"),
  "ARFC selector and candidates", "Fixed and literature comparators"
)
zero_labels <- unname(config$method_labels[structural_plot_methods])
zero_data$label <- factor(zero_data$label, levels = zero_labels, ordered = TRUE)
zero_colors <- c(
  "ARFC--Select" = "#000000", "ASR--MCP" = "#0072B2",
  "AMR--MCP" = "#56B4E9", "ASR--SCAD" = "#D55E00",
  "AMR--SCAD" = "#E69F00", "fMCP" = "#6A3D9A",
  "fSCAD" = "#009E73", "FLiRTI" = "#CC79A7", "Oracle" = "#666666"
)
zero_types <- setNames(
  rep(c("solid", "longdash", "dotdash", "dashed", "twodash"),
      length.out = length(zero_labels)),
  zero_labels
)
zero_plot <- ggplot2::ggplot(zero_data) +
  ggplot2::annotate(
    "rect", xmin = 0.3, xmax = 0.7, ymin = -Inf, ymax = Inf,
    fill = "grey70", alpha = 0.20
  ) +
  ggplot2::geom_vline(
    xintercept = c(0.3, 0.7), colour = "grey45", linetype = "dotted",
    linewidth = 0.35
  ) +
  ggplot2::geom_line(
    ggplot2::aes(x = t, y = zero_frequency, colour = label,
                 linetype = label), linewidth = 0.8, na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~panel, ncol = 1) +
  ggplot2::scale_colour_manual(values = zero_colors, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = zero_types, drop = FALSE) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::labs(
    title = "Case II pointwise exact-zero selection",
    subtitle = paste0(
      mode_caption,
      "; nonsparse Smooth, OLS, and PCR are intentionally omitted"
    ),
    x = "t", y = "Selection frequency", colour = NULL, linetype = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey94"),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(size = 9)
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
  )

coefficient_base <- paste0(
  "case_ii_coefficient_recovery_", report_version, mode_suffix
)
zero_base <- paste0("case_ii_zero_selection_", report_version, mode_suffix)
figure_dir <- file.path(reporting_staging, "figures")
ggplot2::ggsave(
  file.path(figure_dir, paste0(coefficient_base, ".pdf")),
  coefficient_plot, width = 8.2, height = 5.5, units = "in",
  device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  file.path(figure_dir, paste0(coefficient_base, ".png")),
  coefficient_plot, width = 8.2, height = 5.5, units = "in", dpi = 300
)
ggplot2::ggsave(
  file.path(figure_dir, paste0(zero_base, ".pdf")),
  zero_plot, width = 8.2, height = 7.2, units = "in",
  device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  file.path(figure_dir, paste0(zero_base, ".png")),
  zero_plot, width = 8.2, height = 7.2, units = "in", dpi = 300
)

figure_relative_root <- file.path(
  "generated", prefix, paste0(mode, "_reporting_", report_version),
  "figures"
)
figure_fragment <- c(
  "\\begin{figure}[!htbp]",
  "\\centering",
  paste0("\\includegraphics[width=0.92\\linewidth]{",
         file.path(figure_relative_root, paste0(coefficient_base, ".pdf")),
         "}"),
  paste0("\\caption{", mode_caption,
         ". Case-II mean coefficient recovery. The dashed black curve is ",
         "the truth; the grey region is the true null interval.}"),
  "\\label{fig:arfc-legacy-coefficient}",
  "\\end{figure}",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  paste0("\\includegraphics[width=0.92\\linewidth]{",
         file.path(figure_relative_root, paste0(zero_base, ".pdf")), "}"),
  paste0("\\caption{", mode_caption,
         ". Pointwise exact-zero selection frequencies for structurally ",
         "sparse estimators. Smooth, OLS, and PCR do not produce structural ",
         "local zeros and are omitted from this panel.}"),
  "\\label{fig:arfc-legacy-zero}",
  "\\end{figure}"
)
writeLines(
  figure_fragment,
  file.path(reporting_staging, paste0("figure_includes_", report_version,
                                      mode_suffix,
                                      ".tex")),
  useBytes = TRUE
)

methods_fragment <- c(
  "\\paragraph{Secondary classical-comparator benchmark.}",
  paste0(
    "The frozen ARFC--Select estimator and its four candidates were compared ",
    "on the identical task-matched samples with fixed-resolution fMCP and ",
    "fSCAD, Smooth, OLS, PCR, FLiRTI, and a truth-aware diagnostic Oracle. ",
    "The comparator library was not eligible for ARFC selection. fMCP and ",
    "fSCAD used the same $M=32$ basis, task-specific smooth-selected ",
    "$\\kappa$, coefficient cutoff, active-set refit, and lambda path, with ",
    "$\\gamma=2.7$ and $a=3.7$ to match scalar maximum concavity. No separate ",
    "SLoS row was introduced. The two-stage estimator was not assigned a ",
    "numerical row because no source-verified implementation was available."
  )
)
writeLines(
  methods_fragment,
  file.path(reporting_staging, paste0("methods_benchmark_", report_version,
                                      mode_suffix,
                                      ".tex")),
  useBytes = TRUE
)

if (smoke) {
  results_fragment <- paste(
    "\\paragraph{Technical smoke.}",
    "The four-task smoke verified the complete comparator and reporting path.",
    "Its numerical values are noninferential and are not scientific evidence."
  )
} else {
  case_ii_arfc <- case_ii[as.character(case_ii$method) == "arfc_select", ]
  case_ii_fmcp <- case_ii[as.character(case_ii$method) == "fixed_fmcp", ]
  case_ii_fscad <- case_ii[as.character(case_ii$method) == "fixed_fscad", ]
  results_fragment <- paste0(
    "\\paragraph{Classical-comparator benchmark.} In Case II, ARFC--Select ",
    "had response PMSE ",
    format_number(case_ii_arfc$response_prediction_mse, pmse_digits),
    ", compared with ",
    format_number(case_ii_fmcp$response_prediction_mse, pmse_digits),
    " for fMCP and ",
    format_number(case_ii_fscad$response_prediction_mse, pmse_digits),
    " for fSCAD. Its true-zero rate was ",
    format_percent(case_ii_arfc$true_zero_rate),
    ". All four cases were summarized separately; no cross-case overall ",
    "average was used."
  )
}
writeLines(
  results_fragment,
  file.path(reporting_staging, paste0("results_benchmark_", report_version,
                                      mode_suffix,
                                      ".tex")),
  useBytes = TRUE
)

figure_files <- file.path(
  figure_dir,
  c(paste0(coefficient_base, c(".pdf", ".png")),
    paste0(zero_base, c(".pdf", ".png")))
)
png_files <- figure_files[grepl("[.]png$", figure_files)]
pdf_files <- figure_files[grepl("[.]pdf$", figure_files)]
png_checks <- vapply(png_files, function(path) {
  image <- png::readPNG(path)
  length(dim(image)) >= 2L && all(dim(image)[1:2] > 500L) &&
    is.finite(stats::sd(as.numeric(image))) && stats::sd(as.numeric(image)) > 0
}, logical(1L))
pdfinfo <- Sys.which("pdfinfo")
pdf_checks <- if (nzchar(pdfinfo)) {
  vapply(pdf_files, function(path) {
    output <- suppressWarnings(system2(pdfinfo, path, stdout = TRUE,
                                       stderr = TRUE))
    any(grepl("^Pages:[[:space:]]+1$", output))
  }, logical(1L))
} else {
  rep(FALSE, length(pdf_files))
}
if (!all(png_checks) || !all(pdf_checks)) {
  stop("Figure readability or one-page PDF gate failed.", call. = FALSE)
}

writeLines(
  capture.output(sessionInfo()),
  file.path(reporting_staging,
            paste0(prefix, "_reporting_sessionInfo_", report_version,
                   mode_suffix, ".txt"))
)
relative_files <- list.files(reporting_staging, recursive = TRUE,
                             full.names = FALSE)
relative_files <- relative_files[!grepl("manifest", relative_files)]
payload_paths <- file.path(reporting_staging, relative_files)
manifest <- data.frame(
  scientific_signature = analysis$scientific_signature,
  role = ifelse(
    grepl("[.]pdf$", relative_files), "figure_pdf",
    ifelse(grepl("[.]png$", relative_files), "figure_png",
           ifelse(grepl("[.]tex$", relative_files), "latex", "other"))
  ),
  file = relative_files,
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(
  reporting_staging, paste0(prefix, "_reporting_manifest_", report_version,
                            mode_suffix, ".csv")
)
utils::write.csv(manifest, manifest_path, row.names = FALSE)

for (path in pdf_files) {
  destination <- file.path(pdf_staging, basename(path))
  if (!file.copy(path, destination, overwrite = FALSE, copy.mode = TRUE)) {
    stop("PDF mirror staging failed: ", basename(path), call. = FALSE)
  }
}
mirror_paths <- list.files(pdf_staging, full.names = TRUE)
mirror_manifest <- data.frame(
  scientific_signature = analysis$scientific_signature,
  file = basename(mirror_paths),
  bytes = as.numeric(file.info(mirror_paths)$size),
  md5 = unname(as.character(tools::md5sum(mirror_paths))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  mirror_manifest,
  file.path(pdf_staging, paste0(prefix, "_pdf_manifest_", report_version,
                                mode_suffix, ".csv")),
  row.names = FALSE
)

if (!file.rename(reporting_staging, reporting_dir)) {
  stop("Atomic reporting bundle install failed.", call. = FALSE)
}
if (!file.rename(pdf_staging, pdf_dir)) {
  stop("Atomic PDF mirror install failed.", call. = FALSE)
}
cat("ARFC classical-comparator reporting complete.\n")
cat("Tables: 3\n")
cat("Figures: 2 PDF + 2 PNG\n")
cat("Reporting bundle:", reporting_dir, "\n")
cat("PDF mirror:", pdf_dir, "\n")
