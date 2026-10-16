#!/usr/bin/env Rscript

# Write the compact ARFC-Select real-data table, figure, and LaTeX fragments.
# This script is reporting-only and never refits a statistical model.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("real_data/28_write_arfc_real_data_v1.R", mustWork = TRUE)
}
real_dir <- dirname(script_file)
project_dir <- dirname(real_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
report_version_argument <- grep("^--report-version=", arguments, value = TRUE)
if (length(report_version_argument) > 1L) {
  stop("Use at most one --report-version argument.", call. = FALSE)
}
report_version <- if (length(report_version_argument)) {
  sub("^--report-version=", "", report_version_argument)
} else {
  "v1"
}
if (!report_version %in% c("v1", "v2", "v3", "v4")) {
  stop("--report-version must be v1, v2, v3, or v4.", call. = FALSE)
}
unknown <- setdiff(arguments, c("--smoke", report_version_argument))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

if (!requireNamespace("knitr", quietly = TRUE) ||
    !requireNamespace("kableExtra", quietly = TRUE)) {
  stop("Reporting requires the declared knitr and kableExtra dependencies.",
       call. = FALSE)
}
helper_file <- file.path(real_dir, "R", "arfc_real_data_v1.R")
source(helper_file, local = FALSE)
config <- arfc_real_config(project_dir, smoke)
mode_prefix <- if (smoke) "smoke" else "production"
prefix <- paste0("arfc_real_data_v1_", mode_prefix)
analysis_dir <- file.path(
  config$results_root, paste0(mode_prefix, "_analysis_v1")
)
report_dir <- file.path(
  config$output_root, paste0(mode_prefix, "_reporting_", report_version)
)
figure_dir <- file.path(
  config$figures_root, paste0(mode_prefix, "_figures_", report_version)
)
if (!dir.exists(analysis_dir)) {
  stop("Completed analysis bundle is missing: ", analysis_dir,
       call. = FALSE)
}
if (dir.exists(report_dir) || dir.exists(figure_dir)) {
  stop("Reporting or figure bundle already exists; refusing overwrite.",
       call. = FALSE)
}
manifest_path <- file.path(
  analysis_dir, paste0(prefix, "_analysis_manifest.csv")
)
analysis_path <- file.path(analysis_dir, paste0(prefix, "_analysis.rds"))
manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
paths <- file.path(analysis_dir, manifest$file)
if (!all(file.exists(paths)) ||
    any(as.numeric(file.info(paths)$size) != manifest$bytes) ||
    any(arfc_real_md5(paths) != manifest$md5)) {
  stop("Analysis manifest verification failed.", call. = FALSE)
}
analysis <- readRDS(analysis_path)
if (!identical(analysis$mode, config$mode) ||
    !identical(analysis$scientific_signature,
               unique(manifest$scientific_signature)) ||
    !all(analysis$validation$passed)) {
  stop("Analysis identity/validation guard failed.", call. = FALSE)
}
config <- analysis$config

label_lookup <- c(
  `ARFC-Select` = "ARFC-Select",
  asr_mcp = "ASR--MCP", asr_scad = "ASR--SCAD",
  amr_mcp = "AMR--MCP", amr_scad = "AMR--SCAD",
  Smooth = "Smooth", PCR = "PCR"
)
table_data <- analysis$prediction_metrics
table_data$Method <- unname(label_lookup[table_data$method])
selection_lookup <- analysis$selection_frequency[c(
  "dataset", "method", "frequency"
)]
names(selection_lookup)[3L] <- "nested_selection_frequency"
bootstrap_lookup <- analysis$bootstrap_selection_frequency[c(
  "dataset", "method", "frequency"
)]
names(bootstrap_lookup)[3L] <- "bootstrap_selection_frequency"
table_data <- merge(
  table_data, selection_lookup, by = c("dataset", "method"),
  all.x = TRUE, sort = FALSE
)
table_data <- merge(
  table_data, bootstrap_lookup, by = c("dataset", "method"),
  all.x = TRUE, sort = FALSE
)
table_data$Full <- ""
for (index in seq_len(nrow(analysis$full_selection))) {
  selected <- analysis$full_selection[index, ]
  table_data$Full[
    table_data$dataset == selected$dataset &
      table_data$method == selected$selected_method
  ] <- "selected"
}
method_order <- c(config$selection_label, config$method_order,
                  config$secondary_methods)
table_data <- table_data[
  order(match(table_data$dataset, config$dataset_order),
        match(table_data$method, method_order)), , drop = FALSE
]
if (identical(report_version, "v4")) {
  publication_table <- data.frame(
    Dataset = table_data$dataset,
    Method = table_data$Method,
    RMSE = sprintf("%.4f", table_data$RMSE),
    MAE = sprintf("%.4f", table_data$MAE),
    `$R^2$` = sprintf("%.3f", table_data$R2),
    `Nested selection` = ifelse(
      is.finite(table_data$nested_selection_frequency),
      sprintf("%.1f\\%%", 100 * table_data$nested_selection_frequency), "--"
    ),
    `Bootstrap selection` = ifelse(
      is.finite(table_data$bootstrap_selection_frequency),
      sprintf("%.1f\\%%", 100 * table_data$bootstrap_selection_frequency),
      "--"
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )
} else {
  publication_table <- data.frame(
    Dataset = table_data$dataset,
    Method = table_data$Method,
    RMSE = sprintf("%.4f", table_data$RMSE),
    MAE = sprintf("%.4f", table_data$MAE),
    R2 = sprintf("%.3f", table_data$R2),
    `Nested selection` = ifelse(
      is.finite(table_data$nested_selection_frequency),
      sprintf("%.1f%%", 100 * table_data$nested_selection_frequency), "--"
    ),
    `Bootstrap selection` = ifelse(
      is.finite(table_data$bootstrap_selection_frequency),
      sprintf("%.1f%%", 100 * table_data$bootstrap_selection_frequency), "--"
    ),
    `Full-data choice` = ifelse(nzchar(table_data$Full), "yes", ""),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

report_parent <- dirname(report_dir)
figure_parent <- dirname(figure_dir)
dir.create(report_parent, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_parent, recursive = TRUE, showWarnings = FALSE)
report_staging <- file.path(
  report_parent, paste0(".", basename(report_dir), "_staging_", Sys.getpid())
)
figure_staging <- file.path(
  figure_parent, paste0(".", basename(figure_dir), "_staging_", Sys.getpid())
)
if (dir.exists(report_staging) || dir.exists(figure_staging)) {
  stop("Reporting staging directory already exists.", call. = FALSE)
}
dir.create(report_staging, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed) {
  if (dir.exists(report_staging)) unlink(report_staging, recursive = TRUE)
  if (dir.exists(figure_staging)) unlink(figure_staging, recursive = TRUE)
}, add = TRUE)

table_csv <- paste0(prefix, "_table_prediction_selection.csv")
table_tex <- paste0(prefix, "_table_prediction_selection.tex")
utils::write.csv(publication_table, file.path(report_staging, table_csv),
                 row.names = FALSE)
kable <- knitr::kable(
  publication_table, format = "latex", booktabs = TRUE,
  escape = !identical(report_version, "v4"),
  caption = paste(
    "Dataset-specific nested-CV prediction and candidate selection.",
    "Smooth and PCR are secondary benchmarks and are not selectable",
    "ARFC library members."
  ),
  label = "arfc-real-prediction",
  align = if (identical(report_version, "v4")) "llrrrcc" else "llrrrccc"
)
kable <- kableExtra::kable_styling(
  kable, latex_options = c("hold_position", "scale_down"),
  font_size = 8, position = "center"
)
if (identical(report_version, "v4")) {
  kable <- kableExtra::collapse_rows(
    kable, columns = 1L, valign = "top", latex_hline = "major"
  )
}
table_latex <- as.character(kable)
if (report_version %in% c("v3", "v4")) {
  table_latex <- sub(
    "\\begin{table}[!h]", "\\begin{table}[!htbp]",
    table_latex, fixed = TRUE
  )
}
writeLines(table_latex, file.path(report_staging, table_tex),
           useBytes = TRUE)

figure_pdf <- paste0(prefix, "_coefficient_stability.pdf")
figure_png <- paste0(prefix, "_coefficient_stability.png")
draw_figure <- function() {
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(
    mfrow = c(2, 2), mar = c(4.6, 4.2, 3.2, 1.1),
    oma = c(1.5, 0.5, 2.0, 0.5), las = 1
  )
  for (dataset_name in config$dataset_order) {
    stability <- analysis$stability_summary[
      analysis$stability_summary$dataset == dataset_name, , drop = FALSE
    ]
    full <- analysis$full_selected_curves[
      analysis$full_selected_curves$dataset == dataset_name, , drop = FALSE
    ]
    selected_method <- unique(full$selected_method)
    x <- stability$physical_argument
    y_range <- range(c(stability$beta_q10, stability$beta_q90,
                       full$beta_physical), finite = TRUE)
    padding <- 0.05 * max(diff(y_range), .Machine$double.eps)
    graphics::plot(
      x, stability$beta_mean, type = "n",
      ylim = y_range + c(-padding, padding),
      xlab = "",
      ylab = expression(hat(beta)(t)),
      main = paste0(dataset_name, ": full choice ",
                    label_lookup[selected_method])
    )
    graphics::polygon(
      c(x, rev(x)),
      c(stability$beta_q10, rev(stability$beta_q90)),
      col = grDevices::adjustcolor("#56B4E9", alpha.f = 0.25),
      border = NA
    )
    graphics::abline(h = 0, col = "grey60", lty = 3)
    graphics::lines(x, stability$beta_mean, col = "#0072B2", lwd = 2)
    graphics::lines(full$physical_argument, full$beta_physical,
                    col = "#D55E00", lwd = 1.8, lty = 2)
    graphics::legend(
      "topright", legend = c("Bootstrap mean", "10%--90% band",
                              "Full-data fit"),
      col = c("#0072B2", "#56B4E9", "#D55E00"),
      lty = c(1, NA, 2), lwd = c(2, NA, 1.8),
      pch = c(NA, 15, NA), pt.cex = c(NA, 1.4, NA),
      bty = "n", cex = 0.78
    )
  }
  for (dataset_name in config$dataset_order) {
    stability <- analysis$stability_summary[
      analysis$stability_summary$dataset == dataset_name, , drop = FALSE
    ]
    x <- stability$physical_argument
    graphics::plot(
      x, stability$zero_selection_frequency, type = "l",
      ylim = c(0, 1), lwd = 2.2, col = "#009E73",
      xlab = if (dataset_name == "Truck")
        "Acceleration measurement (seconds)" else "Wavelength (nm)",
      ylab = "Exact-zero selection frequency",
      main = paste0(dataset_name, ": full-reselection bootstrap")
    )
    graphics::abline(h = 0.80, lty = 2, col = "grey45")
    segments <- analysis$stable_zero_segments[
      analysis$stable_zero_segments$dataset == dataset_name, , drop = FALSE
    ]
    if (nrow(segments)) {
      for (index in seq_len(nrow(segments))) {
        graphics::rect(
          segments$physical_start[index], 0,
          segments$physical_end[index], 1,
          col = grDevices::adjustcolor("grey60", alpha.f = 0.15),
          border = NA
        )
      }
      graphics::lines(x, stability$zero_selection_frequency,
                      lwd = 2.2, col = "#009E73")
    } else {
      graphics::mtext("No region reaches 80% stability",
                      side = 3, line = 0.2, cex = 0.75, col = "grey35")
    }
  }
  graphics::mtext(
    if (smoke) "ARFC real-data technical smoke (noninferential)" else
      "ARFC-Select real-data coefficient and local-zero stability",
    side = 3, outer = TRUE, line = 0.4, font = 2, cex = 1.05
  )
}
if (report_version %in% c("v3", "v4")) {
  if (!capabilities("cairo")) {
    stop("Cairo graphics support is required for embedded-font v3 PDF.",
         call. = FALSE)
  }
  grDevices::cairo_pdf(
    file.path(report_staging, figure_pdf),
    width = 10.5, height = 8.2, onefile = TRUE, family = "sans"
  )
} else {
  grDevices::pdf(file.path(report_staging, figure_pdf),
                 width = 10.5, height = 8.2, useDingbats = FALSE)
}
draw_figure()
grDevices::dev.off()
grDevices::png(file.path(report_staging, figure_png),
               width = 2520, height = 1968, res = 240,
               type = "cairo")
draw_figure()
grDevices::dev.off()
copied <- file.copy(
  file.path(report_staging, c(figure_pdf, figure_png)),
  figure_staging, overwrite = FALSE
)
if (!all(copied)) stop("Figure-copy staging failed.", call. = FALSE)

figure_tex <- paste0(prefix, "_figure_include.tex")
figure_relative <- file.path(
  "..", "output", "pdf", "arfc_real_data_v1",
  paste0(mode_prefix, "_reporting_", report_version), figure_pdf
)
writeLines(c(
  if (report_version %in% c("v2", "v3", "v4"))
    "\\begin{figure}[!htbp]" else
    "\\begin{figure}[H]",
  "\\centering",
  paste0("\\includegraphics[width=",
         if (report_version %in% c("v2", "v3", "v4")) "0.93" else "0.98",
         "\\linewidth]{",
         figure_relative, "}"),
  paste0(
    "\\caption{Dataset-specific ARFC-Select coefficient estimates and ",
    "sample-level full-reselection bootstrap stability. Bands are the ",
    "pointwise 10th--90th percentiles. The lower panels show exact-zero ",
    "selection frequency; the dashed horizontal line is 80\\%. Truck rows ",
    "are the distributed sampling units, while NeoSpectra rows are unique ",
    "soil-sample IDs.}"
  ),
  "\\label{fig:arfc-real-stability}",
  "\\end{figure}"
), file.path(report_staging, figure_tex), useBytes = TRUE)

results_tex <- paste0(prefix, "_results_text.tex")
result_lines <- c()
if (smoke) {
  result_lines <- c(
    "\\paragraph{Technical smoke only.}",
    paste(
      "The real-data smoke uses two outer tasks and three bootstrap",
      "reselections per dataset. It validates computation and reporting and",
      "has no inferential interpretation."
    )
  )
} else {
  for (dataset_name in config$dataset_order) {
    metric <- analysis$prediction_metrics[
      analysis$prediction_metrics$dataset == dataset_name &
        analysis$prediction_metrics$method == config$selection_label, ]
    full <- analysis$full_selection[
      analysis$full_selection$dataset == dataset_name, ]
    segments <- analysis$stable_zero_segments[
      analysis$stable_zero_segments$dataset == dataset_name, ]
    segment_text <- if (nrow(segments)) {
      paste0(
        " The 80\\% full-reselection stability rule identified ",
        nrow(segments), " segment(s); the longest ran from ",
        sprintf("%.1f", segments$physical_start[
          which.max(segments$physical_length)
        ]), " to ",
        sprintf("%.1f", segments$physical_end[
          which.max(segments$physical_length)
        ]), "."
      )
    } else {
      " No argument region reached 80\\% exact-zero stability."
    }
    result_lines <- c(
      result_lines,
      paste0(
        "For \\textbf{", dataset_name, "}, ARFC-Select achieved outer-CV ",
        "RMSE ", sprintf("%.4f", metric$RMSE), ", MAE ",
        sprintf("%.4f", metric$MAE), ", and $R^2$ ",
        sprintf("%.3f", metric$R2), ". The full-data candidate was ",
        label_lookup[full$selected_method], ".", segment_text
      )
    )
  }
}
writeLines(result_lines, file.path(report_staging, results_tex),
           useBytes = TRUE)

limitations_tex <- paste0(prefix, "_limitations.tex")
writeLines(c(
  paste0(
    "Truck lacks truck/test-cycle identifiers; its resampling unit is the ",
    "distributed row, and unobserved within-vehicle dependence could make ",
    "uncertainty optimistic. NeoSpectra has unique sample-ID resampling but ",
    "no scanner or geographic metadata. Real data have no known coefficient ",
    "truth, so exact-zero stability is not true-zero recovery and cannot ",
    "establish causal irrelevance. The post-selection outer-training refit is ",
    "a practical nested-CV step outside the literal fixed-candidate statement ",
    "of the finite-library proposition."
  )
), file.path(report_staging, limitations_tex), useBytes = TRUE)

session_file <- paste0(prefix, "_reporting_sessionInfo.txt")
writeLines(capture.output(sessionInfo()),
           file.path(report_staging, session_file), useBytes = TRUE)
payloads <- c(
  table_csv, table_tex, figure_pdf, figure_png, figure_tex,
  results_tex, limitations_tex, session_file
)
manifest_name <- paste0(prefix, "_reporting_manifest.csv")
invisible(arfc_real_manifest(
  report_staging, analysis$scientific_signature, payloads, manifest_name
))
if (!file.rename(figure_staging, figure_dir)) {
  stop("Atomic figure installation failed: ", figure_dir, call. = FALSE)
}
if (!file.rename(report_staging, report_dir)) {
  stop("Atomic reporting installation failed: ", report_dir,
       call. = FALSE)
}
committed <- TRUE

cat("ARFC-Select real-data reporting complete.\n")
cat("Mode:", analysis$mode, "\n")
cat("Table rows:", nrow(publication_table), "\n")
cat("Figure:", file.path(report_dir, figure_pdf), "\n")
cat("Bundle:", report_dir, "\n")
