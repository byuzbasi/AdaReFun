#!/usr/bin/env Rscript

# Reporting-only v8 refinement of the frozen production Case-II figures.
# No model fitting or Monte Carlo simulation is performed here.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath(
    "simulation/59_refine_arfc_legacy_figures_v8.R",
    mustWork = TRUE
  )
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)
prefix <- "arfc_legacy_benchmark_v1"

analysis_dir <- file.path(
  simulation_dir, "results", prefix, "production_analysis_v1"
)
analysis_file <- file.path(
  analysis_dir, paste0(prefix, "_analysis_v1.rds")
)
analysis_manifest_file <- file.path(
  analysis_dir, paste0(prefix, "_analysis_manifest_v1.csv")
)
generated_parent <- file.path(
  project_dir, "config", "generated", prefix
)
reporting_dir <- file.path(
  generated_parent, "production_figure_refinement_v8"
)
pdf_parent <- file.path(project_dir, "output", "pdf", prefix)
pdf_dir <- file.path(pdf_parent, "production_figure_refinement_v8")

if (!file.exists(analysis_file) || !file.exists(analysis_manifest_file)) {
  stop("Validated production analysis is missing.", call. = FALSE)
}
if (dir.exists(reporting_dir) || dir.exists(pdf_dir)) {
  stop("Versioned v8 output already exists; refusing overwrite.",
       call. = FALSE)
}

required_packages <- c("ggplot2", "cowplot", "png")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE,
  FUN.VALUE = logical(1L)
)]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
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
signature <- unique(analysis_manifest$scientific_signature)
if (length(signature) != 1L ||
    !identical(signature, analysis$scientific_signature) ||
    !identical(analysis$study, prefix) ||
    !identical(analysis$mode, "production") ||
    !all(analysis$validation$passed)) {
  stop("Frozen production-analysis identity or validation failed.",
       call. = FALSE)
}

case_curves <- analysis$curve_summary[
  analysis$curve_summary$scenario == "case_ii_sparse", , drop = FALSE
]
if (nrow(case_curves) != 201L * 12L ||
    !all(case_curves$replications == 100L)) {
  stop("Case-II curve replication contract failed.", call. = FALSE)
}

coefficient_ids <- c(
  "arfc_select", "fixed_fmcp", "fixed_fscad",
  "smooth", "flirti", "oracle"
)
coefficient_labels <- unname(
  analysis$config$method_labels[coefficient_ids]
)
coefficient_data <- case_curves[
  case_curves$method %in% coefficient_ids, , drop = FALSE
]
coefficient_data$label <- factor(
  coefficient_data$label,
  levels = coefficient_labels,
  ordered = TRUE
)
truth <- unique(coefficient_data[c("t", "truth")])
truth <- truth[order(truth$t), , drop = FALSE]

case_ii <- analysis$case_ii_performance
rate_ids <- c("arfc_select", "fixed_fmcp", "fixed_fscad", "flirti")
rates <- setNames(
  case_ii$true_zero_rate[match(rate_ids, case_ii$method)],
  rate_ids
)
expected_rates <- c(0.539400, 0.276925, 0.342550, 0.196000)
if (any(!is.finite(rates)) ||
    !isTRUE(all.equal(unname(rates), expected_rates, tolerance = 5e-7))) {
  stop("Exact-zero callouts differ from the frozen production analysis.",
       call. = FALSE)
}

coefficient_colours <- c(
  "ARFC--Select" = "#0072B2", "fMCP" = "#D55E00",
  "fSCAD" = "#009E73", "Smooth" = "#6B6B6B",
  "FLiRTI" = "#CC79A7", "Oracle" = "#E69F00"
)
coefficient_linetypes <- c(
  "ARFC--Select" = "solid", "fMCP" = "longdash",
  "fSCAD" = "dotdash", "Smooth" = "dotted",
  "FLiRTI" = "twodash", "Oracle" = "dashed"
)
band_colours <- c(
  "ARFC--Select" = "#CFE8F3", "fMCP" = "#F5DCCF",
  "fSCAD" = "#D2EDE3"
)
band_data <- coefficient_data[
  coefficient_data$method %in%
    c("arfc_select", "fixed_fmcp", "fixed_fscad"),
  , drop = FALSE
]

theme_coefficient <- ggplot2::theme_bw(
  base_size = 10.5, base_family = "serif"
) + ggplot2::theme(
  legend.position = "bottom",
  legend.box = "horizontal",
  legend.key.width = grid::unit(1.25, "cm"),
  legend.spacing.x = grid::unit(0.12, "cm"),
  legend.text = ggplot2::element_text(size = 8.4),
  panel.grid.minor = ggplot2::element_blank(),
  panel.grid.major = ggplot2::element_line(
    colour = "grey91", linewidth = 0.25
  ),
  plot.title.position = "plot",
  plot.title = ggplot2::element_text(face = "bold", size = 13),
  plot.subtitle = ggplot2::element_text(size = 9.2),
  plot.margin = ggplot2::margin(8, 10, 4, 8)
)

coefficient_base_plot <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect", xmin = 0.30, xmax = 0.70, ymin = -Inf, ymax = Inf,
    fill = "#EFEFEF", colour = NA
  ) +
  ggplot2::geom_ribbon(
    data = band_data,
    ggplot2::aes(
      x = t, ymin = q05, ymax = q95, fill = label, group = label
    ),
    alpha = 0.26, colour = NA
  ) +
  ggplot2::geom_hline(
    yintercept = 0, colour = "grey45", linetype = "dashed",
    linewidth = 0.42
  ) +
  ggplot2::geom_line(
    data = coefficient_data,
    ggplot2::aes(
      x = t, y = mean, colour = label, linetype = label, group = label
    ),
    linewidth = 0.86, na.rm = TRUE
  ) +
  ggplot2::geom_line(
    data = truth, ggplot2::aes(x = t, y = truth),
    colour = "black", linetype = "longdash", linewidth = 0.95
  ) +
  ggplot2::annotate(
    "rect", xmin = 0.30, xmax = 0.70, ymin = -0.25, ymax = 0.25,
    fill = "#F7D8D6", colour = "#B2182B", alpha = 0.26,
    linetype = "longdash", linewidth = 0.95
  ) +
  ggplot2::scale_colour_manual(
    values = coefficient_colours, drop = FALSE
  ) +
  ggplot2::scale_linetype_manual(
    values = coefficient_linetypes, drop = FALSE
  ) +
  ggplot2::scale_fill_manual(values = band_colours, guide = "none") +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.2), limits = c(0, 1),
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::coord_cartesian(ylim = c(-2.7, 2.7), clip = "on") +
  ggplot2::labs(
    title = "Case II coefficient recovery",
    subtitle = paste(
      "Production benchmark (100 replications); grey marks the true-null",
      "interval and light bands are pointwise 90% Monte Carlo envelopes"
    ),
    x = "Functional argument t",
    y = expression(hat(beta)(t)),
    colour = NULL,
    linetype = NULL
  ) +
  theme_coefficient +
  ggplot2::guides(
    colour = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  )

inset_subtitle <- paste0(
  "t=[0.30,0.70]; beta(t)=[-0.25,0.25]; ",
  "exact zeros: ARFC ",
  sprintf("%.1f", 100 * rates[["arfc_select"]]),
  "%, fMCP ", sprintf("%.1f", 100 * rates[["fixed_fmcp"]]),
  "%, fSCAD ", sprintf("%.1f", 100 * rates[["fixed_fscad"]]),
  "%, FLiRTI ", sprintf("%.1f", 100 * rates[["flirti"]]), "%"
)

inset_plot <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect", xmin = 0.30, xmax = 0.70, ymin = -Inf, ymax = Inf,
    fill = "#F5F5F5", colour = NA
  ) +
  ggplot2::geom_ribbon(
    data = band_data,
    ggplot2::aes(
      x = t, ymin = q05, ymax = q95, fill = label, group = label
    ),
    alpha = 0.28, colour = NA
  ) +
  ggplot2::geom_hline(
    yintercept = 0, colour = "black", linetype = "longdash",
    linewidth = 0.50
  ) +
  ggplot2::geom_line(
    data = coefficient_data,
    ggplot2::aes(
      x = t, y = mean, colour = label, linetype = label, group = label
    ),
    linewidth = 0.74, na.rm = TRUE
  ) +
  ggplot2::geom_line(
    data = truth, ggplot2::aes(x = t, y = truth),
    colour = "black", linetype = "longdash", linewidth = 0.66
  ) +
  ggplot2::scale_colour_manual(
    values = coefficient_colours, guide = "none"
  ) +
  ggplot2::scale_linetype_manual(
    values = coefficient_linetypes, guide = "none"
  ) +
  ggplot2::scale_fill_manual(values = band_colours, guide = "none") +
  ggplot2::scale_x_continuous(breaks = c(0.30, 0.50, 0.70)) +
  ggplot2::scale_y_continuous(breaks = c(-0.2, 0, 0.2)) +
  ggplot2::coord_cartesian(
    xlim = c(0.30, 0.70), ylim = c(-0.25, 0.25),
    expand = FALSE, clip = "on"
  ) +
  ggplot2::labs(
    title = "True-null interval zoom",
    subtitle = inset_subtitle,
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_bw(base_size = 7.8, base_family = "serif") +
  ggplot2::theme(
    legend.position = "none",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      colour = "grey90", linewidth = 0.2
    ),
    panel.border = ggplot2::element_rect(
      colour = "#B2182B", linewidth = 0.82
    ),
    plot.background = ggplot2::element_rect(
      fill = "white", colour = "#B2182B", linewidth = 0.82
    ),
    plot.title = ggplot2::element_text(
      face = "bold", size = 8.8, colour = "#B2182B"
    ),
    plot.subtitle = ggplot2::element_text(size = 6.4, lineheight = 0.98),
    axis.text = ggplot2::element_text(size = 7.0),
    plot.margin = ggplot2::margin(5, 6, 5, 6)
  )

# Data-coordinate leaders start exactly at the source rectangle's lower
# corners. The inset uses the same coordinate system, so the joins remain
# aligned after resizing or manuscript inclusion. Zero expansion anchors
# the inset flush to the panel's top-right corner.
inset_grob <- ggplot2::ggplotGrob(inset_plot)
coefficient_plot <- coefficient_base_plot +
  ggplot2::annotate(
    "segment", x = 0.30, y = -0.25, xend = 0.44, yend = 0.10,
    colour = "#B2182B", linewidth = 0.82, linetype = "longdash"
  ) +
  ggplot2::annotate(
    "segment", x = 0.70, y = -0.25, xend = 1.00, yend = 0.10,
    colour = "#B2182B", linewidth = 0.82, linetype = "longdash"
  ) +
  ggplot2::annotation_custom(
    grob = inset_grob,
    xmin = 0.44, xmax = 1.00, ymin = 0.10, ymax = 2.70
  )

zero_ids <- c(
  "arfc_select", "asr_mcp", "amr_mcp", "asr_scad", "amr_scad",
  "fixed_fmcp", "fixed_fscad", "flirti", "oracle"
)
zero_labels <- unname(analysis$config$method_labels[zero_ids])
zero_data <- case_curves[
  case_curves$method %in% zero_ids &
    is.finite(case_curves$zero_frequency), , drop = FALSE
]
zero_data$panel <- ifelse(
  zero_data$method %in%
    c("arfc_select", "asr_mcp", "amr_mcp", "asr_scad", "amr_scad"),
  "ARFC selector and candidates",
  "Fixed and literature comparators"
)
zero_data$panel <- factor(
  zero_data$panel,
  levels = c(
    "ARFC selector and candidates",
    "Fixed and literature comparators"
  )
)
zero_data$label <- factor(
  zero_data$label,
  levels = zero_labels,
  ordered = TRUE
)
if (!setequal(unique(as.character(zero_data$label)), zero_labels)) {
  stop("Structural zero-selection series are incomplete.", call. = FALSE)
}

zero_colours <- c(
  "ARFC--Select" = "#000000", "ASR--MCP" = "#0072B2",
  "AMR--MCP" = "#56B4E9", "ASR--SCAD" = "#D55E00",
  "AMR--SCAD" = "#E69F00", "fMCP" = "#6A3D9A",
  "fSCAD" = "#009E73", "FLiRTI" = "#CC79A7",
  "Oracle" = "#666666"
)
zero_linetypes <- setNames(
  rep(
    c("solid", "longdash", "dotdash", "dashed", "twodash"),
    length.out = length(zero_labels)
  ),
  zero_labels
)

zero_plot <- ggplot2::ggplot(zero_data) +
  ggplot2::annotate(
    "rect", xmin = 0.30, xmax = 0.70, ymin = -Inf, ymax = Inf,
    fill = "#EFEFEF", colour = NA
  ) +
  ggplot2::geom_vline(
    xintercept = c(0.30, 0.70), colour = "grey45",
    linetype = "dotted", linewidth = 0.35
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      x = t, y = zero_frequency, colour = label, linetype = label
    ),
    linewidth = 0.82, na.rm = TRUE
  ) +
  ggplot2::facet_wrap(~panel, ncol = 1) +
  ggplot2::scale_colour_manual(values = zero_colours, drop = FALSE) +
  ggplot2::scale_linetype_manual(values = zero_linetypes, drop = FALSE) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.25), limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::labs(
    title = "Case II pointwise exact-zero selection",
    subtitle = paste(
      "Production benchmark; grey marks the true-null interval; nonsparse",
      "Smooth, OLS, and PCR are omitted"
    ),
    x = "Functional argument t",
    y = "Selection frequency",
    colour = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_bw(base_size = 10.5, base_family = "serif") +
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = grid::unit(0.45, "cm"),
    legend.spacing.x = grid::unit(0.015, "cm"),
    legend.spacing.y = grid::unit(0, "cm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.text = ggplot2::element_text(size = 8.4),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      colour = "grey91", linewidth = 0.25
    ),
    strip.background = ggplot2::element_rect(fill = "grey94"),
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    plot.title.position = "plot",
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    plot.subtitle = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(7, 8, 2, 7)
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(linewidth = 0.72)
    ),
    linetype = ggplot2::guide_legend(
      nrow = 1, byrow = TRUE,
      override.aes = list(linewidth = 0.72)
    )
  )

dir.create(generated_parent, recursive = TRUE, showWarnings = FALSE)
dir.create(pdf_parent, recursive = TRUE, showWarnings = FALSE)
reporting_staging <- file.path(
  generated_parent,
  paste0(".", basename(reporting_dir), "_staging_", Sys.getpid())
)
pdf_staging <- file.path(
  pdf_parent,
  paste0(".", basename(pdf_dir), "_staging_", Sys.getpid())
)
if (dir.exists(reporting_staging) || dir.exists(pdf_staging)) {
  stop("Figure staging directory already exists.", call. = FALSE)
}
dir.create(file.path(reporting_staging, "figures"), recursive = TRUE)
dir.create(pdf_staging, recursive = TRUE)
on.exit({
  if (dir.exists(reporting_staging)) {
    unlink(reporting_staging, recursive = TRUE)
  }
  if (dir.exists(pdf_staging)) {
    unlink(pdf_staging, recursive = TRUE)
  }
}, add = TRUE)

save_formats <- function(plot, base, width, height) {
  paths <- file.path(
    reporting_staging,
    "figures",
    paste0(base, c(".eps", ".pdf", ".png"))
  )
  ggplot2::ggsave(
    paths[2L], plot, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf, bg = "white"
  )
  ggplot2::ggsave(
    paths[3L], plot, width = width, height = height, units = "in",
    device = "png", dpi = 320, bg = "white"
  )
  grDevices::cairo_ps(
    filename = paths[1L], width = width, height = height,
    onefile = FALSE, family = "serif", bg = "white",
    fallback_resolution = 600
  )
  print(plot)
  grDevices::dev.off()
  paths
}

coefficient_base <- "case_ii_coefficient_recovery_v8"
zero_base <- "case_ii_zero_selection_v8"
coefficient_paths <- save_formats(
  coefficient_plot, coefficient_base, width = 9.6, height = 5.55
)
zero_paths <- save_formats(
  zero_plot, zero_base, width = 8.4, height = 5.85
)

figure_root <- file.path(
  "generated", prefix, "production_figure_refinement_v8", "figures"
)
figure_fragment <- c(
  "\\begin{figure}[!htbp]",
  "\\centering",
  paste0(
    "\\includegraphics[width=0.98\\linewidth]{",
    file.path(figure_root, paste0(coefficient_base, ".eps")),
    "}"
  ),
  paste0(
    "\\caption{Production Case-II coefficient recovery over 100 paired ",
    "replications. The dashed black curve is the truth, the grey field is ",
    "the complete true-null interval, and light bands are pointwise 90\\% ",
    "Monte Carlo envelopes for ARFC--Select, fMCP, and fSCAD. The ",
    "dark-red dashed rectangle and aligned leaders identify the source ",
    "of the top-right-aligned enlarged inset, ",
    "which displays $t\\in[0.30,0.70]$ and ",
    "$\\beta(t)\\in[-0.25,0.25]$. Exact-zero rates summarize discrete ",
    "structural selection and cannot be inferred from nearly overlapping ",
    "mean curves; Oracle is omitted from the callout because it is ",
    "truth-aware.}"
  ),
  "\\label{fig:arfc-legacy-coefficient}",
  "\\end{figure}",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  paste0(
    "\\includegraphics[width=0.98\\linewidth]{",
    file.path(figure_root, paste0(zero_base, ".eps")),
    "}"
  ),
  paste0(
    "\\caption{Production Case-II pointwise exact-zero selection ",
    "frequencies for structurally sparse estimators. A single-row legend ",
    "retains the full names of all nine structural methods. Smooth, OLS, ",
    "and PCR do not produce structural local zeros and are omitted.}"
  ),
  "\\label{fig:arfc-legacy-zero}",
  "\\end{figure}"
)
writeLines(
  figure_fragment,
  file.path(reporting_staging, "figure_includes_v8.tex"),
  useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(
    reporting_staging,
    paste0(prefix, "_figure_refinement_sessionInfo_v8.txt")
  )
)

all_paths <- c(coefficient_paths, zero_paths)
png_paths <- all_paths[grepl("[.]png$", all_paths)]
pdf_paths <- all_paths[grepl("[.]pdf$", all_paths)]
eps_paths <- all_paths[grepl("[.]eps$", all_paths)]
png_checks <- vapply(png_paths, function(path) {
  image <- png::readPNG(path)
  length(dim(image)) >= 2L && all(dim(image)[1:2] > 1000L) &&
    is.finite(stats::sd(as.numeric(image))) &&
    stats::sd(as.numeric(image)) > 0
}, logical(1L))
pdfinfo <- Sys.which("pdfinfo")
if (!nzchar(pdfinfo)) {
  stop("`pdfinfo` is required for figure validation.", call. = FALSE)
}
pdf_checks <- vapply(pdf_paths, function(path) {
  output <- suppressWarnings(system2(
    pdfinfo, path, stdout = TRUE, stderr = TRUE
  ))
  any(grepl("^Pages:[[:space:]]+1$", output))
}, logical(1L))
eps_checks <- vapply(eps_paths, function(path) {
  header <- readLines(path, n = 30L, warn = FALSE)
  any(grepl("^%!PS-Adobe-.*EPSF", header)) &&
    any(grepl("^%%(HiRes)?BoundingBox:", header))
}, logical(1L))
if (!all(png_checks) || !all(pdf_checks) || !all(eps_checks)) {
  stop("PNG, PDF, or EPS readability gate failed.", call. = FALSE)
}

relative_files <- list.files(
  reporting_staging, recursive = TRUE, full.names = FALSE
)
payload_paths <- file.path(reporting_staging, relative_files)
manifest <- data.frame(
  scientific_signature = analysis$scientific_signature,
  parent_analysis_md5 = unname(as.character(tools::md5sum(analysis_file))),
  reporting_script_md5 = unname(as.character(tools::md5sum(script_file))),
  role = ifelse(
    grepl("[.]eps$", relative_files), "figure_eps",
    ifelse(
      grepl("[.]pdf$", relative_files), "figure_pdf",
      ifelse(
        grepl("[.]png$", relative_files), "figure_png",
        ifelse(grepl("[.]tex$", relative_files), "latex", "other")
      )
    )
  ),
  file = relative_files,
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest,
  file.path(
    reporting_staging,
    paste0(prefix, "_figure_refinement_manifest_v8.csv")
  ),
  row.names = FALSE
)

for (path in all_paths) {
  destination <- file.path(pdf_staging, basename(path))
  if (!file.copy(path, destination, overwrite = FALSE, copy.mode = TRUE)) {
    stop("Figure mirror staging failed: ", basename(path), call. = FALSE)
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
  file.path(
    pdf_staging,
    paste0(prefix, "_figure_refinement_pdf_manifest_v8.csv")
  ),
  row.names = FALSE
)

if (!file.rename(reporting_staging, reporting_dir) ||
    !file.rename(pdf_staging, pdf_dir)) {
  stop("Atomic v8 figure-bundle install failed.", call. = FALSE)
}

cat("ARFC Case-II v8 figure refinement complete.\n")
cat("Scientific signature:", analysis$scientific_signature, "\n")
cat("Coefficient inset: flush top-right with data-aligned leaders.\n")
cat("Zero-selection legend: nine full labels in one row.\n")
cat("Reporting bundle:", reporting_dir, "\n")
