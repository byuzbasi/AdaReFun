script_dir <- local({
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argument)) {
    dirname(normalizePath(sub("^--file=", "", argument[1L])))
  } else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "neospectra_v3.R"))

mode <- real_mode()
config <- neospectra_v3_config(mode, project_root)
analysis_path <- neospectra_v3_result_path(config, "analysis", "rds")
nested_path <- neospectra_v3_result_path(config, "nested_cv", "rds")
simulation_path <- neospectra_v3_result_path(
  config, "real_design_simulation", "rds"
)
validation_path <- neospectra_v3_result_path(
  config, "validation_checks", "csv"
)
validation_report_path <- neospectra_v3_result_path(config, "validation", "rds")
checkpoint_paths <- c(
  nested_cv = file.path(
    config$results_dir, "checkpoints",
    neospectra_v3_name("nested_cv_checkpoint", config, "rds")
  ),
  real_design_simulation = file.path(
    config$results_dir, "checkpoints",
    neospectra_v3_name("real_design_simulation_checkpoint", config, "rds")
  ),
  stability = file.path(
    config$results_dir, "checkpoints",
    neospectra_v3_name("stability_checkpoint", config, "rds")
  )
)
required_inputs <- c(
  analysis_path, nested_path, simulation_path, validation_path,
  validation_report_path, checkpoint_paths
)
if (any(!file.exists(required_inputs))) {
  stop(
    "NeoSpectra v3 manuscript inputs are missing: ",
    paste(basename(required_inputs[!file.exists(required_inputs)]),
          collapse = ", "), call. = FALSE
  )
}
analysis <- readRDS(analysis_path)
nested <- readRDS(nested_path)
simulation <- readRDS(simulation_path)
validation_report <- readRDS(validation_report_path)
validation <- utils::read.csv(
  validation_path, stringsAsFactors = FALSE, check.names = FALSE
)
processed <- load_neospectra_v3_processed(config)
if (!nrow(validation) || !all(validation$passed)) {
  stop("NeoSpectra v3 validation must pass before manuscript output.",
       call. = FALSE)
}
configuration_keys <- c(
  "version", "mode", "master_seed", "methods", "primary_sparse_methods",
  "M", "degree", "penalty_order", "lambda", "kappa", "mcp_gamma",
  "mcp_sensitivity", "scad_a", "max_iter", "tolerance", "norm_floor",
  "relative_cutoff_multiplier", "zero_tolerance", "integration",
  "pcr_max_components",
  "outer_repeats", "outer_folds", "inner_folds", "bootstrap_replicates",
  "simulation_replicates", "simulation_n", "simulation_snr",
  "simulation_snr_definition", "simulation_scenarios"
)
configuration_matches <-
  identical(analysis$version, config$version) &&
  identical(analysis$mode, mode) &&
  identical(nested$version, config$version) &&
  identical(nested$mode, mode) &&
  identical(simulation$version, config$version) &&
  identical(simulation$mode, mode) &&
  identical(validation_report$version, config$version) &&
  identical(validation_report$mode, mode) &&
  all(vapply(
    list(
      analysis$configuration, nested$configuration,
      simulation$configuration, validation_report$configuration
    ),
    function(value) identical(
      value[configuration_keys], config[configuration_keys]
    ), FUN.VALUE = logical(1L)
  ))
if (!configuration_matches) {
  stop("NeoSpectra analysis and manuscript configuration do not match.",
       call. = FALSE)
}
checkpoint_components <- names(checkpoint_paths)
checkpoint_signatures_match <- all(vapply(
  checkpoint_components, function(component) {
    checkpoint <- readRDS(checkpoint_paths[[component]])
    identical(
      checkpoint$signature,
      neospectra_v3_checkpoint_signature(config, processed, component)
    )
  }, FUN.VALUE = logical(1L)
))
if (!checkpoint_signatures_match) {
  stop("NeoSpectra checkpoints do not match the current code/data signature.",
       call. = FALSE)
}
if (!identical(validation_report$checks, validation)) {
  stop("NeoSpectra validation CSV and report object do not match.",
       call. = FALSE)
}
provenance_paths <- normalizePath(c(
  analysis_path, nested_path, simulation_path,
  neospectra_v3_processed_path(config)
))
inventory_index <- match(provenance_paths, validation_report$inventory$path)
provenance_matches <- all(!is.na(inventory_index)) && identical(
  unname(tools::md5sum(provenance_paths)),
  validation_report$inventory$md5[inventory_index]
)
if (!provenance_matches) {
  stop("Validated NeoSpectra inputs changed after validation.", call. = FALSE)
}

output_suffix <- if (identical(mode, "smoke")) {
  "real_data_v3_smoke"
} else "real_data_v3"
generated_dir <- file.path(
  project_root, "config", "generated", output_suffix
)
paths <- list(
  prediction = file.path(generated_dir, "table_neospectra_prediction_v3.tex"),
  semisim = file.path(generated_dir, "table_neospectra_semisim_v3.tex"),
  full = file.path(generated_dir, "table_neospectra_full_fit_v3.tex"),
  gamma = file.path(generated_dir, "table_neospectra_gamma_v3.tex"),
  results = file.path(generated_dir, "neospectra_results_v3.tex"),
  figures = file.path(generated_dir, "neospectra_figures_v3.tex"),
  session = file.path(generated_dir, "sessionInfo_neospectra_outputs_v3.txt"),
  prediction_bootstrap = neospectra_v3_result_path(
    config, "prediction_cluster_bootstrap", "csv"
  ),
  prediction_bootstrap_summary = neospectra_v3_result_path(
    config, "prediction_cluster_bootstrap_summary", "csv"
  )
)
ensure_new_output(unlist(paths), overwrite = real_overwrite())
dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)

fmt <- function(value, digits = 3L) {
  ifelse(
    is.finite(value), formatC(value, digits = digits, format = "fg",
                              flag = "#"), "--"
  )
}
fmt_sci <- function(value, digits = 3L) {
  ifelse(
    is.finite(value), formatC(value, digits = digits, format = "e"), "--"
  )
}
paired_ci <- function(value, confidence = 0.95) {
  value <- value[is.finite(value)]
  n <- length(value)
  estimate <- if (n) mean(value) else NA_real_
  se <- if (n > 1L) stats::sd(value) / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt((1 + confidence) / 2, n - 1L) else
    NA_real_
  c(n = n, estimate = estimate, lower = estimate - critical * se,
    upper = estimate + critical * se, se = se)
}
write_tabular <- function(data, path, align) {
  escaped <- data
  escaped[] <- lapply(escaped, function(column) {
    latex_escape(as.character(column))
  })
  header <- paste(latex_escape(names(escaped)), collapse = " & ")
  body <- apply(escaped, 1L, function(row) {
    paste0(paste(row, collapse = " & "), " \\\\")
  })
  atomic_write_lines(c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule", paste0(header, " \\\\"), "\\midrule", body,
    "\\bottomrule", "\\end{tabular}"
  ), path)
}

prediction <- analysis$prediction$aggregate
prediction <- prediction[prediction$evaluation == "outer_cv", , drop = FALSE]
prediction <- prediction[match(config$methods, prediction$method), , drop = FALSE]
prediction_table <- data.frame(
  Method = prediction$method,
  RMSE = fmt(prediction$RMSE, 4L),
  MAE = fmt(prediction$MAE, 4L),
  `R-squared` = fmt(prediction$R2, 4L),
  CCC = fmt(prediction$CCC, 4L), check.names = FALSE
)
write_tabular(prediction_table, paths$prediction, "lrrrr")

outer_predictions <- nested$predictions[
  nested$predictions$evaluation == "outer_cv" &
    nested$predictions$method %in% config$primary_sparse_methods,
  , drop = FALSE
]
fmcp_predictions <- outer_predictions[
  outer_predictions$method == "FMCP",
  c("repeat_index", "sample_id", "observed", "predicted")
]
fscad_predictions <- outer_predictions[
  outer_predictions$method == "fSCAD",
  c("repeat_index", "sample_id", "observed", "predicted")
]
paired_predictions <- merge(
  fmcp_predictions, fscad_predictions,
  by = c("repeat_index", "sample_id", "observed"),
  suffixes = c("_FMCP", "_fSCAD"), all = TRUE
)
expected_prediction_rows <- nrow(analysis$prediction$aggregate[
  analysis$prediction$aggregate$evaluation == "outer_cv" &
    analysis$prediction$aggregate$method == "FMCP", , drop = FALSE
])
if (nrow(paired_predictions) !=
    length(unique(outer_predictions$sample_id)) * config$outer_repeats ||
    any(!is.finite(paired_predictions$predicted_FMCP)) ||
    any(!is.finite(paired_predictions$predicted_fSCAD)) ||
    expected_prediction_rows != 1L) {
  stop("Sparse-method outer predictions are not completely paired.",
       call. = FALSE)
}
paired_predictions$squared_error_FMCP <-
  (paired_predictions$observed - paired_predictions$predicted_FMCP)^2
paired_predictions$squared_error_fSCAD <-
  (paired_predictions$observed - paired_predictions$predicted_fSCAD)^2
subject_error <- stats::aggregate(
  cbind(squared_error_FMCP, squared_error_fSCAD) ~ sample_id,
  paired_predictions, mean
)
bootstrap_replicates <- 2000L
set.seed(config$master_seed + 700001L)
bootstrap_difference <- vapply(seq_len(bootstrap_replicates), function(index) {
  sampled <- sample.int(nrow(subject_error), nrow(subject_error), replace = TRUE)
  sqrt(mean(subject_error$squared_error_fSCAD[sampled])) -
    sqrt(mean(subject_error$squared_error_FMCP[sampled]))
}, FUN.VALUE = numeric(1L))
prediction_difference <-
  sqrt(mean(subject_error$squared_error_fSCAD)) -
  sqrt(mean(subject_error$squared_error_FMCP))
prediction_ci <- unname(stats::quantile(
  bootstrap_difference, c(0.025, 0.975), type = 8L
))
prediction_bootstrap <- data.frame(
  replication = seq_len(bootstrap_replicates),
  rmse_fscad_minus_fmcp = bootstrap_difference,
  seed = config$master_seed + 700001L,
  resampling_unit = "soil sample ID; both outer repeats kept together",
  stringsAsFactors = FALSE
)
prediction_bootstrap_summary <- data.frame(
  contrast = "RMSE(fSCAD)-RMSE(FMCP)", estimate = prediction_difference,
  lower_95 = prediction_ci[1L], upper_95 = prediction_ci[2L],
  bootstrap_replicates = bootstrap_replicates,
  unique_samples = nrow(subject_error), seed = config$master_seed + 700001L,
  stringsAsFactors = FALSE
)
atomic_write_csv(prediction_bootstrap, paths$prediction_bootstrap)
atomic_write_csv(
  prediction_bootstrap_summary, paths$prediction_bootstrap_summary
)

paired <- analysis$simulation_paired
expected_cells <- expand.grid(
  scenario = config$simulation_scenarios,
  n_train = config$simulation_n,
  snr_label = as.character(config$simulation_snr),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
expected_cells <- expected_cells[
  expected_cells$scenario != "all_null" |
    expected_cells$snr_label == as.character(config$simulation_snr[1L]),
  , drop = FALSE
]
expected_cells$snr_label[expected_cells$scenario == "all_null"] <- "all-null"
cell_id <- function(data) paste(
  data$scenario, data$n_train, data$snr_label, sep = "\r"
)
actual_counts <- table(cell_id(paired))
expected_ids <- cell_id(expected_cells)
if (!setequal(names(actual_counts), expected_ids) ||
    any(actual_counts != config$simulation_replicates) ||
    nrow(paired) != nrow(expected_cells) * config$simulation_replicates ||
    anyDuplicated(paired$task_id)) {
  stop("Real-design output does not contain every prespecified paired task.",
       call. = FALSE)
}
group_key <- interaction(
  paired[c("scenario", "n_train", "snr_label")],
  drop = TRUE, lex.order = TRUE
)
semisim <- bind_rows_base(lapply(split(paired, group_key), function(rows) {
  usable <- rows$success_FMCP %in% TRUE & rows$success_fSCAD %in% TRUE
  successful <- rows[usable, , drop = FALSE]
  ise <- paired_ci(successful$ise_advantage_FMCP)
  null <- paired_ci(successful$null_recovery_advantage_FMCP)
  false_null <- paired_ci(successful$false_null_advantage_FMCP)
  runtime <- successful$runtime_ratio_FMCP_to_fSCAD
  runtime <- runtime[is.finite(runtime)]
  data.frame(
    Scenario = rows$scenario[1L], n = rows$n_train[1L],
    SNR = rows$snr_label[1L],
    Pairs = paste0(nrow(successful), "/", nrow(rows)),
    `ISE fSCAD-FMCP (95% CI)` = paste0(
      fmt_sci(ise["estimate"], 2L), " [", fmt_sci(ise["lower"], 2L),
      ", ", fmt_sci(ise["upper"], 2L), "]"
    ),
    `FMCP ISE win (%)` = fmt(
      100 * mean(successful$ise_advantage_FMCP > 0), 3L
    ),
    `Null recovery FMCP-fSCAD` = fmt(null["estimate"], 3L),
    `False-null fSCAD-FMCP` = fmt(false_null["estimate"], 3L),
    `Median runtime ratio` = fmt(
      if (length(runtime)) stats::median(runtime) else NA_real_, 3L
    ), check.names = FALSE
  )
}))
semisim <- semisim[order(
  match(semisim$Scenario, config$simulation_scenarios),
  semisim$n, semisim$SNR
), , drop = FALSE]
if (nrow(semisim) != nrow(expected_cells)) {
  stop("Real-design table lost a prespecified cell.", call. = FALSE)
}
write_tabular(semisim, paths$semisim, "lrrrlrrrr")

full_selection <- analysis$full$selection
full_curves <- analysis$full$curves
stability <- analysis$stability$summary
stability <- stability[stability$source == "bootstrap", , drop = FALSE]
full_table <- bind_rows_base(lapply(config$primary_sparse_methods, function(method) {
  selected <- full_selection[full_selection$method == method, , drop = FALSE]
  curve <- full_curves[full_curves$method == method, , drop = FALSE]
  stable <- stability[stability$method == method, , drop = FALSE]
  successful_range <- range(stable$successful_fits, na.rm = TRUE)
  successful_label <- if (successful_range[1L] == successful_range[2L]) {
    as.character(successful_range[1L])
  } else paste(successful_range, collapse = "--")
  data.frame(
    Method = method, M = selected$M, Lambda = fmt_sci(selected$lambda, 2L),
    Kappa = fmt_sci(selected$kappa, 2L),
    `Full zero fraction` = fmt(mean(curve$estimated_null), 3L),
    `Stable-zero fraction among successful fits` = fmt(
      mean(stable$zero_selection_frequency >= 0.8), 3L
    ),
    `Successful bootstraps` = successful_label,
    `Elapsed seconds` = fmt(selected$elapsed_seconds, 3L),
    check.names = FALSE
  )
}))
write_tabular(full_table, paths$full, "lrrrrrrr")

gamma_selection <- analysis$mcp_gamma_sensitivity$selection
gamma_curves <- analysis$mcp_gamma_sensitivity$curves
gamma_table <- bind_rows_base(lapply(gamma_selection$mcp_gamma, function(gamma) {
  selected <- gamma_selection[gamma_selection$mcp_gamma == gamma, , drop = FALSE]
  curve <- gamma_curves[gamma_curves$mcp_gamma == gamma, , drop = FALSE]
  data.frame(
    Gamma = fmt(gamma, 2L), Primary = ifelse(selected$primary, "yes", "no"),
    `Fixed candidate` = selected$candidate_id,
    `Zero fraction` = fmt(mean(curve$estimated_null), 3L),
    `Elapsed seconds` = fmt(selected$elapsed_seconds, 3L),
    check.names = FALSE
  )
}))
write_tabular(gamma_table, paths$gamma, "rllrr")

cell_summaries <- bind_rows_base(lapply(split(paired, group_key), function(rows) {
  successful <- rows[
    rows$success_FMCP %in% TRUE & rows$success_fSCAD %in% TRUE,
    , drop = FALSE
  ]
  ise <- paired_ci(successful$ise_advantage_FMCP)
  data.frame(
    positive = is.finite(ise["lower"]) && ise["lower"] > 0,
    negative = is.finite(ise["upper"]) && ise["upper"] < 0,
    successful_pairs = nrow(successful), planned_pairs = nrow(rows)
  )
}))
fmcp_prediction <- prediction[prediction$method == "FMCP", , drop = FALSE]
fscad_prediction <- prediction[prediction$method == "fSCAD", , drop = FALSE]
fmcp_full <- full_table[full_table$Method == "FMCP", , drop = FALSE]
fscad_full <- full_table[full_table$Method == "fSCAD", , drop = FALSE]
atomic_write_lines(c(
  paste0(
    "In grouped repeated nested cross-validation, pooled RMSE was ",
    fmt(fmcp_prediction$RMSE, 4L), " for FMCP and ",
    fmt(fscad_prediction$RMSE, 4L), " for fSCAD.  The pooled ",
    "fSCAD-minus-FMCP RMSE difference was ",
    fmt_sci(prediction_difference, 3L), " (95\\% soil-sample cluster-",
    "bootstrap interval ", fmt_sci(prediction_ci[1L], 3L), " to ",
    fmt_sci(prediction_ci[2L], 3L), "; both outer repeats were kept ",
    "together)."
  ),
  paste0(
    "Across all ", nrow(cell_summaries),
    " prespecified real-design cells, the paired 95\\% interval for ",
    "ISE(fSCAD)-ISE(FMCP) lay wholly above zero in ",
    sum(cell_summaries$positive), " cells and wholly below zero in ",
    sum(cell_summaries$negative),
    " cells.  There were ", sum(cell_summaries$successful_pairs), " successful ",
    "pairs among ", sum(cell_summaries$planned_pairs),
    " planned pairs; every cell remains in ",
    "Table~\\ref{tab:neospectra-semisim}."
  ),
  paste0(
    "The full-data estimated-zero fractions were ",
    fmcp_full$`Full zero fraction`, " for FMCP and ",
    fscad_full$`Full zero fraction`,
    " for fSCAD.  The corresponding fractions selected as zero in at least ",
    "80\\% of successful bootstrap fits, conditional on fixed CV tuning, ",
    "were ", fmcp_full$`Stable-zero fraction among successful fits`, " and ",
    fscad_full$`Stable-zero fraction among successful fits`, ", respectively ",
    "(", fmcp_full$`Successful bootstraps`, " and ",
    fscad_full$`Successful bootstraps`, " successful fits per wavelength)."
  )
), paths$results)

figure_root <- paste0("../real_data/figures/v3/neospectra/")
figure_version <- if (identical(mode, "smoke")) "_v3_smoke.pdf" else
  "_v3.pdf"
atomic_write_lines(c(
  "\\begin{figure}[p]",
  "\\centering",
  paste0("\\includegraphics[width=0.90\\linewidth]{", figure_root,
         "nested_cv_rmse", figure_version, "}"),
  "\\caption{Grouped repeated nested-CV error on log1p soil organic carbon. FMCP and fSCAD share outer and inner folds and candidate grids; Smooth and PCR are training-only BIC references.}",
  "\\label{fig:neospectra-cv}",
  "\\end{figure}",
  "\\begin{figure}[p]",
  "\\centering",
  paste0("\\includegraphics[width=0.94\\linewidth]{", figure_root,
         "coefficient_bootstrap", figure_version, "}"),
  "\\caption{Full-data NeoSpectra coefficient functions and 80\\% quantile bands among successful sample-bootstrap fits. The bands condition on the globally CV-selected tuning candidate.}",
  "\\label{fig:neospectra-coefficient}",
  "\\end{figure}",
  "\\begin{figure}[p]",
  "\\centering",
  paste0("\\includegraphics[width=0.94\\linewidth]{", figure_root,
         "zero_selection_stability", figure_version, "}"),
  "\\caption{Zero-selection frequency among successful bootstrap fits by wavelength. Tuning is fixed rather than repeated within each bootstrap sample.}",
  "\\label{fig:neospectra-zero}",
  "\\end{figure}",
  "\\begin{figure}[p]",
  "\\centering",
  paste0("\\includegraphics[width=0.98\\linewidth]{", figure_root,
         "simulation_ise_advantage", figure_version, "}"),
  "\\caption{Real-design simulation contrasts. Positive ISE(fSCAD)-ISE(FMCP) values favour FMCP; all prespecified cells and negative controls are retained.}",
  "\\label{fig:neospectra-semisim}",
  "\\end{figure}"
), paths$figures)
capture_session_info(paths$session)
cat("NeoSpectra manuscript tables and LaTeX fragments:", generated_dir, "\n")
