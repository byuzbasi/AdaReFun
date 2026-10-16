script_dir <- local({
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argument)) {
    dirname(normalizePath(sub("^--file=", "", argument[1L])))
  } else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
for (file in c("utils.R", "methods.R", "neospectra_v3.R")) {
  source(file.path(script_dir, "R", file))
}
activate_fmcp_libraries()
assert_packages(c("FMCP", "Rcpp", "RcppArmadillo", "ggplot2"))

mode <- real_mode()
resume <- real_resume()
config <- neospectra_v3_config(mode, project_root)
processed <- load_neospectra_v3_processed(config)
dataset <- processed$sample_dataset
checkpoint_signature <- neospectra_v3_checkpoint_signature(
  config, processed, "stability"
)
input_paths <- list(
  nested = neospectra_v3_result_path(config, "nested_cv", "rds"),
  simulation = neospectra_v3_result_path(
    config, "real_design_simulation", "rds"
  )
)
missing_inputs <- unlist(input_paths)[!file.exists(unlist(input_paths))]
if (length(missing_inputs)) {
  stop(
    "Required NeoSpectra v3 results are missing: ",
    paste(basename(missing_inputs), collapse = ", "), call. = FALSE
  )
}
nested <- readRDS(input_paths$nested)
simulation <- readRDS(input_paths$simulation)

paths <- list(
  result = neospectra_v3_result_path(config, "analysis", "rds"),
  prediction_task = neospectra_v3_result_path(
    config, "prediction_metrics_by_task", "csv"
  ),
  prediction_aggregate = neospectra_v3_result_path(
    config, "prediction_metrics_aggregate", "csv"
  ),
  prediction_paired = neospectra_v3_result_path(
    config, "prediction_fmcp_minus_fscad_paired", "csv"
  ),
  simulation_summary = neospectra_v3_result_path(
    config, "real_design_simulation_summary", "csv"
  ),
  simulation_paired = neospectra_v3_result_path(
    config, "real_design_simulation_paired", "csv"
  ),
  simulation_paired_summary = neospectra_v3_result_path(
    config, "real_design_simulation_paired_summary", "csv"
  ),
  candidates = neospectra_v3_result_path(
    config, "selected_global_candidates", "csv"
  ),
  full_selection = neospectra_v3_result_path(
    config, "full_data_selection", "csv"
  ),
  full_curves = neospectra_v3_result_path(
    config, "full_data_coefficient_curves", "csv"
  ),
  gamma_selection = neospectra_v3_result_path(
    config, "mcp_gamma_sensitivity_selection", "csv"
  ),
  gamma_curves = neospectra_v3_result_path(
    config, "mcp_gamma_sensitivity_curves", "csv"
  ),
  stability_curves = neospectra_v3_result_path(
    config, "stability_curves", "csv"
  ),
  stability_summary = neospectra_v3_result_path(
    config, "stability_summary", "csv"
  ),
  stability_runtime = neospectra_v3_result_path(
    config, "stability_runtime", "csv"
  ),
  stability_draws = neospectra_v3_result_path(
    config, "bootstrap_sample_draws", "csv"
  ),
  runtime_detail = neospectra_v3_result_path(
    config, "runtime_detail", "csv"
  ),
  runtime_summary = neospectra_v3_result_path(
    config, "runtime_summary", "csv"
  ),
  stability_checkpoint = file.path(
    config$results_dir, "checkpoints",
    neospectra_v3_name("stability_checkpoint", config, "rds")
  ),
  figure_coefficients_pdf = neospectra_v3_figure_path(
    config, "coefficient_bootstrap", "pdf"
  ),
  figure_coefficients_png = neospectra_v3_figure_path(
    config, "coefficient_bootstrap", "png"
  ),
  figure_zero_pdf = neospectra_v3_figure_path(
    config, "zero_selection_stability", "pdf"
  ),
  figure_zero_png = neospectra_v3_figure_path(
    config, "zero_selection_stability", "png"
  ),
  figure_cv_pdf = neospectra_v3_figure_path(
    config, "nested_cv_rmse", "pdf"
  ),
  figure_cv_png = neospectra_v3_figure_path(
    config, "nested_cv_rmse", "png"
  ),
  figure_simulation_pdf = neospectra_v3_figure_path(
    config, "simulation_ise_advantage", "pdf"
  ),
  figure_simulation_png = neospectra_v3_figure_path(
    config, "simulation_ise_advantage", "png"
  ),
  session = neospectra_v3_result_path(
    config, "analysis_sessionInfo", "txt"
  )
)
ensure_new_output(
  unlist(paths[setdiff(names(paths), "stability_checkpoint")]),
  resume = resume, overwrite = real_overwrite()
)
if (!resume && file.exists(paths$stability_checkpoint) &&
    !real_overwrite()) {
  stop("A NeoSpectra v3 stability checkpoint exists; use --resume.",
       call. = FALSE)
}

safe_mean_v3 <- function(value) {
  value <- value[is.finite(value)]
  if (length(value)) mean(value) else NA_real_
}
safe_sd_v3 <- function(value) {
  value <- value[is.finite(value)]
  if (length(value) > 1L) stats::sd(value) else NA_real_
}
safe_se_v3 <- function(value) {
  value <- value[is.finite(value)]
  if (length(value) > 1L) stats::sd(value) / sqrt(length(value)) else NA_real_
}
safe_quantile_v3 <- function(value, probability) {
  value <- value[is.finite(value)]
  if (length(value)) {
    unname(stats::quantile(value, probability, names = FALSE, type = 8L))
  } else NA_real_
}

prediction <- summarize_prediction_results_v3(nested$predictions)
fmcp_prediction <- prediction$task[
  prediction$task$method == "FMCP", , drop = FALSE
]
fscad_prediction <- prediction$task[
  prediction$task$method == "fSCAD", , drop = FALSE
]
prediction_paired <- merge(
  fmcp_prediction, fscad_prediction,
  by = c("evaluation", "task_id"), suffixes = c("_FMCP", "_fSCAD"),
  all = TRUE, sort = TRUE
)
prediction_paired$rmse_advantage_FMCP <-
  prediction_paired$RMSE_fSCAD - prediction_paired$RMSE_FMCP
prediction_paired$mae_advantage_FMCP <-
  prediction_paired$MAE_fSCAD - prediction_paired$MAE_FMCP
prediction_paired$r2_advantage_FMCP <-
  prediction_paired$R2_FMCP - prediction_paired$R2_fSCAD
prediction_paired$ccc_advantage_FMCP <-
  prediction_paired$CCC_FMCP - prediction_paired$CCC_fSCAD

simulation_metrics <- simulation$metrics
simulation_metrics$snr_label <- ifelse(
  is.na(simulation_metrics$snr), "all-null",
  format(simulation_metrics$snr, trim = TRUE)
)
simulation_variables <- c(
  "ISE", "true_null_recovery", "false_null_rate", "signal_pmse",
  "noisy_rmse", "empirical_variance_snr", "elapsed_seconds"
)
simulation_groups <- interaction(
  simulation_metrics[c("scenario", "n_train", "snr_label", "method")],
  drop = TRUE, lex.order = TRUE
)
simulation_summary <- bind_rows_base(lapply(
  split(simulation_metrics, simulation_groups), function(rows) {
    output <- data.frame(
      scenario = rows$scenario[1L], n_train = rows$n_train[1L],
      snr_label = rows$snr_label[1L], method = rows$method[1L],
      replications = length(unique(rows$replication)),
      success_rate = mean(rows$success), stringsAsFactors = FALSE
    )
    for (variable in simulation_variables) {
      output[[paste0(variable, "_mean")]] <- safe_mean_v3(rows[[variable]])
      output[[paste0(variable, "_se")]] <- safe_se_v3(rows[[variable]])
    }
    output
  }
))

fmcp_simulation <- simulation_metrics[
  simulation_metrics$method == "FMCP", , drop = FALSE
]
fscad_simulation <- simulation_metrics[
  simulation_metrics$method == "fSCAD", , drop = FALSE
]
simulation_paired <- merge(
  fmcp_simulation, fscad_simulation, by = c(
    "task_id", "task_index", "scenario", "n_train", "n_test", "snr",
    "snr_label", "replication"
  ), suffixes = c("_FMCP", "_fSCAD"), all = TRUE, sort = TRUE
)
simulation_paired$ise_advantage_FMCP <-
  simulation_paired$ISE_fSCAD - simulation_paired$ISE_FMCP
simulation_paired$null_recovery_advantage_FMCP <-
  simulation_paired$true_null_recovery_FMCP -
  simulation_paired$true_null_recovery_fSCAD
simulation_paired$false_null_advantage_FMCP <-
  simulation_paired$false_null_rate_fSCAD -
  simulation_paired$false_null_rate_FMCP
simulation_paired$signal_pmse_advantage_FMCP <-
  simulation_paired$signal_pmse_fSCAD - simulation_paired$signal_pmse_FMCP
simulation_paired$runtime_ratio_FMCP_to_fSCAD <-
  simulation_paired$elapsed_seconds_FMCP /
  simulation_paired$elapsed_seconds_fSCAD
paired_groups <- interaction(
  simulation_paired[c("scenario", "n_train", "snr_label")],
  drop = TRUE, lex.order = TRUE
)
simulation_paired_summary <- bind_rows_base(lapply(
  split(simulation_paired, paired_groups), function(rows) {
    data.frame(
      scenario = rows$scenario[1L], n_train = rows$n_train[1L],
      snr_label = rows$snr_label[1L], pairs = nrow(rows),
      ise_advantage_mean = safe_mean_v3(rows$ise_advantage_FMCP),
      ise_advantage_se = safe_se_v3(rows$ise_advantage_FMCP),
      ise_win_rate_FMCP = safe_mean_v3(rows$ise_advantage_FMCP > 0),
      null_recovery_advantage_mean = safe_mean_v3(
        rows$null_recovery_advantage_FMCP
      ),
      false_null_advantage_mean = safe_mean_v3(
        rows$false_null_advantage_FMCP
      ),
      signal_pmse_advantage_mean = safe_mean_v3(
        rows$signal_pmse_advantage_FMCP
      ),
      runtime_ratio_median = stats::median(
        rows$runtime_ratio_FMCP_to_fSCAD[is.finite(
          rows$runtime_ratio_FMCP_to_fSCAD
        )], na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }
))

initialize_real_method_environment(project_root)
selected_candidates <- select_global_sparse_candidates_v3(
  nested$tuning, config
)
full <- fit_full_selected_v3(dataset, selected_candidates, config)
gamma_sensitivity <- fit_mcp_gamma_sensitivity_v3(
  dataset, selected_candidates, config
)

checkpoint <- if (resume && file.exists(paths$stability_checkpoint)) {
  readRDS(paths$stability_checkpoint)
} else list(completed = list())
if (!is.null(checkpoint$sample_id) &&
    !identical(checkpoint$sample_id, dataset$sample_id)) {
  stop("Stability checkpoint sample IDs do not match processed v3 data.",
       call. = FALSE)
}
if (resume && file.exists(paths$stability_checkpoint) &&
    !identical(checkpoint$signature, checkpoint_signature)) {
  stop("Stability checkpoint configuration/source signature is incompatible.",
       call. = FALSE)
}
if (resume && file.exists(paths$stability_checkpoint) &&
    !identical(checkpoint$selected_candidates, selected_candidates)) {
  stop("Stability checkpoint uses different fixed CV-selected candidates.",
       call. = FALSE)
}
completed <- checkpoint$completed
bootstrap_tasks <- as.list(seq_len(config$bootstrap_replicates))
remaining <- Filter(
  function(index) !as.character(index) %in% names(completed),
  bootstrap_tasks
)
cat(length(bootstrap_tasks), "bootstrap tasks planned;", length(remaining),
    "remaining.\n")
batch_size <- max(1L, 2L * config$workers)
if (length(remaining)) {
  batches <- split(remaining, ceiling(seq_along(remaining) / batch_size))
  for (batch in batches) {
    worker <- function(index) tryCatch(
      run_bootstrap_stability_task_v3(
        index, dataset, selected_candidates, config
      ),
      error = function(condition) list(
        key = as.character(index), fatal_error = conditionMessage(condition)
      )
    )
    batch_results <- run_task_batch(batch, worker, config$workers)
    for (result in batch_results) completed[[result$key]] <- result
    atomic_save_rds(
      list(
        completed = completed, sample_id = dataset$sample_id,
        signature = checkpoint_signature,
        selected_candidates = selected_candidates, configuration = config
      ), paths$stability_checkpoint
    )
    fatal <- vapply(
      batch_results, function(result) !is.null(result$fatal_error),
      FUN.VALUE = logical(1L)
    )
    if (any(fatal)) {
      messages <- vapply(batch_results[fatal], function(result) {
        paste(result$key, result$fatal_error, sep = ": ")
      }, FUN.VALUE = character(1L))
      stop("Bootstrap failure after checkpoint: ",
           paste(messages, collapse = " | "), call. = FALSE)
    }
    cat("Bootstrap progress:", length(completed), "/",
        length(bootstrap_tasks), "tasks\n")
  }
}
bootstrap_keys <- as.character(seq_len(config$bootstrap_replicates))
bootstrap <- completed[bootstrap_keys]
if (any(vapply(bootstrap, is.null, logical(1L)))) {
  stop("Stability checkpoint is incomplete.", call. = FALSE)
}
bootstrap_curves <- bind_rows_base(lapply(bootstrap, `[[`, "curves"))
bootstrap_runtime <- bind_rows_base(lapply(bootstrap, `[[`, "runtime"))
bootstrap_draws <- bind_rows_base(lapply(bootstrap, `[[`, "draws"))
scanner <- if (is.null(processed$scan_dataset)) {
  list(curves = data.frame(), runtime = data.frame())
} else run_scanner_stability_v3(
  processed$scan_dataset, selected_candidates, config
)
stability_curves <- bind_rows_base(list(
  bootstrap_curves, scanner$curves
))
stability_runtime <- bind_rows_base(list(
  bootstrap_runtime, scanner$runtime
))
stability_curves$scanner_label <- ifelse(
  is.na(stability_curves$scanner_id), "pooled", stability_curves$scanner_id
)
valid_stability <- stability_curves[
  stability_curves$success & is.finite(stability_curves$beta) &
    !is.na(stability_curves$estimated_null), , drop = FALSE
]
stability_groups <- interaction(
  valid_stability[c("source", "scanner_label", "method", "wavelength")],
  drop = TRUE, lex.order = TRUE
)
stability_summary <- bind_rows_base(lapply(
  split(valid_stability, stability_groups), function(rows) {
    data.frame(
      source = rows$source[1L], scanner_id = rows$scanner_label[1L],
      method = rows$method[1L], wavelength = rows$wavelength[1L],
      successful_fits = nrow(rows), beta_mean = safe_mean_v3(rows$beta),
      beta_sd = safe_sd_v3(rows$beta),
      beta_q10 = safe_quantile_v3(rows$beta, 0.10),
      beta_q90 = safe_quantile_v3(rows$beta, 0.90),
      zero_selection_frequency = mean(rows$estimated_null),
      stringsAsFactors = FALSE
    )
  }
))

nested_runtime <- nested$selection
nested_runtime_detail <- data.frame(
  source = ifelse(
    nested_runtime$evaluation == "outer_cv", "nested_outer_cv",
    "geographic_external"
  ),
  task_id = nested_runtime$task_id, method = nested_runtime$method,
  elapsed_seconds = ifelse(
    is.finite(nested_runtime$tuning_elapsed_seconds),
    nested_runtime$tuning_elapsed_seconds +
      nested_runtime$final_elapsed_seconds,
    nested_runtime$final_elapsed_seconds
  ),
  success = nested_runtime$fit_success,
  error = nested_runtime$error, stringsAsFactors = FALSE
)
simulation_runtime_detail <- data.frame(
  source = "real_design_simulation", task_id = simulation$runtime$task_id,
  method = simulation$runtime$method,
  elapsed_seconds = simulation$runtime$elapsed_seconds,
  success = simulation$runtime$success, error = simulation$runtime$error,
  stringsAsFactors = FALSE
)
full_runtime_detail <- data.frame(
  source = "full_data", task_id = "full",
  method = full$selection$method,
  elapsed_seconds = full$selection$elapsed_seconds,
  success = full$selection$success, error = full$selection$error,
  stringsAsFactors = FALSE
)
gamma_runtime_detail <- data.frame(
  source = "mcp_gamma_sensitivity",
  task_id = paste0("gamma_", gamma_sensitivity$selection$mcp_gamma),
  method = "FMCP",
  elapsed_seconds = gamma_sensitivity$selection$elapsed_seconds,
  success = gamma_sensitivity$selection$success,
  error = gamma_sensitivity$selection$error,
  stringsAsFactors = FALSE
)
stability_runtime_detail <- data.frame(
  source = stability_runtime$source,
  task_id = ifelse(
    stability_runtime$source == "bootstrap",
    paste0("bootstrap_", stability_runtime$replication),
    paste0("scanner_", stability_runtime$scanner_id)
  ),
  method = stability_runtime$method,
  elapsed_seconds = stability_runtime$elapsed_seconds,
  success = stability_runtime$success, error = stability_runtime$error,
  stringsAsFactors = FALSE
)
runtime_detail <- bind_rows_base(list(
  nested_runtime_detail, simulation_runtime_detail,
  full_runtime_detail, gamma_runtime_detail, stability_runtime_detail
))
runtime_groups <- interaction(
  runtime_detail[c("source", "method")], drop = TRUE, lex.order = TRUE
)
runtime_summary <- bind_rows_base(lapply(
  split(runtime_detail, runtime_groups), function(rows) {
    data.frame(
      source = rows$source[1L], method = rows$method[1L],
      tasks = nrow(rows), success_rate = mean(rows$success),
      elapsed_mean = safe_mean_v3(rows$elapsed_seconds),
      elapsed_median = stats::median(
        rows$elapsed_seconds[is.finite(rows$elapsed_seconds)], na.rm = TRUE
      ),
      elapsed_p90 = safe_quantile_v3(rows$elapsed_seconds, 0.90),
      elapsed_total = sum(rows$elapsed_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))

bootstrap_summary <- stability_summary[
  stability_summary$source == "bootstrap", , drop = FALSE
]
coefficient_plot_data <- merge(
  full$curves,
  bootstrap_summary[c(
    "method", "wavelength", "beta_q10", "beta_q90",
    "zero_selection_frequency"
  )],
  by = c("method", "wavelength"), all.x = TRUE, sort = TRUE
)
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
    title = paste(
      "Full-data coefficient functions with 80% bootstrap bands",
      "conditional on fixed CV tuning"
    )
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "top")

zero_plot <- ggplot2::ggplot(
  bootstrap_summary,
  ggplot2::aes(
    x = wavelength, y = zero_selection_frequency, colour = method
  )
) +
  ggplot2::geom_hline(yintercept = c(0.5, 0.8), linetype = 3,
                      colour = "grey65") +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::scale_colour_manual(values = method_colours) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = "NIR wavelength (nm)",
    y = "Conditional bootstrap zero-selection frequency",
    colour = NULL,
    title = "Local-zero stability conditional on fixed CV tuning"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "top")

cv_plot_data <- prediction$task[
  prediction$task$evaluation == "outer_cv", , drop = FALSE
]
cv_plot <- ggplot2::ggplot(
  cv_plot_data,
  ggplot2::aes(x = method, y = RMSE, colour = method)
) +
  ggplot2::geom_boxplot(outlier.shape = NA, width = 0.62) +
  ggplot2::geom_point(
    ggplot2::aes(group = task_id), position = ggplot2::position_jitter(
      width = 0.08, height = 0, seed = config$master_seed
    ), alpha = 0.65, size = 1.5
  ) +
  ggplot2::labs(
    x = NULL, y = "Held-out RMSE (log1p(SOC))", colour = NULL,
    title = "Grouped nested-CV prediction error by outer fold"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "none")

simulation_plot <- ggplot2::ggplot(
  simulation_paired_summary,
  ggplot2::aes(x = n_train, y = ise_advantage_mean, group = snr_label,
               colour = snr_label)
) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = ise_advantage_mean - 1.96 * ise_advantage_se,
      ymax = ise_advantage_mean + 1.96 * ise_advantage_se
    ), width = 0
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

dir.create(config$figures_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(paths$figure_coefficients_pdf, coefficient_plot,
                width = 7.2, height = 4.6)
ggplot2::ggsave(paths$figure_coefficients_png, coefficient_plot,
                width = 7.2, height = 4.6, dpi = 320)
ggplot2::ggsave(paths$figure_zero_pdf, zero_plot, width = 7.2, height = 4.2)
ggplot2::ggsave(paths$figure_zero_png, zero_plot, width = 7.2, height = 4.2,
                dpi = 320)
ggplot2::ggsave(paths$figure_cv_pdf, cv_plot, width = 6.6, height = 4.4)
ggplot2::ggsave(paths$figure_cv_png, cv_plot, width = 6.6, height = 4.4,
                dpi = 320)
ggplot2::ggsave(paths$figure_simulation_pdf, simulation_plot,
                width = 8.2, height = 5.4)
ggplot2::ggsave(paths$figure_simulation_png, simulation_plot,
                width = 8.2, height = 5.4, dpi = 320)

analysis <- list(
  version = config$version, mode = mode, configuration = config,
  analysis_availability = processed$availability,
  prediction = prediction, prediction_paired = prediction_paired,
  simulation_summary = simulation_summary,
  simulation_paired = simulation_paired,
  simulation_paired_summary = simulation_paired_summary,
  selected_candidates = selected_candidates, full = full,
  mcp_gamma_sensitivity = gamma_sensitivity,
  stability = list(
    curves = stability_curves, summary = stability_summary,
    runtime = stability_runtime, draws = bootstrap_draws
  ),
  runtime_detail = runtime_detail, runtime_summary = runtime_summary
)
atomic_save_rds(analysis, paths$result)
atomic_write_csv(prediction$task, paths$prediction_task)
atomic_write_csv(prediction$aggregate, paths$prediction_aggregate)
atomic_write_csv(prediction_paired, paths$prediction_paired)
atomic_write_csv(simulation_summary, paths$simulation_summary)
atomic_write_csv(simulation_paired, paths$simulation_paired)
atomic_write_csv(simulation_paired_summary, paths$simulation_paired_summary)
atomic_write_csv(selected_candidates, paths$candidates)
atomic_write_csv(full$selection, paths$full_selection)
atomic_write_csv(full$curves, paths$full_curves)
atomic_write_csv(gamma_sensitivity$selection, paths$gamma_selection)
atomic_write_csv(gamma_sensitivity$curves, paths$gamma_curves)
atomic_write_csv(stability_curves, paths$stability_curves)
atomic_write_csv(stability_summary, paths$stability_summary)
atomic_write_csv(stability_runtime, paths$stability_runtime)
atomic_write_csv(bootstrap_draws, paths$stability_draws)
atomic_write_csv(runtime_detail, paths$runtime_detail)
atomic_write_csv(runtime_summary, paths$runtime_summary)
capture_session_info(paths$session)
cat("NeoSpectra v3 analysis, stability summaries, and figures completed.\n")
