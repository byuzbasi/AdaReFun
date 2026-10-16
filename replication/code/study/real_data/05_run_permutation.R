script_dir <- local({
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) dirname(normalizePath(sub("^--file=", "", file_argument[1L]))) else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
for (file in c("utils.R", "config.R", "data.R", "methods.R", "resampling.R")) {
  source(file.path(script_dir, "R", file))
}
activate_fmcp_libraries()
assert_packages(c(
  "FMCP", "PFLR", "Rcpp", "RcppArmadillo", "fda", "glmnet", "flare"
))

mode <- real_mode()
resume <- real_resume()
config <- apply_resolved_grid(real_data_config(mode, project_root))
datasets <- load_processed_datasets(config)
paths <- list(
  result = file.path(config$results_dir,
                     versioned_name("permutation", mode, "rds")),
  statistics = file.path(config$results_dir,
                         versioned_name("permutation_statistics", mode, "csv")),
  summary = file.path(config$results_dir,
                      versioned_name("permutation_summary", mode, "csv")),
  selection = file.path(config$results_dir,
                        versioned_name("permutation_selection", mode, "csv")),
  checkpoint = file.path(config$results_dir, "checkpoints",
                         versioned_name("permutation_checkpoint", mode, "rds")),
  session = file.path(config$results_dir,
                      versioned_name("permutation_sessionInfo", mode, "txt"))
)
ensure_new_output(
  unlist(paths[c("result", "statistics", "summary", "selection", "session")]),
  resume = resume, overwrite = real_overwrite()
)
if (!resume && file.exists(paths$checkpoint) && !real_overwrite()) {
  stop("A permutation checkpoint exists; use --resume.", call. = FALSE)
}

initialize_real_method_environment(project_root)

cat("Fitting observed sparse models for permutation statistics...\n")
observed_rows <- list()
observed_index <- 1L
for (dataset_index in seq_along(datasets)) {
  dataset <- datasets[[dataset_index]]
  fitted <- fit_real_method_set(
    dataset, seq_along(dataset$y), NULL, config,
    methods = config$primary_sparse_methods
  )
  for (method in config$primary_sparse_methods) {
    fit <- fitted$fits[[method]]
    statistics <- permutation_statistics(
      fit, dataset$y / fitted$scaler$y_scale,
      fitted$scaler$weights, config
    )
    observed_rows[[observed_index]] <- data.frame(
      task_key = paste(dataset$name, "observed", sep = "::"),
      dataset = dataset$name, permutation = 0L, method = method,
      statistics, success = fit$success,
      seed = task_seed(config, "permutation", dataset_index),
      error = fit$error, warnings = fit$warnings,
      stringsAsFactors = FALSE
    )
    observed_index <- observed_index + 1L
  }
}
observed <- bind_rows_base(observed_rows)

tasks <- list()
index <- 1L
for (dataset_index in seq_along(datasets)) {
  dataset_name <- names(datasets)[dataset_index]
  for (replication in seq_len(config$permutation_replicates)) {
    tasks[[index]] <- list(
      key = permutation_task_key(dataset_name, replication),
      dataset = dataset_name, replication = replication,
      seed = task_seed(config, "permutation", dataset_index, replication)
    )
    index <- index + 1L
  }
}
completed <- if (resume && file.exists(paths$checkpoint)) {
  readRDS(paths$checkpoint)$completed
} else list()
remaining <- Filter(function(task) !task$key %in% names(completed), tasks)
cat(length(tasks), "permutation tasks planned;", length(remaining), "remaining.\n")

batch_size <- max(1L, 2L * config$workers)
if (length(remaining)) {
  batches <- split(remaining, ceiling(seq_along(remaining) / batch_size))
  for (batch_index in seq_along(batches)) {
    worker <- function(task) tryCatch(
      run_permutation_task(task, datasets, config),
      error = function(condition) list(
        key = task$key, fatal_error = conditionMessage(condition)
      )
    )
    batch_results <- run_task_batch(
      batches[[batch_index]], worker, config$workers
    )
    for (result in batch_results) completed[[result$key]] <- result
    atomic_save_rds(list(completed = completed, config = config), paths$checkpoint)
    fatal <- vapply(batch_results, function(x) !is.null(x$fatal_error), logical(1L))
    if (any(fatal)) {
      stop(
        "Permutation task failure after checkpoint: ",
        paste(vapply(batch_results[fatal], function(x) {
          paste(x$key, x$fatal_error, sep = ": ")
        }, character(1L)), collapse = " | "), call. = FALSE
      )
    }
    cat("Permutation progress:", length(completed), "/", length(tasks), "tasks\n")
  }
}

ordered <- completed[vapply(tasks, `[[`, character(1L), "key")]
null_statistics <- bind_rows_base(lapply(ordered, `[[`, "statistics"))
selection <- bind_rows_base(lapply(ordered, `[[`, "selection"))
statistics <- bind_rows_base(list(observed, null_statistics))

summary_rows <- list()
summary_index <- 1L
for (dataset_name in names(datasets)) {
  for (method in config$primary_sparse_methods) {
    current_observed <- observed[
      observed$dataset == dataset_name & observed$method == method,
      , drop = FALSE
    ]
    current_null <- null_statistics[
      null_statistics$dataset == dataset_name &
        null_statistics$method == method & null_statistics$success,
      , drop = FALSE
    ]
    for (statistic in c(
      "bic_gain_over_intercept", "beta_l2_scaled", "active_grid_fraction"
    )) {
      observed_value <- current_observed[[statistic]][1L]
      null_values <- current_null[[statistic]]
      null_values <- null_values[is.finite(null_values)]
      summary_rows[[summary_index]] <- data.frame(
        dataset = dataset_name, method = method, statistic = statistic,
        observed = observed_value,
        null_mean = mean(null_values),
        null_sd = stats::sd(null_values),
        null_q95 = unname(stats::quantile(null_values, 0.95, names = FALSE)),
        p_value_upper = (1 + sum(null_values >= observed_value)) /
          (1 + length(null_values)),
        successful_permutations = length(null_values),
        planned_permutations = config$permutation_replicates,
        primary_statistic = identical(statistic, "bic_gain_over_intercept"),
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L
    }
  }
}
permutation_summary <- bind_rows_base(summary_rows)
result <- list(
  version = config$version, mode = mode, configuration = config,
  statistics = statistics, summary = permutation_summary,
  selection = selection
)
atomic_save_rds(result, paths$result)
atomic_write_csv(statistics, paths$statistics)
atomic_write_csv(permutation_summary, paths$summary)
atomic_write_csv(selection, paths$selection)
capture_session_info(paths$session)
cat("Permutation analysis completed with", nrow(null_statistics),
    "method-specific null fits.\n")

