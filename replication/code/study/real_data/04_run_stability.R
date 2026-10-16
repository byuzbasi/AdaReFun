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
                     versioned_name("stability", mode, "rds")),
  probabilities = file.path(config$results_dir,
                            versioned_name("stability_null_probabilities", mode, "csv")),
  segments = file.path(config$results_dir,
                       versioned_name("stable_null_segments", mode, "csv")),
  selection = file.path(config$results_dir,
                        versioned_name("stability_selection", mode, "csv")),
  checkpoint = file.path(config$results_dir, "checkpoints",
                         versioned_name("stability_checkpoint", mode, "rds")),
  session = file.path(config$results_dir,
                      versioned_name("stability_sessionInfo", mode, "txt"))
)
ensure_new_output(
  unlist(paths[c("result", "probabilities", "segments", "selection", "session")]),
  resume = resume, overwrite = real_overwrite()
)
if (!resume && file.exists(paths$checkpoint) && !real_overwrite()) {
  stop("A stability checkpoint exists; use --resume.", call. = FALSE)
}

initialize_real_method_environment(project_root)
tasks <- list()
index <- 1L
for (dataset_index in seq_along(datasets)) {
  dataset_name <- names(datasets)[dataset_index]
  for (replication in seq_len(config$stability_replicates)) {
    tasks[[index]] <- list(
      key = stability_task_key(dataset_name, replication),
      dataset = dataset_name, replication = replication,
      seed = task_seed(config, "stability", dataset_index, replication)
    )
    index <- index + 1L
  }
}
completed <- if (resume && file.exists(paths$checkpoint)) {
  readRDS(paths$checkpoint)$completed
} else list()
remaining <- Filter(function(task) !task$key %in% names(completed), tasks)
cat(length(tasks), "stability tasks planned;", length(remaining), "remaining.\n")

batch_size <- max(1L, 2L * config$workers)
if (length(remaining)) {
  batches <- split(remaining, ceiling(seq_along(remaining) / batch_size))
  for (batch_index in seq_along(batches)) {
    worker <- function(task) tryCatch(
      run_stability_task(task, datasets, config),
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
        "Stability task failure after checkpoint: ",
        paste(vapply(batch_results[fatal], function(x) {
          paste(x$key, x$fatal_error, sep = ": ")
        }, character(1L)), collapse = " | "), call. = FALSE
      )
    }
    cat("Stability progress:", length(completed), "/", length(tasks), "tasks\n")
  }
}

ordered <- completed[vapply(tasks, `[[`, character(1L), "key")]
curves <- bind_rows_base(lapply(ordered, `[[`, "curves"))
selection <- bind_rows_base(lapply(ordered, `[[`, "selection"))
sampled_subjects <- bind_rows_base(lapply(ordered, `[[`, "sampled_subjects"))

groups <- split(curves, interaction(
  curves$dataset, curves$method, curves$grid_index, drop = TRUE
))
probabilities <- bind_rows_base(lapply(groups, function(group) {
  successful <- group$success & !is.na(group$selected_null)
  data.frame(
    dataset = group$dataset[1L], method = group$method[1L],
    grid_index = group$grid_index[1L], argval_u = group$argval_u[1L],
    argval_physical = group$argval_physical[1L],
    null_selection_probability = if (any(successful)) {
      mean(group$selected_null[successful])
    } else NA_real_,
    successful_replicates = sum(successful),
    planned_replicates = config$stability_replicates,
    stringsAsFactors = FALSE
  )
}))
probabilities <- probabilities[order(
  probabilities$dataset, probabilities$method, probabilities$grid_index
), ]

segments <- list()
segment_index <- 1L
for (dataset_name in names(datasets)) {
  for (method in config$primary_sparse_methods) {
    current <- probabilities[
      probabilities$dataset == dataset_name & probabilities$method == method,
      , drop = FALSE
    ]
    for (threshold in config$stability_thresholds) {
      segments[[segment_index]] <- contiguous_segments(
        current$argval_physical,
        current$null_selection_probability >= threshold,
        dataset_name, method, threshold
      )
      segment_index <- segment_index + 1L
    }
  }
}
segments <- bind_rows_base(segments)

result <- list(
  version = config$version, mode = mode, configuration = config,
  replicate_curves = curves, sampled_subjects = sampled_subjects,
  selection = selection, probabilities = probabilities, segments = segments
)
atomic_save_rds(result, paths$result)
atomic_write_csv(probabilities, paths$probabilities)
atomic_write_csv(segments, paths$segments)
atomic_write_csv(selection, paths$selection)
capture_session_info(paths$session)
cat("Stability analysis completed with", nrow(curves), "curve-grid rows.\n")

