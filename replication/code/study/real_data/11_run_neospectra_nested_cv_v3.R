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
assert_packages(c("FMCP", "Rcpp", "RcppArmadillo"))

mode <- real_mode()
resume <- real_resume()
config <- neospectra_v3_config(mode, project_root)
processed <- load_neospectra_v3_processed(config)
dataset <- processed$sample_dataset
checkpoint_signature <- neospectra_v3_checkpoint_signature(
  config, processed, "nested_cv"
)
paths <- list(
  result = neospectra_v3_result_path(config, "nested_cv", "rds"),
  predictions = neospectra_v3_result_path(
    config, "nested_cv_predictions", "csv"
  ),
  tuning = neospectra_v3_result_path(
    config, "nested_cv_tuning", "csv"
  ),
  inner_folds = neospectra_v3_result_path(
    config, "nested_cv_inner_folds", "csv"
  ),
  selection = neospectra_v3_result_path(
    config, "nested_cv_selection", "csv"
  ),
  folds = neospectra_v3_result_path(config, "nested_cv_folds", "csv"),
  external = neospectra_v3_result_path(
    config, "geographic_external_status", "csv"
  ),
  validation = neospectra_v3_result_path(
    config, "nested_cv_validation", "csv"
  ),
  checkpoint = file.path(
    config$results_dir, "checkpoints",
    neospectra_v3_name("nested_cv_checkpoint", config, "rds")
  ),
  session = neospectra_v3_result_path(
    config, "nested_cv_sessionInfo", "txt"
  )
)
ensure_new_output(
  unlist(paths[setdiff(names(paths), "checkpoint")]),
  resume = resume, overwrite = real_overwrite()
)
if (!resume && file.exists(paths$checkpoint) && !real_overwrite()) {
  stop("A NeoSpectra v3 CV checkpoint exists; use --resume.", call. = FALSE)
}
if (nrow(dataset$x) < config$outer_folds) {
  stop("Outer folds exceed the number of unique soil samples.", call. = FALSE)
}

initialize_real_method_environment(project_root)

tasks <- list()
fold_rows <- list()
task_index <- 1L
for (repeat_index in seq_len(config$outer_repeats)) {
  seed <- neospectra_v3_seed(config, "outer", repeat_index)
  labels <- make_grouped_folds_v3(
    dataset$sample_id, config$outer_folds, seed
  )
  fold_rows[[repeat_index]] <- data.frame(
    repeat_index = repeat_index, sample_id = dataset$sample_id,
    fold = labels, seed = seed, stringsAsFactors = FALSE
  )
  for (fold in seq_len(config$outer_folds)) {
    key <- paste0("outer_r", repeat_index, "_f", fold)
    tasks[[task_index]] <- list(
      key = key, repeat_index = repeat_index, fold = fold,
      outer_labels = labels
    )
    task_index <- task_index + 1L
  }
}
fold_table <- bind_rows_base(fold_rows)

checkpoint <- if (resume && file.exists(paths$checkpoint)) {
  readRDS(paths$checkpoint)
} else list(completed = list())
if (!is.null(checkpoint$sample_id) &&
    !identical(checkpoint$sample_id, dataset$sample_id)) {
  stop("Checkpoint sample IDs do not match the processed v3 data.",
       call. = FALSE)
}
if (resume && file.exists(paths$checkpoint) &&
    !identical(checkpoint$signature, checkpoint_signature)) {
  stop("Checkpoint configuration/source signature does not match this run.",
       call. = FALSE)
}
completed <- checkpoint$completed
remaining <- Filter(function(task) !task$key %in% names(completed), tasks)
cat(length(tasks), "nested-CV tasks planned;", length(remaining), "remaining.\n")

for (task in remaining) {
  result <- tryCatch(
    outer_task_v3(
      dataset, task$repeat_index, task$fold, task$outer_labels, config
    ),
    error = function(condition) list(
      key = task$key, fatal_error = conditionMessage(condition)
    )
  )
  completed[[result$key]] <- result
  atomic_save_rds(
    list(
      completed = completed, sample_id = dataset$sample_id,
      signature = checkpoint_signature, configuration = config
    ),
    paths$checkpoint
  )
  if (!is.null(result$fatal_error)) {
    stop(
      "Nested-CV task failure after checkpoint: ", result$key, ": ",
      result$fatal_error, call. = FALSE
    )
  }
  cat("Nested-CV progress:", length(completed), "/", length(tasks), "tasks\n")
}

task_keys <- vapply(tasks, `[[`, character(1L), "key")
ordered <- completed[task_keys]
if (any(vapply(ordered, is.null, logical(1L)))) {
  stop("Nested-CV checkpoint is incomplete.", call. = FALSE)
}
predictions <- bind_rows_base(lapply(ordered, `[[`, "predictions"))
tuning <- bind_rows_base(lapply(ordered, `[[`, "tuning"))
inner_folds <- bind_rows_base(lapply(ordered, `[[`, "inner_folds"))
selection <- bind_rows_base(lapply(ordered, `[[`, "selection"))

external <- run_geographic_external_v3(dataset, config)
if (external$available) {
  predictions <- bind_rows_base(list(predictions, external$predictions))
  tuning <- bind_rows_base(list(tuning, external$tuning))
  inner_folds <- bind_rows_base(list(inner_folds, external$inner_folds))
  selection <- bind_rows_base(list(selection, external$selection))
}
external_status <- data.frame(
  evaluation = "US_train_Africa_test", available = external$available,
  reason = external$reason,
  n_us = sum(neospectra_v3_country_class(dataset$country) == "US"),
  n_africa = sum(neospectra_v3_country_class(dataset$country) == "Africa"),
  stringsAsFactors = FALSE
)
result <- list(
  version = config$version, mode = mode, configuration = config,
  folds = fold_table, predictions = predictions, tuning = tuning,
  inner_folds = inner_folds, selection = selection,
  external_status = external_status
)
validation <- validate_grouped_nested_cv_v3(result, processed, config)
if (!all(validation$passed)) {
  failed <- validation$check[!validation$passed]
  stop("Nested-CV validation failed: ", paste(failed, collapse = ", "),
       call. = FALSE)
}

atomic_save_rds(result, paths$result)
atomic_write_csv(predictions, paths$predictions)
atomic_write_csv(tuning, paths$tuning)
atomic_write_csv(inner_folds, paths$inner_folds)
atomic_write_csv(selection, paths$selection)
atomic_write_csv(fold_table, paths$folds)
atomic_write_csv(external_status, paths$external)
atomic_write_csv(validation, paths$validation)
capture_session_info(paths$session)
cat("NeoSpectra v3 nested CV completed with", nrow(predictions),
    "prediction rows.\n")
