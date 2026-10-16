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
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
paths <- list(
  result = file.path(config$results_dir,
                     versioned_name("repeated_cv", mode, "rds")),
  predictions = file.path(config$results_dir,
                          versioned_name("repeated_cv_predictions", mode, "csv")),
  selection = file.path(config$results_dir,
                        versioned_name("repeated_cv_selection", mode, "csv")),
  folds = file.path(config$results_dir,
                    versioned_name("repeated_cv_folds", mode, "csv")),
  checkpoint = file.path(config$results_dir, "checkpoints",
                         versioned_name("repeated_cv_checkpoint", mode, "rds")),
  session = file.path(config$results_dir,
                      versioned_name("repeated_cv_sessionInfo", mode, "txt"))
)
ensure_new_output(
  unlist(paths[c("result", "predictions", "selection", "folds", "session")]),
  resume = resume, overwrite = real_overwrite()
)
if (!resume && file.exists(paths$checkpoint) && !real_overwrite()) {
  stop("A CV checkpoint exists; use --resume.", call. = FALSE)
}

initialize_real_method_environment(project_root)

tasks <- list()
fold_rows <- list()
task_index <- 1L
fold_index <- 1L
for (dataset_index in seq_along(datasets)) {
  dataset_name <- names(datasets)[dataset_index]
  dataset <- datasets[[dataset_index]]
  folds <- make_repeated_folds(
    length(dataset$y), config$cv_repeats, config$cv_folds,
    task_seed(config, "cv", dataset_index)
  )
  for (repeat_index in seq_len(config$cv_repeats)) {
    fold_rows[[fold_index]] <- data.frame(
      dataset = dataset_name, repeat_index = repeat_index,
      subject_id = dataset$subject_id,
      fold = folds[[repeat_index]],
      seed = task_seed(config, "cv", dataset_index, repeat_index),
      stringsAsFactors = FALSE
    )
    fold_index <- fold_index + 1L
    for (fold in seq_len(config$cv_folds)) {
      key <- cv_task_key(dataset_name, repeat_index, fold)
      tasks[[task_index]] <- list(
        key = key, dataset = dataset_name,
        repeat_index = repeat_index, fold = fold,
        fold_labels = folds[[repeat_index]],
        seed = task_seed(config, "cv", dataset_index, repeat_index, fold)
      )
      task_index <- task_index + 1L
    }
  }
}
fold_table <- bind_rows_base(fold_rows)

completed <- if (resume && file.exists(paths$checkpoint)) {
  readRDS(paths$checkpoint)$completed
} else list()
remaining <- Filter(function(task) !task$key %in% names(completed), tasks)
cat(length(tasks), "CV tasks planned;", length(remaining), "remaining.\n")

batch_size <- max(1L, 2L * config$workers)
if (length(remaining)) {
  batches <- split(remaining, ceiling(seq_along(remaining) / batch_size))
  for (batch_index in seq_along(batches)) {
    worker <- function(task) tryCatch(
      run_cv_task(task, datasets, config),
      error = function(condition) list(
        key = task$key, fatal_error = conditionMessage(condition)
      )
    )
    batch_results <- run_task_batch(
      batches[[batch_index]], worker, config$workers
    )
    for (result in batch_results) completed[[result$key]] <- result
    atomic_save_rds(list(completed = completed, config = config), paths$checkpoint)
    fatal <- vapply(
      batch_results, function(result) !is.null(result$fatal_error),
      FUN.VALUE = logical(1L)
    )
    if (any(fatal)) {
      messages <- vapply(
        batch_results[fatal], function(result) {
          paste(result$key, result$fatal_error, sep = ": ")
        }, FUN.VALUE = character(1L)
      )
      stop("CV task failure after checkpoint: ", paste(messages, collapse = " | "),
           call. = FALSE)
    }
    cat("CV progress:", length(completed), "/", length(tasks), "tasks\n")
  }
}

ordered <- completed[vapply(tasks, `[[`, character(1L), "key")]
predictions <- bind_rows_base(lapply(ordered, `[[`, "predictions"))
selection <- bind_rows_base(lapply(ordered, `[[`, "selection"))
expected <- sum(vapply(datasets, function(dataset) length(dataset$y), integer(1L))) *
  config$cv_repeats * length(config$methods)
if (nrow(predictions) != expected) {
  stop("CV row-count validation failed: expected ", expected,
       ", observed ", nrow(predictions), call. = FALSE)
}
result <- list(
  version = config$version, mode = mode, configuration = config,
  folds = fold_table, predictions = predictions, selection = selection
)
atomic_save_rds(result, paths$result)
atomic_write_csv(predictions, paths$predictions)
atomic_write_csv(selection, paths$selection)
atomic_write_csv(fold_table, paths$folds)
capture_session_info(paths$session)
cat("Repeated CV completed with", nrow(predictions), "paired predictions.\n")

