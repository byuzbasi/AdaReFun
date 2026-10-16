make_repeated_folds <- function(n, repeats, folds, seed) {
  if (folds > n) stop("Number of folds exceeds sample size.", call. = FALSE)
  result <- vector("list", repeats)
  for (repeat_index in seq_len(repeats)) {
    set.seed(seed + repeat_index)
    labels <- sample(rep(seq_len(folds), length.out = n))
    result[[repeat_index]] <- labels
  }
  result
}

cv_task_key <- function(dataset, repeat_index, fold) {
  paste(dataset, repeat_index, fold, sep = "::")
}

run_cv_task <- function(task, datasets, config) {
  dataset <- datasets[[task$dataset]]
  test_index <- which(task$fold_labels == task$fold)
  train_index <- setdiff(seq_along(dataset$y), test_index)
  set.seed(task$seed)
  fitted <- fit_real_method_set(
    dataset, train_index, test_index, config, methods = config$methods
  )
  rows <- lapply(fitted$fits, function(fit) {
    data.frame(
      task_key = task$key, dataset = dataset$name,
      repeat_index = task$repeat_index, fold = task$fold,
      subject_id = dataset$subject_id[test_index],
      observed = dataset$y[test_index], predicted = fit$prediction,
      residual = dataset$y[test_index] - fit$prediction,
      method = fit$method, success = fit$success,
      seed = task$seed, cutoff_scaled = fitted$cutoff_scaled,
      y_scale = fitted$scaler$y_scale,
      x_functional_rms_scale = fitted$scaler$x_scale,
      error = fit$error, warnings = fit$warnings,
      stringsAsFactors = FALSE
    )
  })
  list(
    key = task$key, predictions = bind_rows_base(rows),
    selection = selection_rows(
      fitted, dataset$name,
      context = paste0("cv_r", task$repeat_index, "_f", task$fold)
    )
  )
}

stability_task_key <- function(dataset, replication) {
  paste(dataset, replication, sep = "::")
}

run_stability_task <- function(task, datasets, config) {
  dataset <- datasets[[task$dataset]]
  set.seed(task$seed)
  bootstrap_index <- sample(
    seq_along(dataset$y), length(dataset$y), replace = TRUE
  )
  fitted <- fit_real_method_set(
    dataset, bootstrap_index, test_index = NULL, config,
    methods = config$primary_sparse_methods
  )
  curves <- bind_rows_base(lapply(fitted$fits, function(fit) {
    data.frame(
      task_key = task$key, dataset = dataset$name,
      replication = task$replication, method = fit$method,
      grid_index = seq_along(dataset$argvals_u),
      argval_u = dataset$argvals_u,
      argval_physical = dataset$argvals_physical,
      beta_scaled = fit$beta_scaled,
      selected_null = abs(fit$beta_scaled) <= config$zero_tolerance,
      success = fit$success, seed = task$seed,
      stringsAsFactors = FALSE
    )
  }))
  list(
    key = task$key, curves = curves,
    selection = selection_rows(
      fitted, dataset$name,
      context = paste0("stability_", task$replication)
    ),
    sampled_subjects = data.frame(
      task_key = task$key, draw = seq_along(bootstrap_index),
      subject_id = dataset$subject_id[bootstrap_index],
      stringsAsFactors = FALSE
    )
  )
}

permutation_task_key <- function(dataset, replication) {
  paste(dataset, replication, sep = "::")
}

permutation_statistics <- function(fit, y_scaled, weights, config) {
  selected <- fit$selected
  null_rss <- sum((y_scaled - mean(y_scaled))^2)
  null_bic <- length(y_scaled) * log(max(
    null_rss / length(y_scaled), .Machine$double.xmin
  ))
  model_bic <- if (!is.null(selected) && is.finite(selected$bic[1L])) {
    selected$bic[1L]
  } else NA_real_
  data.frame(
    bic_gain_over_intercept = null_bic - model_bic,
    beta_l2_scaled = integral_l2(fit$beta_scaled, weights),
    active_grid_fraction = mean(abs(fit$beta_scaled) > config$zero_tolerance),
    null_bic = null_bic, model_bic = model_bic,
    stringsAsFactors = FALSE
  )
}

run_permutation_task <- function(task, datasets, config) {
  original <- datasets[[task$dataset]]
  set.seed(task$seed)
  dataset <- original
  dataset$y <- sample(original$y, replace = FALSE)
  fitted <- fit_real_method_set(
    dataset, seq_along(dataset$y), test_index = NULL, config,
    methods = config$primary_sparse_methods
  )
  rows <- lapply(fitted$fits, function(fit) {
    statistics <- permutation_statistics(
      fit, dataset$y / fitted$scaler$y_scale,
      fitted$scaler$weights, config
    )
    data.frame(
      task_key = task$key, dataset = dataset$name,
      permutation = task$replication, method = fit$method,
      statistics, success = fit$success, seed = task$seed,
      error = fit$error, warnings = fit$warnings,
      stringsAsFactors = FALSE
    )
  })
  list(
    key = task$key, statistics = bind_rows_base(rows),
    selection = selection_rows(
      fitted, dataset$name,
      context = paste0("permutation_", task$replication)
    )
  )
}

contiguous_segments <- function(values, selected, dataset, method,
                                threshold) {
  runs <- rle(selected)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  keep <- which(runs$values)
  if (!length(keep)) {
    return(data.frame(
      dataset = character(), method = character(), threshold = numeric(),
      segment = integer(), grid_start = integer(), grid_end = integer(),
      physical_start = numeric(), physical_end = numeric(),
      physical_length = numeric(), grid_points = integer(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    dataset = dataset, method = method, threshold = threshold,
    segment = seq_along(keep), grid_start = starts[keep],
    grid_end = ends[keep], physical_start = values[starts[keep]],
    physical_end = values[ends[keep]],
    physical_length = values[ends[keep]] - values[starts[keep]],
    grid_points = runs$lengths[keep], stringsAsFactors = FALSE
  )
}
