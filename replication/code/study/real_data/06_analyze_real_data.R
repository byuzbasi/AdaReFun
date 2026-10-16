script_dir <- local({
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) dirname(normalizePath(sub("^--file=", "", file_argument[1L]))) else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
for (file in c("utils.R", "config.R", "data.R")) {
  source(file.path(script_dir, "R", file))
}
activate_fmcp_libraries()

mode <- real_mode()
config <- apply_resolved_grid(real_data_config(mode, project_root))
datasets <- load_processed_datasets(config)
input_paths <- list(
  full = file.path(config$results_dir,
                   versioned_name("full_models", mode, "rds")),
  cv = file.path(config$results_dir,
                 versioned_name("repeated_cv", mode, "rds")),
  stability = file.path(config$results_dir,
                        versioned_name("stability", mode, "rds")),
  permutation = file.path(config$results_dir,
                          versioned_name("permutation", mode, "rds"))
)
missing_inputs <- unlist(input_paths)[!file.exists(unlist(input_paths))]
if (length(missing_inputs)) {
  stop("Missing analysis inputs: ", paste(basename(missing_inputs), collapse = ", "),
       call. = FALSE)
}

paths <- list(
  analysis = file.path(config$results_dir,
                       versioned_name("analysis", mode, "rds")),
  performance = file.path(config$results_dir,
                          versioned_name("cv_performance_summary", mode, "csv")),
  repeat_metrics = file.path(config$results_dir,
                             versioned_name("cv_repeat_metrics", mode, "csv")),
  paired = file.path(config$results_dir,
                     versioned_name("paired_bootstrap_comparisons", mode, "csv")),
  ranks = file.path(config$results_dir,
                    versioned_name("cv_method_ranks", mode, "csv")),
  full_summary = file.path(config$results_dir,
                           versioned_name("full_model_summary", mode, "csv")),
  gamma = file.path(config$results_dir,
                    versioned_name("gamma_sensitivity_summary", mode, "csv")),
  validation = file.path(config$results_dir,
                         versioned_name("analysis_validation", mode, "csv")),
  session = file.path(config$results_dir,
                      versioned_name("analysis_sessionInfo", mode, "txt"))
)
ensure_new_output(unlist(paths), overwrite = real_overwrite())

full <- readRDS(input_paths$full)
cv <- readRDS(input_paths$cv)
stability <- readRDS(input_paths$stability)
permutation <- readRDS(input_paths$permutation)
predictions <- cv$predictions

metric_rows <- function(data, grouping) {
  groups <- split(data, interaction(data[grouping], drop = TRUE))
  bind_rows_base(lapply(groups, function(group) {
    metrics <- metric_vector(group$observed, group$predicted)
    row <- group[1L, grouping, drop = FALSE]
    row$RMSE <- metrics[["RMSE"]]
    row$MAE <- metrics[["MAE"]]
    row$R2 <- metrics[["R2"]]
    row$n_predictions <- nrow(group)
    row$n_subjects <- length(unique(group$subject_id))
    row$all_predictions_finite <- all(is.finite(group$predicted))
    row
  }))
}

repeat_metrics <- metric_rows(
  predictions, c("dataset", "method", "repeat_index")
)
point_metrics <- metric_rows(predictions, c("dataset", "method"))

cat("Running", config$subject_bootstrap_replicates,
    "paired subject-level bootstrap replicates per dataset...\n")
bootstrap_rows <- list()
bootstrap_index <- 1L
for (dataset_index in seq_along(datasets)) {
  dataset_name <- names(datasets)[dataset_index]
  current <- predictions[predictions$dataset == dataset_name, , drop = FALSE]
  subjects <- sort(unique(current$subject_id))
  method_data <- lapply(config$methods, function(method) {
    data <- current[current$method == method, , drop = FALSE]
    list(data = data, rows = split(seq_len(nrow(data)), data$subject_id))
  })
  names(method_data) <- config$methods
  set.seed(task_seed(config, "analysis", dataset_index))
  bootstrap_subjects <- replicate(
    config$subject_bootstrap_replicates,
    sample(subjects, length(subjects), replace = TRUE),
    simplify = FALSE
  )
  for (replication in seq_along(bootstrap_subjects)) {
    sampled <- bootstrap_subjects[[replication]]
    for (method in config$methods) {
      holder <- method_data[[method]]
      indices <- unlist(holder$rows[as.character(sampled)], use.names = FALSE)
      data <- holder$data[indices, , drop = FALSE]
      metrics <- metric_vector(data$observed, data$predicted)
      bootstrap_rows[[bootstrap_index]] <- data.frame(
        dataset = dataset_name, bootstrap = replication,
        method = method, RMSE = metrics[["RMSE"]],
        MAE = metrics[["MAE"]], R2 = metrics[["R2"]],
        sampled_subjects = length(sampled),
        stringsAsFactors = FALSE
      )
      bootstrap_index <- bootstrap_index + 1L
    }
  }
}
bootstrap_metrics <- bind_rows_base(bootstrap_rows)

performance_rows <- list()
performance_index <- 1L
for (dataset_name in names(datasets)) {
  for (method in config$methods) {
    point <- point_metrics[
      point_metrics$dataset == dataset_name & point_metrics$method == method,
      , drop = FALSE
    ]
    boot <- bootstrap_metrics[
      bootstrap_metrics$dataset == dataset_name &
        bootstrap_metrics$method == method, , drop = FALSE
    ]
    row <- data.frame(
      dataset = dataset_name, method = method,
      n_subjects = point$n_subjects,
      n_predictions = point$n_predictions,
      prediction_success_rate = mean(predictions$success[
        predictions$dataset == dataset_name & predictions$method == method
      ]),
      stringsAsFactors = FALSE
    )
    for (metric in c("RMSE", "MAE", "R2")) {
      interval <- stats::quantile(
        boot[[metric]], c(0.025, 0.975), na.rm = TRUE, names = FALSE
      )
      row[[metric]] <- point[[metric]]
      row[[paste0(metric, "_lower")]] <- interval[1L]
      row[[paste0(metric, "_upper")]] <- interval[2L]
    }
    performance_rows[[performance_index]] <- row
    performance_index <- performance_index + 1L
  }
}
performance <- bind_rows_base(performance_rows)

paired_rows <- list()
paired_index <- 1L
for (dataset_name in names(datasets)) {
  point_fmcp <- performance[
    performance$dataset == dataset_name & performance$method == "FMCP",
    , drop = FALSE
  ]
  bootstrap_fmcp <- bootstrap_metrics[
    bootstrap_metrics$dataset == dataset_name &
      bootstrap_metrics$method == "FMCP", , drop = FALSE
  ]
  bootstrap_fmcp <- bootstrap_fmcp[order(bootstrap_fmcp$bootstrap), ]
  for (comparator in setdiff(config$methods, "FMCP")) {
    point_other <- performance[
      performance$dataset == dataset_name & performance$method == comparator,
      , drop = FALSE
    ]
    bootstrap_other <- bootstrap_metrics[
      bootstrap_metrics$dataset == dataset_name &
        bootstrap_metrics$method == comparator, , drop = FALSE
    ]
    bootstrap_other <- bootstrap_other[order(bootstrap_other$bootstrap), ]
    for (metric in c("RMSE", "MAE", "R2")) {
      differences <- bootstrap_fmcp[[metric]] - bootstrap_other[[metric]]
      comparison_tolerance <- sqrt(.Machine$double.eps) * max(
        1, abs(point_fmcp[[metric]]), abs(point_other[[metric]])
      )
      differences[abs(differences) <= comparison_tolerance] <- 0
      interval <- stats::quantile(
        differences, c(0.025, 0.975), na.rm = TRUE, names = FALSE
      )
      lower_tail <- mean(differences <= 0, na.rm = TRUE)
      upper_tail <- mean(differences >= 0, na.rm = TRUE)
      fmcp_better <- if (metric %in% c("RMSE", "MAE")) {
        mean(differences < 0, na.rm = TRUE)
      } else mean(differences > 0, na.rm = TRUE)
      point_difference <- point_fmcp[[metric]] - point_other[[metric]]
      if (abs(point_difference) <= comparison_tolerance) point_difference <- 0
      paired_rows[[paired_index]] <- data.frame(
        dataset = dataset_name, comparator = comparator, metric = metric,
        contrast = "FMCP minus comparator",
        estimate = point_difference,
        lower = interval[1L], upper = interval[2L],
        probability_fmcp_better = fmcp_better,
        probability_numerical_tie = mean(differences == 0, na.rm = TRUE),
        bootstrap_two_sided_p = min(1, 2 * min(lower_tail, upper_tail)),
        bootstrap_replicates = nrow(bootstrap_fmcp),
        comparison_tolerance = comparison_tolerance,
        primary_comparison = identical(comparator, "fSCAD"),
        stringsAsFactors = FALSE
      )
      paired_index <- paired_index + 1L
    }
  }
}
paired <- bind_rows_base(paired_rows)

rank_rows <- list()
rank_index <- 1L
for (dataset_name in names(datasets)) {
  for (repeat_index in seq_len(config$cv_repeats)) {
    current <- repeat_metrics[
      repeat_metrics$dataset == dataset_name &
        repeat_metrics$repeat_index == repeat_index, , drop = FALSE
    ]
    for (metric in c("RMSE", "MAE", "R2")) {
      decreasing <- identical(metric, "R2")
      ranks <- rank(if (decreasing) -current[[metric]] else current[[metric]],
                    ties.method = "average")
      rank_rows[[rank_index]] <- data.frame(
        dataset = dataset_name, repeat_index = repeat_index,
        method = current$method, metric = metric, rank = ranks,
        stringsAsFactors = FALSE
      )
      rank_index <- rank_index + 1L
    }
  }
}
rank_replicates <- bind_rows_base(rank_rows)
rank_groups <- split(
  rank_replicates,
  interaction(rank_replicates$dataset, rank_replicates$method,
              rank_replicates$metric, drop = TRUE),
  drop = TRUE
)
ranks <- bind_rows_base(lapply(rank_groups, function(group) {
  data.frame(
    dataset = group$dataset[1L], method = group$method[1L],
    metric = group$metric[1L], mean_rank = mean(group$rank),
    median_rank = stats::median(group$rank),
    best_rank_rate = mean(group$rank == 1), stringsAsFactors = FALSE
  )
}))

null_fraction <- aggregate(
  estimated_null ~ dataset + method,
  full$coefficients[full$coefficients$context == "full", ], mean
)
names(null_fraction)[names(null_fraction) == "estimated_null"] <-
  "estimated_null_grid_fraction"
full_summary <- merge(
  full$selection, null_fraction, by = c("dataset", "method"),
  all.x = TRUE, sort = FALSE
)

gamma_groups <- split(
  full$gamma_curves,
  interaction(full$gamma_curves$dataset, full$gamma_curves$context, drop = TRUE),
  drop = TRUE
)
gamma_summary_alternatives <- bind_rows_base(lapply(gamma_groups, function(group) {
  weights <- trapezoid_weights(group$argval_u)
  data.frame(
    dataset = group$dataset[1L], context = group$context[1L],
    gamma = group$shape[1L], M = group$M[1L],
    lambda = group$lambda[1L], kappa = group$kappa[1L],
    bic = group$bic[1L],
    beta_l2_scaled = integral_l2(group$beta_scaled, weights),
    estimated_null_grid_fraction = mean(group$estimated_null),
    stringsAsFactors = FALSE
  )
}))
primary_gamma_rows <- list()
for (dataset_name in names(datasets)) {
  curve <- full$coefficients[
    full$coefficients$dataset == dataset_name &
      full$coefficients$method == "FMCP" &
      full$coefficients$context == "full", , drop = FALSE
  ]
  selected <- full$selection[
    full$selection$dataset == dataset_name &
      full$selection$method == "FMCP", , drop = FALSE
  ]
  primary_gamma_rows[[dataset_name]] <- data.frame(
    dataset = dataset_name, context = "gamma_3_primary", gamma = 3,
    M = selected$M[1L], lambda = selected$lambda[1L],
    kappa = selected$kappa[1L], bic = selected$bic[1L],
    beta_l2_scaled = integral_l2(
      curve$beta_scaled, trapezoid_weights(curve$argval_u)
    ),
    estimated_null_grid_fraction = mean(curve$estimated_null),
    stringsAsFactors = FALSE
  )
}
gamma_summary <- bind_rows_base(list(
  gamma_summary_alternatives, bind_rows_base(primary_gamma_rows)
))
gamma_summary <- gamma_summary[order(gamma_summary$dataset, gamma_summary$gamma), ]

validation <- data.frame(
  check = c(
    "CV prediction row count", "CV predictions finite",
    "CV method success", "Stability planned fits complete",
    "Stability fits successful", "Permutation planned fits complete",
    "Permutation fits successful", "One grid expansion maximum respected"
  ),
  expected = c(
    sum(vapply(datasets, function(x) length(x$y), integer(1L))) *
      config$cv_repeats * length(config$methods),
    1, 1,
    length(datasets) * length(config$primary_sparse_methods) *
      config$stability_replicates,
    1,
    length(datasets) * length(config$primary_sparse_methods) *
      config$permutation_replicates,
    1, 1
  ),
  observed = c(
    nrow(predictions), all(is.finite(predictions$predicted)),
    all(predictions$success),
    nrow(unique(stability$replicate_curves[
      c("dataset", "method", "replication")
    ])),
    all(stability$replicate_curves$success),
    nrow(permutation$statistics[permutation$statistics$permutation > 0, ]),
    all(permutation$statistics$success),
    full$resolved_grid$maximum_grid_expansions == 1L
  ),
  stringsAsFactors = FALSE
)
validation$passed <- as.character(validation$expected) ==
  as.character(validation$observed)

analysis <- list(
  version = config$version, mode = mode, configuration = config,
  performance = performance, repeat_metrics = repeat_metrics,
  bootstrap_metrics = bootstrap_metrics, paired = paired,
  ranks = ranks, rank_replicates = rank_replicates,
  full_summary = full_summary, gamma_summary = gamma_summary,
  stability_probabilities = stability$probabilities,
  stable_null_segments = stability$segments,
  permutation_summary = permutation$summary,
  permutation_statistics = permutation$statistics,
  validation = validation
)
atomic_save_rds(analysis, paths$analysis)
atomic_write_csv(performance, paths$performance)
atomic_write_csv(repeat_metrics, paths$repeat_metrics)
atomic_write_csv(paired, paths$paired)
atomic_write_csv(ranks, paths$ranks)
atomic_write_csv(full_summary, paths$full_summary)
atomic_write_csv(gamma_summary, paths$gamma)
atomic_write_csv(validation, paths$validation)
capture_session_info(paths$session)

cat("Analysis completed. Validation:\n")
print(validation)
