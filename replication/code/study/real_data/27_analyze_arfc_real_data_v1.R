#!/usr/bin/env Rscript

# Analyze a completed ARFC-Select real-data run without refitting any model.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("real_data/27_analyze_arfc_real_data_v1.R", mustWork = TRUE)
}
real_dir <- dirname(script_file)
project_dir <- dirname(real_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
unknown <- setdiff(arguments, "--smoke")
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

helper_file <- file.path(real_dir, "R", "arfc_real_data_v1.R")
source(helper_file, local = FALSE)
config <- arfc_real_config(project_dir, smoke)
mode_prefix <- if (smoke) "smoke" else "production"
prefix <- paste0("arfc_real_data_v1_", mode_prefix)
run_dir <- file.path(config$results_root, paste0(mode_prefix, "_run_v1"))
analysis_dir <- file.path(
  config$results_root, paste0(mode_prefix, "_analysis_v1")
)
if (!dir.exists(run_dir)) {
  stop("Completed run bundle is missing: ", run_dir, call. = FALSE)
}
if (dir.exists(analysis_dir)) {
  stop("Analysis bundle already exists: ", analysis_dir, call. = FALSE)
}
manifest_path <- file.path(run_dir, paste0(prefix, "_manifest.csv"))
run_path <- file.path(run_dir, paste0(prefix, "_run.rds"))
if (!file.exists(manifest_path) || !file.exists(run_path)) {
  stop("Run manifest/RDS is incomplete.", call. = FALSE)
}
run_manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
manifest_files <- file.path(run_dir, run_manifest$file)
if (anyDuplicated(run_manifest$file) || !all(file.exists(manifest_files)) ||
    any(as.numeric(file.info(manifest_files)$size) != run_manifest$bytes) ||
    any(arfc_real_md5(manifest_files) != run_manifest$md5)) {
  stop("Run manifest byte/MD5 verification failed.", call. = FALSE)
}
run <- readRDS(run_path)
if (!identical(run$mode, config$mode) ||
    length(unique(run_manifest$scientific_signature)) != 1L ||
    !identical(unique(run_manifest$scientific_signature),
               run$scientific_signature)) {
  stop("Run mode/scientific-signature mismatch.", call. = FALSE)
}
config <- run$config

metric_row <- function(value) {
  residual <- value$observed - value$predicted
  tss <- sum((value$observed - mean(value$observed))^2)
  data.frame(
    dataset = value$dataset[1L], method = value$method[1L],
    n_ids = length(unique(value$id)), n_predictions = nrow(value),
    RMSE = sqrt(mean(residual^2)), MAE = mean(abs(residual)),
    R2 = if (tss > 0) 1 - sum(residual^2) / tss else NA_real_,
    stringsAsFactors = FALSE
  )
}
prediction_split <- split(
  run$predictions,
  interaction(run$predictions$dataset, run$predictions$method, drop = TRUE)
)
prediction_metrics <- arfc_real_bind_rows(lapply(
  prediction_split, metric_row
))
prediction_metrics <- prediction_metrics[
  order(match(prediction_metrics$dataset, config$dataset_order),
        match(prediction_metrics$method,
              c(config$selection_label, config$method_order,
                config$secondary_methods))), , drop = FALSE
]

selection_frequency <- arfc_real_bind_rows(lapply(
  config$dataset_order, function(dataset_name) {
    selected <- run$outer_selection$selected_method[
      run$outer_selection$dataset == dataset_name
    ]
    counts <- table(factor(selected, levels = config$method_order))
    data.frame(
      dataset = dataset_name, method = config$method_order,
      count = as.integer(counts), total_tasks = length(selected),
      frequency = as.integer(counts) / length(selected),
      source = "nested_outer_tasks", stringsAsFactors = FALSE
    )
  }
))
bootstrap_selection_frequency <- arfc_real_bind_rows(lapply(
  config$dataset_order, function(dataset_name) {
    selected <- run$bootstrap_selection$selected_method[
      run$bootstrap_selection$dataset == dataset_name
    ]
    counts <- table(factor(selected, levels = config$method_order))
    data.frame(
      dataset = dataset_name, method = config$method_order,
      count = as.integer(counts), total_tasks = length(selected),
      frequency = as.integer(counts) / length(selected),
      source = "full_reselection_bootstrap", stringsAsFactors = FALSE
    )
  }
))

cluster_bootstrap_dataset <- function(dataset_name) {
  values <- run$predictions[run$predictions$dataset == dataset_name, ]
  methods <- c(config$selection_label, config$method_order,
               config$secondary_methods)
  values <- values[values$method %in% methods, , drop = FALSE]
  ids <- sort(unique(values$id))
  sse <- matrix(0, nrow = length(ids), ncol = length(methods),
                dimnames = list(ids, methods))
  counts <- matrix(0L, nrow = length(ids), ncol = length(methods),
                   dimnames = list(ids, methods))
  for (method in methods) {
    current <- values[values$method == method, , drop = FALSE]
    current$residual_square <- (current$observed - current$predicted)^2
    summed <- tapply(current$residual_square, current$id, sum)
    counted <- table(current$id)
    sse[names(summed), method] <- summed
    counts[names(counted), method] <- as.integer(counted)
  }
  if (any(counts == 0L)) {
    stop("Incomplete matched prediction matrix for ", dataset_name,
         call. = FALSE)
  }
  point_rmse <- sqrt(colSums(sse) / colSums(counts))
  set.seed(as.integer(
    config$master_seed + 900000L +
      match(dataset_name, config$dataset_order) * 10000L
  ))
  differences <- matrix(
    NA_real_, nrow = config$uncertainty_replicates,
    ncol = length(methods) - 1L,
    dimnames = list(NULL, methods[-1L])
  )
  for (replication in seq_len(config$uncertainty_replicates)) {
    draw <- sample.int(length(ids), length(ids), replace = TRUE)
    rmse <- sqrt(colSums(sse[draw, , drop = FALSE]) /
                   colSums(counts[draw, , drop = FALSE]))
    differences[replication, ] <- rmse[config$selection_label] -
      rmse[methods[-1L]]
  }
  arfc_real_bind_rows(lapply(methods[-1L], function(method) {
    quantiles <- stats::quantile(
      differences[, method], c(0.025, 0.975), names = FALSE,
      type = 7
    )
    data.frame(
      dataset = dataset_name,
      contrast = paste(config$selection_label, "minus", method),
      reference_method = method,
      rmse_difference = point_rmse[config$selection_label] -
        point_rmse[method],
      lower_95 = quantiles[1L], upper_95 = quantiles[2L],
      bootstrap_replicates = config$uncertainty_replicates,
      cluster_unit = if (dataset_name == "Truck")
        "distributed row" else "unique soil-sample ID",
      stringsAsFactors = FALSE
    )
  }))
}
prediction_contrasts <- arfc_real_bind_rows(lapply(
  config$dataset_order, cluster_bootstrap_dataset
))

stability_summary <- arfc_real_bind_rows(lapply(
  config$dataset_order, function(dataset_name) {
    curves <- run$bootstrap_curves[
      run$bootstrap_curves$dataset == dataset_name, , drop = FALSE
    ]
    replications <- sort(unique(curves$replication))
    points <- sort(unique(curves$argument_u))
    if (length(replications) != config$bootstrap_replicates ||
        length(points) != config$support_grid_points ||
        nrow(curves) != length(replications) * length(points)) {
      stop("Bootstrap curve grid is incomplete for ", dataset_name,
           call. = FALSE)
    }
    curves <- curves[order(curves$replication, curves$argument_u), ]
    beta <- matrix(
      curves$beta_physical, nrow = length(replications), byrow = TRUE
    )
    zeros <- matrix(
      curves$is_exact_zero, nrow = length(replications), byrow = TRUE
    )
    first <- curves[curves$replication == replications[1L], ]
    data.frame(
      dataset = dataset_name, argument_u = points,
      physical_argument = first$physical_argument,
      beta_mean = colMeans(beta),
      beta_q10 = apply(beta, 2L, stats::quantile, probs = 0.10,
                       names = FALSE),
      beta_q90 = apply(beta, 2L, stats::quantile, probs = 0.90,
                       names = FALSE),
      zero_selection_frequency = colMeans(zeros),
      successful_reselections = length(replications),
      stringsAsFactors = FALSE
    )
  }
))

stable_zero_segments <- arfc_real_bind_rows(lapply(
  config$dataset_order, function(dataset_name) {
    values <- stability_summary[
      stability_summary$dataset == dataset_name, , drop = FALSE
    ]
    selected <- values$zero_selection_frequency >= 0.80
    runs <- rle(selected)
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1L
    indices <- which(runs$values)
    if (!length(indices)) return(NULL)
    arfc_real_bind_rows(lapply(seq_along(indices), function(segment) {
      run_index <- indices[segment]
      left_index <- starts[run_index]
      right_index <- ends[run_index]
      data.frame(
        dataset = dataset_name, threshold = 0.80,
        segment = segment,
        argument_u_start = values$argument_u[left_index],
        argument_u_end = values$argument_u[right_index],
        physical_start = values$physical_argument[left_index],
        physical_end = values$physical_argument[right_index],
        physical_length = values$physical_argument[right_index] -
          values$physical_argument[left_index],
        grid_points = right_index - left_index + 1L,
        stringsAsFactors = FALSE
      )
    }))
  }
))
if (!nrow(stable_zero_segments)) {
  stable_zero_segments <- data.frame(
    dataset = character(), threshold = numeric(), segment = integer(),
    argument_u_start = numeric(), argument_u_end = numeric(),
    physical_start = numeric(), physical_end = numeric(),
    physical_length = numeric(), grid_points = integer(),
    stringsAsFactors = FALSE
  )
}

full_selected_curves <- merge(
  run$full_curves, run$full_selection[c("dataset", "selected_method")],
  by = "dataset", all.x = TRUE, sort = FALSE
)
full_selected_curves <- full_selected_curves[
  full_selected_curves$method == full_selected_curves$selected_method,
  , drop = FALSE
]

lambda_diagnostics <- arfc_real_bind_rows(lapply(
  split(run$outer_candidate_metrics,
        interaction(run$outer_candidate_metrics$dataset,
                    run$outer_candidate_metrics$method, drop = TRUE)),
  function(value) {
    data.frame(
      dataset = value$dataset[1L], method = value$method[1L],
      fits = nrow(value), strict_fit_rate = mean(
        value$success & value$converged & value$refit_success &
          value$fit_usable
      ),
      lower_lambda_boundary_rate = mean(
        value$selected_lambda == min(config$lambda)
      ),
      upper_lambda_boundary_rate = mean(
        value$selected_lambda == max(config$lambda)
      ),
      stringsAsFactors = FALSE
    )
  }
))

checks <- list()
add_check <- function(check, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check, passed = isTRUE(passed), detail = detail,
    stringsAsFactors = FALSE
  )
}
add_check("run_manifest_md5", TRUE,
          paste(nrow(run_manifest), "payloads verified"))
add_check("run_signature", identical(run$scientific_signature,
                                      unique(run_manifest$scientific_signature)),
          run$scientific_signature)
add_check("datasets_fixed", identical(run$dataset_metadata$dataset,
                                        config$dataset_order),
          paste(run$dataset_metadata$dataset, collapse = "|"))
add_check("outer_tasks_complete", nrow(run$outer_selection) ==
            if (smoke) 4L else 110L,
          paste(nrow(run$outer_selection), "tasks"))
add_check("selection_test_blind",
          !any(run$outer_selection$outer_test_used_for_selection),
          "All outer selections used inner validation only.")
add_check("prediction_methods_complete",
          setequal(unique(run$predictions$method),
                   c(config$method_order, config$selection_label,
                     config$secondary_methods)),
          paste(sort(unique(run$predictions$method)), collapse = "|"))
add_check("prediction_values_finite",
          all(is.finite(run$predictions$observed)) &&
            all(is.finite(run$predictions$predicted)),
          paste(nrow(run$predictions), "predictions"))
add_check("full_selection_complete",
          nrow(run$full_selection) == 2L &&
            all(run$full_selection$selected_method %in% config$method_order),
          paste(run$full_selection$selected_method, collapse = "|"))
add_check("bootstrap_reselection_complete",
          nrow(run$bootstrap_selection) ==
            2L * config$bootstrap_replicates &&
            all(run$bootstrap_selection$full_reselection),
          paste(nrow(run$bootstrap_selection), "full reselections"))
add_check("stability_grid_complete",
          nrow(stability_summary) ==
            2L * config$support_grid_points,
          paste(nrow(stability_summary), "point summaries"))
add_check("uncertainty_contrasts_complete",
          nrow(prediction_contrasts) ==
            2L * (length(config$method_order) +
                    length(config$secondary_methods)),
          paste(nrow(prediction_contrasts), "contrasts"))
analysis_validation <- arfc_real_bind_rows(checks)
if (!all(analysis_validation$passed)) {
  stop("ARFC real-data analysis validation failed: ",
       paste(analysis_validation$check[!analysis_validation$passed],
             collapse = ", "), call. = FALSE)
}

analysis <- list(
  study = run$study, mode = run$mode,
  scientific_signature = run$scientific_signature,
  config = config, run_manifest = run_manifest,
  dataset_metadata = run$dataset_metadata,
  prediction_metrics = prediction_metrics,
  prediction_contrasts = prediction_contrasts,
  selection_frequency = selection_frequency,
  bootstrap_selection_frequency = bootstrap_selection_frequency,
  full_selection = run$full_selection,
  full_curves = run$full_curves,
  full_selected_curves = full_selected_curves,
  stability_summary = stability_summary,
  stable_zero_segments = stable_zero_segments,
  lambda_diagnostics = lambda_diagnostics,
  validation = analysis_validation,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

staging <- file.path(
  config$results_root,
  paste0(".", basename(analysis_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) stop("Analysis staging path exists.",
                              call. = FALSE)
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed && dir.exists(staging)) {
  unlink(staging, recursive = TRUE, force = TRUE)
}, add = TRUE)
payloads <- c(
  analysis = paste0(prefix, "_analysis.rds"),
  metrics = paste0(prefix, "_prediction_metrics.csv"),
  contrasts = paste0(prefix, "_prediction_contrasts.csv"),
  selection = paste0(prefix, "_selection_frequency.csv"),
  bootstrap_selection = paste0(
    prefix, "_bootstrap_selection_frequency.csv"
  ),
  full_selection = paste0(prefix, "_full_selection.csv"),
  full_curves = paste0(prefix, "_full_curves.csv"),
  stability = paste0(prefix, "_stability_summary.csv"),
  segments = paste0(prefix, "_stable_zero_segments.csv"),
  diagnostics = paste0(prefix, "_lambda_diagnostics.csv"),
  validation = paste0(prefix, "_analysis_validation.csv"),
  session = paste0(prefix, "_analysis_sessionInfo.txt")
)
saveRDS(analysis, file.path(staging, payloads["analysis"]),
        compress = "xz", version = 3)
write_csv <- function(object, role) utils::write.csv(
  object, file.path(staging, payloads[role]), row.names = FALSE
)
write_csv(prediction_metrics, "metrics")
write_csv(prediction_contrasts, "contrasts")
write_csv(selection_frequency, "selection")
write_csv(bootstrap_selection_frequency, "bootstrap_selection")
write_csv(run$full_selection, "full_selection")
write_csv(run$full_curves, "full_curves")
write_csv(stability_summary, "stability")
write_csv(stable_zero_segments, "segments")
write_csv(lambda_diagnostics, "diagnostics")
write_csv(analysis_validation, "validation")
writeLines(capture.output(sessionInfo()),
           file.path(staging, payloads["session"]), useBytes = TRUE)
manifest_name <- paste0(prefix, "_analysis_manifest.csv")
invisible(arfc_real_manifest(
  staging, run$scientific_signature, unname(payloads), manifest_name
))
if (!file.rename(staging, analysis_dir)) {
  stop("Atomic analysis installation failed: ", analysis_dir,
       call. = FALSE)
}
committed <- TRUE

cat("ARFC-Select real-data analysis complete.\n")
cat("Mode:", run$mode, "\n")
cat("Prediction rows:", nrow(run$predictions), "\n")
cat("Prediction contrasts:", nrow(prediction_contrasts), "\n")
cat("Stable-zero segments:", nrow(stable_zero_segments), "\n")
cat("Bundle:", analysis_dir, "\n")
