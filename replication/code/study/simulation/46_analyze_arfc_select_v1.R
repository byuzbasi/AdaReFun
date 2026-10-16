#!/usr/bin/env Rscript

# Analysis layer for the frozen ARFC four-candidate selection experiment.
# The four cases are summarized separately; no cross-case overall is formed.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/46_analyze_arfc_select_v1.R", mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
unknown <- setdiff(arguments, "--smoke")
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

study <- "arfc_select_v1"
mode <- if (smoke) "smoke" else "production"
results_root <- file.path(simulation_dir, "results", study)
run_dir <- file.path(results_root, paste0(mode, "_run"))
analysis_dir <- file.path(results_root, paste0(mode, "_analysis_v1"))
if (!dir.exists(run_dir)) {
  stop("Complete ARFC selection run bundle not found: ", run_dir,
       call. = FALSE)
}
if (dir.exists(analysis_dir) || file.exists(analysis_dir)) {
  stop("Refusing to overwrite ARFC analysis bundle: ", analysis_dir,
       call. = FALSE)
}

prefix <- study
run_file <- file.path(run_dir, paste0(prefix, "_run.rds"))
run_manifest_file <- file.path(run_dir, paste0(prefix, "_manifest.csv"))
if (!file.exists(run_file) || !file.exists(run_manifest_file)) {
  stop("Run RDS or run manifest is missing.", call. = FALSE)
}

verify_manifest <- function(directory, manifest_file) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required <- c("scientific_signature", "file", "bytes", "md5")
  if (!all(required %in% names(manifest)) || !nrow(manifest) ||
      anyDuplicated(manifest$file)) {
    stop("Manifest schema is invalid: ", manifest_file, call. = FALSE)
  }
  paths <- file.path(directory, manifest$file)
  if (!all(file.exists(paths))) {
    stop("Manifest references missing payloads.", call. = FALSE)
  }
  observed_md5 <- unname(as.character(tools::md5sum(paths)))
  observed_bytes <- as.numeric(file.info(paths)$size)
  if (!identical(observed_md5, as.character(manifest$md5)) ||
      any(observed_bytes != as.numeric(manifest$bytes))) {
    stop("Run manifest integrity verification failed.", call. = FALSE)
  }
  manifest
}

run_manifest <- verify_manifest(run_dir, run_manifest_file)
run <- readRDS(run_file)
manifest_signature <- unique(as.character(run_manifest$scientific_signature))
if (length(manifest_signature) != 1L ||
    !identical(run$study, study) || !identical(run$mode, mode) ||
    !identical(manifest_signature, run$scientific_signature)) {
  stop("Run identity, mode, or scientific signature mismatch.",
       call. = FALSE)
}

config <- run$config
tasks <- run$tasks
methods <- as.character(config$method_order)
method_metadata <- run$method_metadata
candidates <- run$candidate_metrics
selected <- run$selected_arfc_metrics
selections <- run$selection_outcomes
diagnostics <- run$diagnostics
coefficients <- run$selected_coefficients

required_candidate_columns <- c(
  "task_id", "scenario", "case_label", "role", "replication",
  "method", "penalty_family", "stage", "fit_usable", "bic_selected",
  "success", "converged", "refit_success", "validation_response_mse",
  "signal_prediction_mse", "response_prediction_mse"
)
required_selection_columns <- c(
  "task_id", "scenario", "case_label", "role", "replication",
  "selected_method", "selected_penalty_family", "selected_stage",
  "test_oracle_method", "prediction_regret", "test_oracle_agreement",
  "selection_used_test"
)
strict <- with(
  candidates,
  fit_usable & bic_selected & success & converged & refit_success &
    is.finite(bic)
)
if (!all(required_candidate_columns %in% names(candidates)) ||
    !all(required_selection_columns %in% names(selections)) ||
    nrow(tasks) != 4L * config$replications ||
    nrow(candidates) != nrow(tasks) * length(methods) ||
    nrow(selected) != nrow(tasks) || nrow(selections) != nrow(tasks) ||
    nrow(coefficients) != nrow(tasks) * length(methods) ||
    anyDuplicated(candidates[c("task_id", "method")]) ||
    anyDuplicated(coefficients[c("task_id", "method")]) ||
    anyDuplicated(selections$task_id) || !all(strict) ||
    !all(selections$selection_used_test == FALSE) ||
    any(selections$prediction_regret < -1e-12) ||
    !identical(sort(unique(as.character(candidates$method))),
               sort(methods)) ||
    any(table(tasks$scenario) != config$replications)) {
  stop("Run completeness or selected-fit audit failed.", call. = FALSE)
}

selected_keys <- paste(selected$task_id, selected$method, sep = "::")
outcome_keys <- paste(selections$task_id, selections$selected_method,
                      sep = "::")
if (!identical(selected_keys[order(selected$task_id)],
               outcome_keys[order(selections$task_id)])) {
  stop("Selected-metric rows do not match validation outcomes.",
       call. = FALSE)
}

mean_ci <- function(values, confidence = 0.95) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  count <- length(values)
  estimate <- if (count) mean(values) else NA_real_
  standard_deviation <- if (count >= 2L) stats::sd(values) else NA_real_
  standard_error <- if (count >= 2L) {
    standard_deviation / sqrt(count)
  } else NA_real_
  critical <- if (count >= 2L) {
    stats::qt((1 + confidence) / 2, df = count - 1L)
  } else NA_real_
  c(
    n = count, mean = estimate, sd = standard_deviation,
    se = standard_error,
    lower = estimate - critical * standard_error,
    upper = estimate + critical * standard_error
  )
}

wilson_ci <- function(successes, total, confidence = 0.95) {
  if (!is.finite(total) || total <= 0) return(c(lower = NA, upper = NA))
  z <- stats::qnorm((1 + confidence) / 2)
  proportion <- successes / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  radius <- z * sqrt(
    proportion * (1 - proportion) / total + z^2 / (4 * total^2)
  ) / denominator
  c(lower = max(0, center - radius), upper = min(1, center + radius))
}

scenario_order <- as.character(run$scenarios$scenario)
scenario_lookup <- run$scenarios[
  match(scenario_order, run$scenarios$scenario), , drop = FALSE
]
method_lookup <- method_metadata[
  match(methods, method_metadata$method), , drop = FALSE
]

selection_frequency <- do.call(rbind, lapply(seq_along(scenario_order),
  function(scenario_index) {
    scenario <- scenario_order[scenario_index]
    block <- selections[selections$scenario == scenario, , drop = FALSE]
    do.call(rbind, lapply(seq_along(methods), function(method_index) {
      method <- methods[method_index]
      count <- sum(block$selected_method == method)
      interval <- wilson_ci(count, nrow(block))
      data.frame(
        scenario = scenario,
        case_label = scenario_lookup$case_label[scenario_index],
        role = scenario_lookup$role[scenario_index],
        method = method,
        label = method_lookup$label[method_index],
        penalty_family = method_lookup$penalty_family[method_index],
        stage = method_lookup$stage[method_index],
        selected_count = count,
        n_tasks = nrow(block),
        selection_frequency = count / nrow(block),
        lower_95 = unname(interval["lower"]),
        upper_95 = unname(interval["upper"]),
        stringsAsFactors = FALSE
      )
    }))
  }
))
rownames(selection_frequency) <- NULL

selection_margins <- do.call(rbind, lapply(seq_along(scenario_order),
  function(scenario_index) {
    scenario <- scenario_order[scenario_index]
    block <- selections[selections$scenario == scenario, , drop = FALSE]
    dimensions <- list(
      penalty_family = c("mcp", "scad"),
      stage = c("ASR", "AMR")
    )
    do.call(rbind, lapply(names(dimensions), function(dimension) {
      column <- if (dimension == "penalty_family") {
        "selected_penalty_family"
      } else "selected_stage"
      do.call(rbind, lapply(dimensions[[dimension]], function(level) {
        count <- sum(block[[column]] == level)
        interval <- wilson_ci(count, nrow(block))
        data.frame(
          scenario = scenario,
          case_label = scenario_lookup$case_label[scenario_index],
          role = scenario_lookup$role[scenario_index],
          dimension = dimension,
          level = level,
          selected_count = count,
          n_tasks = nrow(block),
          selection_frequency = count / nrow(block),
          lower_95 = unname(interval["lower"]),
          upper_95 = unname(interval["upper"]),
          stringsAsFactors = FALSE
        )
      }))
    }))
  }
))
rownames(selection_margins) <- NULL

performance_metrics <- c(
  "validation_response_mse", "signal_prediction_mse",
  "response_prediction_mse", "ise_total", "ise_null", "ise_nonnull",
  "true_zero_rate", "false_zero_rate", "false_active_proportion",
  "null_iou", "null_symmetric_difference", "boundary_error",
  "estimated_zero_length"
)
if (!all(performance_metrics %in% names(selected))) {
  stop("Selected-performance schema is incomplete.", call. = FALSE)
}

summarize_metric_block <- function(block, metrics, extra = list()) {
  do.call(rbind, lapply(metrics, function(metric) {
    summary <- mean_ci(block[[metric]])
    data.frame(
      extra,
      metric = metric,
      n = unname(summary["n"]),
      mean = unname(summary["mean"]),
      sd = unname(summary["sd"]),
      se = unname(summary["se"]),
      lower_95 = unname(summary["lower"]),
      upper_95 = unname(summary["upper"]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
}

selected_performance <- do.call(rbind, lapply(scenario_order, function(scenario) {
  block <- selected[selected$scenario == scenario, , drop = FALSE]
  summarize_metric_block(block, performance_metrics, list(
    scenario = scenario,
    case_label = block$case_label[1L],
    role = block$role[1L],
    estimator = "ARFC-Select"
  ))
}))
rownames(selected_performance) <- NULL

candidate_performance <- do.call(rbind, lapply(
  split(candidates, interaction(candidates$scenario, candidates$method,
                                drop = TRUE)),
  function(block) {
    summarize_metric_block(block, performance_metrics, list(
      scenario = block$scenario[1L],
      case_label = block$case_label[1L],
      role = block$role[1L],
      method = block$method[1L],
      penalty_family = block$penalty_family[1L],
      stage = block$stage[1L]
    ))
  }
))
rownames(candidate_performance) <- NULL

oracle_summary <- do.call(rbind, lapply(scenario_order, function(scenario) {
  block <- selections[selections$scenario == scenario, , drop = FALSE]
  regret <- mean_ci(block$prediction_regret)
  agreement_count <- sum(block$test_oracle_agreement)
  agreement_interval <- wilson_ci(agreement_count, nrow(block))
  data.frame(
    scenario = scenario,
    case_label = block$case_label[1L],
    role = block$role[1L],
    n_tasks = nrow(block),
    oracle_agreement_count = agreement_count,
    oracle_agreement_rate = agreement_count / nrow(block),
    agreement_lower_95 = unname(agreement_interval["lower"]),
    agreement_upper_95 = unname(agreement_interval["upper"]),
    mean_prediction_regret = unname(regret["mean"]),
    regret_se = unname(regret["se"]),
    regret_lower_95 = unname(regret["lower"]),
    regret_upper_95 = unname(regret["upper"]),
    median_prediction_regret = stats::median(block$prediction_regret),
    validation_tie_rate = mean(block$validation_tie_count > 1L),
    stringsAsFactors = FALSE
  )
}))
rownames(oracle_summary) <- NULL

candidate_match <- match(candidates$task_id, selected$task_id)
if (anyNA(candidate_match)) {
  stop("Candidate-to-selected task matching failed.", call. = FALSE)
}
contrast_metrics <- c("signal_prediction_mse", "response_prediction_mse")
selected_candidate_contrasts <- do.call(rbind, lapply(
  contrast_metrics, function(metric) {
    data.frame(
      task_id = candidates$task_id,
      scenario = candidates$scenario,
      case_label = candidates$case_label,
      role = candidates$role,
      replication = candidates$replication,
      selected_method = selected$method[candidate_match],
      candidate_method = candidates$method,
      metric = metric,
      selected_value = selected[[metric]][candidate_match],
      candidate_value = candidates[[metric]],
      difference_selected_minus_candidate =
        selected[[metric]][candidate_match] - candidates[[metric]],
      stringsAsFactors = FALSE
    )
  }
))
rownames(selected_candidate_contrasts) <- NULL

selected_candidate_case_summaries <- do.call(rbind, lapply(
  split(
    selected_candidate_contrasts,
    interaction(
      selected_candidate_contrasts$scenario,
      selected_candidate_contrasts$candidate_method,
      selected_candidate_contrasts$metric,
      drop = TRUE
    )
  ),
  function(block) {
    summary <- mean_ci(block$difference_selected_minus_candidate)
    data.frame(
      scenario = block$scenario[1L],
      case_label = block$case_label[1L],
      role = block$role[1L],
      candidate_method = block$candidate_method[1L],
      metric = block$metric[1L],
      n_pairs = unname(summary["n"]),
      selected_mean = mean(block$selected_value),
      candidate_mean = mean(block$candidate_value),
      mean_difference = unname(summary["mean"]),
      se = unname(summary["se"]),
      lower_95 = unname(summary["lower"]),
      upper_95 = unname(summary["upper"]),
      stringsAsFactors = FALSE
    )
  }
))
rownames(selected_candidate_case_summaries) <- NULL

# Reconstruct all BIC-selected coefficient curves from saved coefficients.
library_override <- Sys.getenv(
  "ARFC_FMCP_LIBRARY", unset = ".replication-library"
)
if (!dir.exists(library_override)) {
  stop("Frozen FMCP library unavailable for curve reconstruction: ",
       library_override, call. = FALSE)
}
.libPaths(c(normalizePath(library_override), .libPaths()))
if (!requireNamespace("FMCP", quietly = TRUE) ||
    utils::packageVersion("FMCP") != numeric_version("0.2.0")) {
  stop("Curve reconstruction requires frozen FMCP 0.2.0.",
       call. = FALSE)
}
argvals <- seq(0, 1, length.out = config$observed_grid_points)
support_grid <- seq(0, 1, length.out = config$support_grid_points)
basis <- FMCP:::.fmcp_basis_setup(
  argvals, config$finest_M, config$degree, config$penalty_order
)
dense_basis <- unclass(splines::splineDesign(
  knots = basis$full_knots,
  x = support_grid,
  ord = basis$order,
  derivs = rep.int(0L, length(support_grid)),
  outer.ok = TRUE
))
coefficient_names <- grep("^b_[0-9]+$", names(coefficients), value = TRUE)
if (length(coefficient_names) != ncol(dense_basis)) {
  stop("Saved coefficient dimension is incompatible with the frozen basis.",
       call. = FALSE)
}

truth_function <- function(scenario, values) {
  if (identical(scenario, "case_i_null")) return(numeric(length(values)))
  if (identical(scenario, "case_ii_sparse")) {
    beta <- numeric(length(values))
    left <- values <= 0.3
    right <- values >= 0.7
    beta[left] <- 2 * (1 - values[left]) *
      sin(2 * pi * (values[left] + 0.2))
    beta[right] <- 2 * values[right] *
      sin(2 * pi * (values[right] - 0.2))
    return(beta)
  }
  if (identical(scenario, "case_iii_crossing")) {
    return(7 * values^3 + 2 * sin(4 * pi * values + 0.2))
  }
  if (identical(scenario, "case_iv_away")) {
    return(4 * sqrt(values) + exp(values^2) * cos(3 * pi * values))
  }
  stop("Unknown scenario: ", scenario, call. = FALSE)
}
true_null_mask <- function(scenario, values) {
  if (identical(scenario, "case_i_null")) return(rep(TRUE, length(values)))
  if (identical(scenario, "case_ii_sparse")) {
    return(values >= 0.3 & values <= 0.7)
  }
  rep(FALSE, length(values))
}

coefficient_matrix <- as.matrix(coefficients[coefficient_names])
storage.mode(coefficient_matrix) <- "double"
all_curves <- coefficient_matrix %*% t(dense_basis)

summarize_curve_block <- function(curve_matrix, scenario, method,
                                  label, task_count) {
  truth <- truth_function(scenario, support_grid)
  data.frame(
    scenario = scenario,
    case_label = scenario_lookup$case_label[
      match(scenario, scenario_lookup$scenario)
    ],
    role = scenario_lookup$role[match(scenario, scenario_lookup$scenario)],
    method = method,
    label = label,
    t = support_grid,
    truth = truth,
    mean_estimate = colMeans(curve_matrix),
    lower_10 = apply(curve_matrix, 2L, stats::quantile, probs = 0.10,
                     names = FALSE),
    upper_90 = apply(curve_matrix, 2L, stats::quantile, probs = 0.90,
                     names = FALSE),
    zero_selection_frequency = colMeans(
      abs(curve_matrix) <= config$zero_tolerance
    ),
    true_null = true_null_mask(scenario, support_grid),
    n = task_count,
    stringsAsFactors = FALSE
  )
}

candidate_curve_summary <- do.call(rbind, lapply(scenario_order,
  function(scenario) {
    do.call(rbind, lapply(methods, function(method) {
      rows <- which(
        coefficients$scenario == scenario & coefficients$method == method
      )
      if (length(rows) != config$replications) {
        stop("Candidate curve block is incomplete.", call. = FALSE)
      }
      summarize_curve_block(
        all_curves[rows, , drop = FALSE], scenario, method,
        method_lookup$label[match(method, method_lookup$method)], length(rows)
      )
    }))
  }
))
rownames(candidate_curve_summary) <- NULL

coefficient_keys <- paste(coefficients$task_id, coefficients$method,
                          sep = "::")
selected_coefficient_rows <- match(outcome_keys, coefficient_keys)
if (anyNA(selected_coefficient_rows)) {
  stop("Selected coefficient reconstruction failed.", call. = FALSE)
}
selected_coefficient_rows <- selected_coefficient_rows[order(selections$task_id)]
selection_order <- order(selections$task_id)
selected_curve_summary <- do.call(rbind, lapply(scenario_order,
  function(scenario) {
    selection_rows <- which(selections$scenario[selection_order] == scenario)
    rows <- selected_coefficient_rows[selection_rows]
    if (length(rows) != config$replications) {
      stop("Selected curve block is incomplete.", call. = FALSE)
    }
    summarize_curve_block(
      all_curves[rows, , drop = FALSE], scenario, "arfc_select",
      "ARFC-Select", length(rows)
    )
  }
))
rownames(selected_curve_summary) <- NULL

lambda_diagnostics <- do.call(rbind, lapply(
  split(diagnostics, interaction(diagnostics$scenario, diagnostics$method,
                                 drop = TRUE)),
  function(block) data.frame(
    scenario = block$scenario[1L],
    case_label = block$case_label[1L],
    method = block$method[1L],
    n = nrow(block),
    lower_boundary_rate = mean(block$selected_lambda_lower_boundary),
    upper_boundary_rate = mean(block$selected_lambda_upper_boundary),
    mean_usable_candidates = mean(block$usable_candidates),
    convergence_rate = sum(block$converged_candidates) /
      (nrow(block) * length(config$lambda)),
    objective_increase_candidate_rate =
      sum(block$objective_increase_candidates) /
      (nrow(block) * length(config$lambda)),
    median_elapsed_seconds = stats::median(block$path_elapsed_seconds),
    runtime_is_scientific = FALSE,
    stringsAsFactors = FALSE
  )
))
rownames(lambda_diagnostics) <- NULL

analysis_status <- data.frame(
  study = study,
  mode = mode,
  scientific_decision_allowed = !smoke,
  cross_case_overall_computed = FALSE,
  task_count = nrow(tasks),
  candidate_fit_count = nrow(candidates),
  selection_count = nrow(selections),
  interpretation = if (smoke) {
    paste(
      "Technical smoke only; all summaries are noninferential and no",
      "scientific comparison is permitted."
    )
  } else {
    paste(
      "Confirmatory per-case summaries are available; the four cases are",
      "not pooled into an overall result."
    )
  },
  stringsAsFactors = FALSE
)

analysis <- list(
  study = study,
  mode = mode,
  scientific_signature = run$scientific_signature,
  config = config,
  scenarios = run$scenarios,
  method_metadata = method_metadata,
  analysis_script_md5 = unname(as.character(tools::md5sum(script_file))),
  run_manifest_md5 = unname(as.character(tools::md5sum(run_manifest_file))),
  curve_envelope_probabilities = c(0.10, 0.90),
  selection_frequency = selection_frequency,
  selection_margins = selection_margins,
  selected_performance = selected_performance,
  candidate_performance = candidate_performance,
  oracle_summary = oracle_summary,
  selected_candidate_contrasts = selected_candidate_contrasts,
  selected_candidate_case_summaries = selected_candidate_case_summaries,
  selected_curve_summary = selected_curve_summary,
  candidate_curve_summary = candidate_curve_summary,
  lambda_diagnostics = lambda_diagnostics,
  analysis_status = analysis_status,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
staging <- file.path(
  results_root, paste0(".", basename(analysis_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) {
  stop("Analysis staging directory already exists: ", staging,
       call. = FALSE)
}
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed && dir.exists(staging)) {
  unlink(staging, recursive = TRUE, force = TRUE)
}, add = TRUE)

files <- c(
  analysis = paste0(prefix, "_analysis_v1.rds"),
  selection = paste0(prefix, "_selection_frequency_v1.csv"),
  margins = paste0(prefix, "_selection_margins_v1.csv"),
  selected_performance = paste0(prefix, "_selected_performance_v1.csv"),
  candidate_performance = paste0(prefix, "_candidate_performance_v1.csv"),
  oracle = paste0(prefix, "_oracle_summary_v1.csv"),
  contrasts = paste0(prefix, "_selected_candidate_contrasts_v1.csv"),
  contrast_summary = paste0(
    prefix, "_selected_candidate_case_summaries_v1.csv"
  ),
  selected_curves = paste0(prefix, "_selected_curve_summary_v1.csv"),
  candidate_curves = paste0(prefix, "_candidate_curve_summary_v1.csv"),
  lambda = paste0(prefix, "_lambda_diagnostics_v1.csv"),
  status = paste0(prefix, "_analysis_status_v1.csv"),
  session = paste0(prefix, "_analysis_sessionInfo_v1.txt")
)
saveRDS(analysis, file.path(staging, files[["analysis"]]),
        compress = "xz", version = 3)
table_objects <- list(
  selection_frequency, selection_margins, selected_performance,
  candidate_performance, oracle_summary, selected_candidate_contrasts,
  selected_candidate_case_summaries, selected_curve_summary,
  candidate_curve_summary, lambda_diagnostics, analysis_status
)
table_files <- unname(files[c(
  "selection", "margins", "selected_performance",
  "candidate_performance", "oracle", "contrasts", "contrast_summary",
  "selected_curves", "candidate_curves", "lambda", "status"
)])
for (index in seq_along(table_objects)) {
  utils::write.csv(
    table_objects[[index]], file.path(staging, table_files[index]),
    row.names = FALSE
  )
}
writeLines(capture.output(sessionInfo()),
           file.path(staging, files[["session"]]), useBytes = TRUE)

payload_paths <- file.path(staging, unname(files))
analysis_manifest <- data.frame(
  scientific_signature = run$scientific_signature,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  analysis_manifest,
  file.path(staging, paste0(prefix, "_analysis_manifest_v1.csv")),
  row.names = FALSE
)
if (!file.rename(staging, analysis_dir)) {
  stop("Could not atomically install ARFC analysis bundle.", call. = FALSE)
}
committed <- TRUE

cat("ARFC selection analysis complete.\n")
cat("Mode:", mode, "\n")
cat("Tasks:", nrow(tasks), "\n")
cat("Selection-frequency rows:", nrow(selection_frequency), "\n")
cat("Cross-case overall computed: FALSE\n")
cat("Scientific interpretation allowed:", !smoke, "\n")
cat("Bundle:", analysis_dir, "\n")
