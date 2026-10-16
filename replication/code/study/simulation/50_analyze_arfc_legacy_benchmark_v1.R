#!/usr/bin/env Rscript

# Analysis for the secondary ARFC classical-comparator benchmark.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/50_analyze_arfc_legacy_benchmark_v1.R",
                mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
unknown <- setdiff(arguments, "--smoke")
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}
mode <- if (smoke) "smoke" else "production"
prefix <- "arfc_legacy_benchmark_v1"
results_root <- file.path(simulation_dir, "results", prefix)
run_dir <- file.path(results_root, paste0(mode, "_run_v1"))
run_file <- file.path(run_dir, paste0(prefix, "_run_v1.rds"))
run_manifest_file <- file.path(run_dir, paste0(prefix, "_manifest_v1.csv"))
analysis_dir <- file.path(results_root, paste0(mode, "_analysis_v1"))

parent_root <- file.path(simulation_dir, "results", "arfc_select_v1")
parent_run_dir <- file.path(parent_root, paste0(mode, "_run"))
parent_run_file <- file.path(parent_run_dir, "arfc_select_v1_run.rds")
parent_manifest_file <- file.path(parent_run_dir, "arfc_select_v1_manifest.csv")
required <- c(
  run_file, run_manifest_file, parent_run_file, parent_manifest_file,
  file.path(project_dir, "config", "protocol_arfc_legacy_benchmark_v1.md"),
  script_file
)
if (!all(file.exists(required))) {
  stop("Required benchmark analysis input is missing.", call. = FALSE)
}
if (dir.exists(analysis_dir)) {
  stop("Analysis bundle already exists: ", analysis_dir, call. = FALSE)
}

verify_manifest <- function(directory, manifest_file) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required_columns <- c("scientific_signature", "file", "bytes", "md5")
  if (!all(required_columns %in% names(manifest)) ||
      anyDuplicated(manifest$file)) {
    stop("Manifest schema is invalid: ", manifest_file, call. = FALSE)
  }
  paths <- file.path(directory, manifest$file)
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == manifest$bytes) ||
      !all(unname(as.character(tools::md5sum(paths))) == manifest$md5)) {
    stop("Manifest integrity failure: ", manifest_file, call. = FALSE)
  }
  manifest
}

run_manifest <- verify_manifest(run_dir, run_manifest_file)
parent_manifest <- verify_manifest(parent_run_dir, parent_manifest_file)
run <- readRDS(run_file)
parent <- readRDS(parent_run_file)
if (!identical(run$study, prefix) || !identical(run$mode, mode) ||
    length(unique(run_manifest$scientific_signature)) != 1L ||
    !identical(unique(run_manifest$scientific_signature),
               run$scientific_signature) ||
    !identical(parent$study, "arfc_select_v1") ||
    !identical(parent$mode, mode) ||
    !identical(run$parent_signature, parent$scientific_signature) ||
    !identical(run$parent_manifest_md5,
               unname(as.character(tools::md5sum(parent_manifest_file))))) {
  stop("Run/parent identity or provenance mismatch.", call. = FALSE)
}

config <- run$config
tasks <- run$tasks
scenario_order <- as.character(run$scenarios$scenario)
reporting_order <- config$reporting_order
metric_names <- c(
  "response_prediction_mse", "signal_prediction_mse", "ise_total",
  "ise_null", "ise_nonnull", "true_zero_rate", "false_zero_rate",
  "false_active_proportion", "null_iou", "null_symmetric_difference",
  "estimated_zero_length", "boundary_error"
)
identity_names <- c(
  "task_id", "scenario", "case_label", "role", "replication",
  "seed_train", "seed_validation", "seed_test"
)

standardize_parent <- function(data, method, label, selected_component = NULL) {
  missing <- setdiff(c(identity_names, metric_names), names(data))
  if (length(missing)) {
    stop("Parent metric schema is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  result <- data[c(identity_names, metric_names)]
  result$method <- method
  result$label <- label
  result$applicable <- TRUE
  result$fit_success <- TRUE
  result$structural_sparsity <- TRUE
  result$truth_aware <- FALSE
  result$selected_component <- if (is.null(selected_component)) {
    method
  } else {
    selected_component
  }
  result
}

selected_component <- as.character(parent$selected_arfc_metrics$method)
selected_metrics <- standardize_parent(
  parent$selected_arfc_metrics,
  method = "arfc_select",
  label = unname(config$method_labels["arfc_select"]),
  selected_component = selected_component
)

candidate_order <- c("asr_mcp", "amr_mcp", "asr_scad", "amr_scad")
candidate_metrics <- do.call(rbind, lapply(candidate_order, function(method) {
  rows <- parent$candidate_metrics$method == method
  standardize_parent(
    parent$candidate_metrics[rows, , drop = FALSE], method,
    unname(config$method_labels[method])
  )
}))

new_metrics <- run$metrics
missing_new <- setdiff(
  c(identity_names, metric_names, "method", "label", "applicable",
    "fit_success", "structural_sparsity", "truth_aware"),
  names(new_metrics)
)
if (length(missing_new)) {
  stop("Benchmark metric schema is missing: ",
       paste(missing_new, collapse = ", "), call. = FALSE)
}
new_metrics <- new_metrics[c(
  identity_names, metric_names, "method", "label", "applicable",
  "fit_success", "structural_sparsity", "truth_aware"
)]
new_metrics$selected_component <- new_metrics$method

combined_metrics <- rbind(selected_metrics, candidate_metrics, new_metrics)
combined_metrics$method <- factor(
  combined_metrics$method, levels = reporting_order, ordered = TRUE
)
combined_metrics <- combined_metrics[
  order(combined_metrics$task_id, combined_metrics$method), , drop = FALSE
]
combined_metrics$method <- as.character(combined_metrics$method)
rownames(combined_metrics) <- NULL

expected_rows <- nrow(tasks) * length(reporting_order)
expected_applicable <- nrow(tasks) * (length(reporting_order) - 1L) +
  sum(tasks$scenario %in% c("case_i_null", "case_ii_sparse"))
applicable <- combined_metrics$applicable %in% TRUE
if (nrow(combined_metrics) != expected_rows ||
    anyDuplicated(combined_metrics[c("task_id", "method")]) ||
    !identical(sort(unique(combined_metrics$method)), sort(reporting_order)) ||
    sum(applicable) != expected_applicable ||
    !all(combined_metrics$fit_success[applicable] %in% TRUE) ||
    any(grepl("overall", combined_metrics$scenario, ignore.case = TRUE))) {
  stop("Combined metric completeness audit failed.", call. = FALSE)
}

fmcp_library <- Sys.getenv(
  "ARFC_FMCP_LIBRARY", unset = ".replication-library"
)
if (!dir.exists(fmcp_library)) {
  stop("Frozen FMCP library is unavailable: ", fmcp_library, call. = FALSE)
}
.libPaths(c(normalizePath(fmcp_library), .libPaths()))
if (!requireNamespace("FMCP", quietly = TRUE) ||
    utils::packageVersion("FMCP") != numeric_version("0.2.0")) {
  stop("Curve reconstruction requires frozen FMCP 0.2.0.", call. = FALSE)
}
argvals <- seq(0, 1, length.out = config$observed_grid_points)
basis <- FMCP:::.fmcp_basis_setup(
  argvals, config$M, config$degree, config$penalty_order
)
coefficient_names <- grep("^b_[0-9]+$", names(parent$selected_coefficients),
                          value = TRUE)
if (length(coefficient_names) != ncol(basis$basis_observed)) {
  stop("Parent coefficient/basis dimension mismatch.", call. = FALSE)
}

parent_curve_rows <- lapply(seq_len(nrow(parent$selected_coefficients)),
                            function(index) {
  row <- parent$selected_coefficients[index, , drop = FALSE]
  coefficient <- as.numeric(row[coefficient_names])
  beta <- as.numeric(basis$basis_observed %*% coefficient)
  identity <- row[c("task_id", "scenario", "case_label", "replication",
                    "method")]
  payload <- as.data.frame(t(beta))
  names(payload) <- sprintf("beta_%04d", seq_along(beta))
  cbind(identity, payload)
})
parent_curves <- do.call(rbind, parent_curve_rows)

selection_lookup <- parent$selection_outcomes[c("task_id", "selected_method")]
selected_curve_rows <- lapply(seq_len(nrow(selection_lookup)), function(index) {
  chosen <- parent_curves$task_id == selection_lookup$task_id[index] &
    parent_curves$method == selection_lookup$selected_method[index]
  if (sum(chosen) != 1L) {
    stop("Selected ARFC curve lookup failed for task ",
         selection_lookup$task_id[index], call. = FALSE)
  }
  row <- parent_curves[chosen, , drop = FALSE]
  row$method <- "arfc_select"
  row
})
selected_curves <- do.call(rbind, selected_curve_rows)
parent_curves <- parent_curves[
  parent_curves$method %in% candidate_order, , drop = FALSE
]

new_curves <- run$curve_estimates
curve_names <- grep("^beta_[0-9]+$", names(new_curves), value = TRUE)
if (!identical(curve_names,
               sprintf("beta_%04d", seq_along(argvals)))) {
  stop("New curve grid schema is invalid.", call. = FALSE)
}
combined_curves <- rbind(
  selected_curves[c("task_id", "scenario", "case_label", "replication",
                    "method", curve_names)],
  parent_curves[c("task_id", "scenario", "case_label", "replication",
                  "method", curve_names)],
  new_curves[c("task_id", "scenario", "case_label", "replication",
               "method", curve_names)]
)
combined_curves$method <- factor(
  combined_curves$method, levels = reporting_order, ordered = TRUE
)
combined_curves <- combined_curves[
  order(combined_curves$task_id, combined_curves$method), , drop = FALSE
]
combined_curves$method <- as.character(combined_curves$method)
rownames(combined_curves) <- NULL
if (nrow(combined_curves) != expected_rows ||
    anyDuplicated(combined_curves[c("task_id", "method")])) {
  stop("Combined curve completeness audit failed.", call. = FALSE)
}

mean_ci <- function(values, level = 0.95) {
  values <- values[is.finite(values)]
  n <- length(values)
  if (!n) {
    return(c(n = 0, mean = NA, sd = NA, se = NA, lower = NA, upper = NA))
  }
  estimate <- mean(values)
  standard_deviation <- if (n > 1L) stats::sd(values) else NA_real_
  standard_error <- if (n > 1L) standard_deviation / sqrt(n) else NA_real_
  critical <- if (n > 1L) stats::qt((1 + level) / 2, df = n - 1L) else NA_real_
  c(
    n = n, mean = estimate, sd = standard_deviation, se = standard_error,
    lower = if (n > 1L) estimate - critical * standard_error else NA_real_,
    upper = if (n > 1L) estimate + critical * standard_error else NA_real_
  )
}

summary_rows <- list()
cursor <- 1L
for (scenario in scenario_order) {
  for (method in reporting_order) {
    rows <- combined_metrics$scenario == scenario &
      combined_metrics$method == method &
      combined_metrics$applicable %in% TRUE
    for (metric in metric_names) {
      values <- combined_metrics[[metric]][rows]
      summary_rows[[cursor]] <- data.frame(
        scenario = scenario,
        case_label = unique(combined_metrics$case_label[
          combined_metrics$scenario == scenario
        ])[1L],
        method = method,
        label = unname(config$method_labels[method]),
        metric = metric,
        as.list(mean_ci(values)),
        stringsAsFactors = FALSE
      )
      cursor <- cursor + 1L
    }
  }
}
method_case_summary <- do.call(rbind, summary_rows)

summary_value <- function(scenario, method, metric, column = "mean") {
  rows <- method_case_summary$scenario == scenario &
    method_case_summary$method == method &
    method_case_summary$metric == metric
  if (sum(rows) != 1L) return(NA_real_)
  method_case_summary[[column]][rows]
}

case_ii <- "case_ii_sparse"
case_ii_performance <- do.call(rbind, lapply(reporting_order, function(method) {
  data.frame(
    method = method,
    label = unname(config$method_labels[method]),
    n = summary_value(case_ii, method, "response_prediction_mse", "n"),
    response_prediction_mse = summary_value(
      case_ii, method, "response_prediction_mse"
    ),
    ise_null = summary_value(case_ii, method, "ise_null"),
    ise_nonnull = summary_value(case_ii, method, "ise_nonnull"),
    true_zero_rate = summary_value(case_ii, method, "true_zero_rate"),
    false_zero_rate = summary_value(case_ii, method, "false_zero_rate"),
    null_iou = summary_value(case_ii, method, "null_iou"),
    boundary_error = summary_value(case_ii, method, "boundary_error"),
    stringsAsFactors = FALSE
  )
}))

safety_summary <- do.call(rbind, lapply(scenario_order, function(scenario) {
  do.call(rbind, lapply(reporting_order, function(method) {
    safety_metric <- if (identical(scenario, "case_i_null")) {
      "false_active_proportion"
    } else {
      "false_zero_rate"
    }
    data.frame(
      scenario = scenario,
      case_label = unique(combined_metrics$case_label[
        combined_metrics$scenario == scenario
      ])[1L],
      method = method,
      label = unname(config$method_labels[method]),
      n = summary_value(scenario, method, "response_prediction_mse", "n"),
      response_prediction_mse = summary_value(
        scenario, method, "response_prediction_mse"
      ),
      safety_metric = safety_metric,
      safety_value = summary_value(scenario, method, safety_metric),
      stringsAsFactors = FALSE
    )
  }))
}))

contrast_pairs <- data.frame(
  contrast = c(
    "arfc_select_minus_fixed_fmcp", "arfc_select_minus_fixed_fscad",
    "arfc_select_minus_smooth", "arfc_select_minus_flirti",
    "fixed_fmcp_minus_fixed_fscad"
  ),
  left = c(
    "arfc_select", "arfc_select", "arfc_select", "arfc_select",
    "fixed_fmcp"
  ),
  right = c("fixed_fmcp", "fixed_fscad", "smooth", "flirti", "fixed_fscad"),
  stringsAsFactors = FALSE
)
contrast_metrics <- c(
  "response_prediction_mse", "ise_total", "ise_null", "ise_nonnull",
  "true_zero_rate", "false_zero_rate", "false_active_proportion"
)
paired_rows <- list()
cursor <- 1L
for (scenario in scenario_order) {
  for (pair_index in seq_len(nrow(contrast_pairs))) {
    pair <- contrast_pairs[pair_index, ]
    left_rows <- combined_metrics$scenario == scenario &
      combined_metrics$method == pair$left &
      combined_metrics$applicable %in% TRUE
    right_rows <- combined_metrics$scenario == scenario &
      combined_metrics$method == pair$right &
      combined_metrics$applicable %in% TRUE
    left <- combined_metrics[left_rows, c("task_id", contrast_metrics)]
    right <- combined_metrics[right_rows, c("task_id", contrast_metrics)]
    matched <- merge(left, right, by = "task_id", suffixes = c("_left", "_right"))
    for (metric in contrast_metrics) {
      difference <- matched[[paste0(metric, "_left")]] -
        matched[[paste0(metric, "_right")]]
      summary <- mean_ci(difference)
      paired_rows[[cursor]] <- data.frame(
        scenario = scenario,
        case_label = unique(combined_metrics$case_label[
          combined_metrics$scenario == scenario
        ])[1L],
        contrast = pair$contrast, left = pair$left, right = pair$right,
        direction = paste0(pair$left, " - ", pair$right),
        metric = metric, as.list(summary), stringsAsFactors = FALSE
      )
      cursor <- cursor + 1L
    }
  }
}
paired_contrasts <- do.call(rbind, paired_rows)

case_ii_curves <- combined_curves[
  combined_curves$scenario == case_ii, , drop = FALSE
]
structural_lookup <- setNames(
  vapply(reporting_order, function(method) {
    values <- unique(combined_metrics$structural_sparsity[
      combined_metrics$method == method
    ])
    any(values %in% TRUE)
  }, logical(1L)), reporting_order
)
truth_case_ii <- {
  beta <- numeric(length(argvals))
  left <- argvals <= 0.3
  right <- argvals >= 0.7
  beta[left] <- 2 * (1 - argvals[left]) *
    sin(2 * pi * (argvals[left] + 0.2))
  beta[right] <- 2 * argvals[right] *
    sin(2 * pi * (argvals[right] - 0.2))
  beta
}
curve_summary <- do.call(rbind, lapply(reporting_order, function(method) {
  rows <- case_ii_curves$method == method
  matrix_values <- as.matrix(case_ii_curves[rows, curve_names, drop = FALSE])
  finite_rows <- apply(matrix_values, 1L, function(values) all(is.finite(values)))
  matrix_values <- matrix_values[finite_rows, , drop = FALSE]
  if (!nrow(matrix_values)) {
    means <- lower <- upper <- zero_frequency <- rep(NA_real_, length(argvals))
  } else {
    means <- colMeans(matrix_values)
    lower <- apply(matrix_values, 2L, stats::quantile, probs = 0.05,
                   names = FALSE, type = 8)
    upper <- apply(matrix_values, 2L, stats::quantile, probs = 0.95,
                   names = FALSE, type = 8)
    zero_frequency <- if (structural_lookup[method]) {
      colMeans(abs(matrix_values) <= config$zero_tolerance)
    } else {
      rep(NA_real_, length(argvals))
    }
  }
  data.frame(
    scenario = case_ii, case_label = "Case II",
    method = method, label = unname(config$method_labels[method]),
    t = argvals, truth = truth_case_ii,
    mean = means, q05 = lower, q95 = upper,
    zero_frequency = zero_frequency,
    structural_sparsity = structural_lookup[method],
    replications = nrow(matrix_values), stringsAsFactors = FALSE
  )
}))

analysis_validation <- data.frame(
  check = c(
    "run_manifest_integrity", "parent_manifest_integrity",
    "parent_signature_link", "combined_method_order",
    "combined_task_method_completeness", "applicable_fit_completeness",
    "curve_task_method_completeness", "case_specific_only",
    "case_ii_curve_replication_count", "smoke_noninferential_label"
  ),
  passed = c(
    TRUE, TRUE,
    identical(run$parent_signature, parent$scientific_signature),
    identical(sort(unique(combined_metrics$method)), sort(reporting_order)),
    nrow(combined_metrics) == expected_rows &&
      !anyDuplicated(combined_metrics[c("task_id", "method")]),
    all(combined_metrics$fit_success[applicable] %in% TRUE),
    nrow(combined_curves) == expected_rows &&
      !anyDuplicated(combined_curves[c("task_id", "method")]),
    !any(grepl("overall", c(
      method_case_summary$scenario, safety_summary$scenario,
      paired_contrasts$scenario
    ), ignore.case = TRUE)),
    all(curve_summary$replications == config$replications),
    if (smoke) config$mode == "smoke" else TRUE
  ),
  details = c(
    paste(nrow(run_manifest), "payloads"),
    paste(nrow(parent_manifest), "parent payloads"),
    run$parent_signature,
    paste(reporting_order, collapse = ","),
    paste(nrow(combined_metrics), "rows"),
    paste(sum(applicable), "applicable rows"),
    paste(nrow(combined_curves), "curve rows"),
    "no cross-case pooling",
    paste(config$replications, "Case-II curves per method"),
    if (smoke) "noninferential smoke" else "production"
  ),
  stringsAsFactors = FALSE
)
if (!all(analysis_validation$passed)) {
  stop("Analysis validation failed: ", paste(
    analysis_validation$check[!analysis_validation$passed], collapse = ", "
  ), call. = FALSE)
}

analysis_status <- data.frame(
  mode = mode,
  scientific_signature = run$scientific_signature,
  parent_signature = parent$scientific_signature,
  tasks = nrow(tasks),
  reporting_methods = length(reporting_order),
  applicable_task_method_rows = sum(applicable),
  inferential = !smoke,
  interpretation = if (smoke) {
    "Technical smoke only; no scientific comparison is permitted."
  } else {
    paste(
      "Secondary case-specific benchmark; ARFC-Select candidate library and",
      "primary decision remain unchanged."
    )
  },
  stringsAsFactors = FALSE
)

analysis <- list(
  study = prefix, mode = mode,
  scientific_signature = run$scientific_signature,
  parent_signature = parent$scientific_signature,
  config = config, tasks = tasks, scenarios = run$scenarios,
  combined_metrics = combined_metrics,
  combined_curves = combined_curves,
  method_case_summary = method_case_summary,
  case_ii_performance = case_ii_performance,
  safety_summary = safety_summary,
  paired_contrasts = paired_contrasts,
  curve_summary = curve_summary,
  status = analysis_status,
  validation = analysis_validation,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

staging <- file.path(
  results_root, paste0(".", basename(analysis_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) stop("Staging directory exists.", call. = FALSE)
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
payloads <- c(
  analysis = paste0(prefix, "_analysis_v1.rds"),
  combined_metrics = paste0(prefix, "_combined_metrics_v1.csv"),
  method_case = paste0(prefix, "_method_case_summary_v1.csv"),
  case_ii = paste0(prefix, "_case_ii_performance_v1.csv"),
  safety = paste0(prefix, "_safety_summary_v1.csv"),
  contrasts = paste0(prefix, "_paired_contrasts_v1.csv"),
  curves = paste0(prefix, "_curve_summary_v1.csv"),
  status = paste0(prefix, "_analysis_status_v1.csv"),
  validation = paste0(prefix, "_analysis_validation_v1.csv"),
  session = paste0(prefix, "_analysis_sessionInfo_v1.txt")
)
saveRDS(analysis, file.path(staging, payloads["analysis"]),
        compress = "xz", version = 3)
utils::write.csv(combined_metrics,
                 file.path(staging, payloads["combined_metrics"]),
                 row.names = FALSE)
utils::write.csv(method_case_summary,
                 file.path(staging, payloads["method_case"]), row.names = FALSE)
utils::write.csv(case_ii_performance,
                 file.path(staging, payloads["case_ii"]), row.names = FALSE)
utils::write.csv(safety_summary,
                 file.path(staging, payloads["safety"]), row.names = FALSE)
utils::write.csv(paired_contrasts,
                 file.path(staging, payloads["contrasts"]), row.names = FALSE)
utils::write.csv(curve_summary,
                 file.path(staging, payloads["curves"]), row.names = FALSE)
utils::write.csv(analysis_status,
                 file.path(staging, payloads["status"]), row.names = FALSE)
utils::write.csv(analysis_validation,
                 file.path(staging, payloads["validation"]), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(staging, payloads["session"]))
payload_paths <- file.path(staging, unname(payloads))
manifest <- data.frame(
  scientific_signature = run$scientific_signature,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest,
  file.path(staging, paste0(prefix, "_analysis_manifest_v1.csv")),
  row.names = FALSE
)
if (!file.rename(staging, analysis_dir)) {
  stop("Atomic analysis bundle install failed: ", analysis_dir,
       call. = FALSE)
}
cat("ARFC classical-comparator analysis complete.\n")
cat("Combined task-method rows:", nrow(combined_metrics), "\n")
cat("Applicable rows:", sum(applicable), "\n")
cat("Analysis checks:", sum(analysis_validation$passed), "/",
    nrow(analysis_validation), "PASS\n")
cat("Bundle:", analysis_dir, "\n")
