#!/usr/bin/env Rscript

# Independent validation for the ARFC classical-comparator benchmark.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/52_validate_arfc_legacy_benchmark_v1.R",
                mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
core_only <- "--core-only" %in% arguments
report_argument <- grep("^--report-version=", arguments, value = TRUE)
if (length(report_argument) > 1L) {
  stop("Specify --report-version at most once.", call. = FALSE)
}
report_version <- if (length(report_argument)) {
  sub("^--report-version=", "", report_argument)
} else {
  "v1"
}
if (!report_version %in% c("v1", "v2")) {
  stop("--report-version must be v1 or v2.", call. = FALSE)
}
unknown <- setdiff(arguments, c("--smoke", "--core-only", report_argument))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}
mode <- if (smoke) "smoke" else "production"
prefix <- "arfc_legacy_benchmark_v1"
results_root <- file.path(simulation_dir, "results", prefix)
run_dir <- file.path(results_root, paste0(mode, "_run_v1"))
analysis_dir <- file.path(results_root, paste0(mode, "_analysis_v1"))
validation_dir <- file.path(
  results_root,
  paste0(mode, if (core_only) "_validation_core_" else "_validation_",
         report_version)
)
run_file <- file.path(run_dir, paste0(prefix, "_run_v1.rds"))
run_manifest_file <- file.path(run_dir, paste0(prefix, "_manifest_v1.csv"))
analysis_file <- file.path(analysis_dir, paste0(prefix, "_analysis_v1.rds"))
analysis_manifest_file <- file.path(
  analysis_dir, paste0(prefix, "_analysis_manifest_v1.csv")
)
parent_dir <- file.path(
  simulation_dir, "results", "arfc_select_v1", paste0(mode, "_run")
)
parent_file <- file.path(parent_dir, "arfc_select_v1_run.rds")
parent_manifest_file <- file.path(parent_dir, "arfc_select_v1_manifest.csv")
reporting_dir <- file.path(
  project_dir, "config", "generated", prefix,
  paste0(mode, "_reporting_", report_version)
)
report_manifest_file <- file.path(
  reporting_dir,
  paste0(prefix, "_reporting_manifest_", report_version,
         if (smoke) "_smoke" else "", ".csv")
)
pdf_dir <- file.path(
  project_dir, "output", "pdf", prefix,
  paste0(mode, "_reporting_", report_version)
)
pdf_manifest_file <- file.path(
  pdf_dir,
  paste0(prefix, "_pdf_manifest_", report_version,
         if (smoke) "_smoke" else "", ".csv")
)

required_core <- c(
  run_file, run_manifest_file, analysis_file, analysis_manifest_file,
  parent_file, parent_manifest_file,
  file.path(project_dir, "config", "protocol_arfc_legacy_benchmark_v1.md"),
  script_file
)
if (!all(file.exists(required_core))) {
  stop("A core validation input is missing.", call. = FALSE)
}
if (!core_only && (!file.exists(report_manifest_file) ||
                   !file.exists(pdf_manifest_file))) {
  stop("Publication artifacts are missing; use --core-only only for an ",
       "intentional pre-publication validation.", call. = FALSE)
}
if (dir.exists(validation_dir)) {
  stop("Validation bundle already exists: ", validation_dir, call. = FALSE)
}

verify_manifest <- function(directory, manifest_file) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required <- c("scientific_signature", "file", "bytes", "md5")
  if (!all(required %in% names(manifest)) || anyDuplicated(manifest$file)) {
    return(list(passed = FALSE, manifest = manifest,
                details = "invalid manifest schema"))
  }
  paths <- file.path(directory, manifest$file)
  passed <- all(file.exists(paths)) &&
    all(as.numeric(file.info(paths)$size) == manifest$bytes) &&
    all(unname(as.character(tools::md5sum(paths))) == manifest$md5)
  list(
    passed = passed, manifest = manifest,
    details = paste(nrow(manifest), "payloads")
  )
}

run_manifest_check <- verify_manifest(run_dir, run_manifest_file)
analysis_manifest_check <- verify_manifest(
  analysis_dir, analysis_manifest_file
)
parent_manifest_check <- verify_manifest(parent_dir, parent_manifest_file)
if (!run_manifest_check$passed || !analysis_manifest_check$passed ||
    !parent_manifest_check$passed) {
  stop("A required input manifest failed before validation.", call. = FALSE)
}
run <- readRDS(run_file)
analysis <- readRDS(analysis_file)
parent <- readRDS(parent_file)
config <- run$config
tasks <- run$tasks

checks <- list()
add_check <- function(check, passed, details, applicable = TRUE) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    applicable = applicable,
    passed = if (applicable) isTRUE(passed) else NA,
    details = as.character(details),
    stringsAsFactors = FALSE
  )
}

add_check("run_manifest_integrity", run_manifest_check$passed,
          run_manifest_check$details)
add_check("analysis_manifest_integrity", analysis_manifest_check$passed,
          analysis_manifest_check$details)
add_check("parent_manifest_integrity", parent_manifest_check$passed,
          parent_manifest_check$details)
add_check(
  "run_identity_signature",
  identical(run$study, prefix) && identical(run$mode, mode) &&
    identical(unique(run_manifest_check$manifest$scientific_signature),
              run$scientific_signature),
  run$scientific_signature
)
add_check(
  "analysis_identity_signature",
  identical(analysis$study, prefix) && identical(analysis$mode, mode) &&
    identical(analysis$scientific_signature, run$scientific_signature) &&
    identical(unique(analysis_manifest_check$manifest$scientific_signature),
              run$scientific_signature),
  analysis$scientific_signature
)
add_check(
  "parent_identity_signature",
  identical(parent$study, "arfc_select_v1") &&
    identical(parent$mode, mode) &&
    identical(run$parent_signature, parent$scientific_signature) &&
    identical(analysis$parent_signature, parent$scientific_signature),
  parent$scientific_signature
)
add_check(
  "parent_manifest_hash_link",
  identical(
    run$parent_manifest_md5,
    unname(as.character(tools::md5sum(parent_manifest_file)))
  ),
  run$parent_manifest_md5
)

expected_methods <- c(
  "arfc_select", "asr_mcp", "amr_mcp", "asr_scad", "amr_scad",
  "fixed_fmcp", "fixed_fscad", "smooth", "ols", "pcr", "flirti",
  "oracle"
)
expected_new <- c(
  "fixed_fmcp", "fixed_fscad", "smooth", "ols", "pcr", "flirti",
  "oracle"
)
expected_tasks <- 4L * config$replications
add_check(
  "task_count_and_unique_seeds",
  nrow(tasks) == expected_tasks && !anyDuplicated(tasks$task_id) &&
    !anyDuplicated(unlist(tasks[c(
      "seed_train", "seed_validation", "seed_test"
    )], use.names = FALSE)),
  paste(nrow(tasks), "tasks")
)
add_check(
  "frozen_reporting_method_order",
  identical(config$reporting_order, expected_methods),
  paste(config$reporting_order, collapse = ",")
)
add_check(
  "frozen_new_method_order",
  identical(config$new_method_order, expected_new),
  paste(config$new_method_order, collapse = ",")
)
add_check(
  "equal_scalar_maximum_concavity",
  abs(1 / config$gamma - 1 / (config$scad_a - 1)) <= 1e-14,
  paste("MCP", 1 / config$gamma, "SCAD", 1 / (config$scad_a - 1))
)
add_check(
  "frozen_package_versions",
  identical(unname(config$package_versions[c("FMCP", "PFLR")]),
            c("0.2.0", "1.1.0")),
  paste(names(config$package_versions), config$package_versions,
        collapse = "; ")
)

expected_new_rows <- nrow(tasks) * length(expected_new)
expected_combined_rows <- nrow(tasks) * length(expected_methods)
expected_applicable_new <- nrow(tasks) * 6L +
  sum(tasks$scenario %in% c("case_i_null", "case_ii_sparse"))
expected_applicable_combined <- nrow(tasks) * 11L +
  sum(tasks$scenario %in% c("case_i_null", "case_ii_sparse"))
run_applicable <- run$metrics$applicable %in% TRUE
combined_applicable <- analysis$combined_metrics$applicable %in% TRUE
add_check(
  "run_task_method_completeness",
  nrow(run$metrics) == expected_new_rows &&
    !anyDuplicated(run$metrics[c("task_id", "method")]) &&
    sum(run_applicable) == expected_applicable_new,
  paste(nrow(run$metrics), "rows;", sum(run_applicable), "applicable")
)
add_check(
  "combined_task_method_completeness",
  nrow(analysis$combined_metrics) == expected_combined_rows &&
    !anyDuplicated(analysis$combined_metrics[c("task_id", "method")]) &&
    sum(combined_applicable) == expected_applicable_combined,
  paste(nrow(analysis$combined_metrics), "rows;",
        sum(combined_applicable), "applicable")
)
add_check(
  "applicable_fit_success",
  all(run$metrics$fit_success[run_applicable] %in% TRUE) &&
    all(analysis$combined_metrics$fit_success[combined_applicable] %in% TRUE),
  paste(sum(run_applicable), "new applicable fits")
)
oracle_expected <- run$metrics$scenario %in%
  c("case_i_null", "case_ii_sparse")
oracle_rows <- run$metrics$method == "oracle"
add_check(
  "oracle_applicability_contract",
  identical(run$metrics$applicable[oracle_rows], oracle_expected[oracle_rows]) &&
    all(is.na(run$metrics$fit_success[
      oracle_rows & !oracle_expected
    ])),
  "Oracle only in Cases I-II"
)
add_check(
  "no_cross_case_overall_pooling",
  !any(grepl("overall", c(
    run$metrics$scenario, analysis$method_case_summary$scenario,
    analysis$safety_summary$scenario, analysis$paired_contrasts$scenario
  ), ignore.case = TRUE)),
  "case-specific summaries only"
)

strict_methods <- c("fixed_fmcp", "fixed_fscad", "smooth")
strict_failures <- character()
for (task_id in tasks$task_id) {
  for (method in strict_methods) {
    rows <- run$tuning$task_id == task_id & run$tuning$method == method
    grid <- run$tuning[rows, , drop = FALSE]
    eligible <- which(
      grid$success %in% TRUE & grid$converged %in% TRUE &
        grid$refit_success %in% TRUE & is.finite(grid$bic)
    )
    selected <- which(grid$selected %in% TRUE)
    metric <- run$metrics[
      run$metrics$task_id == task_id & run$metrics$method == method,
      , drop = FALSE
    ]
    expected <- if (length(eligible)) eligible[which.min(grid$bic[eligible])] else {
      NA_integer_
    }
    passed <- nrow(metric) == 1L && length(selected) == 1L &&
      is.finite(expected) && selected == expected &&
      isTRUE(all.equal(metric$lambda, grid$lambda[selected],
                       tolerance = 0)) &&
      isTRUE(all.equal(metric$kappa, grid$kappa[selected], tolerance = 0))
    if (!passed) strict_failures <- c(strict_failures, paste(task_id, method))
  }
}
add_check(
  "strict_grid_selection_reconstruction",
  !length(strict_failures),
  if (length(strict_failures)) paste(strict_failures, collapse = "; ") else {
    paste(length(tasks$task_id) * length(strict_methods), "selections")
  }
)

path_mismatch <- numeric(nrow(tasks))
kappa_mismatch <- numeric(nrow(tasks))
for (index in seq_len(nrow(tasks))) {
  task_id <- tasks$task_id[index]
  fmcp_path <- run$tuning[
    run$tuning$task_id == task_id & run$tuning$method == "fixed_fmcp",
    c("lambda", "kappa")
  ]
  fscad_path <- run$tuning[
    run$tuning$task_id == task_id & run$tuning$method == "fixed_fscad",
    c("lambda", "kappa")
  ]
  path_mismatch[index] <- if (isTRUE(all.equal(
    unname(as.matrix(fmcp_path)), unname(as.matrix(fscad_path)),
    tolerance = 0
  ))) 0 else Inf
  selected <- run$metrics[
    run$metrics$task_id == task_id &
      run$metrics$method %in% c("fixed_fmcp", "fixed_fscad"),
    , drop = FALSE
  ]
  kappa_mismatch[index] <- diff(range(selected$kappa))
}
add_check(
  "fmcp_fscad_common_grid_and_kappa",
  max(path_mismatch) == 0 && max(kappa_mismatch) <= 0,
  paste("max kappa difference", max(kappa_mismatch))
)
add_check(
  "parent_data_tuning_provenance",
  max(run$diagnostics$data_provenance_max_difference) <= 1e-10 &&
    max(abs(run$diagnostics$cutoff -
              run$diagnostics$reconstructed_cutoff)) <= 1e-10,
  paste("max difference", max(run$diagnostics$data_provenance_max_difference))
)

fmcp_library <- Sys.getenv(
  "ARFC_FMCP_LIBRARY", unset = ".replication-library"
)
if (!dir.exists(fmcp_library)) {
  stop("Frozen FMCP library unavailable for independent reconstruction.",
       call. = FALSE)
}
.libPaths(c(normalizePath(fmcp_library), .libPaths()))
if (!requireNamespace("FMCP", quietly = TRUE) ||
    utils::packageVersion("FMCP") != numeric_version("0.2.0")) {
  stop("Independent reconstruction requires FMCP 0.2.0.", call. = FALSE)
}
argvals <- seq(0, 1, length.out = config$observed_grid_points)
support_grid <- seq(0, 1, length.out = config$support_grid_points)
integration_weights <- FMCP:::.fmcp_integration_weights(
  argvals, config$integration
)
support_gap <- diff(support_grid)
support_weights <- c(
  support_gap[1L] / 2,
  (support_gap[-length(support_gap)] + support_gap[-1L]) / 2,
  support_gap[length(support_gap)] / 2
)
coefficient_basis <- FMCP:::.fmcp_basis_setup(
  argvals, config$M, config$degree, config$penalty_order
)
dense_basis <- unclass(splines::splineDesign(
  knots = coefficient_basis$full_knots, x = support_grid,
  ord = coefficient_basis$order,
  derivs = rep.int(0L, length(support_grid)), outer.ok = TRUE
))
predictor_basis <- local({
  knots <- seq(0, 1, length.out = config$predictor_knots)
  unclass(splines::bs(
    argvals, knots = knots[-c(1L, length(knots))],
    degree = config$predictor_order - 1L, intercept = TRUE,
    Boundary.knots = range(argvals)
  ))
})
generate_predictors <- function(sample_size) {
  scores <- matrix(
    stats::rnorm(sample_size * ncol(predictor_basis)),
    nrow = sample_size, ncol = ncol(predictor_basis)
  )
  scores %*% t(predictor_basis)
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
  stop("Unknown scenario.", call. = FALSE)
}
true_null_mask <- function(scenario, values) {
  if (identical(scenario, "case_i_null")) return(rep(TRUE, length(values)))
  if (identical(scenario, "case_ii_sparse")) {
    return(values >= 0.3 & values <= 0.7)
  }
  rep(FALSE, length(values))
}
generate_sample <- function(sample_size, seed, truth_observed, noise_sd) {
  set.seed(as.integer(seed))
  x <- generate_predictors(sample_size)
  linear <- as.numeric(x %*% (integration_weights * truth_observed))
  signal <- config$intercept + linear
  noise <- stats::rnorm(sample_size, sd = noise_sd)
  list(y = signal + noise, signal = signal)
}
generate_task_response <- function(task) {
  truth_observed <- truth_function(task$scenario, argvals)
  set.seed(as.integer(task$seed_train))
  x_train <- generate_predictors(task$n_train)
  train_linear <- as.numeric(
    x_train %*% (integration_weights * truth_observed)
  )
  signal_variance <- stats::var(train_linear)
  noise_sd <- if (identical(task$scenario, "case_i_null")) {
    config$case_i_noise_sd
  } else {
    sqrt(signal_variance / config$snr)
  }
  invisible(stats::rnorm(task$n_train, sd = noise_sd))
  list(
    validation = generate_sample(
      task$n_validation, task$seed_validation, truth_observed, noise_sd
    ),
    test = generate_sample(
      task$n_test, task$seed_test, truth_observed, noise_sd
    ),
    truth_dense = truth_function(task$scenario, support_grid)
  )
}
weighted_integral <- function(values, mask = rep(TRUE, length(values))) {
  if (!any(mask) || any(!is.finite(values[mask]))) return(NA_real_)
  sum(support_weights[mask] * values[mask])
}
support_metrics <- function(beta, truth, scenario, structural) {
  truth_null <- true_null_mask(scenario, support_grid)
  truth_nonnull <- !truth_null & abs(truth) > 1e-10
  squared <- (beta - truth)^2
  result <- c(
    ise_total = weighted_integral(squared),
    ise_null = weighted_integral(squared, truth_null),
    ise_nonnull = weighted_integral(squared, truth_nonnull)
  )
  if (!structural) return(result)
  zero <- abs(beta) <= config$zero_tolerance
  null_length <- sum(support_weights[truth_null])
  nonnull_length <- sum(support_weights[truth_nonnull])
  intersection <- sum(support_weights[zero & truth_null])
  false_length <- sum(support_weights[zero & truth_nonnull])
  c(
    result,
    true_zero_rate = if (null_length > 0) intersection / null_length else NA,
    false_zero_rate = if (nonnull_length > 0) false_length / nonnull_length else NA,
    false_active_proportion = if (identical(scenario, "case_i_null")) {
      sum(support_weights[!zero])
    } else NA_real_
  )
}

validation_columns <- grep("^validation_[0-9]+$",
                           names(run$validation_predictions), value = TRUE)
test_columns <- grep("^test_[0-9]+$", names(run$test_predictions),
                     value = TRUE)
curve_columns <- grep("^beta_[0-9]+$", names(run$curve_estimates),
                      value = TRUE)
add_check(
  "stored_vector_dimensions",
  length(validation_columns) == config$n_validation &&
    length(test_columns) == config$n_test &&
    length(curve_columns) == config$observed_grid_points,
  paste(length(validation_columns), length(test_columns), length(curve_columns))
)

prediction_differences <- numeric()
structural_differences <- numeric()
na_contract_passed <- TRUE
coefficient_curve_methods <- c(
  "fixed_fmcp", "fixed_fscad", "smooth", "ols"
)
for (task_index in seq_len(nrow(tasks))) {
  task <- as.list(tasks[task_index, , drop = FALSE])
  response <- generate_task_response(task)
  for (method in expected_new) {
    metric <- run$metrics[
      run$metrics$task_id == task$task_id & run$metrics$method == method,
      , drop = FALSE
    ]
    if (nrow(metric) != 1L || !isTRUE(metric$applicable)) next
    validation_row <- run$validation_predictions[
      run$validation_predictions$task_id == task$task_id &
        run$validation_predictions$method == method, , drop = FALSE
    ]
    test_row <- run$test_predictions[
      run$test_predictions$task_id == task$task_id &
        run$test_predictions$method == method, , drop = FALSE
    ]
    validation_prediction <- as.numeric(validation_row[validation_columns])
    test_prediction <- as.numeric(test_row[test_columns])
    recomputed_prediction <- c(
      validation_signal_mse = mean(
        (validation_prediction - response$validation$signal)^2
      ),
      validation_response_mse = mean(
        (validation_prediction - response$validation$y)^2
      ),
      signal_prediction_mse = mean(
        (test_prediction - response$test$signal)^2
      ),
      response_prediction_mse = mean(
        (test_prediction - response$test$y)^2
      )
    )
    prediction_differences <- c(
      prediction_differences,
      abs(recomputed_prediction - as.numeric(metric[names(recomputed_prediction)]))
    )

    curve_row <- run$curve_estimates[
      run$curve_estimates$task_id == task$task_id &
        run$curve_estimates$method == method, , drop = FALSE
    ]
    beta_observed <- as.numeric(curve_row[curve_columns])
    beta_dense <- if (method %in% coefficient_curve_methods) {
      coefficient <- qr.solve(
        coefficient_basis$basis_observed, beta_observed, tol = 1e-12
      )
      as.numeric(dense_basis %*% coefficient)
    } else {
      as.numeric(stats::approx(
        argvals, beta_observed, xout = support_grid,
        rule = 2, ties = "ordered"
      )$y)
    }
    recomputed_structure <- support_metrics(
      beta_dense, response$truth_dense, task$scenario,
      structural = isTRUE(metric$structural_sparsity)
    )
    observed_structure <- as.numeric(metric[names(recomputed_structure)])
    finite <- is.finite(recomputed_structure) & is.finite(observed_structure)
    structural_differences <- c(
      structural_differences,
      abs(recomputed_structure[finite] - observed_structure[finite])
    )
    na_contract_passed <- na_contract_passed && identical(
      unname(is.na(recomputed_structure)), unname(is.na(observed_structure))
    )
  }
}
max_prediction_difference <- if (length(prediction_differences)) {
  max(prediction_differences)
} else Inf
max_structural_difference <- if (length(structural_differences)) {
  max(structural_differences)
} else Inf
add_check(
  "prediction_metrics_independent_reconstruction",
  max_prediction_difference <= 1e-12,
  paste("max absolute difference", format(max_prediction_difference,
                                          scientific = TRUE))
)
add_check(
  "structural_metrics_independent_reconstruction",
  max_structural_difference <= 1e-8 && na_contract_passed,
  paste("max absolute difference", format(max_structural_difference,
                                          scientific = TRUE),
        "; NA contract", na_contract_passed)
)

metric_names <- c(
  "response_prediction_mse", "signal_prediction_mse", "ise_total",
  "ise_null", "ise_nonnull", "true_zero_rate", "false_zero_rate",
  "false_active_proportion", "null_iou", "null_symmetric_difference",
  "estimated_zero_length", "boundary_error"
)
summary_differences <- numeric()
for (index in seq_len(nrow(analysis$method_case_summary))) {
  row <- analysis$method_case_summary[index, ]
  values <- analysis$combined_metrics[
    analysis$combined_metrics$scenario == row$scenario &
      analysis$combined_metrics$method == row$method &
      analysis$combined_metrics$applicable %in% TRUE,
    row$metric
  ]
  values <- values[is.finite(values)]
  expected_mean <- if (length(values)) mean(values) else NA_real_
  if (is.finite(expected_mean) && is.finite(row$mean)) {
    summary_differences <- c(summary_differences, abs(expected_mean - row$mean))
  } else if (!identical(is.na(expected_mean), is.na(row$mean))) {
    summary_differences <- c(summary_differences, Inf)
  }
}
max_summary_difference <- if (length(summary_differences)) {
  max(summary_differences)
} else 0
add_check(
  "case_summary_independent_reconstruction",
  max_summary_difference <= 1e-14,
  paste("max absolute difference", format(max_summary_difference,
                                          scientific = TRUE))
)

curve_summary_differences <- numeric()
case_ii_curves <- analysis$combined_curves[
  analysis$combined_curves$scenario == "case_ii_sparse", , drop = FALSE
]
analysis_curve_columns <- grep("^beta_[0-9]+$",
                               names(case_ii_curves), value = TRUE)
for (method in expected_methods) {
  matrix_values <- as.matrix(case_ii_curves[
    case_ii_curves$method == method, analysis_curve_columns, drop = FALSE
  ])
  finite_rows <- apply(matrix_values, 1L, function(values) all(is.finite(values)))
  matrix_values <- matrix_values[finite_rows, , drop = FALSE]
  reported <- analysis$curve_summary[
    analysis$curve_summary$method == method, , drop = FALSE
  ]
  if (!nrow(matrix_values)) next
  expected_mean <- colMeans(matrix_values)
  curve_summary_differences <- c(
    curve_summary_differences, abs(expected_mean - reported$mean)
  )
}
max_curve_summary_difference <- if (length(curve_summary_differences)) {
  max(curve_summary_differences)
} else Inf
add_check(
  "case_ii_curve_summary_reconstruction",
  max_curve_summary_difference <= 1e-14,
  paste("max absolute difference", format(max_curve_summary_difference,
                                          scientific = TRUE))
)
add_check(
  "analysis_internal_gates",
  all(analysis$validation$passed),
  paste(sum(analysis$validation$passed), "/",
        nrow(analysis$validation), "PASS")
)
add_check(
  "smoke_noninferential_contract",
  if (smoke) !isTRUE(analysis$status$inferential) else {
    isTRUE(analysis$status$inferential)
  },
  analysis$status$interpretation
)

if (!core_only) {
  report_manifest_check <- verify_manifest(reporting_dir, report_manifest_file)
  pdf_manifest_check <- verify_manifest(pdf_dir, pdf_manifest_file)
  add_check("reporting_manifest_integrity", report_manifest_check$passed,
            report_manifest_check$details)
  add_check("pdf_mirror_manifest_integrity", pdf_manifest_check$passed,
            pdf_manifest_check$details)
  add_check(
    "reporting_signature_link",
    report_manifest_check$passed &&
      all(report_manifest_check$manifest$scientific_signature ==
            run$scientific_signature) &&
      all(pdf_manifest_check$manifest$scientific_signature ==
            run$scientific_signature),
    run$scientific_signature
  )
  report_files <- report_manifest_check$manifest$file
  pdf_report <- report_files[grepl("[.]pdf$", report_files)]
  png_report <- report_files[grepl("[.]png$", report_files)]
  tex_report <- report_files[grepl("[.]tex$", report_files)]
  add_check(
    "publication_payload_counts",
    length(pdf_report) == 2L && length(png_report) == 2L &&
      length(tex_report) == 6L &&
      any(grepl("sessionInfo", report_files)),
    paste(length(pdf_report), "PDF;", length(png_report), "PNG;",
          length(tex_report), "TeX")
  )
  report_pdf_paths <- file.path(reporting_dir, pdf_report)
  mirror_pdf_paths <- file.path(pdf_dir, basename(pdf_report))
  add_check(
    "pdf_mirror_identity",
    all(file.exists(mirror_pdf_paths)) &&
      identical(
        unname(as.character(tools::md5sum(report_pdf_paths))),
        unname(as.character(tools::md5sum(mirror_pdf_paths)))
      ),
    paste(basename(report_pdf_paths), collapse = ",")
  )
  pdfinfo <- Sys.which("pdfinfo")
  pdftotext <- Sys.which("pdftotext")
  pdf_readable <- nzchar(pdfinfo) && nzchar(pdftotext) && all(vapply(
    report_pdf_paths, function(path) {
      info <- suppressWarnings(system2(pdfinfo, path, stdout = TRUE,
                                       stderr = TRUE))
      text <- suppressWarnings(system2(pdftotext, c(path, "-"),
                                       stdout = TRUE, stderr = TRUE))
      any(grepl("^Pages:[[:space:]]+1$", info)) &&
        nchar(paste(text, collapse = "")) > 50L
    }, logical(1L)
  ))
  add_check(
    "pdf_one_page_and_text_readable", pdf_readable,
    if (nzchar(pdfinfo) && nzchar(pdftotext)) "Poppler checks" else {
      "Poppler unavailable"
    }
  )
  png_available <- requireNamespace("png", quietly = TRUE)
  png_paths <- file.path(reporting_dir, png_report)
  png_readable <- png_available && all(vapply(png_paths, function(path) {
    image <- png::readPNG(path)
    length(dim(image)) >= 2L && all(dim(image)[1:2] > 500L) &&
      stats::sd(as.numeric(image)) > 0
  }, logical(1L)))
  add_check("png_readability", png_readable,
            paste(basename(png_paths), collapse = ","))
  table_paths <- file.path(
    reporting_dir,
    report_files[grepl("table_", report_files) & grepl("[.]tex$", report_files)]
  )
  table_text <- paste(unlist(lapply(table_paths, readLines, warn = FALSE)),
                      collapse = "\n")
  labels_present <- all(vapply(
    unname(config$method_labels[expected_methods]),
    function(label) grepl(label, table_text, fixed = TRUE), logical(1L)
  ))
  add_check(
    "table_method_labels_and_no_slos_row",
    labels_present && !grepl("SLoS", table_text, fixed = TRUE),
    "12 frozen labels; no SLoS table row"
  )
  add_check(
    "table_latex_labels_valid",
    !grepl("tab:tab:", table_text, fixed = TRUE),
    "knitr table labels contain exactly one tab: prefix"
  )
  add_check(
    "mode_specific_publication_labels",
    if (smoke) {
      grepl("Technical smoke", table_text, fixed = TRUE) &&
        all(grepl("_smoke", c(pdf_report, png_report)))
    } else {
      !grepl("Technical smoke", table_text, fixed = TRUE) &&
        !any(grepl("smoke", c(pdf_report, png_report), ignore.case = TRUE))
    },
    mode
  )
} else {
  for (name in c(
    "reporting_manifest_integrity", "pdf_mirror_manifest_integrity",
    "reporting_signature_link", "publication_payload_counts",
    "pdf_mirror_identity", "pdf_one_page_and_text_readable",
    "png_readability", "table_method_labels_and_no_slos_row",
    "table_latex_labels_valid",
    "mode_specific_publication_labels"
  )) {
    add_check(name, FALSE, "not applicable under --core-only", applicable = FALSE)
  }
}

validation <- do.call(rbind, checks)
failed <- validation$applicable %in% TRUE & !(validation$passed %in% TRUE)
if (any(failed)) {
  stop("Validation failed: ",
       paste(validation$check[failed], collapse = ", "), call. = FALSE)
}

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
staging <- file.path(
  results_root, paste0(".", basename(validation_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) stop("Validation staging exists.", call. = FALSE)
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
mode_suffix <- if (smoke) "_smoke" else ""
core_suffix <- if (core_only) "_core" else ""
base <- paste0(prefix, "_validation_", report_version, mode_suffix,
               core_suffix)
validation_file <- file.path(staging, paste0(base, ".csv"))
utils::write.csv(validation, validation_file, row.names = FALSE)
validation_object <- list(
  study = prefix, mode = mode, core_only = core_only,
  scientific_signature = run$scientific_signature,
  checks = validation,
  reconstruction = list(
    max_prediction_difference = max_prediction_difference,
    max_structural_difference = max_structural_difference,
    max_summary_difference = max_summary_difference,
    max_curve_summary_difference = max_curve_summary_difference
  ),
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(validation_object, file.path(staging, paste0(base, ".rds")),
        compress = "xz", version = 3)
display <- validation
display$Status <- ifelse(
  !display$applicable, "N/A", ifelse(display$passed, "PASS", "FAIL")
)
display <- display[c("check", "Status", "details")]
names(display) <- c("Check", "Status", "Details")
if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("knitr is required for the validation table.", call. = FALSE)
}
table <- knitr::kable(
  display, format = "latex", booktabs = TRUE, escape = TRUE,
  caption = paste0(
    if (smoke) "Technical smoke: " else "Production: ",
    "ARFC classical-comparator validation gates."
  ),
  label = "tab:arfc-legacy-validation", longtable = TRUE, linesep = ""
)
writeLines(as.character(table), file.path(staging, paste0(base, ".tex")))
writeLines(capture.output(sessionInfo()),
           file.path(staging, paste0(base, "_sessionInfo.txt")))
payload_paths <- list.files(staging, full.names = TRUE)
manifest <- data.frame(
  scientific_signature = run$scientific_signature,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  manifest, file.path(staging, paste0(base, "_manifest.csv")),
  row.names = FALSE
)
if (!file.rename(staging, validation_dir)) {
  stop("Atomic validation bundle install failed.", call. = FALSE)
}
cat("ARFC classical-comparator validation complete.\n")
cat("Applicable checks:", sum(validation$applicable), "\n")
cat("PASS:", sum(validation$passed %in% TRUE, na.rm = TRUE), "\n")
cat("N/A:", sum(!validation$applicable), "\n")
cat("Bundle:", validation_dir, "\n")
