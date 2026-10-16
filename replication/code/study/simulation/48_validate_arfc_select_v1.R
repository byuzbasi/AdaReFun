#!/usr/bin/env Rscript

# Independent end-to-end validation for ARFC-Select version 1.
# Reconstructs BIC choices, validation selection, test oracle/regret, analysis
# summaries, coefficient curves, and publication artifacts.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/48_validate_arfc_select_v1.R", mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
version_argument <- grep("^--report-version=", arguments, value = TRUE)
if (length(version_argument) > 1L) {
  stop("Supply at most one --report-version argument.", call. = FALSE)
}
report_version <- if (length(version_argument)) {
  sub("^--report-version=", "", version_argument)
} else {
  "v1"
}
if (!report_version %in% c("v1", "v2")) {
  stop("`--report-version` must be v1 or v2.", call. = FALSE)
}
unknown <- setdiff(arguments, c("--smoke", version_argument))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

study <- "arfc_select_v1"
mode <- if (smoke) "smoke" else "production"
results_root <- file.path(simulation_dir, "results", study)
run_dir <- file.path(results_root, paste0(mode, "_run"))
analysis_dir <- file.path(results_root, paste0(mode, "_analysis_v1"))
report_dir <- file.path(
  project_dir, "output", "pdf", study,
  paste0(mode, "_reporting_", report_version)
)
validation_dir <- file.path(
  results_root, paste0(mode, "_validation_", report_version)
)
required_directories <- c(run_dir, analysis_dir, report_dir)
if (!all(dir.exists(required_directories))) {
  stop(
    "Run, analysis, and reporting bundles must exist before validation: ",
    paste(required_directories[!dir.exists(required_directories)],
          collapse = ", "),
    call. = FALSE
  )
}
if (dir.exists(validation_dir) || file.exists(validation_dir)) {
  stop("Refusing to overwrite ARFC validation bundle: ", validation_dir,
       call. = FALSE)
}

prefix <- study
run_file <- file.path(run_dir, paste0(prefix, "_run.rds"))
run_manifest_file <- file.path(run_dir, paste0(prefix, "_manifest.csv"))
analysis_file <- file.path(analysis_dir, paste0(prefix, "_analysis_v1.rds"))
analysis_manifest_file <- file.path(
  analysis_dir, paste0(prefix, "_analysis_manifest_v1.csv")
)
report_manifest_file <- file.path(
  report_dir, paste0(prefix, "_reporting_manifest_", report_version, ".csv")
)
required_files <- c(
  run_file, run_manifest_file, analysis_file, analysis_manifest_file,
  report_manifest_file
)
if (!all(file.exists(required_files))) {
  stop("One or more canonical input files are missing.", call. = FALSE)
}

verify_manifest <- function(directory, manifest_file, require_role = FALSE) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required <- c("scientific_signature", "file", "bytes", "md5")
  if (require_role) required <- c(required, "role")
  if (!all(required %in% names(manifest)) || !nrow(manifest) ||
      anyDuplicated(manifest$file)) {
    stop("Invalid manifest schema: ", manifest_file, call. = FALSE)
  }
  paths <- file.path(directory, manifest$file)
  observed_md5 <- if (all(file.exists(paths))) {
    unname(as.character(tools::md5sum(paths)))
  } else rep(NA_character_, length(paths))
  observed_bytes <- if (all(file.exists(paths))) {
    as.numeric(file.info(paths)$size)
  } else rep(NA_real_, length(paths))
  valid <- all(file.exists(paths)) &&
    identical(observed_md5, as.character(manifest$md5)) &&
    all(observed_bytes == as.numeric(manifest$bytes))
  list(manifest = manifest, paths = paths, valid = valid)
}

run_manifest <- verify_manifest(run_dir, run_manifest_file)
analysis_manifest <- verify_manifest(analysis_dir, analysis_manifest_file)
report_manifest <- verify_manifest(
  report_dir, report_manifest_file, require_role = TRUE
)
if (!run_manifest$valid || !analysis_manifest$valid ||
    !report_manifest$valid) {
  stop("At least one input manifest failed byte/MD5 verification.",
       call. = FALSE)
}

run <- readRDS(run_file)
analysis <- readRDS(analysis_file)
metadata_path <- report_manifest$paths[
  report_manifest$manifest$role == "reporting_metadata_rds"
]
if (length(metadata_path) != 1L) {
  stop("Reporting metadata role is missing or duplicated.", call. = FALSE)
}
reporting <- readRDS(metadata_path)

checks <- list()
add_check <- function(name, passed, detail, applicable = TRUE) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    applicable = as.logical(applicable),
    passed = if (applicable) isTRUE(passed) else NA,
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

same_numeric <- function(left, right, tolerance = 1e-10) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (length(left) != length(right)) return(FALSE)
  both_na <- is.na(left) & is.na(right)
  if (any(xor(is.na(left), is.na(right)))) return(FALSE)
  finite <- is.finite(left) & is.finite(right)
  if (any(!(both_na | finite))) return(FALSE)
  if (!any(finite)) return(TRUE)
  max(abs(left[finite] - right[finite])) <= tolerance
}

signatures <- c(
  run$scientific_signature,
  analysis$scientific_signature,
  reporting$scientific_signature,
  unique(as.character(run_manifest$manifest$scientific_signature)),
  unique(as.character(analysis_manifest$manifest$scientific_signature)),
  unique(as.character(report_manifest$manifest$scientific_signature))
)
add_check(
  "identity_and_signature_chain",
  identical(run$study, study) && identical(run$mode, mode) &&
    identical(analysis$study, study) && identical(analysis$mode, mode) &&
    identical(reporting$study, study) && identical(reporting$mode, mode) &&
    (identical(report_version, "v1") ||
       identical(reporting$report_version, report_version)) &&
    length(unique(signatures)) == 1L,
  paste("signature", run$scientific_signature)
)
add_check(
  "all_manifests_md5_and_bytes",
  run_manifest$valid && analysis_manifest$valid && report_manifest$valid,
  sprintf(
    "run=%d analysis=%d report=%d payloads",
    nrow(run_manifest$manifest), nrow(analysis_manifest$manifest),
    nrow(report_manifest$manifest)
  )
)

source_paths <- c(
  file.path(project_dir, "config", "protocol_arfc_select_v1.md"),
  file.path(
    project_dir, "config",
    "protocol_adaptive_regional_folded_concave_v1.md"
  ),
  file.path(simulation_dir, "src", "adaptive_regional_core.cpp"),
  file.path(simulation_dir, "45_run_arfc_select_v1.R")
)
source_names <- basename(source_paths)
source_md5_current <- setNames(
  unname(as.character(tools::md5sum(source_paths))), source_names
)
add_check(
  "frozen_scientific_source_hashes",
  all(file.exists(source_paths)) &&
    identical(source_md5_current[names(run$source_md5)], run$source_md5),
  paste(length(source_paths), "protocol/core/runner hashes")
)
package_paths <- c(
  file.path(project_dir, "FMCP-package", "src", "fmcp_core.cpp"),
  file.path(project_dir, "FMCP-package", "src", "fscad_core.cpp"),
  file.path(project_dir, "FMCP-package", "R", "basis.R")
)
package_md5_current <- setNames(
  unname(as.character(tools::md5sum(package_paths))), basename(package_paths)
)
add_check(
  "frozen_package_source_hashes",
  all(file.exists(package_paths)) &&
    identical(
      package_md5_current[names(run$package_source_md5)],
      run$package_source_md5
    ),
  paste(length(package_paths), "package source hashes")
)
add_check(
  "analysis_and_writer_source_hashes",
  identical(
    analysis$analysis_script_md5,
    unname(as.character(tools::md5sum(file.path(
      simulation_dir, "46_analyze_arfc_select_v1.R"
    ))))
  ) && identical(
    reporting$writer_script_md5,
    unname(as.character(tools::md5sum(file.path(
      simulation_dir, "47_write_arfc_select_v1.R"
    ))))
  ),
  "analysis and reporting scripts match their saved provenance"
)
add_check(
  "upstream_manifest_provenance",
  identical(
    analysis$run_manifest_md5,
    unname(as.character(tools::md5sum(run_manifest_file)))
  ) && identical(
    reporting$analysis_manifest_md5,
    unname(as.character(tools::md5sum(analysis_manifest_file)))
  ),
  "run -> analysis -> reporting manifest chain"
)

config <- run$config
tasks <- run$tasks
methods <- as.character(config$method_order)
path <- run$path_tuning
candidates <- run$candidate_metrics
selected <- run$selected_arfc_metrics
selections <- run$selection_outcomes
coefficients <- run$selected_coefficients
expected_tasks <- 4L * config$replications
expected_candidates <- expected_tasks * length(methods)
expected_path <- expected_candidates * length(config$lambda)

add_check(
  "frozen_design_dimensions",
  nrow(tasks) == expected_tasks && nrow(candidates) == expected_candidates &&
    nrow(selected) == expected_tasks && nrow(selections) == expected_tasks &&
    nrow(path) == expected_path && nrow(coefficients) == expected_candidates,
  sprintf(
    "%d tasks, %d candidate fits, %d path rows",
    nrow(tasks), nrow(candidates), nrow(path)
  )
)
add_check(
  "sample_sizes_and_replications",
  all(tasks$n_train == config$n_train) &&
    all(tasks$n_validation == config$n_validation) &&
    all(tasks$n_test == config$n_test) &&
    all(table(tasks$scenario) == config$replications),
  sprintf(
    "train=%d validation=%d test=%d; %d reps/case",
    config$n_train, config$n_validation, config$n_test,
    config$replications
  )
)
all_seeds <- unlist(tasks[c("seed_train", "seed_validation", "seed_test")],
                    use.names = FALSE)
add_check(
  "independent_deterministic_seed_ids",
  !anyDuplicated(all_seeds) && all(is.finite(all_seeds)),
  paste(length(all_seeds), "unique train/validation/test seed IDs")
)
add_check(
  "candidate_library_and_tie_order",
  identical(
    methods, c("asr_mcp", "asr_scad", "amr_mcp", "amr_scad")
  ) && identical(
    run$method_metadata$method[
      order(run$method_metadata$validation_tie_priority)
    ], methods
  ),
  paste(methods, collapse = ", ")
)
add_check(
  "equal_component_maximum_concavity",
  abs(1 / config$gamma - 1 / (config$scad_a - 1)) <= 1e-14,
  sprintf("MCP gamma=%.1f; SCAD a=%.1f", config$gamma, config$scad_a)
)
add_check(
  "canonical_mass_and_envelope_checks",
  nrow(run$mass_checks) == 4L && all(run$mass_checks$passed) &&
    all(run$mass_checks$observed <= run$mass_checks$tolerance),
  sprintf("%d/%d runner checks passed",
          sum(run$mass_checks$passed), nrow(run$mass_checks))
)

path_groups <- split(
  path, interaction(path$task_id, path$method, drop = TRUE)
)
path_grid_complete <- length(path_groups) == expected_candidates &&
  all(vapply(path_groups, function(block) {
    nrow(block) == length(config$lambda) &&
      identical(as.numeric(block$lambda), as.numeric(config$lambda)) &&
      length(unique(block$kappa)) == 1L
  }, logical(1L)))
add_check(
  "common_complete_lambda_paths",
  path_grid_complete,
  sprintf("%d task-method paths x %d lambda values",
          length(path_groups), length(config$lambda))
)

reselected_rows <- vector("list", length(path_groups))
group_names <- names(path_groups)
for (index in seq_along(path_groups)) {
  block <- path_groups[[index]]
  eligible <- which(
    block$fit_usable %in% TRUE & block$success %in% TRUE &
      block$converged %in% TRUE & block$refit_success %in% TRUE &
      is.finite(block$bic)
  )
  chosen <- if (length(eligible)) {
    eligible[which.min(block$bic[eligible])]
  } else NA_integer_
  if (is.finite(chosen)) {
    reselected_rows[[index]] <- block[chosen, c(
      "task_id", "method", "candidate_index", "lambda", "kappa", "bic"
    )]
  } else {
    reselected_rows[[index]] <- data.frame(
      task_id = block$task_id[1L], method = block$method[1L],
      candidate_index = NA, lambda = NA, kappa = NA, bic = NA
    )
  }
}
reselected <- do.call(rbind, reselected_rows)
rownames(reselected) <- NULL
candidate_order <- order(candidates$task_id, match(candidates$method, methods))
reselected_order <- order(reselected$task_id, match(reselected$method, methods))
candidate_sorted <- candidates[candidate_order, , drop = FALSE]
reselected_sorted <- reselected[reselected_order, , drop = FALSE]
strict_reselection_pass <- identical(
  paste(candidate_sorted$task_id, candidate_sorted$method),
  paste(reselected_sorted$task_id, reselected_sorted$method)
) && identical(
  as.integer(candidate_sorted$candidate_index),
  as.integer(reselected_sorted$candidate_index)
) && same_numeric(candidate_sorted$lambda, reselected_sorted$lambda) &&
  same_numeric(candidate_sorted$kappa, reselected_sorted$kappa) &&
  same_numeric(candidate_sorted$bic, reselected_sorted$bic)
add_check(
  "strict_bic_reselection_from_full_paths",
  strict_reselection_pass,
  paste(expected_candidates, "candidate choices reconstructed")
)
add_check(
  "all_selected_candidate_fits_usable",
  all(with(
    candidates,
    fit_usable & bic_selected & success & converged & refit_success &
      is.finite(bic)
  )),
  sprintf("%d/%d strict fits", sum(candidates$fit_usable),
          nrow(candidates))
)

reselected_methods <- character(nrow(tasks))
reoracle_methods <- character(nrow(tasks))
recomputed_regret <- numeric(nrow(tasks))
recomputed_ties <- integer(nrow(tasks))
for (index in seq_len(nrow(tasks))) {
  block <- candidates[candidates$task_id == tasks$task_id[index],
                      , drop = FALSE]
  block <- block[match(methods, block$method), , drop = FALSE]
  validation_min <- min(block$validation_response_mse)
  validation_tolerance <- config$validation_tie_relative_tolerance *
    (1 + abs(validation_min))
  validation_tied <- block$method[
    block$validation_response_mse <= validation_min + validation_tolerance
  ]
  reselected_methods[index] <- methods[methods %in% validation_tied][1L]
  recomputed_ties[index] <- length(validation_tied)

  test_min <- min(block$response_prediction_mse)
  test_tolerance <- config$validation_tie_relative_tolerance *
    (1 + abs(test_min))
  test_tied <- block$method[
    block$response_prediction_mse <= test_min + test_tolerance
  ]
  reoracle_methods[index] <- methods[methods %in% test_tied][1L]
  chosen_row <- block[block$method == reselected_methods[index],
                      , drop = FALSE]
  oracle_row <- block[block$method == reoracle_methods[index], , drop = FALSE]
  recomputed_regret[index] <- chosen_row$response_prediction_mse -
    oracle_row$response_prediction_mse
}
selection_order <- match(tasks$task_id, selections$task_id)
selection_sorted <- selections[selection_order, , drop = FALSE]
add_check(
  "independent_validation_selection_reconstruction",
  identical(as.character(selection_sorted$selected_method),
            reselected_methods) &&
    identical(as.integer(selection_sorted$validation_tie_count),
              recomputed_ties) &&
    all(selection_sorted$selection_used_test == FALSE),
  paste(nrow(tasks), "task selections; test_used=FALSE")
)
add_check(
  "test_oracle_and_regret_reconstruction",
  identical(as.character(selection_sorted$test_oracle_method),
            reoracle_methods) &&
    same_numeric(selection_sorted$prediction_regret, recomputed_regret) &&
    identical(
      as.logical(selection_sorted$test_oracle_agreement),
      reselected_methods == reoracle_methods
    ),
  paste(nrow(tasks), "test-only oracle diagnostics")
)

candidate_selected_flags <- vapply(seq_len(nrow(candidates)), function(index) {
  selection_sorted$selected_method[
    match(candidates$task_id[index], selection_sorted$task_id)
  ] == candidates$method[index]
}, logical(1L))
candidate_oracle_flags <- vapply(seq_len(nrow(candidates)), function(index) {
  selection_sorted$test_oracle_method[
    match(candidates$task_id[index], selection_sorted$task_id)
  ] == candidates$method[index]
}, logical(1L))
add_check(
  "candidate_selection_and_oracle_flags",
  identical(as.logical(candidates$selected_by_validation),
            candidate_selected_flags) &&
    identical(as.logical(candidates$test_oracle), candidate_oracle_flags) &&
    all(table(candidates$task_id)[
      as.character(unique(candidates$task_id))
    ] == length(methods)),
  "one validation selection and one deterministic oracle per task"
)

case_i <- candidates$scenario == "case_i_null"
case_ii <- candidates$scenario == "case_ii_sparse"
no_null <- candidates$scenario %in% c("case_iii_crossing", "case_iv_away")
case_semantics <-
  all(is.finite(candidates$ise_null[case_i])) &&
  all(is.na(candidates$ise_nonnull[case_i])) &&
  all(is.finite(candidates$true_zero_rate[case_i])) &&
  all(is.na(candidates$false_zero_rate[case_i])) &&
  all(is.finite(candidates$false_active_proportion[case_i])) &&
  all(is.finite(candidates$ise_null[case_ii])) &&
  all(is.finite(candidates$ise_nonnull[case_ii])) &&
  all(is.finite(candidates$true_zero_rate[case_ii])) &&
  all(is.finite(candidates$false_zero_rate[case_ii])) &&
  all(is.finite(candidates$null_iou[case_ii])) &&
  all(is.finite(candidates$boundary_error[case_ii])) &&
  all(is.na(candidates$ise_null[no_null])) &&
  all(is.finite(candidates$ise_nonnull[no_null])) &&
  all(is.na(candidates$true_zero_rate[no_null])) &&
  all(is.finite(candidates$false_zero_rate[no_null])) &&
  all(is.na(candidates$boundary_error[no_null]))
add_check(
  "case_specific_metric_semantics",
  case_semantics,
  "all-null, local-null, and no-null NA/finite contracts"
)
add_check(
  "finite_prediction_metrics",
  all(is.finite(candidates$validation_response_mse)) &&
    all(is.finite(candidates$validation_signal_mse)) &&
    all(is.finite(candidates$response_prediction_mse)) &&
    all(is.finite(candidates$signal_prediction_mse)),
  paste(nrow(candidates), "candidate train/validation/test records")
)

# Analysis reconstruction helpers.
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

selection_expected <- analysis$selection_frequency
selection_frequency_pass <- TRUE
for (index in seq_len(nrow(selection_expected))) {
  row <- selection_expected[index, , drop = FALSE]
  block <- selections[selections$scenario == row$scenario, , drop = FALSE]
  count <- sum(block$selected_method == row$method)
  selection_frequency_pass <- selection_frequency_pass &&
    count == row$selected_count && nrow(block) == row$n_tasks &&
    same_numeric(count / nrow(block), row$selection_frequency)
}
selection_mass <- tapply(
  selection_expected$selection_frequency,
  selection_expected$scenario, sum
)
add_check(
  "selection_frequency_reconstruction",
  selection_frequency_pass && all(abs(selection_mass - 1) <= 1e-12),
  paste(nrow(selection_expected), "per-case candidate-frequency rows")
)

margin_pass <- TRUE
for (index in seq_len(nrow(analysis$selection_margins))) {
  row <- analysis$selection_margins[index, , drop = FALSE]
  block <- selections[selections$scenario == row$scenario, , drop = FALSE]
  column <- if (row$dimension == "penalty_family") {
    "selected_penalty_family"
  } else "selected_stage"
  count <- sum(block[[column]] == row$level)
  margin_pass <- margin_pass && count == row$selected_count &&
    nrow(block) == row$n_tasks &&
    same_numeric(count / nrow(block), row$selection_frequency)
}
add_check(
  "marginal_family_and_resolution_frequencies",
  margin_pass,
  paste(nrow(analysis$selection_margins), "marginal rows")
)

selected_performance_pass <- TRUE
for (index in seq_len(nrow(analysis$selected_performance))) {
  row <- analysis$selected_performance[index, , drop = FALSE]
  block <- selected[selected$scenario == row$scenario, , drop = FALSE]
  summary <- mean_ci(block[[row$metric]])
  selected_performance_pass <- selected_performance_pass &&
    same_numeric(
      unlist(row[c("n", "mean", "sd", "se", "lower_95", "upper_95")]),
      unname(summary[c("n", "mean", "sd", "se", "lower", "upper")])
    )
}
add_check(
  "selected_performance_reconstruction",
  selected_performance_pass,
  paste(nrow(analysis$selected_performance), "per-case metric summaries")
)

candidate_performance_pass <- TRUE
for (index in seq_len(nrow(analysis$candidate_performance))) {
  row <- analysis$candidate_performance[index, , drop = FALSE]
  block <- candidates[
    candidates$scenario == row$scenario & candidates$method == row$method,
    , drop = FALSE
  ]
  summary <- mean_ci(block[[row$metric]])
  candidate_performance_pass <- candidate_performance_pass &&
    same_numeric(
      unlist(row[c("n", "mean", "sd", "se", "lower_95", "upper_95")]),
      unname(summary[c("n", "mean", "sd", "se", "lower", "upper")])
    )
}
add_check(
  "candidate_performance_reconstruction",
  candidate_performance_pass,
  paste(nrow(analysis$candidate_performance), "candidate metric summaries")
)

oracle_pass <- TRUE
for (index in seq_len(nrow(analysis$oracle_summary))) {
  row <- analysis$oracle_summary[index, , drop = FALSE]
  block <- selections[selections$scenario == row$scenario, , drop = FALSE]
  regret <- mean_ci(block$prediction_regret)
  oracle_pass <- oracle_pass && row$n_tasks == nrow(block) &&
    row$oracle_agreement_count == sum(block$test_oracle_agreement) &&
    same_numeric(
      c(row$oracle_agreement_rate, row$mean_prediction_regret,
        row$regret_se, row$regret_lower_95, row$regret_upper_95,
        row$median_prediction_regret, row$validation_tie_rate),
      c(mean(block$test_oracle_agreement), regret["mean"], regret["se"],
        regret["lower"], regret["upper"],
        stats::median(block$prediction_regret),
        mean(block$validation_tie_count > 1L))
    )
}
add_check(
  "oracle_summary_reconstruction",
  oracle_pass,
  paste(nrow(analysis$oracle_summary), "case-wise oracle rows")
)

contrast <- analysis$selected_candidate_contrasts
selected_map <- match(contrast$task_id, selected$task_id)
candidate_keys <- paste(candidates$task_id, candidates$method, sep = "::")
contrast_candidate_map <- match(
  paste(contrast$task_id, contrast$candidate_method, sep = "::"),
  candidate_keys
)
contrast_values_pass <- !anyNA(selected_map) && !anyNA(contrast_candidate_map)
if (contrast_values_pass) {
  expected_difference <- vapply(seq_len(nrow(contrast)), function(index) {
    metric <- contrast$metric[index]
    selected[[metric]][selected_map[index]] -
      candidates[[metric]][contrast_candidate_map[index]]
  }, numeric(1L))
  contrast_values_pass <- same_numeric(
    contrast$difference_selected_minus_candidate, expected_difference
  )
}
add_check(
  "selected_minus_candidate_contrasts",
  contrast_values_pass &&
    nrow(contrast) == expected_tasks * length(methods) * 2L,
  paste(nrow(contrast), "task-level paired contrast rows")
)

contrast_summary_pass <- TRUE
for (index in seq_len(nrow(analysis$selected_candidate_case_summaries))) {
  row <- analysis$selected_candidate_case_summaries[index, , drop = FALSE]
  block <- contrast[
    contrast$scenario == row$scenario &
      contrast$candidate_method == row$candidate_method &
      contrast$metric == row$metric, , drop = FALSE
  ]
  summary <- mean_ci(block$difference_selected_minus_candidate)
  contrast_summary_pass <- contrast_summary_pass &&
    same_numeric(
      c(row$n_pairs, row$selected_mean, row$candidate_mean,
        row$mean_difference, row$se, row$lower_95, row$upper_95),
      c(summary["n"], mean(block$selected_value),
        mean(block$candidate_value), summary["mean"], summary["se"],
        summary["lower"], summary["upper"])
    )
}
add_check(
  "paired_case_summary_reconstruction",
  contrast_summary_pass,
  paste(nrow(analysis$selected_candidate_case_summaries),
        "case-candidate-metric summaries")
)

# Independent coefficient-curve reconstruction.
library_override <- Sys.getenv(
  "ARFC_FMCP_LIBRARY", unset = ".replication-library"
)
if (!dir.exists(library_override)) {
  stop("Frozen FMCP library unavailable for curve validation: ",
       library_override, call. = FALSE)
}
.libPaths(c(normalizePath(library_override), .libPaths()))
if (!requireNamespace("FMCP", quietly = TRUE) ||
    utils::packageVersion("FMCP") != numeric_version("0.2.0")) {
  stop("Curve validation requires frozen FMCP 0.2.0.", call. = FALSE)
}
argvals <- seq(0, 1, length.out = config$observed_grid_points)
support_grid <- seq(0, 1, length.out = config$support_grid_points)
basis <- FMCP:::.fmcp_basis_setup(
  argvals, config$finest_M, config$degree, config$penalty_order
)
dense_basis <- unclass(splines::splineDesign(
  knots = basis$full_knots, x = support_grid, ord = basis$order,
  derivs = rep.int(0L, length(support_grid)), outer.ok = TRUE
))
coefficient_names <- grep("^b_[0-9]+$", names(coefficients), value = TRUE)
coefficient_matrix <- as.matrix(coefficients[coefficient_names])
storage.mode(coefficient_matrix) <- "double"
all_curves <- coefficient_matrix %*% t(dense_basis)

candidate_curve_pass <- TRUE
for (scenario in as.character(run$scenarios$scenario)) {
  for (method in methods) {
    rows <- which(
      coefficients$scenario == scenario & coefficients$method == method
    )
    observed <- analysis$candidate_curve_summary[
      analysis$candidate_curve_summary$scenario == scenario &
        analysis$candidate_curve_summary$method == method,
      , drop = FALSE
    ]
    observed <- observed[order(observed$t), , drop = FALSE]
    matrix_block <- all_curves[rows, , drop = FALSE]
    candidate_curve_pass <- candidate_curve_pass &&
      length(rows) == config$replications &&
      nrow(observed) == length(support_grid) &&
      same_numeric(observed$t, support_grid) &&
      same_numeric(observed$mean_estimate, colMeans(matrix_block)) &&
      same_numeric(
        observed$zero_selection_frequency,
        colMeans(abs(matrix_block) <= config$zero_tolerance)
      )
  }
}
add_check(
  "candidate_curve_reconstruction",
  candidate_curve_pass,
  paste(length(methods) * 4L, "case-method curve summaries")
)

coefficient_keys <- paste(coefficients$task_id, coefficients$method,
                          sep = "::")
selected_coefficient_rows <- match(
  paste(selection_sorted$task_id, selection_sorted$selected_method,
        sep = "::"),
  coefficient_keys
)
selected_curve_pass <- !anyNA(selected_coefficient_rows)
if (selected_curve_pass) {
  for (scenario in as.character(run$scenarios$scenario)) {
    task_rows <- which(selection_sorted$scenario == scenario)
    matrix_rows <- selected_coefficient_rows[task_rows]
    observed <- analysis$selected_curve_summary[
      analysis$selected_curve_summary$scenario == scenario,
      , drop = FALSE
    ]
    observed <- observed[order(observed$t), , drop = FALSE]
    matrix_block <- all_curves[matrix_rows, , drop = FALSE]
    selected_curve_pass <- selected_curve_pass &&
      nrow(observed) == length(support_grid) &&
      same_numeric(observed$mean_estimate, colMeans(matrix_block)) &&
      same_numeric(
        observed$lower_10,
        apply(matrix_block, 2L, stats::quantile, probs = 0.10,
              names = FALSE)
      ) && same_numeric(
        observed$upper_90,
        apply(matrix_block, 2L, stats::quantile, probs = 0.90,
              names = FALSE)
      ) && same_numeric(
        observed$zero_selection_frequency,
        colMeans(abs(matrix_block) <= config$zero_tolerance)
      )
  }
}
add_check(
  "selected_curve_and_envelope_reconstruction",
  selected_curve_pass,
  "four case-wise selected curves; Case-II 10%-90% envelope"
)

analysis_manifest_rds_row <- analysis_manifest$manifest$file ==
  basename(analysis_file)
analysis_csv_match <- sum(analysis_manifest_rds_row) == 1L &&
  identical(
    analysis$analysis_status$cross_case_overall_computed,
    FALSE
  ) && !any(grepl("overall", names(analysis), ignore.case = TRUE))
add_check(
  "no_cross_case_overall",
  analysis_csv_match &&
    identical(reporting$cross_case_overall_computed, FALSE),
  "all scientific summaries remain case-specific"
)
add_check(
  "smoke_inference_gate",
  identical(analysis$analysis_status$scientific_decision_allowed, !smoke) &&
    identical(reporting$scientific_interpretation_allowed, !smoke),
  if (smoke) "scientific interpretation disabled" else
    "production case-wise interpretation enabled"
)

expected_roles <- c(
  "selection_figure_pdf", "selection_figure_png",
  "case_ii_recovery_pdf", "case_ii_recovery_png",
  "selection_table_tex", "performance_table_tex", "contrast_table_tex",
  "results_text_tex", "figure_includes_tex", "reporting_metadata_rds",
  "session_info_txt"
)
add_check(
  "reporting_role_completeness",
  identical(sort(as.character(report_manifest$manifest$role)),
            sort(expected_roles)),
  paste(length(expected_roles), "versioned reporting roles")
)

performance_tex_path <- report_manifest$paths[
  report_manifest$manifest$role == "performance_table_tex"
]
contrast_tex_path <- report_manifest$paths[
  report_manifest$manifest$role == "contrast_table_tex"
]
v2_table_semantics <- if (identical(report_version, "v2")) {
  length(performance_tex_path) == 1L && length(contrast_tex_path) == 1L &&
    !grepl(
      "Tasks", paste(readLines(performance_tex_path, warn = FALSE),
                     collapse = "\n"), fixed = TRUE
    ) && !grepl(
      "ARFC-Select MSE",
      paste(readLines(contrast_tex_path, warn = FALSE), collapse = "\n"),
      fixed = TRUE
    ) && grepl(
      "common to all four rows",
      paste(readLines(contrast_tex_path, warn = FALSE), collapse = "\n"),
      fixed = TRUE
    )
} else {
  TRUE
}
add_check(
  "report_version_table_semantics",
  v2_table_semantics,
  if (identical(report_version, "v2")) {
    "Tasks removed; repeated ARFC-Select MSE collapsed into caption"
  } else {
    "legacy v1 reporting schema"
  }
)

pdf_paths <- report_manifest$paths[grepl(
  "_pdf$", report_manifest$manifest$role
)]
pdfinfo <- Sys.which("pdfinfo")
pdf_pass <- length(pdf_paths) == 2L && nzchar(pdfinfo) &&
  all(vapply(pdf_paths, function(path_value) {
    output <- suppressWarnings(system2(
      pdfinfo, path_value, stdout = TRUE, stderr = TRUE
    ))
    pages <- grep("^Pages:", output, value = TRUE)
    length(pages) == 1L && grepl("Pages:[[:space:]]+1$", pages) &&
      file.info(path_value)$size > 10000
  }, logical(1L)))
add_check(
  "pdf_readability_and_single_page_layout",
  pdf_pass,
  paste(length(pdf_paths), "single-page readable PDFs")
)

png_paths <- report_manifest$paths[grepl(
  "_png$", report_manifest$manifest$role
)]
png_pass <- requireNamespace("png", quietly = TRUE) &&
  length(png_paths) == 2L &&
  all(vapply(png_paths, function(path_value) {
    image <- png::readPNG(path_value)
    dimensions <- dim(image)
    length(dimensions) >= 2L && all(dimensions[1:2] >= 1000L) &&
      stats::sd(as.numeric(image)) > 0.01
  }, logical(1L)))
add_check(
  "png_dimensions_and_nonblank_content",
  png_pass,
  paste(length(png_paths), "readable nonblank PNGs")
)

tex_paths <- report_manifest$paths[grepl(
  "_tex$", report_manifest$manifest$role
)]
tex_text <- paste(vapply(tex_paths, function(path_value) {
  paste(readLines(path_value, warn = FALSE), collapse = "\n")
}, character(1L)), collapse = "\n")
add_check(
  "latex_fragment_semantics",
  length(tex_paths) == 5L &&
    grepl("tab:arfc-select-frequency", tex_text, fixed = TRUE) &&
    grepl("tab:arfc-select-performance", tex_text, fixed = TRUE) &&
    grepl("tab:arfc-select-caseii-contrast", tex_text, fixed = TRUE) &&
    grepl("fig:arfc-select-frequency", tex_text, fixed = TRUE) &&
    grepl("fig:arfc-select-caseii", tex_text, fixed = TRUE),
  paste(length(tex_paths), "LaTeX fragments with unique labels")
)
smoke_language_pass <- if (smoke) {
  grepl("noninferential", tex_text, ignore.case = TRUE) &&
    grepl("technical smoke", tex_text, ignore.case = TRUE)
} else {
  !grepl("technical smoke", tex_text, ignore.case = TRUE) &&
    !any(grepl("smoke", basename(report_manifest$paths), ignore.case = TRUE))
}
add_check(
  "mode_appropriate_reporting_language",
  smoke_language_pass,
  if (smoke) "smoke is explicitly noninferential" else
    "no smoke language or file leakage"
)

validation <- do.call(rbind, checks)
rownames(validation) <- NULL
applicable <- validation$applicable
all_passed <- all(validation$passed[applicable])

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
target_dir <- if (all_passed) {
  validation_dir
} else {
  file.path(
    results_root,
    paste0(
      "failed_", mode, "_validation_", report_version, "_",
      format(Sys.time(), "%Y%m%d_%H%M%S")
    )
  )
}
staging <- file.path(
  results_root, paste0(".", basename(target_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging) || dir.exists(target_dir)) {
  stop("Validation staging or target path already exists.", call. = FALSE)
}
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed && dir.exists(staging)) {
  unlink(staging, recursive = TRUE, force = TRUE)
}, add = TRUE)

validation_file <- paste0(prefix, "_validation_", report_version, ".csv")
validation_rds_file <- paste0(prefix, "_validation_", report_version, ".rds")
validation_tex_file <- paste0(prefix, "_validation_", report_version, ".tex")
session_file <- paste0(
  prefix, "_validation_sessionInfo_", report_version, ".txt"
)
utils::write.csv(validation, file.path(staging, validation_file),
                 row.names = FALSE)
saveRDS(list(
  study = study,
  mode = mode,
  report_version = report_version,
  scientific_signature = run$scientific_signature,
  all_passed = all_passed,
  validation = validation,
  validator_script_md5 = unname(as.character(tools::md5sum(script_file))),
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(staging, validation_rds_file), compress = "xz", version = 3)

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("knitr is required for the validation LaTeX table.", call. = FALSE)
}
validation_table <- validation[c("check", "applicable", "passed", "detail")]
validation_table$applicable <- ifelse(
  validation_table$applicable, "Yes", "No"
)
validation_table$passed <- ifelse(
  is.na(validation_table$passed), "N/A",
  ifelse(validation_table$passed, "PASS", "FAIL")
)
names(validation_table) <- c("Check", "Applicable", "Status", "Detail")
latex <- knitr::kable(
  validation_table,
  format = "latex",
  booktabs = TRUE,
  longtable = TRUE,
  escape = TRUE,
  caption = paste(
    "ARFC-Select end-to-end validation.",
    if (smoke) "Technical smoke only." else "Production analysis."
  ),
  label = "arfc-select-validation",
  align = c("l", "c", "c", "p{7.2cm}")
)
writeLines(as.character(latex), file.path(staging, validation_tex_file),
           useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(staging, session_file),
           useBytes = TRUE)

payload_files <- c(
  validation_file, validation_rds_file, validation_tex_file, session_file
)
payload_paths <- file.path(staging, payload_files)
validation_manifest <- data.frame(
  scientific_signature = run$scientific_signature,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
utils::write.csv(
  validation_manifest,
  file.path(
    staging,
    paste0(prefix, "_validation_manifest_", report_version, ".csv")
  ),
  row.names = FALSE
)
if (!file.rename(staging, target_dir)) {
  stop("Could not atomically install ARFC validation bundle.",
       call. = FALSE)
}
committed <- TRUE

cat("ARFC selection validation complete.\n")
cat("Mode:", mode, "\n")
cat("Applicable checks:", sum(applicable), "\n")
cat("PASS:", sum(validation$passed[applicable]), "\n")
cat("FAIL:", sum(!validation$passed[applicable]), "\n")
cat("Bundle:", target_dir, "\n")
if (!all_passed) {
  failed <- validation$check[applicable & !validation$passed]
  stop("Validation failed: ", paste(failed, collapse = ", "),
       call. = FALSE)
}
