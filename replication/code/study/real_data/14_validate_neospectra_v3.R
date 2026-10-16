script_dir <- local({
  argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argument)) {
    dirname(normalizePath(sub("^--file=", "", argument[1L])))
  } else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
for (file in c("utils.R", "neospectra_v3.R")) {
  source(file.path(script_dir, "R", file))
}

arguments <- commandArgs(trailingOnly = TRUE)
mode <- real_mode(arguments)
preparation_only <- "--preparation-only" %in% arguments
config <- neospectra_v3_config(mode, project_root)
processed <- load_neospectra_v3_processed(config)
scope <- if (preparation_only) "preparation_only" else "complete_study"
stem <- if (preparation_only) "preparation_validation" else "validation"
paths <- list(
  report = neospectra_v3_result_path(config, stem, "rds"),
  checks = neospectra_v3_result_path(
    config, paste0(stem, "_checks"), "csv"
  ),
  inventory = neospectra_v3_result_path(
    config, paste0(stem, "_inventory"), "csv"
  ),
  session = neospectra_v3_result_path(
    config, paste0(stem, "_sessionInfo"), "txt"
  )
)
ensure_new_output(
  unlist(paths), resume = FALSE, overwrite = real_overwrite()
)

checks <- list()
record <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    scope = scope, check = name, passed = isTRUE(passed),
    detail = as.character(detail), stringsAsFactors = FALSE
  )
}

source_exists <- file.exists(processed$source_manifest$path)
record(
  "raw_sources_exist", all(source_exists),
  paste(processed$source_manifest$path[!source_exists], collapse = "|")
)
current_md5 <- rep(NA_character_, nrow(processed$source_manifest))
if (all(source_exists)) {
  current_md5 <- unname(tools::md5sum(processed$source_manifest$path))
}
record(
  "raw_sources_md5_unchanged",
  all(source_exists) && identical(current_md5, processed$source_manifest$md5),
  paste("stored/current", paste(
    processed$source_manifest$md5, current_md5, sep = "/", collapse = "|"
  ))
)

dataset <- processed$sample_dataset
record(
  "one_primary_row_per_soil_sample", !anyDuplicated(dataset$sample_id),
  paste(nrow(dataset$x), "rows and", length(unique(dataset$sample_id)),
        "unique IDs")
)
record(
  "processed_values_finite",
  all(is.finite(dataset$x)) && all(is.finite(dataset$y)),
  "Primary spectra and log1p(SOC) must be finite."
)
record(
  "exact_wavelength_grid",
  identical(as.numeric(dataset$argvals_physical), seq(1350, 2550, by = 2)) &&
    ncol(dataset$x) == 601L,
  paste(ncol(dataset$x), "columns;", range(dataset$argvals_physical)[1L],
        "to", range(dataset$argvals_physical)[2L], "nm")
)
included_manifest <- processed$sample_manifest[
  processed$sample_manifest$included, , drop = FALSE
]
response_difference <- included_manifest$response_log1p_soc -
  log1p(included_manifest$raw_soc)
record(
  "primary_response_is_explicit_log1p_soc",
  identical(processed$schema$outcome_transformation, "log1p(raw SOC)") &&
    all(is.finite(response_difference)) &&
    max(abs(response_difference)) <= 1e-12,
  paste(
    "source", processed$schema$outcome,
    "scale", processed$schema$outcome_source_scale,
    "transformation", processed$schema$outcome_transformation
  )
)
record(
  "sample_manifest_matches_primary_data",
  setequal(included_manifest$sample_id, dataset$sample_id) &&
    nrow(included_manifest) == nrow(dataset$x),
  paste(nrow(included_manifest), "included manifest samples")
)
included_rows <- processed$row_manifest[
  processed$row_manifest$included, , drop = FALSE
]
record(
  "included_rows_map_to_included_samples",
  all(included_rows$sample_id %in% dataset$sample_id),
  paste(nrow(included_rows), "included raw rows")
)
approved_exclusions <- c(
  "included", "missing_sample_id", "missing_or_nonfinite_reflectance",
  "sample_missing_or_nonfinite_soc_consensus",
  "conflicting_soc_within_sample"
)
observed_exclusions <- unique(processed$row_manifest$exclusion_reason)
record(
  "only_declared_exclusion_rules_used",
  all(observed_exclusions %in% approved_exclusions),
  paste(sort(observed_exclusions), collapse = "|")
)
record(
  "source_outcome_and_sample_columns_recorded",
  nzchar(processed$schema$outcome) && nzchar(processed$schema$sample_id),
  paste("sample", processed$schema$sample_id,
        "outcome", processed$schema$outcome)
)

is_1976 <- any(
  basename(processed$source_manifest$path) ==
    "1976_NSlibrary_withmetadata.csv"
)
if (is_1976) {
  record(
    "verified_1976_schema_mapping",
    identical(processed$schema$sample_id, "id.layer_uuid_txt") &&
      identical(processed$schema$scan_id, "id.scan_local_c") &&
      identical(processed$schema$outcome, "oc_usda.c729_w.pct") &&
      identical(processed$schema$spectral_mapping,
                "exact numeric wavelength headers"),
    paste(
      processed$schema$sample_id, processed$schema$scan_id,
      processed$schema$outcome, processed$schema$spectral_mapping,
      sep = "|"
    )
  )
  scanner_status <- processed$availability[
    processed$availability$analysis == "scanner_stability", , drop = FALSE
  ]
  country_status <- processed$availability[
    processed$availability$analysis == "geographic_external_metadata",
    , drop = FALSE
  ]
  record(
    "1976_scanner_analysis_explicitly_unavailable",
    is.null(processed$scan_dataset) && nrow(scanner_status) == 1L &&
      !scanner_status$available,
    if (nrow(scanner_status)) scanner_status$reason else "missing status"
  )
  record(
    "1976_geographic_analysis_explicitly_unavailable",
    nrow(country_status) == 1L && !country_status$available,
    if (nrow(country_status)) country_status$reason else "missing status"
  )
}

full_objects <- list()
if (!preparation_only) {
  required <- c(
    nested = neospectra_v3_result_path(config, "nested_cv", "rds"),
    simulation = neospectra_v3_result_path(
      config, "real_design_simulation", "rds"
    ),
    analysis = neospectra_v3_result_path(config, "analysis", "rds")
  )
  present <- file.exists(required)
  record(
    "complete_study_objects_exist", all(present),
    paste(basename(required[!present]), collapse = "|")
  )
  if (all(present)) {
    nested <- readRDS(required[["nested"]])
    simulation <- readRDS(required[["simulation"]])
    analysis <- readRDS(required[["analysis"]])
    full_objects <- list(
      nested = nested, simulation = simulation, analysis = analysis
    )
    grouped_checks <- validate_grouped_nested_cv_v3(
      nested, processed, config
    )
    for (index in seq_len(nrow(grouped_checks))) {
      record(
        paste0("nested_", grouped_checks$check[index]),
        grouped_checks$passed[index], grouped_checks$detail[index]
      )
    }
    outer <- nested$predictions[
      nested$predictions$evaluation == "outer_cv", , drop = FALSE
    ]
    record(
      "outer_predictions_successful_and_finite",
      all(outer$success) && all(is.finite(outer$predicted)),
      paste(sum(!outer$success | !is.finite(outer$predicted)), "failed rows")
    )
    simulation_pairs <- table(
      simulation$metrics$task_id, simulation$metrics$method
    )
    record(
      "simulation_has_one_paired_result_per_method",
      ncol(simulation_pairs) == 2L && all(simulation_pairs == 1L),
      paste(nrow(simulation_pairs), "tasks")
    )
    record(
      "simulation_splits_are_sample_disjoint",
      !any(duplicated(simulation$splits[c("task_id", "sample_id")])) &&
        all(simulation$splits$sample_id %in% dataset$sample_id),
      paste(nrow(simulation$splits), "split-membership rows")
    )
    record(
      "simulation_fits_successful",
      all(simulation$metrics$success) &&
        all(is.finite(simulation$metrics$ISE)) &&
        all(is.finite(simulation$metrics$signal_pmse)),
      paste(sum(!simulation$metrics$success), "failed fits")
    )
    nonnull <- simulation$metrics$scenario != "all_null"
    record(
      "simulation_uses_variance_snr_protocol",
      identical(
        simulation$protocol$snr_definition,
        "Var(signal) / Var(error)"
      ) &&
        all(simulation$metrics$snr_definition ==
              "Var(signal) / Var(error)") &&
        all(abs(
          simulation$metrics$noise_sd[nonnull]^2 -
            1 / simulation$metrics$target_variance_snr[nonnull]
        ) <= 1e-12) &&
        all(abs(
          simulation$metrics$training_signal_variance[nonnull] - 1
        ) <= 1e-10) &&
        all(is.finite(
          simulation$metrics$empirical_variance_snr[nonnull]
        )),
      "Signal is training-SD normalized and noise SD is 1/sqrt(target SNR)."
    )
    no_null <- simulation$metrics$scenario == "no_null"
    all_null <- simulation$metrics$scenario == "all_null"
    local_null <- simulation$metrics$scenario %in%
      c("interior_null", "moderate_shoulder")
    record(
      "simulation_support_metrics_follow_scenario_null_sets",
      all(simulation$metrics$true_null_fraction[no_null] == 0) &&
        all(is.na(simulation$metrics$true_null_recovery[no_null])) &&
        all(is.finite(simulation$metrics$false_null_rate[no_null])) &&
        all(simulation$metrics$true_null_fraction[all_null] == 1) &&
        all(is.finite(simulation$metrics$true_null_recovery[all_null])) &&
        all(is.na(simulation$metrics$false_null_rate[all_null])) &&
        all(simulation$metrics$true_null_fraction[local_null] > 0 &
              simulation$metrics$true_null_fraction[local_null] < 1) &&
        all(is.finite(
          simulation$metrics$true_null_recovery[local_null]
        )) &&
        all(is.finite(simulation$metrics$false_null_rate[local_null])),
      paste(
        "Isolated zero crossings are not null intervals; all-null has no",
        "nonnull region."
      )
    )
    record(
      "simulation_all_null_snr_provenance_is_explicit",
      all(is.na(simulation$metrics$snr[all_null])) &&
        all(is.na(simulation$metrics$target_variance_snr[all_null])) &&
        all(simulation$metrics$noise_sd[all_null] == 1) &&
        all(grepl("snrallnull", simulation$metrics$task_id[all_null],
                  fixed = TRUE)),
      "All-null tasks must use unit noise and an undefined target SNR."
    )
    curve_fmcp <- simulation$curves[
      simulation$curves$method == "FMCP",
      c("task_id", "wavelength", "beta_truth")
    ]
    curve_fscad <- simulation$curves[
      simulation$curves$method == "fSCAD",
      c("task_id", "wavelength", "beta_truth")
    ]
    truth_pair <- merge(
      curve_fmcp, curve_fscad, by = c("task_id", "wavelength"),
      suffixes = c("_FMCP", "_fSCAD"), all = TRUE
    )
    record(
      "simulation_truth_is_method_identical",
      nrow(truth_pair) == nrow(curve_fmcp) &&
        nrow(truth_pair) == nrow(curve_fscad) &&
        all(truth_pair$beta_truth_FMCP == truth_pair$beta_truth_fSCAD),
      paste(nrow(truth_pair), "paired curve points")
    )
    selected_methods <- sort(analysis$selected_candidates$method)
    record(
      "global_candidates_selected_for_both_sparse_methods",
      identical(selected_methods, sort(config$primary_sparse_methods)),
      paste(selected_methods, collapse = "|")
    )
    expected_curve_rows <- length(config$primary_sparse_methods) *
      length(dataset$argvals_physical)
    record(
      "full_coefficient_curves_complete",
      nrow(analysis$full$curves) == expected_curve_rows &&
        all(is.finite(analysis$full$curves$beta)),
      paste(nrow(analysis$full$curves), "rows")
    )
    expected_gammas <- sort(unique(c(
      config$mcp_gamma, config$mcp_sensitivity
    )))
    observed_gammas <- sort(unique(
      analysis$mcp_gamma_sensitivity$selection$mcp_gamma
    ))
    gamma_curves <- analysis$mcp_gamma_sensitivity$curves
    record(
      "mcp_gamma_sensitivity_complete",
      identical(observed_gammas, expected_gammas) &&
        all(analysis$mcp_gamma_sensitivity$selection$success) &&
        nrow(gamma_curves) == length(expected_gammas) *
          length(dataset$argvals_physical) &&
        all(is.finite(gamma_curves$beta)),
      paste("observed gammas", paste(observed_gammas, collapse = "|"))
    )
    draws_per_replication <- table(analysis$stability$draws$replication)
    record(
      "sample_bootstrap_draws_have_primary_sample_size",
      length(draws_per_replication) == config$bootstrap_replicates &&
        all(draws_per_replication == nrow(dataset$x)) &&
        all(analysis$stability$draws$sample_id %in% dataset$sample_id),
      paste(length(draws_per_replication), "replications")
    )
    scanner_curves <- analysis$stability$curves[
      analysis$stability$curves$source == "scanner", , drop = FALSE
    ]
    scanner_available <- processed$availability$available[
      processed$availability$analysis == "scanner_stability"
    ]
    record(
      "scanner_stability_respects_metadata_availability",
      if (length(scanner_available) && scanner_available) {
        nrow(scanner_curves) > 0L
      } else nrow(scanner_curves) == 0L,
      paste(nrow(scanner_curves), "scanner curve rows")
    )
    successful_runtime <- analysis$runtime_detail[
      analysis$runtime_detail$success, , drop = FALSE
    ]
    record(
      "successful_runtime_rows_are_finite",
      nrow(successful_runtime) > 0L &&
        all(is.finite(successful_runtime$elapsed_seconds)),
      paste(nrow(successful_runtime), "successful runtime rows")
    )
  }
}

inventory_paths <- unique(c(
  neospectra_v3_processed_path(config),
  list.files(config$results_dir, recursive = TRUE, full.names = TRUE),
  list.files(config$figures_dir, recursive = TRUE, full.names = TRUE)
))
inventory_paths <- inventory_paths[file.exists(inventory_paths) &
                                     !dir.exists(inventory_paths)]
inventory <- data.frame(
  path = normalizePath(inventory_paths),
  bytes = as.numeric(file.info(inventory_paths)$size),
  md5 = unname(tools::md5sum(inventory_paths)), stringsAsFactors = FALSE
)
v3_roots <- normalizePath(c(
  config$processed_dir, config$results_dir, config$figures_dir
), mustWork = FALSE)
under_v3 <- vapply(inventory$path, function(path) {
  any(startsWith(path, paste0(v3_roots, .Platform$file.sep))) ||
    path %in% v3_roots
}, FUN.VALUE = logical(1L))
record(
  "all_generated_inventory_is_v3_scoped", all(under_v3),
  paste(inventory$path[!under_v3], collapse = "|")
)

check_table <- bind_rows_base(checks)
report <- list(
  version = config$version, mode = mode, scope = scope,
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  checks = check_table, inventory = inventory,
  exclusion_counts = as.data.frame(table(
    processed$row_manifest$exclusion_reason
  ), stringsAsFactors = FALSE),
  analysis_availability = processed$availability,
  configuration = config
)
atomic_save_rds(report, paths$report)
atomic_write_csv(check_table, paths$checks)
atomic_write_csv(inventory, paths$inventory)
capture_session_info(paths$session)

if (!all(check_table$passed)) {
  stop(
    "NeoSpectra v3 validation failed: ",
    paste(check_table$check[!check_table$passed], collapse = ", "),
    call. = FALSE
  )
}
cat("NeoSpectra v3", scope, "validation passed", nrow(check_table),
    "checks.\n")
