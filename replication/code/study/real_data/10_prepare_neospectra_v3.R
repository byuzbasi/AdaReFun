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

mode <- real_mode()
config <- neospectra_v3_config(mode, project_root)
inputs <- resolve_neospectra_v3_inputs(config)
paths <- list(
  processed = neospectra_v3_processed_path(config),
  row_manifest = neospectra_v3_result_path(
    config, "preparation_row_manifest", "csv"
  ),
  sample_manifest = neospectra_v3_result_path(
    config, "preparation_sample_manifest", "csv"
  ),
  source_manifest = neospectra_v3_result_path(
    config, "preparation_source_manifest", "csv"
  ),
  quality = neospectra_v3_result_path(
    config, "preparation_quality", "csv"
  ),
  transformations = neospectra_v3_result_path(
    config, "preparation_transformations", "csv"
  ),
  schema = neospectra_v3_result_path(
    config, "preparation_expected_schema", "csv"
  ),
  schema_audit = neospectra_v3_result_path(
    config, "preparation_schema_audit", "csv"
  ),
  availability = neospectra_v3_result_path(
    config, "preparation_analysis_availability", "csv"
  ),
  session = neospectra_v3_result_path(
    config, "preparation_sessionInfo", "txt"
  )
)
ensure_new_output(
  unlist(paths), resume = FALSE, overwrite = real_overwrite()
)

processed <- prepare_neospectra_v3(inputs, config)
processed$configuration <- config

schema <- rbind(
  data.frame(
    role = c("sample_id", "scan_id", "scanner_id", "outcome"),
    exact_column = c(
      paste(c(config$sample_id_column, config$sample_id_aliases),
            collapse = "|"),
      paste(c(config$scan_id_column, config$scan_id_aliases), collapse = "|"),
      paste(c(config$scanner_column, config$scanner_aliases), collapse = "|"),
      paste(c(config$outcome_column, config$raw_soc_aliases,
              config$log1p_soc_aliases), collapse = "|")
    ),
    required = c(TRUE, FALSE, FALSE, TRUE),
    source_table = c("NIR/combined", "NIR/combined", "NIR/combined",
                     "laboratory/combined"),
    notes = c(
      "Grouping unit; one final row per soil sample",
      "Optional provenance; never used as the independent grouping unit",
      "Optional; required only for scanner-specific stability",
      paste0(
        "Primary response is log1p(SOC); actual source=",
        processed$schema$outcome, "; source scale=",
        processed$schema$outcome_source_scale
      )
    ),
    stringsAsFactors = FALSE
  ),
  data.frame(
    role = "country_candidate", exact_column = config$country_candidates,
    required = FALSE, source_table = "site/combined",
    notes = "First available candidate; needed only for US-to-Africa validation",
    stringsAsFactors = FALSE
  ),
  data.frame(
    role = "reflectance", exact_column = processed$schema$spectral_columns,
    required = TRUE, source_table = "NIR/combined",
    notes = "Raw reflectance retained; no interpolation or transformation",
    stringsAsFactors = FALSE
  )
)

atomic_save_rds(processed, paths$processed)
atomic_write_csv(processed$row_manifest, paths$row_manifest)
atomic_write_csv(processed$sample_manifest, paths$sample_manifest)
atomic_write_csv(processed$source_manifest, paths$source_manifest)
atomic_write_csv(processed$quality, paths$quality)
atomic_write_csv(processed$transformations, paths$transformations)
atomic_write_csv(schema, paths$schema)
atomic_write_csv(processed$schema_audit, paths$schema_audit)
atomic_write_csv(processed$availability, paths$availability)
capture_session_info(paths$session)

cat(
  "NeoSpectra v3 preparation completed:",
  nrow(processed$sample_dataset$x), "unique soil samples,",
  if (is.null(processed$scan_dataset)) 0L else
    nrow(processed$scan_dataset$x), "sample-by-scanner rows,",
  ncol(processed$sample_dataset$x), "wavelengths.\n"
)
