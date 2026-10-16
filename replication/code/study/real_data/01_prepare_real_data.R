script_dir <- real_data_script_directory <- local({
  file_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
  if (length(file_argument)) {
    dirname(normalizePath(sub("^--file=", "", file_argument[1L])))
  } else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
source(file.path(script_dir, "R", "utils.R"))
source(file.path(script_dir, "R", "config.R"))
source(file.path(script_dir, "R", "data.R"))
activate_fmcp_libraries()
assert_packages("caret")

config <- real_data_config("full", project_root)
dir.create(config$processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)

processed_paths <- vapply(
  config$datasets, function(name) processed_dataset_path(config, name),
  FUN.VALUE = character(1L)
)
manifest_path <- file.path(config$results_dir, "data_manifest_v1.csv")
dictionary_path <- file.path(config$results_dir, "data_dictionary_v1.csv")
session_path <- file.path(config$results_dir, "preparation_sessionInfo_v1.txt")
ensure_new_output(
  c(processed_paths, manifest_path, dictionary_path, session_path),
  overwrite = real_overwrite()
)

datasets <- list(
  Truck = load_truck_source(project_root),
  Tecator = load_tecator_source()
)
datasets <- lapply(datasets, validate_real_dataset)

for (name in names(datasets)) {
  dataset <- datasets[[name]]
  dataset$preparation <- list(
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    raw_values_changed = FALSE,
    exclusions = "None",
    missing_value_handling = "None required; no values were missing",
    outlier_handling = "None",
    transformations = paste(
      "No response or predictor transformation at preparation;",
      "a [0,1] analysis grid was stored alongside physical labels"
    ),
    response_choice = if (identical(name, "Tecator")) {
      "Second endpoint column: fat content (%)"
    } else {
      "Distributed Y column: particulate-matter emission"
    }
  )
  atomic_save_rds(dataset, processed_paths[[name]])
}

manifest <- bind_rows_base(Map(
  dataset_manifest_row, datasets, unname(processed_paths)
))
atomic_write_csv(manifest, manifest_path)

dictionary <- data.frame(
  dataset = rep(names(datasets), each = 6L),
  field = rep(c(
    "subject_id", "y", "x", "argvals_u", "argvals_physical", "notes"
  ), times = length(datasets)),
  definition = c(
    "Sequential row identifier; no source identifier was altered",
    datasets$Truck$response_description,
    datasets$Truck$predictor_description,
    "Equally spaced computational domain [0,1]",
    "Measurement indices 1--90 (seconds in the source description)",
    datasets$Truck$notes,
    "Sequential sample identifier; no source identifier was supplied",
    datasets$Tecator$response_description,
    datasets$Tecator$predictor_description,
    "Equally spaced computational domain [0,1]",
    "Nominal linear map over documented 850--1050 nm range",
    datasets$Tecator$notes
  ),
  stringsAsFactors = FALSE
)
atomic_write_csv(dictionary, dictionary_path)
capture_session_info(session_path)

cat("Prepared immutable analysis copies:\n")
print(manifest[, c("dataset", "n", "grid_points", "processed_path")])

