load_truck_source <- function(project_root) {
  source_path <- file.path(project_root, "PFLR-master", "data", "truck.rda")
  source_environment <- new.env(parent = emptyenv())
  load(source_path, envir = source_environment)
  if (!exists("truck", envir = source_environment, inherits = FALSE)) {
    stop("truck.rda did not contain `truck`.", call. = FALSE)
  }
  source_data <- get("truck", envir = source_environment)
  column_names <- colnames(source_data)
  if (is.null(column_names)) column_names <- names(source_data)
  x <- as.matrix(source_data[, grep("^X", column_names), drop = FALSE])
  storage.mode(x) <- "double"
  y <- as.numeric(source_data[, match("Y", column_names)])
  list(
    name = "Truck",
    subject_id = seq_along(y),
    y = y,
    x = x,
    argvals_u = seq(0, 1, length.out = ncol(x)),
    argvals_physical = seq_len(ncol(x)),
    x_axis_label = "Acceleration measurement index (seconds)",
    y_axis_label = "Particulate-matter emission",
    response_description = "Particulate-matter emission",
    predictor_description = "Engine acceleration recorded at each second",
    source = paste(
      "PFLR truck data; CRC E55/E59 heavy-duty vehicle chassis",
      "dynamometer project (Clark et al., 2007)"
    ),
    source_path = normalizePath(source_path),
    source_md5 = unname(tools::md5sum(source_path)),
    notes = paste(
      "No truck/test-cycle identifier is supplied in the distributed data;",
      "rows are therefore treated as observational units."
    )
  )
}

load_tecator_source <- function() {
  source_environment <- new.env(parent = emptyenv())
  utils::data("tecator", package = "caret", envir = source_environment)
  if (!all(c("absorp", "endpoints") %in% ls(source_environment))) {
    stop("caret::tecator was not loaded as expected.", call. = FALSE)
  }
  x <- as.matrix(get("absorp", envir = source_environment))
  storage.mode(x) <- "double"
  endpoints <- as.matrix(get("endpoints", envir = source_environment))
  y <- as.numeric(endpoints[, 2L])
  list(
    name = "Tecator",
    subject_id = seq_along(y),
    y = y,
    x = x,
    argvals_u = seq(0, 1, length.out = ncol(x)),
    argvals_physical = seq(850, 1050, length.out = ncol(x)),
    x_axis_label = "Nominal wavelength (nm)",
    y_axis_label = "Fat content (%)",
    response_description = "Analytic-chemistry fat content (percent)",
    predictor_description = paste(
      "Near-infrared transmission absorbance (-log10 transmittance),",
      "100 channels over the documented 850--1050 nm range"
    ),
    source = "caret::tecator (public-domain Tecator data)",
    source_path = "R package caret, data object tecator",
    source_md5 = NA_character_,
    notes = paste(
      "All 215 samples are used in repeated cross-validation; the documented",
      "first-129 historical training split is not imposed. Nominal wavelengths",
      "are linearly mapped across the documented range because channel-specific",
      "wavelength labels are not included in caret."
    )
  )
}

validate_real_dataset <- function(dataset) {
  required <- c("name", "subject_id", "y", "x", "argvals_u",
                "argvals_physical")
  missing <- setdiff(required, names(dataset))
  if (length(missing)) {
    stop("Dataset is missing fields: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (!is.matrix(dataset$x) || !is.numeric(dataset$x) ||
      length(dataset$y) != nrow(dataset$x) ||
      length(dataset$argvals_u) != ncol(dataset$x) ||
      length(dataset$argvals_physical) != ncol(dataset$x)) {
    stop("Dataset dimensions are inconsistent for ", dataset$name,
         call. = FALSE)
  }
  if (anyNA(dataset$x) || any(!is.finite(dataset$x)) ||
      anyNA(dataset$y) || any(!is.finite(dataset$y))) {
    stop("Non-finite values found in ", dataset$name, call. = FALSE)
  }
  if (any(diff(dataset$argvals_u) <= 0) ||
      any(diff(dataset$argvals_physical) <= 0)) {
    stop("Argument grids must be strictly increasing.", call. = FALSE)
  }
  dataset
}

processed_dataset_path <- function(config, dataset_name) {
  file.path(
    config$processed_dir,
    paste0(tolower(dataset_name), "_processed_v1.rds")
  )
}

load_processed_datasets <- function(config) {
  paths <- vapply(
    config$datasets, function(name) processed_dataset_path(config, name),
    FUN.VALUE = character(1L)
  )
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Processed data are missing: ", paste(basename(missing), collapse = ", "),
      ". Run 01_prepare_real_data.R first.", call. = FALSE
    )
  }
  datasets <- lapply(paths, readRDS)
  names(datasets) <- config$datasets
  lapply(datasets, validate_real_dataset)
}

dataset_manifest_row <- function(dataset, processed_path) {
  data.frame(
    dataset = dataset$name,
    n = length(dataset$y),
    grid_points = ncol(dataset$x),
    response = dataset$response_description,
    predictor = dataset$predictor_description,
    physical_grid_min = min(dataset$argvals_physical),
    physical_grid_max = max(dataset$argvals_physical),
    missing_y = sum(is.na(dataset$y)),
    missing_x = sum(is.na(dataset$x)),
    transformations_at_preparation = "None",
    analysis_domain = "[0,1] linear map; physical labels retained",
    source = dataset$source,
    source_path = dataset$source_path,
    source_md5 = dataset$source_md5,
    processed_path = normalizePath(processed_path, mustWork = FALSE),
    processed_md5 = unname(tools::md5sum(processed_path)),
    notes = dataset$notes,
    stringsAsFactors = FALSE
  )
}
