# NeoSpectra real-data and real-design simulation helpers (version 3).
#
# This file intentionally contains no download code. Raw files are treated as
# immutable inputs and every derived artifact is written under a v3 directory.

neospectra_v3_env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < minimum) {
    stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
  }
  parsed
}

neospectra_v3_config <- function(mode = c("full", "smoke"),
                                 project_root = NULL) {
  mode <- match.arg(mode)
  if (is.null(project_root)) project_root <- find_fmcp_project_root()
  real_dir <- file.path(project_root, "real_data")
  smoke <- identical(mode, "smoke")
  list(
    version = "neospectra_v3",
    mode = mode,
    master_seed = 20260821L,
    raw_dir = file.path(real_dir, "data", "raw", "neospectra_v3"),
    raw_search_dirs = c(
      file.path(real_dir, "data", "raw", "neospectra_v3"),
      file.path(real_dir, "data", "raw")
    ),
    processed_dir = file.path(real_dir, "data", "processed", "v3"),
    results_dir = file.path(real_dir, "results", "v3", "neospectra"),
    figures_dir = file.path(real_dir, "figures", "v3", "neospectra"),
    raw_combined_candidates = c(
      "1976_NSlibrary_withmetadata.csv",
      "Neospectra_WoodwellKSSL_soil+site+NIR.csv",
      "Neospectra_WoodwellKSSL_avg_soil+site+NIR.csv"
    ),
    raw_nir_candidates = c(
      "neospectra_nir_v1.2.csv.gz", "neospectra_nir_v1.2.csv"
    ),
    raw_lab_candidates = c(
      "neospectra_soillab_v1.2.csv.gz", "neospectra_soillab_v1.2.csv"
    ),
    raw_site_candidates = c(
      "neospectra_soilsite_v1.2.csv.gz", "neospectra_soilsite_v1.2.csv"
    ),
    sample_id_column = "id.sample_local_c",
    sample_id_aliases = c(
      "id.layer_uuid_txt", "kssl_id", "sample_id", "sample.id"
    ),
    scan_id_column = "id.scan_local_c",
    scan_id_aliases = c("scan_id", "scan.id"),
    scanner_column = "scan.nir.model.name_utf8_txt",
    scanner_aliases = c("scanner_name", "scanner", "instrument_name"),
    scanner_serial_column = "scan.nir.model.serialnumber_utf8_in",
    outcome_column = "oc_usda.c729_w.pct",
    raw_soc_aliases = c(
      "SOC", "soc", "soc_percent", "SOC_percent",
      "soil_organic_carbon", "soil_organic_carbon_percent"
    ),
    log1p_soc_aliases = c(
      "log1p_SOC", "log1p_soc", "SOC_log1p", "soc_log1p"
    ),
    ambiguous_log_soc_aliases = c(
      "log_SOC", "log_soc", "logSOC", "logsoc", "ln_SOC", "ln_soc"
    ),
    country_candidates = c(
      "loc.country_src_code", "location.country_iso.3166_txt",
      "country_iso", "country"
    ),
    wavelengths = seq(1350, 2550, by = 2),
    spectral_prefix = "scan_nir.",
    spectral_suffix = "_ref",
    predictor_transform = "none; source reflectance fractions retained",
    methods = c("FMCP", "fSCAD", "Smooth", "PCR"),
    primary_sparse_methods = c("FMCP", "fSCAD"),
    M = if (smoke) c(10L) else c(10L, 20L, 30L),
    degree = 3L,
    penalty_order = 2L,
    lambda = if (smoke) c(1e-3, 1e-1) else 10^seq(-4, 0, length.out = 5L),
    kappa = if (smoke) c(1e-8, 1e-4) else 10^seq(-10, -2, length.out = 4L),
    mcp_gamma = 3.7,
    mcp_sensitivity = c(3.0, 2.7),
    scad_a = 3.7,
    max_iter = if (smoke) 40L else 120L,
    tolerance = 1e-8,
    norm_floor = 1e-8,
    relative_cutoff_multiplier = 1e-4,
    zero_tolerance = sqrt(.Machine$double.eps),
    integration = "trapezoid",
    pcr_max_components = if (smoke) 5L else 30L,
    outer_repeats = if (smoke) 1L else 2L,
    outer_folds = if (smoke) 2L else 5L,
    inner_folds = if (smoke) 2L else 3L,
    bootstrap_replicates = neospectra_v3_env_integer(
      "FMCP_NEOSPECTRA_BOOTSTRAPS", if (smoke) 3L else 200L
    ),
    simulation_replicates = neospectra_v3_env_integer(
      "FMCP_NEOSPECTRA_SIM_REPS", if (smoke) 2L else 100L
    ),
    simulation_n = if (smoke) c(80L) else c(150L, 450L, 1000L),
    simulation_snr = if (smoke) c(3) else c(2, 5),
    simulation_snr_definition = "Var(signal) / Var(error)",
    simulation_scenarios = c(
      "interior_null", "moderate_shoulder", "no_null", "all_null"
    ),
    workers = neospectra_v3_env_integer(
      "FMCP_NEOSPECTRA_WORKERS",
      if (.Platform$OS.type == "windows") 1L else 4L
    )
  )
}

neospectra_v3_name <- function(stem, config, extension = "rds") {
  suffix <- if (identical(config$mode, "smoke")) "_smoke" else ""
  paste0(stem, "_v3", suffix, ".", extension)
}

neospectra_v3_result_path <- function(config, stem, extension = "rds") {
  file.path(config$results_dir, neospectra_v3_name(stem, config, extension))
}

neospectra_v3_figure_path <- function(config, stem, extension = "pdf") {
  file.path(config$figures_dir, neospectra_v3_name(stem, config, extension))
}

neospectra_v3_processed_path <- function(config) {
  file.path(config$processed_dir, neospectra_v3_name(
    "neospectra_soc_processed", config, "rds"
  ))
}

neospectra_v3_checkpoint_signature <- function(config, processed,
                                                component) {
  real_dir <- normalizePath(
    file.path(config$raw_dir, "..", "..", ".."), mustWork = TRUE
  )
  code_files <- c(
    "R/utils.R", "R/methods.R", "R/neospectra_v3.R",
    "10_prepare_neospectra_v3.R", "11_run_neospectra_nested_cv_v3.R",
    "12_run_neospectra_real_design_sim_v3.R",
    "13_analyze_neospectra_v3.R", "14_validate_neospectra_v3.R"
  )
  code_paths <- file.path(real_dir, code_files)
  if (any(!file.exists(code_paths))) {
    stop("NeoSpectra v3 checkpoint code inventory is incomplete.",
         call. = FALSE)
  }
  code_md5 <- unname(tools::md5sum(code_paths))
  names(code_md5) <- code_files
  list(
    component = component, version = config$version, mode = config$mode,
    master_seed = config$master_seed,
    source_md5 = processed$source_manifest$md5,
    n_samples = nrow(processed$sample_dataset$x),
    sample_ids = processed$sample_dataset$sample_id,
    wavelengths = processed$sample_dataset$argvals_physical,
    code_md5 = code_md5,
    methods = config$methods,
    primary_sparse_methods = config$primary_sparse_methods,
    M = config$M, lambda = config$lambda, kappa = config$kappa,
    mcp_gamma = config$mcp_gamma,
    mcp_sensitivity = config$mcp_sensitivity, scad_a = config$scad_a,
    degree = config$degree, penalty_order = config$penalty_order,
    max_iter = config$max_iter, tolerance = config$tolerance,
    norm_floor = config$norm_floor,
    relative_cutoff_multiplier = config$relative_cutoff_multiplier,
    zero_tolerance = config$zero_tolerance,
    integration = config$integration,
    pcr_max_components = config$pcr_max_components,
    outer_repeats = config$outer_repeats,
    outer_folds = config$outer_folds, inner_folds = config$inner_folds,
    simulation_n = config$simulation_n,
    simulation_snr = config$simulation_snr,
    simulation_snr_definition = config$simulation_snr_definition,
    simulation_scenarios = config$simulation_scenarios,
    simulation_replicates = config$simulation_replicates,
    bootstrap_replicates = config$bootstrap_replicates
  )
}

neospectra_v3_seed <- function(config, component, index = 0L,
                               subindex = 0L) {
  offsets <- c(
    preparation = 100000L, outer = 200000L, inner = 300000L,
    external = 400000L, simulation = 500000L, stability = 600000L
  )
  if (!component %in% names(offsets)) {
    stop("Unknown NeoSpectra seed component: ", component, call. = FALSE)
  }
  as.integer(config$master_seed + offsets[[component]] +
               as.integer(index) * 100L + as.integer(subindex))
}

neospectra_v3_existing_override <- function(environment_name) {
  value <- Sys.getenv(environment_name, unset = "")
  if (!nzchar(value)) return(NULL)
  path <- normalizePath(value, mustWork = FALSE)
  if (!file.exists(path)) {
    stop(environment_name, " points to a missing file: ", path,
         call. = FALSE)
  }
  path
}

neospectra_v3_unique_candidate <- function(raw_dir, candidates, label,
                                           required = TRUE) {
  paths <- file.path(raw_dir, candidates)
  existing <- normalizePath(paths[file.exists(paths)], mustWork = TRUE)
  if (length(existing) > 1L) {
    stop("Multiple ", label, " candidates found; keep one or use the explicit ",
         "environment override. Files: ", paste(basename(existing), collapse = ", "),
         call. = FALSE)
  }
  if (!length(existing)) {
    if (required) stop("No ", label, " file found in ", raw_dir, call. = FALSE)
    return(NULL)
  }
  existing
}

neospectra_v3_unique_candidate_dirs <- function(raw_dirs, candidates, label,
                                                required = TRUE) {
  paths <- unique(unlist(lapply(raw_dirs, function(directory) {
    file.path(directory, candidates)
  }), use.names = FALSE))
  existing <- normalizePath(paths[file.exists(paths)], mustWork = TRUE)
  if (length(existing) > 1L) {
    stop(
      "Multiple ", label, " candidates found; use an explicit environment ",
      "override. Files: ", paste(existing, collapse = ", "), call. = FALSE
    )
  }
  if (!length(existing)) {
    if (required) {
      stop("No ", label, " file found in the configured raw directories.",
           call. = FALSE)
    }
    return(NULL)
  }
  existing
}

resolve_neospectra_v3_inputs <- function(config) {
  combined_override <- neospectra_v3_existing_override(
    "FMCP_NEOSPECTRA_COMBINED"
  )
  nir_override <- neospectra_v3_existing_override("FMCP_NEOSPECTRA_NIR")
  lab_override <- neospectra_v3_existing_override("FMCP_NEOSPECTRA_LAB")
  site_override <- neospectra_v3_existing_override("FMCP_NEOSPECTRA_SITE")

  if (!is.null(combined_override) &&
      (!is.null(nir_override) || !is.null(lab_override))) {
    stop("Choose either a combined file or separate NIR/lab files, not both.",
         call. = FALSE)
  }
  if (!is.null(combined_override)) {
    return(list(mode = "combined", combined = combined_override,
                nir = NULL, lab = NULL, site = NULL))
  }
  if (!is.null(nir_override) || !is.null(lab_override)) {
    if (is.null(nir_override) || is.null(lab_override)) {
      stop("Both FMCP_NEOSPECTRA_NIR and FMCP_NEOSPECTRA_LAB are required.",
           call. = FALSE)
    }
    return(list(mode = "separate", combined = NULL, nir = nir_override,
                lab = lab_override, site = site_override))
  }

  combined <- neospectra_v3_unique_candidate_dirs(
    config$raw_search_dirs, config$raw_combined_candidates,
    "combined NeoSpectra",
    required = FALSE
  )
  nir <- neospectra_v3_unique_candidate_dirs(
    config$raw_search_dirs, config$raw_nir_candidates, "NeoSpectra NIR",
    required = FALSE
  )
  lab <- neospectra_v3_unique_candidate_dirs(
    config$raw_search_dirs, config$raw_lab_candidates,
    "NeoSpectra laboratory",
    required = FALSE
  )
  site <- neospectra_v3_unique_candidate_dirs(
    config$raw_search_dirs, config$raw_site_candidates, "NeoSpectra site",
    required = FALSE
  )
  if (!is.null(combined) && (!is.null(nir) || !is.null(lab))) {
    stop("Both combined and separate raw inputs were found. Remove ambiguity ",
         "or set explicit environment overrides.", call. = FALSE)
  }
  if (!is.null(combined)) {
    return(list(mode = "combined", combined = combined,
                nir = NULL, lab = NULL, site = NULL))
  }
  if (is.null(nir) || is.null(lab)) {
    stop(
      "NeoSpectra raw data are not present. Place the official combined CSV, ",
      "or the official neospectra_nir_v1.2 and neospectra_soillab_v1.2 ",
      "CSV/CSV.GZ files in one configured raw directory. No downloader is ",
      "provided. Searched: ", paste(config$raw_search_dirs, collapse = ", "),
      call. = FALSE
    )
  }
  list(mode = "separate", combined = NULL, nir = nir, lab = lab, site = site)
}

read_neospectra_v3_csv <- function(path) {
  utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE,
    na.strings = c("", "NA", "NaN", "NULL")
  )
}

neospectra_v3_require_column <- function(data, column, table_name) {
  if (!column %in% names(data)) {
    stop("Expected column `", column, "` is missing from ", table_name, ".",
         call. = FALSE)
  }
  column
}

neospectra_v3_optional_column <- function(data, candidates) {
  found <- candidates[candidates %in% names(data)]
  if (!length(found)) return(NULL)
  found[1L]
}

neospectra_v3_resolve_alias <- function(data, role, preferred, aliases,
                                        override_environment,
                                        required = TRUE) {
  override <- Sys.getenv(override_environment, unset = "")
  candidates <- unique(c(preferred, aliases))
  if (nzchar(override)) {
    if (!override %in% names(data)) {
      stop(
        override_environment, "=`", override, "` but that column is missing ",
        "for ", role, ".", call. = FALSE
      )
    }
    selected <- override
    detection <- paste0("environment:", override_environment)
    matches <- override
  } else {
    matches <- candidates[candidates %in% names(data)]
    if (length(matches) > 1L) {
      stop(
        "Ambiguous ", role, " columns: ", paste(matches, collapse = ", "),
        ". Set ", override_environment, " explicitly.", call. = FALSE
      )
    }
    if (!length(matches)) {
      if (required) {
        stop(
          "No recognized ", role, " column. Expected one of: ",
          paste(candidates, collapse = ", "), ". Or set ",
          override_environment, ".", call. = FALSE
        )
      }
      selected <- NA_character_
      detection <- "unavailable"
    } else {
      selected <- matches
      detection <- "exact_alias"
    }
  }
  list(
    column = if (is.na(selected)) NULL else selected,
    audit = data.frame(
      role = role, required = required,
      selected_column = if (is.na(selected)) "" else selected,
      detection = detection,
      recognized_candidates = paste(candidates, collapse = "|"),
      stringsAsFactors = FALSE
    )
  )
}

neospectra_v3_resolve_outcome <- function(data, config, table_name) {
  override <- Sys.getenv("FMCP_NEOSPECTRA_SOC_COLUMN", unset = "")
  override_scale <- tolower(Sys.getenv(
    "FMCP_NEOSPECTRA_SOC_SCALE", unset = ""
  ))
  raw_candidates <- unique(c(config$outcome_column, config$raw_soc_aliases))
  log1p_candidates <- config$log1p_soc_aliases
  ambiguous_candidates <- config$ambiguous_log_soc_aliases
  if (nzchar(override)) {
    if (!override %in% names(data)) {
      stop(
        "FMCP_NEOSPECTRA_SOC_COLUMN=`", override,
        "` is missing from ", table_name, ".", call. = FALSE
      )
    }
    column <- override
    scale <- override_scale
    if (!nzchar(scale)) {
      if (column %in% raw_candidates) scale <- "raw"
      if (column %in% log1p_candidates) scale <- "log1p"
    }
    if (!nzchar(scale)) {
      stop(
        "The scale of outcome column `", column, "` is ambiguous. Set ",
        "FMCP_NEOSPECTRA_SOC_SCALE to raw, log1p, log, or log10.",
        call. = FALSE
      )
    }
    detection <- "environment"
  } else {
    raw_found <- raw_candidates[raw_candidates %in% names(data)]
    log1p_found <- log1p_candidates[log1p_candidates %in% names(data)]
    ambiguous_found <- ambiguous_candidates[
      ambiguous_candidates %in% names(data)
    ]
    if (length(raw_found) > 1L ||
        (!length(raw_found) && length(log1p_found) > 1L)) {
      found <- if (length(raw_found)) raw_found else log1p_found
      stop(
        "Ambiguous SOC outcome columns in ", table_name, ": ",
        paste(found, collapse = ", "),
        ". Set FMCP_NEOSPECTRA_SOC_COLUMN explicitly.", call. = FALSE
      )
    }
    if (length(raw_found) == 1L) {
      column <- raw_found
      scale <- "raw"
    } else if (length(log1p_found) == 1L) {
      column <- log1p_found
      scale <- "log1p"
    } else if (length(ambiguous_found)) {
      stop(
        "Only ambiguously transformed SOC column(s) were found in ",
        table_name, ": ", paste(ambiguous_found, collapse = ", "),
        ". Set FMCP_NEOSPECTRA_SOC_COLUMN and ",
        "FMCP_NEOSPECTRA_SOC_SCALE explicitly.", call. = FALSE
      )
    } else {
      stop(
        "No recognized SOC column in ", table_name, ". Expected raw: ",
        paste(raw_candidates, collapse = ", "), "; or explicit log1p: ",
        paste(log1p_candidates, collapse = ", "), ".", call. = FALSE
      )
    }
    detection <- "exact_alias"
  }
  allowed <- c("raw", "log1p", "log", "log10")
  if (!scale %in% allowed) {
    stop("FMCP_NEOSPECTRA_SOC_SCALE must be one of: ",
         paste(allowed, collapse = ", "), ".", call. = FALSE)
  }
  source_value <- suppressWarnings(as.numeric(data[[column]]))
  if (identical(scale, "raw")) {
    raw_soc <- source_value
    response <- suppressWarnings(log1p(raw_soc))
  } else if (identical(scale, "log1p")) {
    response <- source_value
    raw_soc <- suppressWarnings(expm1(source_value))
  } else if (identical(scale, "log")) {
    raw_soc <- suppressWarnings(exp(source_value))
    response <- suppressWarnings(log1p(raw_soc))
  } else {
    raw_soc <- suppressWarnings(10^source_value)
    response <- suppressWarnings(log1p(raw_soc))
  }
  status <- rep("usable", length(source_value))
  status[!is.finite(source_value)] <- "missing_or_nonfinite_source_soc"
  status[is.finite(raw_soc) & raw_soc < 0] <- "negative_soc"
  status[status == "usable" &
           (!is.finite(raw_soc) | !is.finite(response))] <-
    "nonfinite_after_soc_transformation"
  list(
    raw_soc = raw_soc, response = response, status = status,
    column = column, source_scale = scale,
    transformation = if (scale == "raw") "log1p(raw SOC)" else
      if (scale == "log1p") "source already log1p(SOC)" else
        paste0("recover raw SOC from ", scale, "; then log1p"),
    audit = data.frame(
      role = "outcome", required = TRUE, selected_column = column,
      detection = paste(detection, scale, sep = ":"),
      recognized_candidates = paste(
        c(raw_candidates, log1p_candidates, ambiguous_candidates),
        collapse = "|"
      ), stringsAsFactors = FALSE
    )
  )
}

neospectra_v3_spectral_columns <- function(data, config) {
  expected <- paste0(
    config$spectral_prefix, config$wavelengths, config$spectral_suffix
  )
  if (all(expected %in% names(data))) {
    attr(expected, "mapping") <- "OSSL scan_nir wavelength headers"
    return(expected)
  }
  numeric_headers <- suppressWarnings(as.numeric(names(data)))
  indices <- match(config$wavelengths, numeric_headers)
  if (all(!is.na(indices)) && !anyDuplicated(indices)) {
    result <- names(data)[indices]
    attr(result, "mapping") <- "exact numeric wavelength headers"
    return(result)
  }
  missing_ossl <- setdiff(expected, names(data))
  missing_numeric <- config$wavelengths[
    !config$wavelengths %in% numeric_headers
  ]
  stop(
    "The NIR grid must be complete (1350--2550 nm by 2 nm). Neither the ",
    "OSSL scan_nir headers nor exact numeric wavelength headers are ",
    "complete. First missing OSSL header: ", missing_ossl[1L],
    "; first missing numeric wavelength: ", missing_numeric[1L],
    ". No interpolation is performed.", call. = FALSE
  )
}

neospectra_v3_consensus <- function(values) {
  values <- unique(as.character(values[!is.na(values) & nzchar(values)]))
  if (!length(values)) return(NA_character_)
  if (length(values) > 1L) return(NA_character_)
  values
}

neospectra_v3_lab_map <- function(lab, config) {
  sample_mapping <- neospectra_v3_resolve_alias(
    lab, "laboratory sample ID", config$sample_id_column,
    config$sample_id_aliases, "FMCP_NEOSPECTRA_SAMPLE_ID_COLUMN", TRUE
  )
  outcome <- neospectra_v3_resolve_outcome(
    lab, config, "laboratory table"
  )
  sample_id <- trimws(as.character(lab[[sample_mapping$column]]))
  ids <- unique(sample_id[!is.na(sample_id) & nzchar(sample_id)])
  rows <- lapply(ids, function(id) {
    index <- which(sample_id == id & outcome$status == "usable")
    values <- unique(outcome$response[index])
    raw_values <- unique(outcome$raw_soc[index])
    source_status <- unique(outcome$status[sample_id == id])
    data.frame(
      sample_id = id,
      response = if (length(values) == 1L) values else NA_real_,
      raw_soc = if (length(raw_values) == 1L) raw_values else NA_real_,
      lab_rows = sum(sample_id == id),
      lab_status = if (!length(values)) {
        paste(sort(source_status), collapse = "|")
      } else
        if (length(values) > 1L) "conflicting_outcome" else "usable",
      stringsAsFactors = FALSE
    )
  })
  result <- bind_rows_base(rows)
  audit <- bind_rows_base(list(sample_mapping$audit, outcome$audit))
  audit$source_table <- "laboratory"
  attr(result, "schema_audit") <- audit
  attr(result, "outcome_mapping") <- outcome[c(
    "column", "source_scale", "transformation"
  )]
  result
}

neospectra_v3_site_map <- function(site, config) {
  if (is.null(site)) return(NULL)
  sample_mapping <- neospectra_v3_resolve_alias(
    site, "site sample ID", config$sample_id_column,
    config$sample_id_aliases, "FMCP_NEOSPECTRA_SAMPLE_ID_COLUMN", TRUE
  )
  country_mapping <- neospectra_v3_resolve_alias(
    site, "country", config$country_candidates[1L],
    config$country_candidates[-1L], "FMCP_NEOSPECTRA_COUNTRY_COLUMN", FALSE
  )
  audit <- bind_rows_base(list(sample_mapping$audit, country_mapping$audit))
  audit$source_table <- "site"
  if (is.null(country_mapping$column)) {
    result <- data.frame(sample_id = character(), country = character())
    attr(result, "schema_audit") <- audit
    return(result)
  }
  sample_id <- trimws(as.character(site[[sample_mapping$column]]))
  ids <- unique(sample_id[!is.na(sample_id) & nzchar(sample_id)])
  rows <- lapply(ids, function(id) {
    data.frame(
      sample_id = id,
      country = neospectra_v3_consensus(
        site[[country_mapping$column]][sample_id == id]
      ),
      stringsAsFactors = FALSE
    )
  })
  result <- bind_rows_base(rows)
  attr(result, "schema_audit") <- audit
  result
}

neospectra_v3_raw_table <- function(inputs, config) {
  if (identical(inputs$mode, "combined")) {
    raw <- read_neospectra_v3_csv(inputs$combined)
    sample_mapping <- neospectra_v3_resolve_alias(
      raw, "sample ID", config$sample_id_column, config$sample_id_aliases,
      "FMCP_NEOSPECTRA_SAMPLE_ID_COLUMN", TRUE
    )
    scan_mapping <- neospectra_v3_resolve_alias(
      raw, "scan ID", config$scan_id_column, config$scan_id_aliases,
      "FMCP_NEOSPECTRA_SCAN_ID_COLUMN", FALSE
    )
    scanner_mapping <- neospectra_v3_resolve_alias(
      raw, "scanner", config$scanner_column, config$scanner_aliases,
      "FMCP_NEOSPECTRA_SCANNER_COLUMN", FALSE
    )
    country_mapping <- neospectra_v3_resolve_alias(
      raw, "country", config$country_candidates[1L],
      config$country_candidates[-1L], "FMCP_NEOSPECTRA_COUNTRY_COLUMN", FALSE
    )
    outcome <- neospectra_v3_resolve_outcome(
      raw, config, "combined table"
    )
    raw$.fmcp_sample_id <- as.character(raw[[sample_mapping$column]])
    raw$.fmcp_scan_id <- if (is.null(scan_mapping$column)) NA_character_ else
      as.character(raw[[scan_mapping$column]])
    raw$.fmcp_scanner_id <- if (is.null(scanner_mapping$column)) {
      NA_character_
    } else as.character(raw[[scanner_mapping$column]])
    raw$.fmcp_country <- if (is.null(country_mapping$column)) {
      NA_character_
    } else as.character(raw[[country_mapping$column]])
    raw$.fmcp_soc_raw <- outcome$raw_soc
    raw$.fmcp_response <- outcome$response
    raw$.fmcp_outcome_status <- outcome$status
    audit <- bind_rows_base(list(
      sample_mapping$audit, scan_mapping$audit, scanner_mapping$audit,
      country_mapping$audit, outcome$audit
    ))
    audit$source_table <- "combined"
    attr(raw, "schema_audit") <- audit
    attr(raw, "outcome_mapping") <- outcome[c(
      "column", "source_scale", "transformation"
    )]
    attr(raw, "column_mapping") <- list(
      sample_id = sample_mapping$column, scan_id = scan_mapping$column,
      scanner = scanner_mapping$column, country = country_mapping$column
    )
    return(raw)
  }

  nir <- read_neospectra_v3_csv(inputs$nir)
  lab <- read_neospectra_v3_csv(inputs$lab)
  site <- if (is.null(inputs$site)) NULL else
    read_neospectra_v3_csv(inputs$site)
  sample_mapping <- neospectra_v3_resolve_alias(
    nir, "NIR sample ID", config$sample_id_column,
    config$sample_id_aliases, "FMCP_NEOSPECTRA_SAMPLE_ID_COLUMN", TRUE
  )
  scan_mapping <- neospectra_v3_resolve_alias(
    nir, "NIR scan ID", config$scan_id_column, config$scan_id_aliases,
    "FMCP_NEOSPECTRA_SCAN_ID_COLUMN", FALSE
  )
  scanner_mapping <- neospectra_v3_resolve_alias(
    nir, "NIR scanner", config$scanner_column, config$scanner_aliases,
    "FMCP_NEOSPECTRA_SCANNER_COLUMN", FALSE
  )
  lab_map <- neospectra_v3_lab_map(lab, config)
  nir_sample_id <- trimws(as.character(nir[[sample_mapping$column]]))
  match_index <- match(
    nir_sample_id, lab_map$sample_id
  )
  nir$.fmcp_sample_id <- nir_sample_id
  nir$.fmcp_scan_id <- if (is.null(scan_mapping$column)) NA_character_ else
    as.character(nir[[scan_mapping$column]])
  nir$.fmcp_scanner_id <- if (is.null(scanner_mapping$column)) {
    NA_character_
  } else as.character(nir[[scanner_mapping$column]])
  nir$.fmcp_response <- lab_map$response[match_index]
  nir$.fmcp_soc_raw <- lab_map$raw_soc[match_index]
  nir$.fmcp_outcome_status <- lab_map$lab_status[match_index]
  nir$.fmcp_outcome_status[is.na(nir$.fmcp_outcome_status)] <-
    "sample_missing_from_laboratory_table"
  site_map <- neospectra_v3_site_map(site, config)
  nir$.fmcp_country <- if (is.null(site_map) || !nrow(site_map)) {
    NA_character_
  } else {
    site_map$country[match(
      nir_sample_id, site_map$sample_id
    )]
  }
  nir_audit <- bind_rows_base(list(
    sample_mapping$audit, scan_mapping$audit, scanner_mapping$audit
  ))
  nir_audit$source_table <- "NIR"
  site_audit <- if (is.null(site_map)) data.frame() else
    attr(site_map, "schema_audit")
  attr(nir, "schema_audit") <- bind_rows_base(list(
    nir_audit, attr(lab_map, "schema_audit"), site_audit
  ))
  attr(nir, "outcome_mapping") <- attr(lab_map, "outcome_mapping")
  attr(nir, "column_mapping") <- list(
    sample_id = sample_mapping$column, scan_id = scan_mapping$column,
    scanner = scanner_mapping$column,
    country = if (is.null(site_map) || !nrow(site_map)) NULL else {
      site_audit$selected_column[site_audit$role == "country"][1L]
    }
  )
  nir
}

prepare_neospectra_v3 <- function(inputs, config) {
  raw <- neospectra_v3_raw_table(inputs, config)
  column_mapping <- attr(raw, "column_mapping")
  outcome_mapping <- attr(raw, "outcome_mapping")
  schema_audit <- attr(raw, "schema_audit")
  spectral_columns <- neospectra_v3_spectral_columns(raw, config)
  spectral_mapping <- attr(spectral_columns, "mapping")

  sample_id <- trimws(as.character(raw$.fmcp_sample_id))
  scan_id <- trimws(as.character(raw$.fmcp_scan_id))
  scanner_id <- trimws(as.character(raw$.fmcp_scanner_id))
  raw_soc <- suppressWarnings(as.numeric(raw$.fmcp_soc_raw))
  response <- suppressWarnings(as.numeric(raw$.fmcp_response))
  outcome_source_status <- as.character(raw$.fmcp_outcome_status)
  country <- trimws(as.character(raw$.fmcp_country))
  x <- as.matrix(raw[, spectral_columns, drop = FALSE])
  suppressWarnings(storage.mode(x) <- "double")

  reason <- rep("included", nrow(raw))
  set_reason <- function(condition, value) {
    replace <- condition & reason == "included"
    reason[replace] <<- value
  }
  set_reason(is.na(sample_id) | !nzchar(sample_id), "missing_sample_id")
  set_reason(rowSums(!is.finite(x)) > 0L, "missing_or_nonfinite_reflectance")

  valid_ids <- unique(sample_id[!is.na(sample_id) & nzchar(sample_id)])
  outcome_rows <- lapply(valid_ids, function(id) {
    index <- which(sample_id == id)
    usable_index <- index[outcome_source_status[index] == "usable"]
    response_values <- unique(response[usable_index][
      is.finite(response[usable_index])
    ])
    raw_values <- unique(raw_soc[usable_index][is.finite(raw_soc[usable_index])])
    status <- if (!length(response_values)) {
      "missing_or_nonfinite_soc_consensus"
    } else if (length(response_values) > 1L || length(raw_values) > 1L) {
      "conflicting_soc_within_sample"
    } else "usable"
    data.frame(
      sample_id = id,
      raw_soc = if (status == "usable") raw_values[1L] else NA_real_,
      response = if (status == "usable") response_values[1L] else NA_real_,
      outcome_status = status,
      finite_outcome_rows = length(usable_index),
      missing_outcome_rows = length(index) - length(usable_index),
      stringsAsFactors = FALSE
    )
  })
  outcome_by_sample <- bind_rows_base(outcome_rows)
  outcome_match <- match(sample_id, outcome_by_sample$sample_id)
  sample_outcome_status <- outcome_by_sample$outcome_status[outcome_match]
  set_reason(
    sample_outcome_status == "missing_or_nonfinite_soc_consensus",
    "sample_missing_or_nonfinite_soc_consensus"
  )
  set_reason(
    sample_outcome_status == "conflicting_soc_within_sample",
    "conflicting_soc_within_sample"
  )

  row_manifest <- data.frame(
    raw_row = seq_len(nrow(raw)), sample_id = sample_id,
    scan_id = scan_id, scanner_id = scanner_id,
    soc_raw = raw_soc, response_log1p_soc = response,
    outcome_source_status = outcome_source_status,
    included = reason == "included", exclusion_reason = reason,
    nonpositive_reflectance = rowSums(x <= 0, na.rm = TRUE),
    above_one_reflectance = rowSums(x > 1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  keep <- reason == "included"
  if (!any(keep)) stop("No usable NeoSpectra rows remain.", call. = FALSE)
  kept_sample_id <- sample_id[keep]
  kept_scan_id <- scan_id[keep]
  kept_scanner_id <- scanner_id[keep]
  kept_country <- country[keep]
  kept_x <- x[keep, , drop = FALSE]
  sample_ids <- unique(kept_sample_id)
  sample_factor_raw <- factor(kept_sample_id, levels = sample_ids)
  raw_rows_per_sample <- as.numeric(table(sample_factor_raw))
  scanner_available <- !is.null(column_mapping$scanner) &&
    all(!is.na(kept_scanner_id) & nzchar(kept_scanner_id))

  if (scanner_available) {
    scan_key <- paste(kept_sample_id, kept_scanner_id, sep = "\r")
    scan_factor <- factor(scan_key, levels = unique(scan_key))
    scan_counts <- as.numeric(table(scan_factor))
    scan_x <- rowsum(kept_x, scan_factor, reorder = FALSE) / scan_counts
    first <- match(levels(scan_factor), scan_key)
    scan_sample_id <- kept_sample_id[first]
    scan_scanner_id <- kept_scanner_id[first]
    scan_country <- vapply(
      split(kept_country, scan_factor), neospectra_v3_consensus,
      FUN.VALUE = character(1L)
    )
    scan_ids <- vapply(split(kept_scan_id, scan_factor), function(values) {
      values <- unique(values[!is.na(values) & nzchar(values)])
      if (length(values)) paste(values, collapse = "|") else NA_character_
    }, FUN.VALUE = character(1L))
    scan_sample_match <- match(scan_sample_id, outcome_by_sample$sample_id)
    scan_y <- outcome_by_sample$response[scan_sample_match]
    sample_factor <- factor(scan_sample_id, levels = sample_ids)
    scanner_counts <- as.numeric(table(sample_factor))
    sample_x <- rowsum(scan_x, sample_factor, reorder = FALSE) /
      scanner_counts
    sample_country <- vapply(
      split(scan_country, sample_factor), neospectra_v3_consensus,
      FUN.VALUE = character(1L)
    )
    scan_dataset <- list(
      name = "NeoSpectra_SOC_scan_level",
      subject_id = scan_sample_id, sample_id = scan_sample_id,
      scan_id = unname(scan_ids), scanner_id = scan_scanner_id,
      replicate_rows = scan_counts, country = unname(scan_country),
      y = unname(scan_y), x = unclass(scan_x)
    )
    aggregation_operation <- paste(
      "Arithmetic mean within sample x scanner, then equal-weight mean",
      "across scanners; no outlier deletion"
    )
    scanner_reason <- "Complete scanner metadata in included rows."
  } else {
    scanner_counts <- rep(NA_integer_, length(sample_ids))
    sample_x <- rowsum(kept_x, sample_factor_raw, reorder = FALSE) /
      raw_rows_per_sample
    sample_country <- vapply(
      split(kept_country, sample_factor_raw), neospectra_v3_consensus,
      FUN.VALUE = character(1L)
    )
    scan_dataset <- NULL
    aggregation_operation <- paste(
      "Equal-weight arithmetic mean of all available anonymous spectra",
      "within soil sample; scanner leave-out unavailable"
    )
    scanner_reason <- if (is.null(column_mapping$scanner)) {
      "No scanner column in the source; anonymous rows averaged within sample."
    } else {
      "Scanner column is incomplete; anonymous rows averaged within sample."
    }
  }
  sample_match <- match(sample_ids, outcome_by_sample$sample_id)
  sample_y <- outcome_by_sample$response[sample_match]
  sample_raw_soc <- outcome_by_sample$raw_soc[sample_match]

  argvals_physical <- as.numeric(config$wavelengths)
  argvals_u <- (argvals_physical - min(argvals_physical)) /
    diff(range(argvals_physical))
  sample_dataset <- list(
    name = "NeoSpectra_SOC",
    subject_id = sample_ids,
    sample_id = sample_ids,
    scanner_count = scanner_counts,
    country = unname(sample_country),
    y = unname(sample_y), x = unclass(sample_x),
    argvals_u = argvals_u, argvals_physical = argvals_physical,
    x_axis_label = "NIR wavelength (nm)",
    y_axis_label = "log1p soil organic carbon",
    response_description = paste0(
      "log1p(SOC) from source column ", outcome_mapping$column,
      " (source scale: ", outcome_mapping$source_scale, ")"
    ),
    predictor_description = paste(
      "NeoSpectra NIR reflectance, 1350--2550 nm by 2 nm;",
      if (scanner_available) "scanner means averaged once per soil sample" else
        "anonymous replicate spectra averaged once per soil sample"
    )
  )
  if (!is.null(scan_dataset)) {
    scan_dataset$argvals_u <- argvals_u
    scan_dataset$argvals_physical <- argvals_physical
  }

  all_ids <- unique(row_manifest$sample_id[
    !is.na(row_manifest$sample_id) & nzchar(row_manifest$sample_id)
  ])
  sample_manifest <- bind_rows_base(lapply(all_ids, function(id) {
    raw_rows <- row_manifest[row_manifest$sample_id == id, , drop = FALSE]
    included <- id %in% sample_ids
    data.frame(
      sample_id = id, included = included,
      raw_rows = nrow(raw_rows), valid_raw_rows = sum(raw_rows$included),
      finite_outcome_rows = sum(raw_rows$outcome_source_status == "usable"),
      missing_outcome_rows = sum(raw_rows$outcome_source_status != "usable"),
      raw_soc = if (included) sample_raw_soc[match(id, sample_ids)] else NA_real_,
      response_log1p_soc = if (included) {
        sample_y[match(id, sample_ids)]
      } else NA_real_,
      scanner_count = if (included) scanner_counts[match(id, sample_ids)] else
        NA_integer_,
      exclusion_reason = if (included) "included" else
        paste(sort(unique(raw_rows$exclusion_reason)), collapse = "|"),
      stringsAsFactors = FALSE
    )
  }))

  source_paths <- unlist(
    Filter(Negate(is.null), inputs[c("combined", "nir", "lab", "site")]),
    use.names = TRUE
  )
  source_manifest <- data.frame(
    role = names(source_paths), path = normalizePath(source_paths),
    bytes = as.numeric(file.info(source_paths)$size),
    md5 = unname(tools::md5sum(source_paths)),
    stringsAsFactors = FALSE
  )
  quality <- data.frame(
    item = c(
      "raw_rows", "included_raw_rows", "excluded_raw_rows",
      "unique_included_samples", "scan_level_rows", "wavelengths",
      "samples_with_country", "samples_without_country",
      "scanner_metadata_available", "sample_aggregation"
    ),
    value = as.character(c(
      nrow(raw), sum(row_manifest$included), sum(!row_manifest$included),
      length(sample_ids), if (is.null(scan_dataset)) 0L else nrow(scan_dataset$x),
      ncol(sample_x),
      sum(!is.na(sample_country) & nzchar(sample_country)),
      sum(is.na(sample_country) | !nzchar(sample_country)),
      scanner_available, aggregation_operation
    )),
    stringsAsFactors = FALSE
  )
  spectral_audit <- data.frame(
    role = "reflectance", required = TRUE,
    selected_column = as.character(spectral_columns),
    detection = spectral_mapping,
    recognized_candidates = paste0("wavelength_nm=", config$wavelengths),
    source_table = if (inputs$mode == "combined") "combined" else "NIR",
    canonical_column = paste0("wavelength_", config$wavelengths),
    stringsAsFactors = FALSE
  )
  schema_audit <- bind_rows_base(list(schema_audit, spectral_audit))
  country_available <- !is.null(column_mapping$country) &&
    any(!is.na(sample_country) & nzchar(sample_country))
  availability <- data.frame(
    analysis = c("scanner_stability", "geographic_external_metadata"),
    available = c(scanner_available, country_available),
    reason = c(
      scanner_reason,
      if (country_available) "Country metadata are present." else
        "No usable country column in the source."
    ), stringsAsFactors = FALSE
  )
  list(
    version = config$version, mode = config$mode,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    sample_dataset = validate_neospectra_v3_dataset(sample_dataset, TRUE),
    scan_dataset = if (is.null(scan_dataset)) NULL else
      validate_neospectra_v3_dataset(scan_dataset, FALSE),
    row_manifest = row_manifest, sample_manifest = sample_manifest,
    source_manifest = source_manifest, quality = quality,
    availability = availability, schema_audit = schema_audit,
    schema = list(
      sample_id = column_mapping$sample_id, scan_id = column_mapping$scan_id,
      scanner = column_mapping$scanner, outcome = outcome_mapping$column,
      outcome_source_scale = outcome_mapping$source_scale,
      outcome_transformation = outcome_mapping$transformation,
      country = column_mapping$country,
      country_candidates = config$country_candidates,
      spectral_columns = as.character(spectral_columns),
      spectral_mapping = spectral_mapping
    ),
    transformations = data.frame(
      stage = c("spectra", "sample aggregation", "outcome", "exclusions"),
      operation = c(
        config$predictor_transform,
        aggregation_operation,
        paste0(outcome_mapping$transformation, "; raw SOC retained in manifests"),
        paste(
          "No imputation/outlier deletion; exclude missing sample IDs,",
          "nonfinite spectra, or samples without a unique finite SOC consensus"
        )
      ),
      stringsAsFactors = FALSE
    )
  )
}

validate_neospectra_v3_dataset <- function(dataset, unique_samples = TRUE) {
  required <- c("subject_id", "sample_id", "y", "x", "argvals_u",
                "argvals_physical")
  missing <- setdiff(required, names(dataset))
  if (length(missing)) stop("NeoSpectra dataset missing: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  if (!is.matrix(dataset$x) || nrow(dataset$x) != length(dataset$y) ||
      ncol(dataset$x) != length(dataset$argvals_u) ||
      length(dataset$sample_id) != length(dataset$y)) {
    stop("NeoSpectra dataset dimensions are inconsistent.", call. = FALSE)
  }
  if (any(!is.finite(dataset$x)) || any(!is.finite(dataset$y))) {
    stop("NeoSpectra processed data contain non-finite values.", call. = FALSE)
  }
  if (anyDuplicated(dataset$sample_id) && unique_samples) {
    stop("Primary NeoSpectra data must have one row per sample ID.",
         call. = FALSE)
  }
  if (any(diff(dataset$argvals_u) <= 0) ||
      any(diff(dataset$argvals_physical) <= 0)) {
    stop("NeoSpectra grids must be strictly increasing.", call. = FALSE)
  }
  dataset
}

load_neospectra_v3_processed <- function(config) {
  path <- neospectra_v3_processed_path(config)
  if (!file.exists(path)) {
    stop("Processed NeoSpectra v3 data are missing. Run ",
         "10_prepare_neospectra_v3.R first.", call. = FALSE)
  }
  object <- readRDS(path)
  object$sample_dataset <- validate_neospectra_v3_dataset(
    object$sample_dataset, TRUE
  )
  if (!is.null(object$scan_dataset)) {
    object$scan_dataset <- validate_neospectra_v3_dataset(
      object$scan_dataset, FALSE
    )
  }
  object
}

make_grouped_folds_v3 <- function(groups, folds, seed) {
  groups <- as.character(groups)
  unique_groups <- unique(groups)
  if (folds > length(unique_groups)) {
    stop("Number of grouped folds exceeds unique sample count.", call. = FALSE)
  }
  set.seed(as.integer(seed))
  shuffled <- sample(unique_groups)
  assignment <- rep(seq_len(folds), length.out = length(shuffled))
  labels <- assignment[match(groups, shuffled)]
  split_check <- tapply(labels, groups, function(x) length(unique(x)))
  if (any(split_check != 1L)) stop("Grouped fold construction failed.",
                                  call. = FALSE)
  as.integer(labels)
}

sparse_candidate_grid_v3 <- function(config, method = NULL) {
  grid <- expand.grid(
    M = as.integer(config$M), lambda = as.numeric(config$lambda),
    kappa = as.numeric(config$kappa), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid <- grid[order(grid$M, grid$lambda, grid$kappa), , drop = FALSE]
  grid$candidate_id <- sprintf(
    "M%d_L%.12g_K%.12g", grid$M, grid$lambda, grid$kappa
  )
  if (!is.null(method)) grid$method <- method
  rownames(grid) <- NULL
  grid
}

candidate_config_v3 <- function(config, candidate) {
  current <- config
  current$M <- as.integer(candidate$M[1L])
  current$lambda <- as.numeric(candidate$lambda[1L])
  current$kappa <- as.numeric(candidate$kappa[1L])
  current
}

fit_neospectra_candidate_v3 <- function(dataset, train_index, test_index,
                                        config, method, candidate, seed) {
  set.seed(as.integer(seed))
  current <- candidate_config_v3(config, candidate)
  fit_set <- fit_real_method_set(
    dataset, train_index, test_index, current, methods = method
  )
  fit <- fit_set$fits[[method]]
  list(fit_set = fit_set, fit = fit, prediction = fit$prediction)
}

inner_tune_sparse_v3 <- function(dataset, train_index, config, method,
                                 seed, task_id, evaluation = "outer_cv") {
  candidate_grid <- sparse_candidate_grid_v3(config, method)
  inner_labels <- make_grouped_folds_v3(
    dataset$sample_id[train_index], config$inner_folds, seed
  )
  tasks <- lapply(seq_len(nrow(candidate_grid)), function(candidate_index) {
    list(index = candidate_index, candidate = candidate_grid[candidate_index, ])
  })
  worker <- function(task) {
    candidate <- task$candidate
    rows <- vector("list", config$inner_folds)
    for (fold in seq_len(config$inner_folds)) {
      validation_index <- train_index[inner_labels == fold]
      fitting_index <- train_index[inner_labels != fold]
      captured <- tryCatch({
        fitted <- fit_neospectra_candidate_v3(
          dataset, fitting_index, validation_index, config, method,
          candidate, seed + task$index * 100L + fold
        )
        metrics <- metric_vector(
          dataset$y[validation_index], fitted$prediction
        )
        list(
          success = fitted$fit$success && all(is.finite(fitted$prediction)),
          metrics = metrics, elapsed = fitted$fit$elapsed,
          error = fitted$fit$error, warnings = fitted$fit$warnings
        )
      }, error = function(condition) list(
        success = FALSE, metrics = c(RMSE = NA, MAE = NA, R2 = NA),
        elapsed = NA_real_, error = conditionMessage(condition), warnings = ""
      ))
      rows[[fold]] <- data.frame(
        evaluation = evaluation, task_id = task_id, method = method,
        candidate_id = candidate$candidate_id, M = candidate$M,
        lambda = candidate$lambda, kappa = candidate$kappa,
        inner_fold = fold, n_validation = length(validation_index),
        rmse = unname(captured$metrics["RMSE"]),
        mae = unname(captured$metrics["MAE"]),
        r2 = unname(captured$metrics["R2"]),
        elapsed_seconds = captured$elapsed, success = captured$success,
        error = captured$error, warnings = captured$warnings,
        stringsAsFactors = FALSE
      )
    }
    bind_rows_base(rows)
  }
  results <- run_task_batch(tasks, worker, config$workers)
  tuning <- bind_rows_base(results)
  split_rows <- split(tuning, tuning$candidate_id)
  summary_rows <- lapply(split_rows, function(rows) {
    all_success <- nrow(rows) == config$inner_folds && all(rows$success) &&
      all(is.finite(rows$rmse))
    data.frame(
      candidate_id = rows$candidate_id[1L], M = rows$M[1L],
      lambda = rows$lambda[1L], kappa = rows$kappa[1L],
      success = all_success,
      pooled_rmse = if (all_success) sqrt(weighted.mean(
        rows$rmse^2, rows$n_validation
      )) else Inf,
      total_elapsed_seconds = sum(rows$elapsed_seconds, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary <- bind_rows_base(summary_rows)
  valid <- which(summary$success & is.finite(summary$pooled_rmse))
  if (!length(valid)) {
    stop("Every inner-CV candidate failed for ", method, " in ", task_id,
         call. = FALSE)
  }
  best <- valid[order(
    summary$pooled_rmse[valid], summary$M[valid],
    summary$lambda[valid], summary$kappa[valid]
  )[1L]]
  selected <- summary[best, , drop = FALSE]
  tuning$selected <- tuning$candidate_id == selected$candidate_id
  inner_folds <- data.frame(
    evaluation = evaluation, task_id = task_id, method = method,
    sample_id = dataset$sample_id[train_index], inner_fold = inner_labels,
    seed = as.integer(seed), stringsAsFactors = FALSE
  )
  list(
    selected = selected, tuning = tuning, inner_labels = inner_labels,
    inner_folds = inner_folds
  )
}

outer_task_v3 <- function(dataset, repeat_index, fold, outer_labels,
                          config) {
  test_index <- which(outer_labels == fold)
  train_index <- which(outer_labels != fold)
  task_id <- paste0("outer_r", repeat_index, "_f", fold)
  prediction_rows <- list()
  tuning_rows <- list()
  inner_fold_rows <- list()
  selection_rows_list <- list()
  method_index <- 1L
  for (method in config$primary_sparse_methods) {
    tuning <- inner_tune_sparse_v3(
      dataset, train_index, config, method,
      neospectra_v3_seed(config, "inner", repeat_index * 10L + fold),
      task_id
    )
    candidate <- tuning$selected
    fitted <- fit_neospectra_candidate_v3(
      dataset, train_index, test_index, config, method, candidate,
      neospectra_v3_seed(config, "outer", repeat_index * 10L + fold,
                         method_index)
    )
    prediction_rows[[method]] <- data.frame(
      evaluation = "outer_cv", task_id = task_id,
      repeat_index = repeat_index, fold = fold,
      sample_id = dataset$sample_id[test_index],
      country = dataset$country[test_index], observed = dataset$y[test_index],
      predicted = fitted$prediction, method = method,
      success = fitted$fit$success, error = fitted$fit$error,
      warnings = fitted$fit$warnings, stringsAsFactors = FALSE
    )
    current_selection <- fitted$fit$selected
    current_selection$evaluation <- "outer_cv"
    current_selection$task_id <- task_id
    current_selection$tuning_elapsed_seconds <- sum(
      tuning$tuning$elapsed_seconds, na.rm = TRUE
    )
    current_selection$final_elapsed_seconds <- fitted$fit$elapsed
    current_selection$selected_candidate_id <- candidate$candidate_id
    selection_rows_list[[method]] <- current_selection
    tuning_rows[[method]] <- tuning$tuning
    inner_fold_rows[[method]] <- tuning$inner_folds
    method_index <- method_index + 1L
  }

  baseline_methods <- setdiff(config$methods, config$primary_sparse_methods)
  if (length(baseline_methods)) {
    set.seed(neospectra_v3_seed(
      config, "outer", repeat_index * 10L + fold, 90L
    ))
    baseline <- fit_real_method_set(
      dataset, train_index, test_index, config, methods = baseline_methods
    )
    for (method in baseline_methods) {
      fit <- baseline$fits[[method]]
      prediction_rows[[method]] <- data.frame(
        evaluation = "outer_cv", task_id = task_id,
        repeat_index = repeat_index, fold = fold,
        sample_id = dataset$sample_id[test_index],
        country = dataset$country[test_index], observed = dataset$y[test_index],
        predicted = fit$prediction, method = method,
        success = fit$success, error = fit$error, warnings = fit$warnings,
        stringsAsFactors = FALSE
      )
      current_selection <- fit$selected
      current_selection$evaluation <- "outer_cv"
      current_selection$task_id <- task_id
      current_selection$tuning_elapsed_seconds <- NA_real_
      current_selection$final_elapsed_seconds <- fit$elapsed
      current_selection$selected_candidate_id <- NA_character_
      selection_rows_list[[method]] <- current_selection
    }
  }
  list(
    key = task_id, predictions = bind_rows_base(prediction_rows),
    tuning = bind_rows_base(tuning_rows),
    inner_folds = bind_rows_base(inner_fold_rows),
    selection = bind_rows_base(selection_rows_list)
  )
}

neospectra_v3_country_class <- function(country) {
  normalized <- toupper(gsub("[^A-Z]", "", as.character(country)))
  result <- rep("unknown", length(normalized))
  result[normalized %in% c("US", "USA", "UNITEDSTATES",
                           "UNITEDSTATESOFAMERICA")] <- "US"
  result[normalized %in% c("GH", "GHA", "GHANA", "KE", "KEN", "KENYA",
                           "NG", "NGA", "NIGERIA")] <- "Africa"
  result
}

run_geographic_external_v3 <- function(dataset, config) {
  country_class <- neospectra_v3_country_class(dataset$country)
  train_index <- which(country_class == "US")
  test_index <- which(country_class == "Africa")
  if (length(train_index) < max(30L, config$inner_folds) ||
      length(test_index) < 10L) {
    return(list(
      available = FALSE,
      reason = paste(
        "External split requires >=30 US and >=10 Ghana/Kenya/Nigeria",
        "samples with country metadata."
      ),
      predictions = data.frame(), tuning = data.frame(),
      inner_folds = data.frame(), selection = data.frame()
    ))
  }
  prediction_rows <- list()
  tuning_rows <- list()
  inner_fold_rows <- list()
  selection_list <- list()
  for (method_index in seq_along(config$primary_sparse_methods)) {
    method <- config$primary_sparse_methods[method_index]
    tuning <- inner_tune_sparse_v3(
      dataset, train_index, config, method,
      neospectra_v3_seed(config, "external"),
      "us_to_africa", evaluation = "geographic_external"
    )
    fitted <- fit_neospectra_candidate_v3(
      dataset, train_index, test_index, config, method, tuning$selected,
      neospectra_v3_seed(config, "external", method_index, 50L)
    )
    prediction_rows[[method]] <- data.frame(
      evaluation = "geographic_external", task_id = "us_to_africa",
      repeat_index = NA_integer_, fold = NA_integer_,
      sample_id = dataset$sample_id[test_index],
      country = dataset$country[test_index], observed = dataset$y[test_index],
      predicted = fitted$prediction, method = method,
      success = fitted$fit$success, error = fitted$fit$error,
      warnings = fitted$fit$warnings, stringsAsFactors = FALSE
    )
    current_selection <- fitted$fit$selected
    current_selection$evaluation <- "geographic_external"
    current_selection$task_id <- "us_to_africa"
    current_selection$tuning_elapsed_seconds <- sum(
      tuning$tuning$elapsed_seconds, na.rm = TRUE
    )
    current_selection$final_elapsed_seconds <- fitted$fit$elapsed
    current_selection$selected_candidate_id <- tuning$selected$candidate_id
    selection_list[[method]] <- current_selection
    tuning_rows[[method]] <- tuning$tuning
    inner_fold_rows[[method]] <- tuning$inner_folds
  }
  baseline_methods <- setdiff(config$methods, config$primary_sparse_methods)
  if (length(baseline_methods)) {
    baseline <- fit_real_method_set(
      dataset, train_index, test_index, config, methods = baseline_methods
    )
    for (method in baseline_methods) {
      fit <- baseline$fits[[method]]
      prediction_rows[[method]] <- data.frame(
        evaluation = "geographic_external", task_id = "us_to_africa",
        repeat_index = NA_integer_, fold = NA_integer_,
        sample_id = dataset$sample_id[test_index],
        country = dataset$country[test_index], observed = dataset$y[test_index],
        predicted = fit$prediction, method = method,
        success = fit$success, error = fit$error, warnings = fit$warnings,
        stringsAsFactors = FALSE
      )
      current_selection <- fit$selected
      current_selection$evaluation <- "geographic_external"
      current_selection$task_id <- "us_to_africa"
      current_selection$tuning_elapsed_seconds <- NA_real_
      current_selection$final_elapsed_seconds <- fit$elapsed
      current_selection$selected_candidate_id <- NA_character_
      selection_list[[method]] <- current_selection
    }
  }
  list(
    available = TRUE, reason = "",
    predictions = bind_rows_base(prediction_rows),
    tuning = bind_rows_base(tuning_rows),
    inner_folds = bind_rows_base(inner_fold_rows),
    selection = bind_rows_base(selection_list)
  )
}

ccc_v3 <- function(observed, predicted) {
  complete <- is.finite(observed) & is.finite(predicted)
  observed <- observed[complete]
  predicted <- predicted[complete]
  if (length(observed) < 2L) return(NA_real_)
  covariance <- stats::cov(observed, predicted)
  denominator <- stats::var(observed) + stats::var(predicted) +
    (mean(observed) - mean(predicted))^2
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  2 * covariance / denominator
}

prediction_metrics_v3 <- function(observed, predicted) {
  metrics <- metric_vector(observed, predicted)
  c(metrics, CCC = ccc_v3(observed, predicted), n = length(observed))
}

neospectra_truth_v3 <- function(argvals, scenario) {
  u <- as.numeric(argvals)
  bump <- function(left, right, amplitude = 1) {
    value <- numeric(length(u))
    inside <- u >= left & u <= right
    value[inside] <- amplitude * sin(pi * (u[inside] - left) /
                                      (right - left))^2
    value
  }
  if (identical(scenario, "interior_null")) {
    return(bump(0.02, 0.30, 1.0) - bump(0.70, 0.98, 0.9))
  }
  if (identical(scenario, "moderate_shoulder")) {
    return(bump(0.05, 0.42, 1.0) + bump(0.55, 0.76, 0.32))
  }
  if (identical(scenario, "no_null")) {
    return(0.35 + 0.8 * sin(2 * pi * u) + 0.35 * cos(5 * pi * u))
  }
  if (identical(scenario, "all_null")) return(numeric(length(u)))
  stop("Unknown NeoSpectra simulation scenario: ", scenario, call. = FALSE)
}

simulation_tasks_v3 <- function(dataset, config) {
  if (nrow(dataset$x) < 40L) {
    stop("Real-design simulation requires at least 40 unique samples.",
         call. = FALSE)
  }
  maximum_training <- max(20L, nrow(dataset$x) - 20L)
  n_values <- config$simulation_n[config$simulation_n <= maximum_training]
  if (!length(n_values)) n_values <- max(20L, floor(0.7 * nrow(dataset$x)))
  grid <- expand.grid(
    scenario = config$simulation_scenarios,
    n_train = as.integer(n_values), snr = config$simulation_snr,
    replication = seq_len(config$simulation_replicates),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  grid <- grid[!(grid$scenario == "all_null" &
                   grid$snr != config$simulation_snr[1L]), , drop = FALSE]
  grid$snr_label <- ifelse(
    grid$scenario == "all_null", "allnull",
    format(grid$snr, trim = TRUE)
  )
  grid$snr[grid$scenario == "all_null"] <- NA_real_
  grid$task_index <- seq_len(nrow(grid))
  grid$task_id <- sprintf(
    "%s_n%d_snr%s_r%d", grid$scenario, grid$n_train,
    grid$snr_label, grid$replication
  )
  grid$snr_label <- NULL
  grid$snr_definition <- config$simulation_snr_definition
  grid
}

truth_null_mask_v3 <- function(beta_truth, scenario, tolerance = 1e-12) {
  if (identical(scenario, "no_null")) {
    return(rep(FALSE, length(beta_truth)))
  }
  if (identical(scenario, "all_null")) {
    return(rep(TRUE, length(beta_truth)))
  }
  abs(beta_truth) <= tolerance
}

support_metrics_v3 <- function(beta_estimate, beta_truth, weights,
                               truth_null = abs(beta_truth) <= 1e-12) {
  if (length(truth_null) != length(beta_truth) || anyNA(truth_null)) {
    stop("Truth-null mask must be a complete coefficient-grid mask.",
         call. = FALSE)
  }
  estimate_null <- abs(beta_estimate) <= sqrt(.Machine$double.eps)
  c(
    ISE = sum(weights * (beta_estimate - beta_truth)^2),
    true_null_recovery = if (any(truth_null)) mean(estimate_null[truth_null]) else NA_real_,
    false_null_rate = if (any(!truth_null)) mean(estimate_null[!truth_null]) else NA_real_,
    estimated_null_fraction = mean(estimate_null),
    true_null_fraction = mean(truth_null)
  )
}

run_real_design_simulation_task_v3 <- function(task, dataset, config) {
  set.seed(neospectra_v3_seed(config, "simulation", task$task_index))
  n_available <- nrow(dataset$x)
  n_train <- min(as.integer(task$n_train), n_available - 20L)
  n_test <- min(max(50L, ceiling(n_train / 4)), n_available - n_train)
  selected <- sample(seq_len(n_available), n_train + n_test, replace = FALSE)
  train_index <- seq_len(n_train)
  test_index <- n_train + seq_len(n_test)
  x <- dataset$x[selected, , drop = FALSE]
  x_mean <- colMeans(x[train_index, , drop = FALSE])
  x_centered <- sweep(x, 2L, x_mean, FUN = "-")
  weights <- trapezoid_weights(dataset$argvals_u)
  beta_base <- neospectra_truth_v3(dataset$argvals_u, task$scenario)
  raw_signal <- as.numeric(x_centered %*% (weights * beta_base))
  if (identical(task$scenario, "all_null")) {
    beta_truth <- beta_base
    signal <- raw_signal
    noise_sd <- 1
    target_variance_snr <- NA_real_
  } else {
    scale <- stats::sd(raw_signal[train_index])
    if (!is.finite(scale) || scale <= sqrt(.Machine$double.eps)) {
      stop("Real-design signal scale is numerically zero.", call. = FALSE)
    }
    beta_truth <- beta_base / scale
    signal <- raw_signal / scale
    target_variance_snr <- as.numeric(task$snr)
    noise_sd <- 1 / sqrt(target_variance_snr)
  }
  truth_null <- truth_null_mask_v3(beta_truth, task$scenario)
  error <- stats::rnorm(length(signal), sd = noise_sd)
  y <- signal + error
  training_signal_variance <- stats::var(signal[train_index])
  training_noise_variance <- stats::var(error[train_index])
  empirical_variance_snr <- if (
    is.finite(training_noise_variance) && training_noise_variance > 0
  ) training_signal_variance / training_noise_variance else NA_real_
  simulated <- dataset
  simulated$name <- paste0("NeoSpectra_real_design_", task$scenario)
  simulated$sample_id <- dataset$sample_id[selected]
  simulated$subject_id <- simulated$sample_id
  simulated$country <- dataset$country[selected]
  simulated$x <- x
  simulated$y <- y
  fitted <- fit_real_method_set(
    simulated, train_index, test_index, config,
    methods = config$primary_sparse_methods
  )
  metric_rows <- list()
  curve_rows <- list()
  for (method in config$primary_sparse_methods) {
    fit <- fitted$fits[[method]]
    beta_u <- fit$beta_scaled * fitted$scaler$y_scale / fitted$scaler$x_scale
    support <- support_metrics_v3(
      beta_u, beta_truth, weights, truth_null = truth_null
    )
    metric_rows[[method]] <- data.frame(
      task_id = task$task_id, task_index = task$task_index,
      scenario = task$scenario, n_train = n_train, n_test = n_test,
      snr = if (task$scenario == "all_null") NA_real_ else task$snr,
      snr_definition = config$simulation_snr_definition,
      target_variance_snr = target_variance_snr,
      noise_sd = noise_sd,
      training_signal_variance = training_signal_variance,
      training_noise_variance = training_noise_variance,
      empirical_variance_snr = empirical_variance_snr,
      replication = task$replication, method = method,
      ISE = unname(support["ISE"]),
      true_null_recovery = unname(support["true_null_recovery"]),
      false_null_rate = unname(support["false_null_rate"]),
      estimated_null_fraction = unname(support["estimated_null_fraction"]),
      true_null_fraction = unname(support["true_null_fraction"]),
      signal_pmse = mean((fit$prediction - signal[test_index])^2),
      noisy_rmse = rmse(y[test_index], fit$prediction),
      elapsed_seconds = fit$elapsed, success = fit$success,
      M = fit$selected$M[1L], lambda = fit$selected$lambda[1L],
      kappa = fit$selected$kappa[1L], error = fit$error,
      warnings = fit$warnings, stringsAsFactors = FALSE
    )
    if (task$replication == 1L) {
      curve_rows[[method]] <- data.frame(
        task_id = task$task_id, scenario = task$scenario,
        n_train = n_train, snr = task$snr, method = method,
        wavelength = dataset$argvals_physical,
        beta_truth = beta_truth, beta_estimate = beta_u,
        coefficient_scale = "normalized wavelength domain [0,1]",
        true_null = truth_null,
        estimated_null = abs(beta_u) <= config$zero_tolerance,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    key = task$task_id, metrics = bind_rows_base(metric_rows),
    curves = bind_rows_base(curve_rows),
    split = data.frame(
      task_id = task$task_id, role = c(rep("train", n_train), rep("test", n_test)),
      sample_id = simulated$sample_id, stringsAsFactors = FALSE
    )
  )
}

select_global_sparse_candidates_v3 <- function(tuning, config) {
  tuning <- tuning[
    tuning$evaluation == "outer_cv" &
      tuning$method %in% config$primary_sparse_methods,
    , drop = FALSE
  ]
  rows <- list()
  for (method in config$primary_sparse_methods) {
    current <- tuning[tuning$method == method, , drop = FALSE]
    candidates <- split(current, current$candidate_id)
    summary <- bind_rows_base(lapply(candidates, function(value) {
      complete <- all(value$success) && all(is.finite(value$rmse))
      data.frame(
        method = method, candidate_id = value$candidate_id[1L],
        M = value$M[1L], lambda = value$lambda[1L],
        kappa = value$kappa[1L],
        pooled_rmse = if (complete) sqrt(weighted.mean(
          value$rmse^2, value$n_validation
        )) else Inf,
        complete = complete,
        evaluations = nrow(value), stringsAsFactors = FALSE
      )
    }))
    valid <- which(summary$complete & is.finite(summary$pooled_rmse))
    if (!length(valid)) stop("No valid global candidate for ", method,
                             call. = FALSE)
    best <- valid[order(
      summary$pooled_rmse[valid], summary$M[valid],
      summary$lambda[valid], summary$kappa[valid]
    )[1L]]
    rows[[method]] <- summary[best, , drop = FALSE]
  }
  bind_rows_base(rows)
}

fit_full_selected_v3 <- function(dataset, selected_candidates, config) {
  rows <- list()
  curves <- list()
  fits <- list()
  for (method in config$primary_sparse_methods) {
    candidate <- selected_candidates[
      selected_candidates$method == method, , drop = FALSE
    ]
    if (nrow(candidate) != 1L) {
      stop("Expected exactly one full-data candidate for ", method, ".",
           call. = FALSE)
    }
    fitted <- fit_neospectra_candidate_v3(
      dataset, seq_along(dataset$y), NULL, config, method, candidate,
      neospectra_v3_seed(config, "stability", match(method, config$methods))
    )
    fit <- fitted$fit
    fits[[method]] <- fitted$fit
    rows[[method]] <- data.frame(
      context = "full", method = method, candidate_id = candidate$candidate_id,
      M = candidate$M, lambda = candidate$lambda, kappa = candidate$kappa,
      elapsed_seconds = fit$elapsed, success = fit$success,
      error = fit$error, warnings = fit$warnings, stringsAsFactors = FALSE
    )
    curves[[method]] <- data.frame(
      context = "full", method = method,
      wavelength = dataset$argvals_physical,
      beta = fit$beta_physical,
      estimated_null = abs(fit$beta_scaled) <= config$zero_tolerance,
      stringsAsFactors = FALSE
    )
  }
  list(fits = fits, selection = bind_rows_base(rows),
       curves = bind_rows_base(curves))
}

fit_mcp_gamma_sensitivity_v3 <- function(dataset, selected_candidates,
                                         config) {
  candidate <- selected_candidates[
    selected_candidates$method == "FMCP", , drop = FALSE
  ]
  if (nrow(candidate) != 1L) {
    stop("MCP-gamma sensitivity requires one selected FMCP candidate.",
         call. = FALSE)
  }
  gammas <- unique(c(config$mcp_gamma, config$mcp_sensitivity))
  curve_rows <- list()
  selection_rows <- list()
  for (index in seq_along(gammas)) {
    gamma <- gammas[index]
    current <- config
    current$mcp_gamma <- gamma
    captured <- tryCatch(list(
      success = TRUE,
      value = fit_neospectra_candidate_v3(
        dataset, seq_along(dataset$y), NULL, current, "FMCP", candidate,
        neospectra_v3_seed(config, "stability", 2000L + index)
      )
    ), error = function(condition) list(
      success = FALSE, error = conditionMessage(condition)
    ))
    fit <- if (captured$success) captured$value$fit else NULL
    fit_success <- captured$success && fit$success
    fit_error <- if (captured$success) fit$error else captured$error
    fit_elapsed <- if (captured$success) fit$elapsed else NA_real_
    beta <- if (fit_success) fit$beta_physical else
      rep(NA_real_, length(dataset$argvals_physical))
    estimated_null <- if (fit_success) {
      abs(fit$beta_scaled) <= config$zero_tolerance
    } else rep(NA, length(dataset$argvals_physical))
    selection_rows[[index]] <- data.frame(
      method = "FMCP", mcp_gamma = gamma,
      primary = gamma == config$mcp_gamma,
      candidate_id = candidate$candidate_id,
      M = candidate$M, lambda = candidate$lambda, kappa = candidate$kappa,
      candidate_strategy = paste(
        "Fixed at the nested-CV candidate selected under primary gamma",
        config$mcp_gamma
      ),
      elapsed_seconds = fit_elapsed, success = fit_success,
      error = fit_error, stringsAsFactors = FALSE
    )
    curve_rows[[index]] <- data.frame(
      method = "FMCP", mcp_gamma = gamma,
      primary = gamma == config$mcp_gamma,
      wavelength = dataset$argvals_physical, beta = beta,
      estimated_null = estimated_null, success = fit_success,
      stringsAsFactors = FALSE
    )
  }
  list(
    selection = bind_rows_base(selection_rows),
    curves = bind_rows_base(curve_rows)
  )
}

run_bootstrap_stability_task_v3 <- function(task_index, dataset,
                                            selected_candidates, config) {
  set.seed(neospectra_v3_seed(config, "stability", task_index, 10L))
  index <- sample(seq_along(dataset$y), length(dataset$y), replace = TRUE)
  rows <- list()
  runtime <- list()
  for (method in config$primary_sparse_methods) {
    candidate <- selected_candidates[
      selected_candidates$method == method, , drop = FALSE
    ]
    fitted <- fit_neospectra_candidate_v3(
      dataset, index, NULL, config, method, candidate,
      neospectra_v3_seed(
        config, "stability", task_index,
        match(method, config$primary_sparse_methods) + 20L
      )
    )
    fit <- fitted$fit
    rows[[method]] <- data.frame(
      source = "bootstrap", replication = task_index, scanner_id = NA_character_,
      method = method, wavelength = dataset$argvals_physical,
      beta = fit$beta_physical,
      estimated_null = abs(fit$beta_scaled) <= config$zero_tolerance,
      success = fit$success, stringsAsFactors = FALSE
    )
    runtime[[method]] <- data.frame(
      source = "bootstrap", replication = task_index,
      scanner_id = NA_character_, method = method,
      elapsed_seconds = fit$elapsed, success = fit$success,
      error = fit$error, stringsAsFactors = FALSE
    )
  }
  list(
    key = as.character(task_index), curves = bind_rows_base(rows),
    runtime = bind_rows_base(runtime),
    draws = data.frame(
      replication = task_index, draw = seq_along(index),
      sample_id = dataset$sample_id[index], stringsAsFactors = FALSE
    )
  )
}

run_scanner_stability_v3 <- function(scan_dataset, selected_candidates,
                                     config) {
  scanners <- sort(unique(as.character(scan_dataset$scanner_id)))
  curve_rows <- list()
  runtime_rows <- list()
  for (scanner_index in seq_along(scanners)) {
    scanner <- scanners[scanner_index]
    keep <- which(scan_dataset$scanner_id == scanner)
    current <- scan_dataset
    current$subject_id <- current$sample_id[keep]
    current$sample_id <- current$sample_id[keep]
    current$country <- current$country[keep]
    current$y <- current$y[keep]
    current$x <- current$x[keep, , drop = FALSE]
    if (anyDuplicated(current$sample_id)) {
      stop("Scan-level aggregation left duplicate sample x scanner rows.",
           call. = FALSE)
    }
    for (method in config$primary_sparse_methods) {
      candidate <- selected_candidates[
        selected_candidates$method == method, , drop = FALSE
      ]
      captured <- tryCatch(list(
        success = TRUE,
        value = fit_neospectra_candidate_v3(
          current, seq_along(current$y), NULL, config, method, candidate,
          neospectra_v3_seed(
            config, "stability", 1000L + scanner_index,
            match(method, config$primary_sparse_methods)
          )
        )
      ), error = function(condition) list(
        success = FALSE, error = conditionMessage(condition)
      ))
      fit <- if (captured$success) captured$value$fit else NULL
      beta <- if (captured$success) fit$beta_physical else
        rep(NA_real_, length(current$argvals_physical))
      estimated_null <- if (captured$success) {
        abs(fit$beta_scaled) <= config$zero_tolerance
      } else rep(NA, length(current$argvals_physical))
      fit_success <- captured$success && fit$success
      fit_error <- if (captured$success) fit$error else captured$error
      fit_elapsed <- if (captured$success) fit$elapsed else NA_real_
      curve_rows[[paste(scanner, method)]] <- data.frame(
        source = "scanner", replication = NA_integer_,
        scanner_id = scanner, method = method,
        wavelength = current$argvals_physical,
        beta = beta, estimated_null = estimated_null,
        success = fit_success, stringsAsFactors = FALSE
      )
      runtime_rows[[paste(scanner, method)]] <- data.frame(
        source = "scanner", replication = NA_integer_,
        scanner_id = scanner, method = method,
        n_samples = length(current$y), elapsed_seconds = fit_elapsed,
        success = fit_success, error = fit_error,
        stringsAsFactors = FALSE
      )
    }
  }
  list(curves = bind_rows_base(curve_rows),
       runtime = bind_rows_base(runtime_rows))
}

summarize_prediction_results_v3 <- function(predictions) {
  keys <- interaction(
    predictions$evaluation, predictions$task_id, predictions$method,
    drop = TRUE, lex.order = TRUE
  )
  task_metrics <- bind_rows_base(lapply(split(predictions, keys), function(rows) {
    metrics <- prediction_metrics_v3(rows$observed, rows$predicted)
    data.frame(
      evaluation = rows$evaluation[1L], task_id = rows$task_id[1L],
      method = rows$method[1L], RMSE = unname(metrics["RMSE"]),
      MAE = unname(metrics["MAE"]), R2 = unname(metrics["R2"]),
      CCC = unname(metrics["CCC"]), n = unname(metrics["n"]),
      stringsAsFactors = FALSE
    )
  }))
  aggregate_keys <- interaction(
    predictions$evaluation, predictions$method, drop = TRUE,
    lex.order = TRUE
  )
  aggregate <- bind_rows_base(lapply(
    split(predictions, aggregate_keys), function(rows) {
      metrics <- prediction_metrics_v3(rows$observed, rows$predicted)
      data.frame(
        evaluation = rows$evaluation[1L], method = rows$method[1L],
        RMSE = unname(metrics["RMSE"]), MAE = unname(metrics["MAE"]),
        R2 = unname(metrics["R2"]), CCC = unname(metrics["CCC"]),
        n_predictions = unname(metrics["n"]), stringsAsFactors = FALSE
      )
    }
  ))
  list(task = task_metrics, aggregate = aggregate)
}

validate_grouped_nested_cv_v3 <- function(result, processed, config) {
  issues <- list()
  add <- function(check, passed, detail) {
    issues[[length(issues) + 1L]] <<- data.frame(
      check = check, passed = isTRUE(passed), detail = detail,
      stringsAsFactors = FALSE
    )
  }
  folds <- result$folds
  add("primary_sample_ids_unique",
      !anyDuplicated(processed$sample_dataset$sample_id),
      "Primary processed rows must be unique by soil sample.")
  group_fold_counts <- aggregate(
    fold ~ repeat_index + sample_id, folds,
    function(x) length(unique(x))
  )
  add("outer_groups_not_split", all(group_fold_counts$fold == 1L),
      "Each sample ID has one outer-fold label per repeat.")
  expected <- nrow(processed$sample_dataset$x) * config$outer_repeats *
    length(config$methods)
  outer_predictions <- result$predictions[
    result$predictions$evaluation == "outer_cv", , drop = FALSE
  ]
  add("outer_prediction_row_count", nrow(outer_predictions) == expected,
      paste("expected", expected, "observed", nrow(outer_predictions)))
  duplicate_predictions <- duplicated(outer_predictions[
    c("repeat_index", "sample_id", "method")
  ])
  add("outer_predictions_unique", !any(duplicate_predictions),
      "No sample may appear twice for a method within an outer repeat.")
  inner <- result$inner_folds[
    result$inner_folds$evaluation == "outer_cv", , drop = FALSE
  ]
  inner_unique <- !any(duplicated(inner[c("task_id", "method", "sample_id")]))
  add("inner_groups_not_duplicated", inner_unique,
      "Each training sample has one inner-fold label per method and task.")
  fmcp_inner <- inner[inner$method == "FMCP",
                      c("task_id", "sample_id", "inner_fold")]
  fscad_inner <- inner[inner$method == "fSCAD",
                       c("task_id", "sample_id", "inner_fold")]
  matched_inner <- merge(
    fmcp_inner, fscad_inner, by = c("task_id", "sample_id"),
    suffixes = c("_FMCP", "_fSCAD"), all = TRUE
  )
  add(
    "matched_sparse_inner_folds",
    nrow(matched_inner) == nrow(fmcp_inner) &&
      nrow(matched_inner) == nrow(fscad_inner) &&
      all(matched_inner$inner_fold_FMCP == matched_inner$inner_fold_fSCAD),
    "FMCP and fSCAD must tune on the same inner validation folds."
  )
  inner$repeat_index <- as.integer(sub(
    "^outer_r([0-9]+)_f[0-9]+$", "\\1", inner$task_id
  ))
  inner$outer_fold <- as.integer(sub(
    "^outer_r[0-9]+_f([0-9]+)$", "\\1", inner$task_id
  ))
  lookup <- paste(folds$repeat_index, folds$sample_id, sep = "\r")
  observed_outer_fold <- folds$fold[match(
    paste(inner$repeat_index, inner$sample_id, sep = "\r"), lookup
  )]
  add(
    "inner_training_excludes_outer_test",
    all(is.finite(observed_outer_fold)) &&
      all(observed_outer_fold != inner$outer_fold),
    "No outer-test sample may enter sparse-method inner tuning."
  )
  sparse_tuning <- result$tuning[
    result$tuning$evaluation == "outer_cv" &
      result$tuning$method %in% config$primary_sparse_methods,
    , drop = FALSE
  ]
  expected_grid <- sparse_candidate_grid_v3(config)
  expected_grid <- expected_grid[
    order(expected_grid$candidate_id),
    c("candidate_id", "M", "lambda", "kappa"), drop = FALSE
  ]
  rownames(expected_grid) <- NULL
  task_ids <- sort(unique(sparse_tuning$task_id))
  task_grid_ok <- vapply(task_ids, function(task_id) {
    method_ok <- vapply(config$primary_sparse_methods, function(method) {
      rows <- sparse_tuning[
        sparse_tuning$task_id == task_id & sparse_tuning$method == method,
        , drop = FALSE
      ]
      no_duplicates <- !any(duplicated(rows[c("candidate_id", "inner_fold")]))
      complete_count <- nrow(rows) == nrow(expected_grid) * config$inner_folds
      grid <- unique(rows[c("candidate_id", "M", "lambda", "kappa")])
      grid <- grid[order(grid$candidate_id), , drop = FALSE]
      rownames(grid) <- NULL
      no_duplicates && complete_count && isTRUE(all.equal(
        grid, expected_grid, check.attributes = FALSE
      ))
    }, FUN.VALUE = logical(1L))
    all(method_ok)
  }, FUN.VALUE = logical(1L))
  add(
    "complete_matched_sparse_tuning_grid_by_task",
    length(task_ids) == config$outer_repeats * config$outer_folds &&
      all(task_grid_ok),
    paste(
      "Every outer task and sparse method must contain each prespecified",
      "M/lambda/kappa candidate in every inner fold."
    )
  )
  sparse_selection <- result$selection[
    result$selection$evaluation == "outer_cv" &
      result$selection$method %in% config$primary_sparse_methods,
    , drop = FALSE
  ]
  selection_columns <- c(
    "task_id", "shape", "cutoff_scaled", "y_scale",
    "x_functional_rms_scale", "integration"
  )
  fmcp_selection <- sparse_selection[
    sparse_selection$method == "FMCP", selection_columns, drop = FALSE
  ]
  fscad_selection <- sparse_selection[
    sparse_selection$method == "fSCAD", selection_columns, drop = FALSE
  ]
  matched_selection <- merge(
    fmcp_selection, fscad_selection, by = "task_id",
    suffixes = c("_FMCP", "_fSCAD"), all = TRUE
  )
  matched_engineering <-
    nrow(matched_selection) == config$outer_repeats * config$outer_folds &&
    all(matched_selection$shape_FMCP == config$mcp_gamma) &&
    all(matched_selection$shape_fSCAD == config$scad_a) &&
    isTRUE(all.equal(
      matched_selection$cutoff_scaled_FMCP,
      matched_selection$cutoff_scaled_fSCAD, tolerance = 1e-12
    )) &&
    isTRUE(all.equal(
      matched_selection$y_scale_FMCP,
      matched_selection$y_scale_fSCAD, tolerance = 1e-12
    )) &&
    isTRUE(all.equal(
      matched_selection$x_functional_rms_scale_FMCP,
      matched_selection$x_functional_rms_scale_fSCAD, tolerance = 1e-12
    )) &&
    all(matched_selection$integration_FMCP == config$integration) &&
    all(matched_selection$integration_fSCAD == config$integration)
  add(
    "matched_sparse_shapes_cutoff_scaling_by_task", matched_engineering,
    "FMCP and fSCAD must share flat-tail shape, cutoff, scaling, and integration."
  )
  bind_rows_base(issues)
}
