# Shared helpers for the versioned ARFC-Select real-data application.
#
# The functions in this file never write raw or processed data.  They expect
# the caller to load the frozen FMCP package and the generic ARFC Rcpp core.

arfc_real_bind_rows <- function(values) {
  values <- values[!vapply(values, is.null, logical(1L))]
  values <- values[vapply(values, nrow, integer(1L)) > 0L]
  if (!length(values)) return(data.frame())
  columns <- unique(unlist(lapply(values, names), use.names = FALSE))
  values <- lapply(values, function(value) {
    missing <- setdiff(columns, names(value))
    for (name in missing) value[[name]] <- NA
    value[columns]
  })
  rownames(values[[1L]]) <- NULL
  do.call(rbind, values)
}

arfc_real_md5 <- function(paths) {
  unname(as.character(tools::md5sum(paths)))
}

arfc_real_object_md5 <- function(object) {
  path <- tempfile("arfc_real_signature_", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3)
  unname(as.character(tools::md5sum(path)))
}

arfc_real_atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(basename(path), tmpdir = dirname(path),
                        fileext = ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "xz", version = 3)
  if (!file.rename(temporary, path)) {
    stop("Atomic RDS installation failed: ", path, call. = FALSE)
  }
  invisible(path)
}

arfc_real_config <- function(project_dir, smoke = FALSE) {
  real_dir <- file.path(project_dir, "real_data")
  list(
    study = "arfc_real_data_v1",
    version = "v1",
    mode = if (smoke) "smoke" else "production",
    master_seed = if (smoke) 1210000L else 1220000L,
    dataset_order = c("Truck", "NeoSpectra"),
    structural_dataset = "Truck",
    external_dataset = "NeoSpectra",
    truck_raw = file.path(
      project_dir, "PFLR-master", "data", "truck.rda"
    ),
    neospectra_raw = file.path(
      real_dir, "data", "raw", "neospectra_v3",
      "1976_NSlibrary_withmetadata.csv"
    ),
    truck_processed = file.path(
      real_dir, "data", "processed", "v1", "truck_processed_v1.rds"
    ),
    neospectra_processed = file.path(
      real_dir, "data", "processed", "v3",
      if (smoke) "neospectra_soc_processed_v3.rds" else
        "neospectra_soc_processed_v3.rds"
    ),
    truck_folds = file.path(
      real_dir, "results", "v1", "repeated_cv_folds_v1.csv"
    ),
    truck_baselines = file.path(
      real_dir, "results", "v1", "repeated_cv_predictions_v1.csv"
    ),
    neospectra_folds = file.path(
      real_dir, "results", "v3", "neospectra",
      "nested_cv_folds_v3.csv"
    ),
    neospectra_baselines = file.path(
      real_dir, "results", "v3", "neospectra",
      "nested_cv_predictions_v3.csv"
    ),
    results_root = file.path(
      real_dir, "results", "v6", "arfc_real_data_v1"
    ),
    figures_root = file.path(
      real_dir, "figures", "v6", "arfc_real_data_v1"
    ),
    output_root = file.path(
      project_dir, "output", "pdf", "arfc_real_data_v1"
    ),
    levels = c(8L, 16L, 32L),
    finest_M = 32L,
    degree = 3L,
    penalty_order = 2L,
    gamma = 2.7,
    scad_a = 3.7,
    level_weights = rep(1 / 3, 3L),
    adaptive_alpha = 1,
    adaptive_epsilon_fraction = 0.05,
    adaptive_clip = c(0.25, 4),
    lambda = if (smoke) c(1e-9, 1e-5, 1e-1, 1e3) else
      10^seq(-9, 3, by = 0.5),
    smooth_kappa = if (smoke) c(1e-10, 1e-6, 1e-2) else
      10^seq(-13, -1, by = 1),
    max_iter = if (smoke) 80L else 120L,
    tolerance = 1e-8,
    norm_floor = 1e-8,
    relative_cutoff_multiplier = 1e-4,
    zero_tolerance = 1e-8,
    refit = TRUE,
    integration = "trapezoid",
    support_grid_points = 1001L,
    inner_folds = if (smoke) 2L else 3L,
    outer_repeats = c(Truck = 20L, NeoSpectra = 2L),
    outer_folds = c(Truck = 5L, NeoSpectra = 5L),
    smoke_outer_repeats = 1L,
    smoke_outer_folds = 2L,
    bootstrap_replicates = if (smoke) 3L else 200L,
    uncertainty_replicates = if (smoke) 100L else 2000L,
    minimum_oob_units = 20L,
    method_order = c("asr_mcp", "asr_scad", "amr_mcp", "amr_scad"),
    method_labels = c(
      asr_mcp = "ASR--MCP", asr_scad = "ASR--SCAD",
      amr_mcp = "AMR--MCP", amr_scad = "AMR--SCAD"
    ),
    secondary_methods = c("Smooth", "PCR"),
    selection_label = "ARFC-Select",
    tie_relative_tolerance = 1e-12,
    workers = {
      requested <- suppressWarnings(as.integer(Sys.getenv(
        "ARFC_REAL_WORKERS", unset = "1"
      )))
      if (!is.finite(requested) || requested < 1L) 1L else requested
    },
    runtime_is_scientific = FALSE
  )
}

arfc_real_method_metadata <- function(config) {
  data.frame(
    method = config$method_order,
    label = unname(config$method_labels[config$method_order]),
    penalty_family = c("mcp", "scad", "mcp", "scad"),
    shape = c(config$gamma, config$scad_a,
              config$gamma, config$scad_a),
    resolution = c("ASR", "ASR", "AMR", "AMR"),
    multi_resolution = c(FALSE, FALSE, TRUE, TRUE),
    adaptive = TRUE,
    tie_priority = seq_along(config$method_order),
    stringsAsFactors = FALSE
  )
}

arfc_real_load_datasets <- function(config) {
  required <- c(
    config$truck_raw, config$neospectra_raw,
    config$truck_processed, config$neospectra_processed
  )
  if (!all(file.exists(required))) {
    stop("One or more frozen processed datasets are missing.", call. = FALSE)
  }
  truck <- readRDS(config$truck_processed)
  neo_object <- readRDS(config$neospectra_processed)
  neo <- neo_object$sample_dataset
  if (!identical(arfc_real_md5(config$truck_raw),
                 as.character(truck$source_md5)) ||
      !identical(arfc_real_md5(config$neospectra_raw),
                 as.character(neo_object$source_manifest$md5[1L]))) {
    stop("Raw-source MD5 no longer matches the frozen processed object.",
         call. = FALSE)
  }
  physical_length_truck <- diff(range(truck$argvals_physical))
  physical_length_neo <- diff(range(neo$argvals_physical))
  if (!identical(dim(truck$x), c(108L, 90L)) ||
      length(truck$y) != 108L ||
      !identical(dim(neo$x), c(1974L, 601L)) ||
      length(neo$y) != 1974L ||
      any(!is.finite(truck$x)) || any(!is.finite(truck$y)) ||
      any(!is.finite(neo$x)) || any(!is.finite(neo$y)) ||
      physical_length_truck <= 0 || physical_length_neo <= 0) {
    stop("Frozen processed-data contract failed.", call. = FALSE)
  }
  values <- list(
    Truck = list(
      dataset = "Truck",
      ids = as.character(truck$subject_id),
      y = as.numeric(truck$y),
      x = unclass(as.matrix(truck$x)),
      argvals_u = as.numeric(truck$argvals_u),
      argvals_physical = as.numeric(truck$argvals_physical),
      x_axis_label = truck$x_axis_label,
      y_axis_label = truck$y_axis_label,
      response_description = truck$response_description,
      source_md5 = as.character(truck$source_md5),
      processed_md5 = arfc_real_md5(config$truck_processed),
      sampling_unit = "distributed row",
      limitation = paste(
        "No truck/test-cycle identifier is supplied; possible within-vehicle",
        "dependence cannot be grouped."
      )
    ),
    NeoSpectra = list(
      dataset = "NeoSpectra",
      ids = as.character(neo$sample_id),
      y = as.numeric(neo$y),
      x = unclass(as.matrix(neo$x)),
      argvals_u = as.numeric(neo$argvals_u),
      argvals_physical = as.numeric(neo$argvals_physical),
      x_axis_label = neo$x_axis_label,
      y_axis_label = "log1p(SOC percent)",
      response_description = neo$response_description,
      source_md5 = paste(neo_object$source_manifest$md5, collapse = "|"),
      processed_md5 = arfc_real_md5(config$neospectra_processed),
      sampling_unit = "unique soil-sample ID",
      limitation = paste(
        "Scanner and country metadata are unavailable; anonymous repeated",
        "spectra were previously averaged within sample."
      )
    )
  )
  for (name in names(values)) {
    value <- values[[name]]
    if (anyDuplicated(value$ids) || nrow(value$x) != length(value$ids) ||
        ncol(value$x) != length(value$argvals_u) ||
        length(value$argvals_physical) != length(value$argvals_u) ||
        max(abs(range(value$argvals_u) - c(0, 1))) > 1e-12) {
      stop("Dataset identity/grid contract failed for ", name,
           call. = FALSE)
    }
  }
  values[config$dataset_order]
}

arfc_real_input_inventory <- function(config) {
  paths <- c(
    truck_raw = config$truck_raw,
    neospectra_raw = config$neospectra_raw,
    truck_processed = config$truck_processed,
    neospectra_processed = config$neospectra_processed,
    truck_folds = config$truck_folds,
    truck_baselines = config$truck_baselines,
    neospectra_folds = config$neospectra_folds,
    neospectra_baselines = config$neospectra_baselines
  )
  if (!all(file.exists(paths))) {
    stop("Frozen real-data input inventory is incomplete.", call. = FALSE)
  }
  data.frame(
    role = names(paths),
    path = normalizePath(paths, mustWork = TRUE),
    bytes = as.numeric(file.info(paths)$size),
    md5 = arfc_real_md5(paths),
    stringsAsFactors = FALSE
  )
}

arfc_real_load_folds <- function(config, datasets) {
  truck <- utils::read.csv(config$truck_folds, stringsAsFactors = FALSE)
  truck <- truck[truck$dataset == "Truck", , drop = FALSE]
  truck <- data.frame(
    dataset = "Truck", repeat_index = as.integer(truck$repeat_index),
    id = as.character(truck$subject_id), fold = as.integer(truck$fold),
    source_seed = as.integer(truck$seed), stringsAsFactors = FALSE
  )
  neo <- utils::read.csv(config$neospectra_folds, stringsAsFactors = FALSE)
  neo <- data.frame(
    dataset = "NeoSpectra", repeat_index = as.integer(neo$repeat_index),
    id = as.character(neo$sample_id), fold = as.integer(neo$fold),
    source_seed = as.integer(neo$seed), stringsAsFactors = FALSE
  )
  folds <- rbind(truck, neo)
  for (dataset_name in config$dataset_order) {
    current <- folds[folds$dataset == dataset_name, , drop = FALSE]
    expected_repeats <- config$outer_repeats[[dataset_name]]
    expected_folds <- config$outer_folds[[dataset_name]]
    expected_ids <- datasets[[dataset_name]]$ids
    if (!identical(sort(unique(current$repeat_index)),
                   seq_len(expected_repeats)) ||
        !identical(sort(unique(current$fold)), seq_len(expected_folds)) ||
        nrow(current) != length(expected_ids) * expected_repeats ||
        anyDuplicated(current[c("repeat_index", "id")]) ||
        !setequal(unique(current$id), expected_ids)) {
      stop("Frozen outer-fold contract failed for ", dataset_name,
           call. = FALSE)
    }
  }
  folds
}

arfc_real_build_tasks <- function(config, folds) {
  values <- unique(folds[c("dataset", "repeat_index", "fold")])
  values <- values[order(match(values$dataset, config$dataset_order),
                         values$repeat_index, values$fold), , drop = FALSE]
  if (identical(config$mode, "smoke")) {
    values <- values[
      values$repeat_index <= config$smoke_outer_repeats &
        values$fold <= config$smoke_outer_folds, , drop = FALSE
    ]
  }
  values$task_id <- sprintf(
    "%s_r%02d_f%02d", values$dataset,
    values$repeat_index, values$fold
  )
  dataset_offset <- match(values$dataset, config$dataset_order) * 100000L
  values$inner_seed <- as.integer(
    config$master_seed + dataset_offset + values$repeat_index * 100L +
      values$fold
  )
  values <- values[c(
    "task_id", "dataset", "repeat_index", "fold", "inner_seed"
  )]
  rownames(values) <- NULL
  if (anyDuplicated(values$task_id)) {
    stop("Outer task IDs are not unique.", call. = FALSE)
  }
  values
}

arfc_real_load_secondary_predictions <- function(config, tasks) {
  truck <- utils::read.csv(config$truck_baselines, stringsAsFactors = FALSE)
  truck <- truck[
    truck$dataset == "Truck" & truck$method %in% config$secondary_methods,
    , drop = FALSE
  ]
  truck_out <- data.frame(
    task_id = sprintf("Truck_r%02d_f%02d", truck$repeat_index, truck$fold),
    dataset = "Truck", repeat_index = as.integer(truck$repeat_index),
    fold = as.integer(truck$fold), id = as.character(truck$subject_id),
    observed = as.numeric(truck$observed),
    predicted = as.numeric(truck$predicted), method = truck$method,
    success = as.logical(truck$success), source = "frozen_v1",
    stringsAsFactors = FALSE
  )
  neo <- utils::read.csv(config$neospectra_baselines,
                         stringsAsFactors = FALSE)
  neo <- neo[
    neo$evaluation == "outer_cv" &
      neo$method %in% config$secondary_methods,
    , drop = FALSE
  ]
  neo_out <- data.frame(
    task_id = sprintf("NeoSpectra_r%02d_f%02d",
                      neo$repeat_index, neo$fold),
    dataset = "NeoSpectra", repeat_index = as.integer(neo$repeat_index),
    fold = as.integer(neo$fold), id = as.character(neo$sample_id),
    observed = as.numeric(neo$observed),
    predicted = as.numeric(neo$predicted), method = neo$method,
    success = as.logical(neo$success), source = "frozen_v3",
    stringsAsFactors = FALSE
  )
  output <- rbind(truck_out, neo_out)
  output <- output[output$task_id %in% tasks$task_id, , drop = FALSE]
  rownames(output) <- NULL
  if (!nrow(output) || any(!output$success) ||
      any(!is.finite(output$observed)) ||
      any(!is.finite(output$predicted)) ||
      anyDuplicated(output[c("task_id", "id", "method")])) {
    stop("Frozen secondary benchmark predictions failed validation.",
         call. = FALSE)
  }
  output
}

arfc_real_make_inner_folds <- function(ids, fold_count, seed) {
  ids <- as.character(ids)
  unique_ids <- unique(ids)
  if (length(unique_ids) < fold_count) {
    stop("Inner-fold count exceeds training units.", call. = FALSE)
  }
  set.seed(as.integer(seed))
  shuffled <- sample(unique_ids, length(unique_ids), replace = FALSE)
  assignment <- rep(seq_len(fold_count), length.out = length(shuffled))
  lookup <- setNames(assignment, shuffled)
  labels <- unname(lookup[ids])
  if (any(!is.finite(labels))) stop("Inner-fold assignment failed.",
                                    call. = FALSE)
  as.integer(labels)
}

arfc_real_build_region_tree <- function(basis, levels, level_weights) {
  basis_count <- ncol(basis$basis_quad)
  region_count <- sum(levels)
  grams <- array(0, dim = c(basis_count, basis_count, region_count))
  metadata <- vector("list", region_count)
  cursor <- 0L
  for (level_index in seq_along(levels)) {
    level_M <- levels[level_index]
    descendants <- basis$M %/% level_M
    for (region_index in seq_len(level_M)) {
      cursor <- cursor + 1L
      first_fine <- (region_index - 1L) * descendants + 1L
      last_fine <- region_index * descendants
      rows <- basis$quad_intervals >= first_fine &
        basis$quad_intervals <= last_fine
      regional_basis <- basis$basis_quad[rows, , drop = FALSE]
      regional_weights <- basis$quad_weights[rows]
      grams[, , cursor] <- crossprod(
        regional_basis, regional_basis * regional_weights
      )
      left <- basis$breaks[first_fine]
      right <- basis$breaks[last_fine + 1L]
      metadata[[cursor]] <- data.frame(
        region_id = cursor, level_index = level_index,
        level_M = level_M, region_index = region_index,
        first_fine_interval = first_fine,
        last_fine_interval = last_fine,
        left = left, right = right, region_length = right - left,
        level_weight = level_weights[level_index],
        canonical_multiplier = level_weights[level_index] * (right - left),
        stringsAsFactors = FALSE
      )
    }
  }
  list(gram = grams, metadata = do.call(rbind, metadata))
}

arfc_real_build_context <- function(dataset, config) {
  argvals <- dataset$argvals_u
  integration_weights <- FMCP:::.fmcp_integration_weights(
    argvals, config$integration
  )
  basis <- FMCP:::.fmcp_basis_setup(
    argvals, config$finest_M, config$degree, config$penalty_order
  )
  roughness <- FMCP:::cpp_weighted_crossprod(
    basis$derivative_quad, basis$quad_weights
  )
  local_gram <- FMCP:::cpp_local_gram(
    basis$basis_quad, basis$quad_weights, basis$quad_intervals,
    config$finest_M, basis$order
  )
  region_tree <- arfc_real_build_region_tree(
    basis, config$levels, config$level_weights
  )
  region_metadata <- region_tree$metadata
  fine_rows <- region_metadata$level_M == config$finest_M
  fine_gram <- region_tree$gram[, , fine_rows, drop = FALSE]
  fine_length <- region_metadata$region_length[fine_rows]
  fine_multiplier <- fine_length
  multi_multiplier <- region_metadata$canonical_multiplier
  support_u <- seq(0, 1, length.out = config$support_grid_points)
  support_physical <- min(dataset$argvals_physical) + support_u *
    diff(range(dataset$argvals_physical))
  dense_basis <- unclass(splines::splineDesign(
    knots = basis$full_knots, x = support_u, ord = basis$order,
    derivs = rep.int(0L, length(support_u)), outer.ok = TRUE
  ))
  asr_mcp <- cpp_arfc_curvature_envelope_v1(
    fine_gram, fine_length, fine_multiplier, "mcp", config$gamma
  )
  asr_scad <- cpp_arfc_curvature_envelope_v1(
    fine_gram, fine_length, fine_multiplier, "scad", config$scad_a
  )
  amr_mcp <- cpp_arfc_curvature_envelope_v1(
    region_tree$gram, region_metadata$region_length, multi_multiplier,
    "mcp", config$gamma
  )
  amr_scad <- cpp_arfc_curvature_envelope_v1(
    region_tree$gram, region_metadata$region_length, multi_multiplier,
    "scad", config$scad_a
  )
  mass_checks <- data.frame(
    dataset = dataset$dataset,
    check = c(
      "asr_amr_mass_matrix_identity",
      "asr_mcp_scad_envelope_identity",
      "amr_mcp_scad_envelope_identity",
      "regional_multiplier_total_mass"
    ),
    observed = c(
      max(abs(asr_mcp$mass_matrix - amr_mcp$mass_matrix)),
      max(abs(asr_mcp$curvature_envelope - asr_scad$curvature_envelope)),
      max(abs(amr_mcp$curvature_envelope - amr_scad$curvature_envelope)),
      max(abs(c(sum(fine_multiplier), sum(multi_multiplier)) - 1))
    ),
    tolerance = c(1e-10, 1e-10, 1e-10, 1e-12),
    stringsAsFactors = FALSE
  )
  mass_checks$passed <- mass_checks$observed <= mass_checks$tolerance
  if (!all(mass_checks$passed)) {
    stop("ARFC mass/envelope check failed for ", dataset$dataset,
         call. = FALSE)
  }
  list(
    integration_weights = integration_weights,
    basis = basis, roughness = roughness, local_gram = local_gram,
    region_gram = region_tree$gram, region_metadata = region_metadata,
    fine_rows = fine_rows, fine_gram = fine_gram,
    fine_length = fine_length, fine_multiplier = fine_multiplier,
    multi_multiplier = multi_multiplier,
    support_u = support_u, support_physical = support_physical,
    dense_basis = dense_basis,
    interval_scale = sqrt(config$finest_M),
    mass_checks = mass_checks
  )
}

arfc_real_prepare_design <- function(dataset, context, train_index,
                                     validation_index) {
  y_train <- dataset$y[train_index]
  x_train <- dataset$x[train_index, , drop = FALSE]
  y_mean <- mean(y_train)
  y_scale <- stats::sd(y_train)
  x_mean <- colMeans(x_train)
  centered_x <- sweep(x_train, 2L, x_mean, FUN = "-")
  integrated_square <- as.numeric(
    centered_x^2 %*% context$integration_weights
  )
  x_scale <- sqrt(mean(integrated_square))
  if (!is.finite(y_scale) || y_scale <= 0 ||
      !is.finite(x_scale) || x_scale <= 0) {
    stop("Training-only response/predictor scale is not positive.",
         call. = FALSE)
  }
  design_for <- function(index) {
    centered <- sweep(dataset$x[index, , drop = FALSE],
                      2L, x_mean, FUN = "-") / x_scale
    FMCP:::cpp_design_matrix(
      centered, context$basis$basis_observed,
      context$integration_weights
    )
  }
  list(
    centered_y = (y_train - y_mean) / y_scale,
    design = design_for(train_index),
    validation_design = design_for(validation_index),
    y_mean = y_mean, y_scale = y_scale,
    x_mean = x_mean, x_scale = x_scale,
    train_index = train_index, validation_index = validation_index
  )
}

arfc_real_strict_candidate <- function(grid) {
  finite_coefficients <- apply(
    grid$coefficients, 2L, function(value) all(is.finite(value))
  )
  eligible <- which(
    grid$success %in% TRUE & grid$converged %in% TRUE &
      grid$refit_success %in% TRUE & is.finite(grid$bic) &
      finite_coefficients
  )
  if (!length(eligible)) return(NA_integer_)
  eligible[which.min(grid$bic[eligible])]
}

arfc_real_adaptive_weights <- function(pilot, context, config) {
  norms <- as.numeric(cpp_arfc_region_norms_v1(
    pilot, context$region_gram,
    context$region_metadata$region_length
  ))
  positive <- norms[is.finite(norms) & norms > 0]
  if (!length(positive)) {
    stop("Smooth pilot has no positive regional norm.", call. = FALSE)
  }
  epsilon <- max(
    1e-8, config$adaptive_epsilon_fraction * stats::median(positive)
  )
  raw <- (norms + epsilon)^(-config$adaptive_alpha)
  normalized <- raw
  for (level_M in config$levels) {
    indices <- which(context$region_metadata$level_M == level_M)
    normalized[indices] <- raw[indices] / stats::median(raw[indices])
  }
  weights <- pmin(
    config$adaptive_clip[2L],
    pmax(config$adaptive_clip[1L], normalized)
  )
  list(norms = norms, raw = raw, weights = weights, epsilon = epsilon)
}

arfc_real_select_method <- function(losses, config) {
  losses <- losses[config$method_order]
  if (any(!is.finite(losses))) {
    stop("Candidate selection losses are incomplete.", call. = FALSE)
  }
  minimum <- min(losses)
  tolerance <- config$tie_relative_tolerance * (1 + abs(minimum))
  tied <- names(losses)[losses <= minimum + tolerance]
  selected <- config$method_order[config$method_order %in% tied][1L]
  list(
    method = selected, minimum = minimum, tolerance = tolerance,
    tie_count = length(tied), tied = tied,
    margin = sort(losses)[2L] - sort(losses)[1L]
  )
}

arfc_real_fit_library <- function(dataset, context, train_index,
                                  validation_index, config, fit_id,
                                  fit_stage) {
  prepared <- arfc_real_prepare_design(
    dataset, context, train_index, validation_index
  )
  smooth <- FMCP:::cpp_fmcp_grid(
    prepared$centered_y, prepared$design, context$roughness,
    context$local_gram, 0, config$smooth_kappa, config$gamma,
    context$interval_scale, config$max_iter, config$tolerance,
    config$norm_floor, 0, config$refit, FALSE
  )
  smooth_selected <- arfc_real_strict_candidate(smooth)
  if (!is.finite(smooth_selected)) {
    stop("No strict smooth reference for ", fit_id, call. = FALSE)
  }
  pilot <- smooth$coefficients[, smooth_selected]
  global_gram <- apply(context$fine_gram, c(1L, 2L), sum)
  pilot_l2 <- sqrt(drop(crossprod(pilot, global_gram %*% pilot)))
  cutoff <- config$relative_cutoff_multiplier * pilot_l2
  if (!is.finite(cutoff) || cutoff <= 0) {
    stop("Positive coefficient cutoff unavailable for ", fit_id,
         call. = FALSE)
  }
  adaptive <- arfc_real_adaptive_weights(pilot, context, config)
  path_kappa <- as.numeric(smooth$kappa[smooth_selected])
  method_metadata <- arfc_real_method_metadata(config)
  path_rows <- metrics_rows <- curve_rows <- coefficient_rows <- list()
  prediction_matrix <- matrix(
    NA_real_, nrow = length(validation_index),
    ncol = length(config$method_order),
    dimnames = list(NULL, config$method_order)
  )
  diagnostic_rows <- list()
  for (method in config$method_order) {
    specification <- method_metadata[
      method_metadata$method == method, , drop = FALSE
    ]
    multi <- specification$multi_resolution
    start <- proc.time()[["elapsed"]]
    grid <- cpp_arfc_grid_v1(
      prepared$centered_y, prepared$design, context$roughness,
      if (multi) context$region_gram else context$fine_gram,
      if (multi) context$region_metadata$region_length else
        context$fine_length,
      if (multi) context$multi_multiplier else context$fine_multiplier,
      if (multi) adaptive$weights else
        adaptive$weights[context$fine_rows],
      config$lambda, path_kappa, specification$penalty_family,
      specification$shape, config$max_iter, config$tolerance,
      config$norm_floor, cutoff, config$refit, FALSE
    )
    elapsed <- proc.time()[["elapsed"]] - start
    finite_coefficients <- apply(
      grid$coefficients, 2L, function(value) all(is.finite(value))
    )
    usable <- grid$success %in% TRUE & grid$converged %in% TRUE &
      grid$refit_success %in% TRUE & is.finite(grid$bic) &
      finite_coefficients
    selected <- arfc_real_strict_candidate(grid)
    if (!is.finite(selected)) {
      stop("No strict candidate for ", method, " in ", fit_id,
           call. = FALSE)
    }
    coefficient <- grid$coefficients[, selected]
    prediction <- prepared$y_mean + prepared$y_scale * as.numeric(
      prepared$validation_design %*% coefficient
    )
    observed <- dataset$y[validation_index]
    residual <- observed - prediction
    prediction_matrix[, method] <- prediction
    beta_standard <- as.numeric(context$dense_basis %*% coefficient)
    physical_length <- diff(range(dataset$argvals_physical))
    beta_physical <- prepared$y_scale /
      (prepared$x_scale * physical_length) * beta_standard
    path_rows[[method]] <- data.frame(
      fit_id = fit_id, fit_stage = fit_stage,
      dataset = dataset$dataset, method = method,
      candidate_index = seq_along(grid$lambda),
      lambda = as.numeric(grid$lambda), kappa = as.numeric(grid$kappa),
      bic = as.numeric(grid$bic), rss = as.numeric(grid$rss),
      df = as.numeric(grid$df), objective = as.numeric(grid$objective),
      iterations = as.integer(grid$iterations),
      active = as.integer(grid$active), success = as.logical(grid$success),
      converged = as.logical(grid$converged),
      refit_success = as.logical(grid$refit_success),
      objective_increased = as.logical(grid$objective_increased),
      fit_usable = usable,
      bic_selected = seq_along(grid$lambda) == selected,
      message = as.character(grid$message), stringsAsFactors = FALSE
    )
    metrics_rows[[method]] <- data.frame(
      fit_id = fit_id, fit_stage = fit_stage,
      dataset = dataset$dataset, method = method,
      penalty_family = specification$penalty_family,
      resolution = specification$resolution,
      shape = specification$shape,
      selected_lambda = grid$lambda[selected], path_kappa = path_kappa,
      bic = grid$bic[selected], df = grid$df[selected],
      n_train = length(train_index), n_validation = length(validation_index),
      validation_sse = sum(residual^2),
      validation_mse = mean(residual^2),
      validation_mae = mean(abs(residual)),
      cutoff = cutoff, y_mean = prepared$y_mean,
      y_scale = prepared$y_scale, x_scale = prepared$x_scale,
      success = grid$success[selected], converged = grid$converged[selected],
      refit_success = grid$refit_success[selected],
      fit_usable = usable[selected], elapsed_seconds = elapsed,
      stringsAsFactors = FALSE
    )
    curve_rows[[method]] <- data.frame(
      fit_id = fit_id, fit_stage = fit_stage,
      dataset = dataset$dataset, method = method,
      argument_u = context$support_u,
      physical_argument = context$support_physical,
      beta_standardized = beta_standard,
      beta_physical = beta_physical,
      is_exact_zero = abs(beta_standard) <= config$zero_tolerance,
      stringsAsFactors = FALSE
    )
    coefficient_frame <- as.data.frame(t(coefficient))
    names(coefficient_frame) <- sprintf(
      "b_%03d", seq_along(coefficient)
    )
    coefficient_rows[[method]] <- cbind(
      data.frame(
        fit_id = fit_id, fit_stage = fit_stage,
        dataset = dataset$dataset, method = method,
        stringsAsFactors = FALSE
      ), coefficient_frame
    )
    diagnostic_rows[[method]] <- data.frame(
      fit_id = fit_id, fit_stage = fit_stage,
      dataset = dataset$dataset, method = method,
      usable_candidates = sum(usable),
      selected_lambda_lower_boundary =
        grid$lambda[selected] == min(config$lambda),
      selected_lambda_upper_boundary =
        grid$lambda[selected] == max(config$lambda),
      objective_increase_candidates =
        sum(grid$objective_increased %in% TRUE),
      runtime_is_scientific = config$runtime_is_scientific,
      stringsAsFactors = FALSE
    )
  }
  list(
    path = arfc_real_bind_rows(path_rows),
    metrics = arfc_real_bind_rows(metrics_rows),
    curves = arfc_real_bind_rows(curve_rows),
    coefficients = arfc_real_bind_rows(coefficient_rows),
    predictions = prediction_matrix,
    diagnostics = arfc_real_bind_rows(diagnostic_rows),
    adaptive = cbind(
      data.frame(
        fit_id = fit_id, fit_stage = fit_stage,
        dataset = dataset$dataset, epsilon = adaptive$epsilon,
        stringsAsFactors = FALSE
      ), context$region_metadata,
      pilot_norm = adaptive$norms,
      raw_adaptive_weight = adaptive$raw,
      adaptive_weight = adaptive$weights
    ),
    smooth = data.frame(
      fit_id = fit_id, fit_stage = fit_stage,
      dataset = dataset$dataset,
      candidate_index = seq_along(smooth$kappa),
      kappa = smooth$kappa, bic = smooth$bic,
      success = smooth$success, converged = smooth$converged,
      refit_success = smooth$refit_success,
      selected = seq_along(smooth$kappa) == smooth_selected,
      stringsAsFactors = FALSE
    )
  )
}

arfc_real_run_outer_task <- function(task, dataset, folds, context, config) {
  task <- as.list(task)
  current_folds <- folds[
    folds$dataset == task$dataset &
      folds$repeat_index == task$repeat_index, , drop = FALSE
  ]
  fold_lookup <- setNames(current_folds$fold, current_folds$id)
  observed_fold <- unname(fold_lookup[dataset$ids])
  if (any(!is.finite(observed_fold))) {
    stop("Outer fold lookup failed for ", task$task_id, call. = FALSE)
  }
  test_index <- which(observed_fold == task$fold)
  train_index <- which(observed_fold != task$fold)
  inner_labels <- arfc_real_make_inner_folds(
    dataset$ids[train_index], config$inner_folds, task$inner_seed
  )
  inner_rows <- path_rows <- diagnostic_rows <- list()
  for (inner_fold in seq_len(config$inner_folds)) {
    inner_train <- train_index[inner_labels != inner_fold]
    inner_validation <- train_index[inner_labels == inner_fold]
    fit_id <- paste0(task$task_id, "_inner", inner_fold)
    fit <- arfc_real_fit_library(
      dataset, context, inner_train, inner_validation, config,
      fit_id, "inner_validation"
    )
    current <- fit$metrics
    current$task_id <- task$task_id
    current$repeat_index <- task$repeat_index
    current$outer_fold <- task$fold
    current$inner_fold <- inner_fold
    inner_rows[[inner_fold]] <- current
    fit$path$task_id <- task$task_id
    fit$path$repeat_index <- task$repeat_index
    fit$path$outer_fold <- task$fold
    fit$path$inner_fold <- inner_fold
    path_rows[[inner_fold]] <- fit$path
    fit$diagnostics$task_id <- task$task_id
    fit$diagnostics$repeat_index <- task$repeat_index
    fit$diagnostics$outer_fold <- task$fold
    fit$diagnostics$inner_fold <- inner_fold
    diagnostic_rows[[inner_fold]] <- fit$diagnostics
  }
  inner <- arfc_real_bind_rows(inner_rows)
  pooled_sse <- tapply(inner$validation_sse, inner$method, sum)
  pooled_n <- tapply(inner$n_validation, inner$method, sum)
  pooled_mse <- pooled_sse / pooled_n
  selection <- arfc_real_select_method(pooled_mse, config)
  outer <- arfc_real_fit_library(
    dataset, context, train_index, test_index, config,
    task$task_id, "outer_refit"
  )
  outer$path$task_id <- task$task_id
  outer$path$repeat_index <- task$repeat_index
  outer$path$outer_fold <- task$fold
  outer$path$inner_fold <- NA_integer_
  path_rows[[length(path_rows) + 1L]] <- outer$path
  outer$diagnostics$task_id <- task$task_id
  outer$diagnostics$repeat_index <- task$repeat_index
  outer$diagnostics$outer_fold <- task$fold
  outer$diagnostics$inner_fold <- NA_integer_
  diagnostic_rows[[length(diagnostic_rows) + 1L]] <- outer$diagnostics
  fixed_predictions <- do.call(rbind, lapply(config$method_order, function(method) {
    data.frame(
      task_id = task$task_id, dataset = task$dataset,
      repeat_index = task$repeat_index, fold = task$fold,
      id = dataset$ids[test_index], observed = dataset$y[test_index],
      predicted = outer$predictions[, method], method = method,
      success = TRUE, source = "arfc_real_v1",
      stringsAsFactors = FALSE
    )
  }))
  selected_prediction <- fixed_predictions[
    fixed_predictions$method == selection$method, , drop = FALSE
  ]
  selected_prediction$method <- config$selection_label
  selected_prediction$source <- "nested_inner_selection"
  predictions <- rbind(fixed_predictions, selected_prediction)
  selection_row <- data.frame(
    task_id = task$task_id, dataset = task$dataset,
    repeat_index = task$repeat_index, fold = task$fold,
    selected_method = selection$method,
    selected_inner_validation_mse = selection$minimum,
    validation_margin_to_second = selection$margin,
    validation_tie_count = selection$tie_count,
    inner_seed = task$inner_seed,
    outer_test_used_for_selection = FALSE,
    practical_outer_refit = TRUE,
    stringsAsFactors = FALSE
  )
  list(
    task_id = task$task_id,
    inner_validation = inner,
    path = arfc_real_bind_rows(path_rows),
    diagnostics = arfc_real_bind_rows(diagnostic_rows),
    predictions = predictions,
    selection = selection_row,
    outer_candidate_metrics = transform(
      outer$metrics, task_id = task$task_id,
      repeat_index = task$repeat_index, outer_fold = task$fold
    )
  )
}

arfc_real_select_full_candidate <- function(predictions, config) {
  candidate <- predictions[
    predictions$method %in% config$method_order, , drop = FALSE
  ]
  mse <- tapply((candidate$observed - candidate$predicted)^2,
                candidate$method, mean)
  arfc_real_select_method(mse, config)
}

arfc_real_run_full_fit <- function(dataset, context, outer_predictions,
                                   config) {
  selection <- arfc_real_select_full_candidate(
    outer_predictions[outer_predictions$dataset == dataset$dataset, ],
    config
  )
  indices <- seq_along(dataset$y)
  fit <- arfc_real_fit_library(
    dataset, context, indices, indices, config,
    paste0(dataset$dataset, "_full"), "full_data"
  )
  fit$curves$selected_full_candidate <-
    fit$curves$method == selection$method
  fit$metrics$selected_full_candidate <-
    fit$metrics$method == selection$method
  list(
    dataset = dataset$dataset,
    selection = data.frame(
      dataset = dataset$dataset, selected_method = selection$method,
      pooled_outer_cv_mse = selection$minimum,
      cv_margin_to_second = selection$margin,
      cv_tie_count = selection$tie_count,
      rule = "minimum pooled fixed-candidate outer-CV response MSE",
      stringsAsFactors = FALSE
    ),
    curves = fit$curves, metrics = fit$metrics,
    path = fit$path, diagnostics = fit$diagnostics,
    adaptive = fit$adaptive, smooth = fit$smooth
  )
}

arfc_real_run_bootstrap <- function(dataset, context, replication, config) {
  seed <- as.integer(
    config$master_seed + match(dataset$dataset, config$dataset_order) *
      200000L + replication
  )
  set.seed(seed)
  n <- length(dataset$y)
  draw <- sample.int(n, n, replace = TRUE)
  oob <- setdiff(seq_len(n), unique(draw))
  if (length(oob) < config$minimum_oob_units) {
    stop("Bootstrap OOB sample is too small for ", dataset$dataset,
         " replication ", replication, call. = FALSE)
  }
  fit_id <- sprintf("%s_bootstrap_%03d", dataset$dataset, replication)
  fit <- arfc_real_fit_library(
    dataset, context, draw, oob, config, fit_id, "bootstrap_oob"
  )
  losses <- setNames(fit$metrics$validation_mse, fit$metrics$method)
  selection <- arfc_real_select_method(losses, config)
  selected_curve <- fit$curves[
    fit$curves$method == selection$method, , drop = FALSE
  ]
  selected_curve$replication <- replication
  selected_curve$seed <- seed
  list(
    dataset = dataset$dataset,
    replication = replication,
    selection = data.frame(
      dataset = dataset$dataset, replication = replication, seed = seed,
      selected_method = selection$method,
      oob_mse = selection$minimum, margin_to_second = selection$margin,
      tie_count = selection$tie_count,
      n_bootstrap = n, n_unique_bootstrap = length(unique(draw)),
      n_oob = length(oob), full_reselection = TRUE,
      stringsAsFactors = FALSE
    ),
    curves = selected_curve,
    metrics = fit$metrics,
    path = fit$path,
    diagnostics = fit$diagnostics
  )
}

arfc_real_manifest <- function(directory, scientific_signature,
                               files, manifest_name) {
  paths <- file.path(directory, files)
  if (!all(file.exists(paths))) {
    stop("Manifest payload is incomplete.", call. = FALSE)
  }
  manifest <- data.frame(
    scientific_signature = scientific_signature,
    file = files,
    bytes = as.numeric(file.info(paths)$size),
    md5 = arfc_real_md5(paths),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(directory, manifest_name),
                   row.names = FALSE)
  manifest
}
