#!/usr/bin/env Rscript

# Fresh confirmatory four-candidate ARFC selection experiment.
# Training tunes each candidate; an independent validation sample selects
# ASR/AMR and MCP/SCAD; the test sample is evaluation-only.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/45_run_arfc_select_v1.R", mustWork = TRUE)
}
simulation_dir <- dirname(script_file)
project_dir <- dirname(simulation_dir)

arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
resume <- "--resume" %in% arguments
validate_design_only <- "--validate-design" %in% arguments
unknown <- setdiff(arguments, c("--smoke", "--resume", "--validate-design"))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}
if (resume && validate_design_only) {
  stop("--resume and --validate-design cannot be combined.", call. = FALSE)
}

library_override <- Sys.getenv(
  "ARFC_FMCP_LIBRARY", unset = ".replication-library"
)
if (!dir.exists(library_override)) {
  stop("Frozen FMCP 0.2.0 library unavailable: ", library_override,
       call. = FALSE)
}
.libPaths(c(normalizePath(library_override), .libPaths()))
required_packages <- c("FMCP", "Rcpp", "RcppArmadillo")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE,
  FUN.VALUE = logical(1L)
)]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
if (utils::packageVersion("FMCP") != numeric_version("0.2.0")) {
  stop("ARFC selection requires frozen FMCP 0.2.0.", call. = FALSE)
}

protocol_file <- file.path(
  project_dir, "config", "protocol_arfc_select_v1.md"
)
framework_protocol_file <- file.path(
  project_dir, "config",
  "protocol_adaptive_regional_folded_concave_v1.md"
)
core_file <- file.path(simulation_dir, "src", "adaptive_regional_core.cpp")
required_sources <- c(
  protocol_file, framework_protocol_file, core_file, script_file
)
if (!all(file.exists(required_sources))) {
  stop("One or more frozen ARFC selection sources are missing.",
       call. = FALSE)
}

Rcpp::sourceCpp(core_file, rebuild = FALSE, verbose = FALSE)
required_exports <- c(
  "cpp_arfc_grid_v1", "cpp_arfc_curvature_envelope_v1",
  "cpp_arfc_region_norms_v1"
)
if (!all(vapply(required_exports, exists, mode = "function",
                FUN.VALUE = logical(1L)))) {
  stop("Generic ARFC RcppArmadillo exports were not loaded.",
       call. = FALSE)
}

config <- list(
  study = "arfc_select_v1",
  protocol_id = "arfc_select_v1",
  mode = if (smoke) "smoke" else "production",
  confirmatory = !smoke,
  master_seed = if (smoke) 1110000L else 1120000L,
  n_train = if (smoke) 60L else 150L,
  n_validation = if (smoke) 60L else 150L,
  n_test = if (smoke) 250L else 2000L,
  replications = if (smoke) 1L else 100L,
  observed_grid_points = 201L,
  support_grid_points = 1001L,
  predictor_order = 5L,
  predictor_knots = 71L,
  intercept = 1,
  snr = 4,
  case_i_noise_sd = 1,
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
  lambda = if (smoke) {
    c(1e-9, 1e-5, 1e-1, 1e3)
  } else {
    10^seq(-9, 3, by = 0.5)
  },
  smooth_kappa = if (smoke) {
    c(1e-10, 1e-6, 1e-2)
  } else {
    10^seq(-13, -1, by = 1)
  },
  max_iter = if (smoke) 80L else 120L,
  tolerance = 1e-8,
  norm_floor = 1e-8,
  relative_cutoff_multiplier = 1e-4,
  zero_tolerance = 1e-8,
  refit = TRUE,
  integration = "simpson",
  method_order = c("asr_mcp", "asr_scad", "amr_mcp", "amr_scad"),
  method_labels = c(
    asr_mcp = "ASR--MCP", asr_scad = "ASR--SCAD",
    amr_mcp = "AMR--MCP", amr_scad = "AMR--SCAD"
  ),
  validation_tie_relative_tolerance = 1e-12,
  selection_rule = paste(
    "minimum independent-validation response MSE; numerical ties use",
    "asr_mcp, asr_scad, amr_mcp, amr_scad"
  ),
  strict_bic_tie_rule =
    "first exact minimum in common ascending-lambda order",
  runtime_is_scientific = FALSE,
  test_used_for_selection = FALSE,
  package_versions = vapply(
    required_packages,
    function(package) as.character(utils::packageVersion(package)),
    FUN.VALUE = character(1L)
  )
)

scenarios <- data.frame(
  scenario_index = 1:4,
  scenario = c(
    "case_i_null", "case_ii_sparse",
    "case_iii_crossing", "case_iv_away"
  ),
  case_label = c("Case I", "Case II", "Case III", "Case IV"),
  role = c(
    "safety_all_null", "primary_local_null",
    "safety_no_null", "safety_no_null"
  ),
  null_type = c("all", "local", "none", "none"),
  null_left = c(0, 0.3, NA, NA),
  null_right = c(1, 0.7, NA, NA),
  stringsAsFactors = FALSE
)

tasks <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(index) {
  scenario <- scenarios[index, , drop = FALSE]
  replications <- seq_len(config$replications)
  seed_base <- config$master_seed + index * 10000L + replications * 10L
  data.frame(
    task_id = (index - 1L) * config$replications + replications,
    scenario_index = scenario$scenario_index,
    scenario = scenario$scenario,
    case_label = scenario$case_label,
    role = scenario$role,
    null_type = scenario$null_type,
    null_left = scenario$null_left,
    null_right = scenario$null_right,
    n_train = config$n_train,
    n_validation = config$n_validation,
    n_test = config$n_test,
    replication = replications,
    seed_train = seed_base + 1L,
    seed_validation = seed_base + 2L,
    seed_test = seed_base + 3L,
    stringsAsFactors = FALSE
  )
}))
rownames(tasks) <- NULL

seed_values <- unlist(tasks[c("seed_train", "seed_validation", "seed_test")],
                      use.names = FALSE)
expected_tasks <- 4L * config$replications
if (nrow(tasks) != expected_tasks || anyDuplicated(tasks$task_id) ||
    anyDuplicated(seed_values) ||
    !identical(config$levels, c(8L, 16L, 32L)) ||
    any(config$finest_M %% config$levels != 0L) ||
    abs(sum(config$level_weights) - 1) > 1e-14 ||
    abs(1 / config$gamma - 1 / (config$scad_a - 1)) > 1e-14 ||
    !identical(config$method_order,
               c("asr_mcp", "asr_scad", "amr_mcp", "amr_scad")) ||
    (!smoke && length(config$lambda) != 25L)) {
  stop("The frozen ARFC selection design failed its static audit.",
       call. = FALSE)
}

md5_character <- function(paths) {
  unname(as.character(tools::md5sum(paths)))
}
object_md5 <- function(object) {
  temporary <- tempfile("arfc_select_signature_", fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  unname(as.character(tools::md5sum(temporary)))
}
source_md5 <- setNames(
  md5_character(required_sources), basename(required_sources)
)
package_source_files <- c(
  file.path(project_dir, "FMCP-package", "src", "fmcp_core.cpp"),
  file.path(project_dir, "FMCP-package", "src", "fscad_core.cpp"),
  file.path(project_dir, "FMCP-package", "R", "basis.R")
)
if (!all(file.exists(package_source_files))) {
  stop("Package provenance sources are missing.", call. = FALSE)
}
package_source_md5 <- setNames(
  md5_character(package_source_files), basename(package_source_files)
)
scientific_signature <- object_md5(list(
  config = config, scenarios = scenarios, tasks = tasks,
  source_md5 = source_md5, package_source_md5 = package_source_md5
))

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
domain_length <- diff(range(argvals))

coefficient_basis <- FMCP:::.fmcp_basis_setup(
  argvals, config$finest_M, config$degree, config$penalty_order
)
dense_basis <- unclass(splines::splineDesign(
  knots = coefficient_basis$full_knots,
  x = support_grid,
  ord = coefficient_basis$order,
  derivs = rep.int(0L, length(support_grid)),
  outer.ok = TRUE
))
local_gram <- FMCP:::cpp_local_gram(
  coefficient_basis$basis_quad,
  coefficient_basis$quad_weights,
  coefficient_basis$quad_intervals,
  config$finest_M,
  coefficient_basis$order
)
roughness <- FMCP:::cpp_weighted_crossprod(
  coefficient_basis$derivative_quad, coefficient_basis$quad_weights
)
interval_scale <- sqrt(config$finest_M / domain_length)

build_region_tree <- function(basis, levels, level_weights) {
  basis_count <- ncol(basis$basis_quad)
  region_count <- sum(levels)
  grams <- array(0, dim = c(basis_count, basis_count, region_count))
  metadata <- vector("list", region_count)
  cursor <- 0L
  for (level_index in seq_along(levels)) {
    level_M <- levels[level_index]
    descendant_count <- basis$M %/% level_M
    for (region_index in seq_len(level_M)) {
      cursor <- cursor + 1L
      first_fine <- (region_index - 1L) * descendant_count + 1L
      last_fine <- region_index * descendant_count
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
        region_id = cursor,
        level_index = level_index,
        level_M = level_M,
        region_index = region_index,
        first_fine_interval = first_fine,
        last_fine_interval = last_fine,
        left = left,
        right = right,
        region_length = right - left,
        level_weight = level_weights[level_index],
        canonical_multiplier = level_weights[level_index] * (right - left),
        stringsAsFactors = FALSE
      )
    }
  }
  list(gram = grams, metadata = do.call(rbind, metadata))
}

region_tree <- build_region_tree(
  coefficient_basis, config$levels, config$level_weights
)
region_metadata <- region_tree$metadata
fine_rows <- region_metadata$level_M == config$finest_M
fine_gram <- region_tree$gram[, , fine_rows, drop = FALSE]
fine_length <- region_metadata$region_length[fine_rows]
fine_multiplier <- fine_length / domain_length
multi_multiplier <- region_metadata$canonical_multiplier / domain_length

mass_asr_mcp <- cpp_arfc_curvature_envelope_v1(
  fine_gram, fine_length, fine_multiplier, "mcp", config$gamma
)
mass_asr_scad <- cpp_arfc_curvature_envelope_v1(
  fine_gram, fine_length, fine_multiplier, "scad", config$scad_a
)
mass_amr_mcp <- cpp_arfc_curvature_envelope_v1(
  region_tree$gram, region_metadata$region_length, multi_multiplier,
  "mcp", config$gamma
)
mass_amr_scad <- cpp_arfc_curvature_envelope_v1(
  region_tree$gram, region_metadata$region_length, multi_multiplier,
  "scad", config$scad_a
)
mass_checks <- data.frame(
  check = c(
    "asr_amr_mass_matrix_identity",
    "asr_mcp_scad_envelope_identity",
    "amr_mcp_scad_envelope_identity",
    "regional_multiplier_total_mass"
  ),
  observed = c(
    max(abs(mass_asr_mcp$mass_matrix - mass_amr_mcp$mass_matrix)),
    max(abs(mass_asr_mcp$curvature_envelope -
              mass_asr_scad$curvature_envelope)),
    max(abs(mass_amr_mcp$curvature_envelope -
              mass_amr_scad$curvature_envelope)),
    max(abs(c(sum(fine_multiplier), sum(multi_multiplier)) - 1))
  ),
  tolerance = c(1e-10, 1e-10, 1e-10, 1e-12),
  stringsAsFactors = FALSE
)
mass_checks$passed <- mass_checks$observed <= mass_checks$tolerance
if (!all(mass_checks$passed)) {
  stop("Canonical ARFC mass/envelope checks failed.", call. = FALSE)
}

method_metadata <- data.frame(
  method = config$method_order,
  label = unname(config$method_labels[config$method_order]),
  penalty_family = c("mcp", "scad", "mcp", "scad"),
  shape = c(config$gamma, config$scad_a,
            config$gamma, config$scad_a),
  stage = c("ASR", "ASR", "AMR", "AMR"),
  multi_resolution = c(FALSE, FALSE, TRUE, TRUE),
  adaptive = TRUE,
  region_count = c(
    config$finest_M, config$finest_M,
    nrow(region_metadata), nrow(region_metadata)
  ),
  canonical_total_mass = 1,
  validation_tie_priority = seq_along(config$method_order),
  stringsAsFactors = FALSE
)

if (validate_design_only) {
  cat("Study:", config$study, "\n")
  cat("Mode:", config$mode, "\n")
  cat("Tasks:", nrow(tasks), "\n")
  cat("Samples: train", config$n_train, "validation", config$n_validation,
      "test", config$n_test, "\n")
  cat("Methods:", paste(config$method_order, collapse = ", "), "\n")
  cat("Lambda candidates:", length(config$lambda), "\n")
  cat("Mass checks:", sum(mass_checks$passed), "/",
      nrow(mass_checks), "PASS\n")
  cat("Scientific signature:", scientific_signature, "\n")
  quit(save = "no", status = 0L)
}

predictor_basis <- local({
  knots <- seq(0, 1, length.out = config$predictor_knots)
  unclass(splines::bs(
    argvals,
    knots = knots[-c(1L, length(knots))],
    degree = config$predictor_order - 1L,
    intercept = TRUE,
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
  if (identical(scenario, "case_i_null")) {
    return(numeric(length(values)))
  }
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
  if (identical(scenario, "case_i_null")) {
    return(rep(TRUE, length(values)))
  }
  if (identical(scenario, "case_ii_sparse")) {
    return(values >= 0.3 & values <= 0.7)
  }
  if (scenario %in% c("case_iii_crossing", "case_iv_away")) {
    return(rep(FALSE, length(values)))
  }
  stop("Unknown scenario: ", scenario, call. = FALSE)
}

generate_sample <- function(sample_size, seed, truth_observed, noise_sd) {
  set.seed(as.integer(seed))
  x <- generate_predictors(sample_size)
  linear <- as.numeric(x %*% (integration_weights * truth_observed))
  signal <- config$intercept + linear
  noise <- stats::rnorm(sample_size, sd = noise_sd)
  list(x = x, y = signal + noise, signal = signal, noise = noise,
       linear = linear)
}

generate_task_data <- function(task) {
  truth_observed <- truth_function(task$scenario, argvals)
  truth_dense <- truth_function(task$scenario, support_grid)

  set.seed(as.integer(task$seed_train))
  x_train <- generate_predictors(task$n_train)
  train_linear <- as.numeric(
    x_train %*% (integration_weights * truth_observed)
  )
  signal_variance <- stats::var(train_linear)
  case_i <- identical(task$scenario, "case_i_null")
  noise_sd <- if (case_i) {
    config$case_i_noise_sd
  } else {
    sqrt(signal_variance / config$snr)
  }
  if (!is.finite(noise_sd) || noise_sd <= 0) {
    stop("Invalid response-noise scale for task ", task$task_id,
         call. = FALSE)
  }
  train_noise <- stats::rnorm(task$n_train, sd = noise_sd)
  signal_train <- config$intercept + train_linear
  train <- list(
    x = x_train, y = signal_train + train_noise,
    signal = signal_train, noise = train_noise, linear = train_linear
  )
  validation <- generate_sample(
    task$n_validation, task$seed_validation, truth_observed, noise_sd
  )
  test <- generate_sample(
    task$n_test, task$seed_test, truth_observed, noise_sd
  )
  list(
    truth_observed = truth_observed,
    truth_dense = truth_dense,
    true_null_dense = true_null_mask(task$scenario, support_grid),
    train = train, validation = validation, test = test,
    signal_variance = signal_variance,
    noise_sd = noise_sd,
    target_variance_snr = if (case_i) NA_real_ else {
      signal_variance / noise_sd^2
    },
    empirical_variance_snr = if (case_i) NA_real_ else {
      signal_variance / stats::var(train_noise)
    }
  )
}

prepare_design <- function(data) {
  means <- list(y = mean(data$train$y), x = colMeans(data$train$x))
  centered_y <- data$train$y - means$y
  centered_x <- sweep(data$train$x, 2L, means$x, FUN = "-")
  training_design <- FMCP:::cpp_design_matrix(
    centered_x, coefficient_basis$basis_observed, integration_weights
  )
  make_design <- function(x) {
    centered <- sweep(x, 2L, means$x, FUN = "-")
    FMCP:::cpp_design_matrix(
      centered, coefficient_basis$basis_observed, integration_weights
    )
  }
  list(
    centered_y = centered_y,
    design = training_design,
    validation_design = make_design(data$validation$x),
    test_design = make_design(data$test$x),
    means = means
  )
}

strict_candidate <- function(grid) {
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

fit_smooth_reference <- function(prepared) {
  grid <- FMCP:::cpp_fmcp_grid(
    prepared$centered_y, prepared$design, roughness, local_gram,
    0, config$smooth_kappa, config$gamma, interval_scale,
    config$max_iter, config$tolerance, config$norm_floor,
    0, config$refit, FALSE
  )
  selected <- strict_candidate(grid)
  if (!is.finite(selected)) {
    stop("No strict smooth-reference candidate was available.",
         call. = FALSE)
  }
  list(grid = grid, selected = selected)
}

adaptive_weights <- function(pilot_coefficients) {
  norms <- as.numeric(cpp_arfc_region_norms_v1(
    pilot_coefficients, region_tree$gram, region_metadata$region_length
  ))
  positive <- norms[is.finite(norms) & norms > 0]
  if (!length(positive)) {
    stop("Smooth pilot produced no positive regional norm.", call. = FALSE)
  }
  epsilon <- max(
    1e-8,
    config$adaptive_epsilon_fraction * stats::median(positive)
  )
  raw <- (norms + epsilon)^(-config$adaptive_alpha)
  normalized <- raw
  for (level_M in config$levels) {
    indices <- which(region_metadata$level_M == level_M)
    normalized[indices] <- raw[indices] / stats::median(raw[indices])
  }
  clipped <- pmin(
    config$adaptive_clip[2L],
    pmax(config$adaptive_clip[1L], normalized)
  )
  list(norms = norms, epsilon = epsilon, raw = raw, weights = clipped)
}

fit_one_method <- function(method, prepared, path_kappa, cutoff, weights) {
  specification <- method_metadata[method_metadata$method == method, ]
  multi <- specification$multi_resolution
  start <- proc.time()[["elapsed"]]
  grid <- cpp_arfc_grid_v1(
    prepared$centered_y, prepared$design, roughness,
    if (multi) region_tree$gram else fine_gram,
    if (multi) region_metadata$region_length else fine_length,
    if (multi) multi_multiplier else fine_multiplier,
    if (multi) weights else weights[fine_rows],
    config$lambda, path_kappa, specification$penalty_family,
    specification$shape, config$max_iter, config$tolerance,
    config$norm_floor, cutoff, config$refit, FALSE
  )
  list(grid = grid, elapsed = proc.time()[["elapsed"]] - start)
}

weighted_integral <- function(values, mask = rep(TRUE, length(values))) {
  if (!any(mask) || any(!is.finite(values[mask]))) return(NA_real_)
  sum(support_weights[mask] * values[mask])
}

support_metric_names <- c(
  "ise_total", "ise_null", "ise_nonnull", "true_zero_rate",
  "false_zero_rate", "false_active_proportion", "null_iou",
  "null_symmetric_difference", "estimated_zero_length",
  "null_component_count", "selected_null_left", "selected_null_right",
  "boundary_error"
)

support_metrics <- function(beta, truth, task) {
  if (length(beta) != length(support_grid) || any(!is.finite(beta))) {
    return(setNames(rep(NA_real_, length(support_metric_names)),
                    support_metric_names))
  }
  truth_null <- true_null_mask(task$scenario, support_grid)
  truth_nonnull <- !truth_null & abs(truth) > 1e-10
  estimated_zero <- abs(beta) <= config$zero_tolerance
  null_length <- sum(support_weights[truth_null])
  nonnull_length <- sum(support_weights[truth_nonnull])
  estimated_length <- sum(support_weights[estimated_zero])
  intersection <- sum(support_weights[estimated_zero & truth_null])
  false_length <- sum(support_weights[estimated_zero & truth_nonnull])
  local_null <- identical(task$scenario, "case_ii_sparse")

  component_count <- NA_integer_
  selected_left <- selected_right <- boundary_error <- NA_real_
  null_iou <- null_symmetric_difference <- NA_real_
  if (local_null) {
    union <- null_length + estimated_length - intersection
    null_iou <- if (union > 0) intersection / union else 1
    null_symmetric_difference <-
      (null_length + estimated_length - 2 * intersection) / domain_length
    runs <- rle(estimated_zero)
    zero_runs <- which(runs$values)
    component_count <- length(zero_runs)
    if (!length(zero_runs)) {
      boundary_error <- 1
    } else {
      ends <- cumsum(runs$lengths)
      starts <- ends - runs$lengths + 1L
      overlap <- vapply(zero_runs, function(run) {
        indices <- starts[run]:ends[run]
        sum(support_weights[indices] * truth_null[indices])
      }, numeric(1L))
      component_length <- vapply(zero_runs, function(run) {
        sum(support_weights[starts[run]:ends[run]])
      }, numeric(1L))
      chosen <- zero_runs[order(-overlap, -component_length)[1L]]
      selected_left <- support_grid[starts[chosen]]
      selected_right <- support_grid[ends[chosen]]
      boundary_error <- mean(abs(c(
        selected_left - task$null_left,
        selected_right - task$null_right
      )))
    }
  }
  c(
    ise_total = weighted_integral((beta - truth)^2),
    ise_null = weighted_integral((beta - truth)^2, truth_null),
    ise_nonnull = weighted_integral((beta - truth)^2, truth_nonnull),
    true_zero_rate = if (null_length > 0) {
      intersection / null_length
    } else NA_real_,
    false_zero_rate = if (nonnull_length > 0) {
      false_length / nonnull_length
    } else NA_real_,
    false_active_proportion = if (
      identical(task$scenario, "case_i_null")
    ) {
      sum(support_weights[!estimated_zero]) / domain_length
    } else NA_real_,
    null_iou = null_iou,
    null_symmetric_difference = null_symmetric_difference,
    estimated_zero_length = estimated_length,
    null_component_count = component_count,
    selected_null_left = selected_left,
    selected_null_right = selected_right,
    boundary_error = boundary_error
  )
}

evaluate_method <- function(fitted, method, task, data, prepared,
                            path_kappa, cutoff) {
  grid <- fitted$grid
  finite_coefficients <- apply(
    grid$coefficients, 2L, function(value) all(is.finite(value))
  )
  usable <- grid$success %in% TRUE & grid$converged %in% TRUE &
    grid$refit_success %in% TRUE & is.finite(grid$bic) &
    finite_coefficients
  selected <- strict_candidate(grid)
  if (!is.finite(selected)) {
    stop("No strict BIC candidate for ", method, " in task ",
         task$task_id, call. = FALSE)
  }
  specification <- method_metadata[method_metadata$method == method, ]
  path <- data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    case_label = task$case_label,
    role = task$role,
    replication = task$replication,
    seed_train = task$seed_train,
    seed_validation = task$seed_validation,
    seed_test = task$seed_test,
    method = method,
    penalty_family = specification$penalty_family,
    stage = specification$stage,
    shape = specification$shape,
    multi_resolution = specification$multi_resolution,
    candidate_index = seq_along(grid$lambda),
    lambda = as.numeric(grid$lambda),
    kappa = as.numeric(grid$kappa),
    bic = as.numeric(grid$bic),
    rss = as.numeric(grid$rss),
    df = as.numeric(grid$df),
    objective = as.numeric(grid$objective),
    iterations = as.integer(grid$iterations),
    active = as.integer(grid$active),
    success = as.logical(grid$success),
    converged = as.logical(grid$converged),
    refit_success = as.logical(grid$refit_success),
    objective_increased = as.logical(grid$objective_increased),
    fit_usable = usable,
    bic_selected = seq_along(grid$lambda) == selected,
    message = as.character(grid$message),
    stringsAsFactors = FALSE
  )
  coefficient <- grid$coefficients[, selected]
  beta <- as.numeric(dense_basis %*% coefficient)
  structural <- support_metrics(beta, data$truth_dense, task)
  validation_prediction <- as.numeric(
    prepared$means$y + prepared$validation_design %*% coefficient
  )
  test_prediction <- as.numeric(
    prepared$means$y + prepared$test_design %*% coefficient
  )
  candidate_metrics <- cbind(
    path[selected, , drop = FALSE],
    as.data.frame(as.list(structural)),
    validation_signal_mse = mean(
      (validation_prediction - data$validation$signal)^2
    ),
    validation_response_mse = mean(
      (validation_prediction - data$validation$y)^2
    ),
    signal_prediction_mse = mean(
      (test_prediction - data$test$signal)^2
    ),
    response_prediction_mse = mean(
      (test_prediction - data$test$y)^2
    )
  )
  diagnostics <- data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    case_label = task$case_label,
    replication = task$replication,
    method = method,
    penalty_family = specification$penalty_family,
    stage = specification$stage,
    path_kappa = path_kappa,
    cutoff = cutoff,
    selected_lambda = grid$lambda[selected],
    selected_lambda_lower_boundary =
      grid$lambda[selected] == min(config$lambda),
    selected_lambda_upper_boundary =
      grid$lambda[selected] == max(config$lambda),
    usable_candidates = sum(usable),
    converged_candidates = sum(grid$converged %in% TRUE),
    objective_increase_candidates =
      sum(grid$objective_increased %in% TRUE),
    path_elapsed_seconds = fitted$elapsed,
    signal_variance = data$signal_variance,
    noise_sd = data$noise_sd,
    target_variance_snr = data$target_variance_snr,
    empirical_variance_snr = data$empirical_variance_snr,
    runtime_is_scientific = config$runtime_is_scientific,
    stringsAsFactors = FALSE
  )
  coefficient_row <- data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    case_label = task$case_label,
    replication = task$replication,
    method = method,
    stringsAsFactors = FALSE
  )
  coefficient_frame <- as.data.frame(t(coefficient))
  names(coefficient_frame) <- sprintf("b_%03d", seq_along(coefficient))
  list(
    path = path,
    candidate_metrics = candidate_metrics,
    diagnostics = diagnostics,
    coefficient = cbind(coefficient_row, coefficient_frame)
  )
}

select_candidate <- function(candidate_metrics, task) {
  losses <- candidate_metrics$validation_response_mse
  if (length(losses) != length(config$method_order) ||
      any(!is.finite(losses)) ||
      !setequal(candidate_metrics$method, config$method_order)) {
    stop("Validation candidate set is incomplete for task ", task$task_id,
         call. = FALSE)
  }
  minimum <- min(losses)
  tolerance <- config$validation_tie_relative_tolerance *
    (1 + abs(minimum))
  tied <- candidate_metrics$method[losses <= minimum + tolerance]
  selected_method <- config$method_order[
    config$method_order %in% tied
  ][1L]
  selected_index <- match(selected_method, candidate_metrics$method)

  test_minimum <- min(candidate_metrics$response_prediction_mse)
  test_tied <- candidate_metrics$method[
    candidate_metrics$response_prediction_mse <=
      test_minimum + config$validation_tie_relative_tolerance *
      (1 + abs(test_minimum))
  ]
  test_oracle_method <- config$method_order[
    config$method_order %in% test_tied
  ][1L]
  test_oracle_index <- match(test_oracle_method, candidate_metrics$method)

  candidate_metrics$selected_by_validation <-
    candidate_metrics$method == selected_method
  candidate_metrics$test_oracle <-
    candidate_metrics$method == test_oracle_method
  selected <- candidate_metrics[selected_index, , drop = FALSE]
  selected$test_oracle_method <- test_oracle_method
  selected$test_oracle_response_mse <-
    candidate_metrics$response_prediction_mse[test_oracle_index]
  selected$prediction_regret <-
    selected$response_prediction_mse - selected$test_oracle_response_mse
  selected$test_oracle_agreement <- selected_method == test_oracle_method
  selected$selection_used_test <- FALSE

  sorted_validation <- sort(losses)
  selection <- data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    case_label = task$case_label,
    role = task$role,
    replication = task$replication,
    seed_train = task$seed_train,
    seed_validation = task$seed_validation,
    seed_test = task$seed_test,
    selected_method = selected_method,
    selected_penalty_family = selected$penalty_family,
    selected_stage = selected$stage,
    selected_validation_response_mse = losses[selected_index],
    minimum_validation_response_mse = minimum,
    validation_margin_to_second = sorted_validation[2L] - sorted_validation[1L],
    validation_tie_count = length(tied),
    test_oracle_method = test_oracle_method,
    selected_test_response_mse = selected$response_prediction_mse,
    test_oracle_response_mse = selected$test_oracle_response_mse,
    prediction_regret = selected$prediction_regret,
    test_oracle_agreement = selected$test_oracle_agreement,
    selection_used_test = FALSE,
    stringsAsFactors = FALSE
  )
  list(
    candidate_metrics = candidate_metrics,
    selected = selected,
    selection = selection
  )
}

run_task <- function(task_row) {
  task <- as.list(task_row)
  integer_names <- c(
    "task_id", "scenario_index", "n_train", "n_validation", "n_test",
    "replication", "seed_train", "seed_validation", "seed_test"
  )
  for (name in integer_names) task[[name]] <- as.integer(task[[name]])
  data <- generate_task_data(task)
  prepared <- prepare_design(data)

  smooth <- fit_smooth_reference(prepared)
  smooth_coefficient <- smooth$grid$coefficients[, smooth$selected]
  global_gram <- apply(fine_gram, c(1L, 2L), sum)
  smooth_l2 <- sqrt(drop(crossprod(
    smooth_coefficient, global_gram %*% smooth_coefficient
  )))
  cutoff <- config$relative_cutoff_multiplier * smooth_l2 / domain_length
  if (!is.finite(cutoff) || cutoff <= 0) {
    stop("Common positive coefficient cutoff unavailable for task ",
         task$task_id, call. = FALSE)
  }
  adaptive <- adaptive_weights(smooth_coefficient)
  path_kappa <- as.numeric(smooth$grid$kappa[smooth$selected])

  fitted <- lapply(config$method_order, function(method) {
    fit_one_method(
      method, prepared, path_kappa, cutoff, adaptive$weights
    )
  })
  names(fitted) <- config$method_order
  evaluated <- lapply(config$method_order, function(method) {
    evaluate_method(
      fitted[[method]], method, task, data, prepared,
      path_kappa, cutoff
    )
  })
  names(evaluated) <- config$method_order

  candidates <- do.call(rbind, lapply(
    evaluated, `[[`, "candidate_metrics"
  ))
  selected <- select_candidate(candidates, task)

  adaptive_table <- cbind(
    data.frame(
      task_id = task$task_id,
      scenario = task$scenario,
      case_label = task$case_label,
      replication = task$replication,
      epsilon = adaptive$epsilon,
      stringsAsFactors = FALSE
    ),
    region_metadata,
    pilot_norm = adaptive$norms,
    raw_adaptive_weight = adaptive$raw,
    adaptive_weight = adaptive$weights
  )
  smooth_tuning <- data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    case_label = task$case_label,
    replication = task$replication,
    candidate_index = seq_along(smooth$grid$kappa),
    lambda = smooth$grid$lambda,
    kappa = smooth$grid$kappa,
    bic = smooth$grid$bic,
    success = smooth$grid$success,
    converged = smooth$grid$converged,
    refit_success = smooth$grid$refit_success,
    selected = seq_along(smooth$grid$kappa) == smooth$selected,
    stringsAsFactors = FALSE
  )
  list(
    task_id = task$task_id,
    path_tuning = do.call(rbind, lapply(evaluated, `[[`, "path")),
    candidate_metrics = selected$candidate_metrics,
    selected_arfc_metrics = selected$selected,
    selection_outcome = selected$selection,
    diagnostics = do.call(rbind, lapply(evaluated, `[[`, "diagnostics")),
    selected_coefficients = do.call(rbind, lapply(
      evaluated, `[[`, "coefficient"
    )),
    adaptive_weights = adaptive_table,
    smooth_tuning = smooth_tuning
  )
}

results_root <- file.path(simulation_dir, "results", config$study)
run_dir <- file.path(
  results_root, if (smoke) "smoke_run" else "production_run"
)
checkpoint_dir <- file.path(
  results_root,
  if (smoke) "smoke_checkpoints" else "production_checkpoints",
  scientific_signature
)
if (dir.exists(run_dir)) {
  stop("Completed ARFC selection bundle already exists: ", run_dir,
       call. = FALSE)
}
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
existing_checkpoints <- if (dir.exists(checkpoint_dir)) {
  list.files(checkpoint_dir, pattern = "^task_[0-9]+[.]rds$",
             full.names = TRUE)
} else character()
if (length(existing_checkpoints) && !resume) {
  stop("Matching checkpoints exist; use --resume.", call. = FALSE)
}
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path,
                              call. = FALSE)
  temporary <- tempfile(
    basename(path), tmpdir = dirname(path), fileext = ".tmp"
  )
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = "xz", version = 3)
  if (!file.rename(temporary, path)) {
    stop("Atomic checkpoint install failed: ", path, call. = FALSE)
  }
}
checkpoint_path <- function(task_id) {
  file.path(checkpoint_dir, sprintf("task_%04d.rds", as.integer(task_id)))
}

for (index in seq_len(nrow(tasks))) {
  task_id <- tasks$task_id[index]
  path <- checkpoint_path(task_id)
  if (file.exists(path)) {
    checkpoint <- readRDS(path)
    if (!identical(checkpoint$scientific_signature, scientific_signature) ||
        !identical(as.integer(checkpoint$task_id), as.integer(task_id))) {
      stop("Checkpoint signature/task mismatch: ", path, call. = FALSE)
    }
    next
  }
  result <- tryCatch(
    run_task(tasks[index, , drop = FALSE]),
    error = function(condition) condition
  )
  if (inherits(result, "condition")) {
    stop("Task ", task_id, " failed without checkpoint: ",
         conditionMessage(result), ". Correct and use --resume.",
         call. = FALSE)
  }
  atomic_save_rds(list(
    scientific_signature = scientific_signature,
    task_id = as.integer(task_id), result = result
  ), path)
  cat(sprintf(
    "[%s] task %d/%d: %s rep %d\n",
    format(Sys.time(), "%H:%M:%S"), index, nrow(tasks),
    tasks$scenario[index], tasks$replication[index]
  ))
}

checkpoint_files <- vapply(tasks$task_id, checkpoint_path, character(1L))
if (!all(file.exists(checkpoint_files))) {
  stop("ARFC selection run incomplete after task loop.", call. = FALSE)
}
checkpoints <- lapply(checkpoint_files, readRDS)
checkpoint_valid <- vapply(seq_along(checkpoints), function(index) {
  value <- checkpoints[[index]]
  identical(value$scientific_signature, scientific_signature) &&
    identical(as.integer(value$task_id), as.integer(tasks$task_id[index])) &&
    is.list(value$result)
}, logical(1L))
if (!all(checkpoint_valid)) {
  stop("Mixed or invalid checkpoints detected.", call. = FALSE)
}
results <- lapply(checkpoints, `[[`, "result")

path_tuning <- do.call(rbind, lapply(results, `[[`, "path_tuning"))
candidate_metrics <- do.call(rbind, lapply(
  results, `[[`, "candidate_metrics"
))
selected_arfc_metrics <- do.call(rbind, lapply(
  results, `[[`, "selected_arfc_metrics"
))
selection_outcomes <- do.call(rbind, lapply(
  results, `[[`, "selection_outcome"
))
diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
selected_coefficients <- do.call(rbind, lapply(
  results, `[[`, "selected_coefficients"
))
adaptive_weight_table <- do.call(rbind, lapply(
  results, `[[`, "adaptive_weights"
))
smooth_tuning <- do.call(rbind, lapply(results, `[[`, "smooth_tuning"))

expected_path <- nrow(tasks) * length(config$method_order) *
  length(config$lambda)
expected_candidates <- nrow(tasks) * length(config$method_order)
strict_selected <- with(
  candidate_metrics,
  fit_usable & success & converged & refit_success & bic_selected &
    is.finite(bic)
)
if (nrow(path_tuning) != expected_path ||
    nrow(candidate_metrics) != expected_candidates ||
    nrow(selected_arfc_metrics) != nrow(tasks) ||
    nrow(selection_outcomes) != nrow(tasks) ||
    nrow(diagnostics) != expected_candidates ||
    nrow(selected_coefficients) != expected_candidates ||
    nrow(adaptive_weight_table) != nrow(tasks) * nrow(region_metadata) ||
    nrow(smooth_tuning) != nrow(tasks) * length(config$smooth_kappa) ||
    !all(strict_selected) ||
    anyDuplicated(candidate_metrics[c("task_id", "method")]) ||
    anyDuplicated(selection_outcomes$task_id) ||
    !all(selection_outcomes$selection_used_test == FALSE) ||
    any(selection_outcomes$prediction_regret < -1e-12)) {
  stop("Final ARFC selection completeness audit failed.", call. = FALSE)
}

run <- list(
  study = config$study,
  mode = config$mode,
  scientific_signature = scientific_signature,
  config = config,
  scenarios = scenarios,
  tasks = tasks,
  source_md5 = source_md5,
  package_source_md5 = package_source_md5,
  mass_checks = mass_checks,
  method_metadata = method_metadata,
  region_metadata = region_metadata,
  path_tuning = path_tuning,
  candidate_metrics = candidate_metrics,
  selected_arfc_metrics = selected_arfc_metrics,
  selection_outcomes = selection_outcomes,
  diagnostics = diagnostics,
  selected_coefficients = selected_coefficients,
  adaptive_weights = adaptive_weight_table,
  smooth_tuning = smooth_tuning,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

staging <- file.path(
  results_root,
  paste0(".", basename(run_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) {
  stop("Staging directory already exists: ", staging, call. = FALSE)
}
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed && dir.exists(staging)) {
  unlink(staging, recursive = TRUE, force = TRUE)
}, add = TRUE)
prefix <- "arfc_select_v1"
payloads <- c(
  run = paste0(prefix, "_run.rds"),
  tasks = paste0(prefix, "_tasks.csv"),
  methods = paste0(prefix, "_method_metadata.csv"),
  regions = paste0(prefix, "_regions.csv"),
  mass = paste0(prefix, "_mass_checks.csv"),
  path = paste0(prefix, "_path_tuning.csv"),
  candidates = paste0(prefix, "_candidate_metrics.csv"),
  selected = paste0(prefix, "_selected_arfc_metrics.csv"),
  selection = paste0(prefix, "_selection_outcomes.csv"),
  diagnostics = paste0(prefix, "_diagnostics.csv"),
  coefficients = paste0(prefix, "_selected_coefficients.csv"),
  weights = paste0(prefix, "_adaptive_weights.csv"),
  smooth = paste0(prefix, "_smooth_tuning.csv"),
  session = paste0(prefix, "_sessionInfo.txt")
)
saveRDS(run, file.path(staging, payloads["run"]),
        compress = "xz", version = 3)
utils::write.csv(tasks, file.path(staging, payloads["tasks"]),
                 row.names = FALSE)
utils::write.csv(method_metadata, file.path(staging, payloads["methods"]),
                 row.names = FALSE)
utils::write.csv(region_metadata, file.path(staging, payloads["regions"]),
                 row.names = FALSE)
utils::write.csv(mass_checks, file.path(staging, payloads["mass"]),
                 row.names = FALSE)
utils::write.csv(path_tuning, file.path(staging, payloads["path"]),
                 row.names = FALSE)
utils::write.csv(candidate_metrics,
                 file.path(staging, payloads["candidates"]),
                 row.names = FALSE)
utils::write.csv(selected_arfc_metrics,
                 file.path(staging, payloads["selected"]),
                 row.names = FALSE)
utils::write.csv(selection_outcomes,
                 file.path(staging, payloads["selection"]),
                 row.names = FALSE)
utils::write.csv(diagnostics,
                 file.path(staging, payloads["diagnostics"]),
                 row.names = FALSE)
utils::write.csv(selected_coefficients,
                 file.path(staging, payloads["coefficients"]),
                 row.names = FALSE)
utils::write.csv(adaptive_weight_table,
                 file.path(staging, payloads["weights"]),
                 row.names = FALSE)
utils::write.csv(smooth_tuning,
                 file.path(staging, payloads["smooth"]),
                 row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(staging, payloads["session"]), useBytes = TRUE)

payload_paths <- file.path(staging, unname(payloads))
manifest <- data.frame(
  scientific_signature = scientific_signature,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
manifest_file <- file.path(staging, paste0(prefix, "_manifest.csv"))
utils::write.csv(manifest, manifest_file, row.names = FALSE)
if (!file.rename(staging, run_dir)) {
  stop("Atomic final bundle installation failed: ", run_dir,
       call. = FALSE)
}
committed <- TRUE

cat("ARFC four-candidate selection run complete.\n")
cat("Mode:", config$mode, "\n")
cat("Tasks:", nrow(tasks), "\n")
cat("Strict candidate fits:", nrow(candidate_metrics), "\n")
cat("Selections:", nrow(selection_outcomes), "\n")
cat("Scientific signature:", scientific_signature, "\n")
cat("Bundle:", run_dir, "\n")
