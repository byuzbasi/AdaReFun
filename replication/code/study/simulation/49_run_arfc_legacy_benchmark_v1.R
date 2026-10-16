#!/usr/bin/env Rscript

# Secondary classical-comparator benchmark for the frozen ARFC-Select study.
# The parent four candidates are never refitted.  This runner fits only fMCP,
# fSCAD, Smooth, OLS, PCR, FLiRTI, and the diagnostic Oracle.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("simulation/49_run_arfc_legacy_benchmark_v1.R",
                mustWork = TRUE)
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

fmcp_library <- Sys.getenv(
  "ARFC_FMCP_LIBRARY", unset = ".replication-library"
)
pflr_library <- Sys.getenv(
  "ARFC_PFLR_LIBRARY", unset = ".replication-library"
)
if (!dir.exists(fmcp_library) || !dir.exists(pflr_library)) {
  stop(
    "Frozen FMCP/PFLR libraries are unavailable. Expected: ",
    fmcp_library, " and ", pflr_library, call. = FALSE
  )
}
.libPaths(c(normalizePath(fmcp_library), normalizePath(pflr_library),
            .libPaths()))
required_packages <- c(
  "FMCP", "Rcpp", "RcppArmadillo", "PFLR", "fda", "glmnet", "flare",
  "psych"
)
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, quietly = TRUE,
  FUN.VALUE = logical(1L)
)]
if (length(missing_packages)) {
  stop("Missing frozen package(s): ", paste(missing_packages, collapse = ", "),
       call. = FALSE)
}
if (utils::packageVersion("FMCP") != numeric_version("0.2.0") ||
    utils::packageVersion("PFLR") != numeric_version("1.1.0")) {
  stop("Benchmark requires frozen FMCP 0.2.0 and PFLR 1.1.0.",
       call. = FALSE)
}

mode <- if (smoke) "smoke" else "production"
prefix <- "arfc_legacy_benchmark_v1"
protocol_file <- file.path(
  project_dir, "config", "protocol_arfc_legacy_benchmark_v1.md"
)
parent_protocol_file <- file.path(
  project_dir, "config", "protocol_arfc_select_v1.md"
)
adapter_file <- file.path(
  simulation_dir, "comparators", "fscad_and_pflr_adapters.R"
)
classic_file <- file.path(
  simulation_dir, "comparators", "classic_methods.R"
)
parent_root <- file.path(simulation_dir, "results", "arfc_select_v1")
parent_run_dir <- file.path(parent_root, paste0(mode, "_run"))
parent_run_file <- file.path(parent_run_dir, "arfc_select_v1_run.rds")
parent_manifest_file <- file.path(
  parent_run_dir, "arfc_select_v1_manifest.csv"
)
required_sources <- c(
  protocol_file, parent_protocol_file, adapter_file, classic_file,
  script_file, parent_run_file, parent_manifest_file
)
if (!all(file.exists(required_sources))) {
  stop("One or more frozen benchmark inputs are missing.", call. = FALSE)
}

verify_manifest <- function(directory, manifest_file) {
  manifest <- utils::read.csv(manifest_file, check.names = FALSE)
  required <- c("scientific_signature", "file", "bytes", "md5")
  if (!all(required %in% names(manifest)) || anyDuplicated(manifest$file)) {
    stop("Parent manifest schema is invalid.", call. = FALSE)
  }
  paths <- file.path(directory, manifest$file)
  if (!all(file.exists(paths)) ||
      !all(as.numeric(file.info(paths)$size) == manifest$bytes) ||
      !all(unname(as.character(tools::md5sum(paths))) == manifest$md5)) {
    stop("Parent run manifest validation failed.", call. = FALSE)
  }
  manifest
}

parent_manifest <- verify_manifest(parent_run_dir, parent_manifest_file)
parent <- readRDS(parent_run_file)
if (!identical(parent$study, "arfc_select_v1") ||
    !identical(parent$mode, mode) ||
    length(unique(parent_manifest$scientific_signature)) != 1L ||
    !identical(unique(parent_manifest$scientific_signature),
               parent$scientific_signature)) {
  stop("Parent identity, mode, or signature mismatch.", call. = FALSE)
}

source(adapter_file, local = FALSE)
source(classic_file, local = FALSE)

parent_config <- parent$config
tasks <- parent$tasks
scenarios <- parent$scenarios
rownames(tasks) <- NULL

config <- list(
  study = prefix,
  protocol_id = prefix,
  mode = mode,
  secondary_benchmark = TRUE,
  parent_study = parent$study,
  parent_signature = parent$scientific_signature,
  observed_grid_points = parent_config$observed_grid_points,
  support_grid_points = parent_config$support_grid_points,
  predictor_order = parent_config$predictor_order,
  predictor_knots = parent_config$predictor_knots,
  intercept = parent_config$intercept,
  snr = parent_config$snr,
  case_i_noise_sd = parent_config$case_i_noise_sd,
  n_train = parent_config$n_train,
  n_validation = parent_config$n_validation,
  n_test = parent_config$n_test,
  replications = parent_config$replications,
  M = parent_config$finest_M,
  degree = parent_config$degree,
  penalty_order = parent_config$penalty_order,
  gamma = parent_config$gamma,
  scad_a = parent_config$scad_a,
  lambda = parent_config$lambda,
  smooth_kappa = parent_config$smooth_kappa,
  max_iter = parent_config$max_iter,
  tolerance = parent_config$tolerance,
  norm_floor = parent_config$norm_floor,
  relative_cutoff_multiplier = parent_config$relative_cutoff_multiplier,
  zero_tolerance = parent_config$zero_tolerance,
  refit = parent_config$refit,
  integration = parent_config$integration,
  pcr_max_components = 30L,
  flirti_M = 40L,
  flirti_tuning = list(
    Mf = c(10L, 20L, 30L),
    lambda = exp(seq(log(5e-4), log(5e-2), length.out = 7L)),
    cons = 4L
  ),
  new_method_order = c(
    "fixed_fmcp", "fixed_fscad", "smooth", "ols", "pcr", "flirti",
    "oracle"
  ),
  reporting_order = c(
    "arfc_select", "asr_mcp", "amr_mcp", "asr_scad", "amr_scad",
    "fixed_fmcp", "fixed_fscad", "smooth", "ols", "pcr", "flirti",
    "oracle"
  ),
  method_labels = c(
    arfc_select = "ARFC--Select",
    asr_mcp = "ASR--MCP", amr_mcp = "AMR--MCP",
    asr_scad = "ASR--SCAD", amr_scad = "AMR--SCAD",
    fixed_fmcp = "fMCP", fixed_fscad = "fSCAD",
    smooth = "Smooth", ols = "OLS", pcr = "PCR",
    flirti = "FLiRTI", oracle = "Oracle"
  ),
  exact_zero_methods = c("fixed_fmcp", "fixed_fscad", "flirti", "oracle"),
  strict_bic_tie_rule =
    "first exact minimum among strict candidates in frozen grid order",
  test_used_for_tuning = FALSE,
  validation_used_for_tuning = FALSE,
  runtime_is_scientific = FALSE,
  package_versions = vapply(
    required_packages,
    function(package) as.character(utils::packageVersion(package)),
    FUN.VALUE = character(1L)
  )
)

expected_tasks <- 4L * config$replications
required_task_columns <- c(
  "task_id", "scenario_index", "scenario", "case_label", "role",
  "null_type", "null_left", "null_right", "n_train", "n_validation",
  "n_test", "replication", "seed_train", "seed_validation", "seed_test"
)
seed_values <- unlist(
  tasks[c("seed_train", "seed_validation", "seed_test")], use.names = FALSE
)
if (!all(required_task_columns %in% names(tasks)) ||
    nrow(tasks) != expected_tasks || anyDuplicated(tasks$task_id) ||
    anyDuplicated(seed_values) || !identical(config$M, 32L) ||
    !identical(sort(config$lambda), config$lambda) ||
    abs(1 / config$gamma - 1 / (config$scad_a - 1)) > 1e-14 ||
    !identical(config$new_method_order,
               c("fixed_fmcp", "fixed_fscad", "smooth", "ols", "pcr",
                 "flirti", "oracle")) ||
    !identical(config$reporting_order,
               c("arfc_select", "asr_mcp", "amr_mcp", "asr_scad",
                 "amr_scad", "fixed_fmcp", "fixed_fscad", "smooth", "ols",
                 "pcr", "flirti", "oracle"))) {
  stop("Frozen benchmark design failed its static audit.", call. = FALSE)
}

md5_character <- function(paths) {
  unname(as.character(tools::md5sum(paths)))
}
object_md5 <- function(object) {
  temporary <- tempfile("arfc_legacy_signature_", fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  unname(as.character(tools::md5sum(temporary)))
}
relative_names <- function(paths) {
  sub(paste0("^", project_dir, "/?"), "", normalizePath(paths))
}
source_md5 <- setNames(md5_character(required_sources),
                       relative_names(required_sources))
fmcp_source_files <- c(
  file.path(project_dir, "FMCP-package", "DESCRIPTION"),
  file.path(project_dir, "FMCP-package", "src", "fmcp_core.cpp"),
  file.path(project_dir, "FMCP-package", "src", "fscad_core.cpp"),
  file.path(project_dir, "FMCP-package", "R", "basis.R")
)
pflr_source_files <- c(
  file.path(project_dir, "PFLR-master", "DESCRIPTION"),
  list.files(file.path(project_dir, "PFLR-master", "R"),
             pattern = "[.]R$", full.names = TRUE)
)
package_source_files <- c(fmcp_source_files, pflr_source_files)
if (!all(file.exists(package_source_files))) {
  stop("Package provenance sources are missing.", call. = FALSE)
}
package_source_md5 <- setNames(
  md5_character(package_source_files), relative_names(package_source_files)
)
parent_manifest_md5 <- unname(as.character(tools::md5sum(parent_manifest_file)))
scientific_signature <- object_md5(list(
  config = config, scenarios = scenarios, tasks = tasks,
  source_md5 = source_md5, package_source_md5 = package_source_md5,
  parent_manifest_md5 = parent_manifest_md5
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
  argvals, config$M, config$degree, config$penalty_order
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
  config$M,
  coefficient_basis$order
)
roughness <- FMCP:::cpp_weighted_crossprod(
  coefficient_basis$derivative_quad, coefficient_basis$quad_weights
)
interval_scale <- sqrt(config$M / domain_length)
global_gram <- FMCP:::cpp_weighted_crossprod(
  coefficient_basis$basis_quad, coefficient_basis$quad_weights
)

method_metadata <- data.frame(
  method = config$new_method_order,
  label = unname(config$method_labels[config$new_method_order]),
  family = c("mcp", "scad", "smooth", "ols", "pcr", "flirti", "oracle"),
  fixed_resolution = c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE),
  adaptive = FALSE,
  multi_resolution = FALSE,
  structural_sparsity = config$new_method_order %in% config$exact_zero_methods,
  truth_aware = config$new_method_order == "oracle",
  selector_eligible = FALSE,
  stringsAsFactors = FALSE
)

if (validate_design_only) {
  cat("Study:", config$study, "\n")
  cat("Mode:", config$mode, "\n")
  cat("Parent signature:", parent$scientific_signature, "\n")
  cat("Tasks:", nrow(tasks), "\n")
  cat("New methods:", paste(config$new_method_order, collapse = ", "), "\n")
  cat("Reporting methods:", paste(config$reporting_order, collapse = ", "),
      "\n")
  cat("Lambda candidates:", length(config$lambda), "\n")
  cat("Equal maximum concavity:", 1 / config$gamma, "\n")
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
  stop("Unknown scenario: ", scenario, call. = FALSE)
}

true_null_mask <- function(scenario, values) {
  if (identical(scenario, "case_i_null")) return(rep(TRUE, length(values)))
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
  noise_sd <- if (case_i) config$case_i_noise_sd else {
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
  design <- FMCP:::cpp_design_matrix(
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
    design = design,
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

grid_tuning <- function(grid, method, selected, task) {
  data.frame(
    task_id = task$task_id,
    scenario = task$scenario,
    replication = task$replication,
    method = method,
    candidate_index = seq_along(grid$lambda),
    M = config$M,
    lambda = as.numeric(grid$lambda),
    kappa = as.numeric(grid$kappa),
    bic = as.numeric(grid$bic),
    rss = as.numeric(grid$rss),
    df = as.numeric(grid$df),
    success = as.logical(grid$success),
    converged = as.logical(grid$converged),
    refit_success = as.logical(grid$refit_success),
    selected = seq_along(grid$lambda) == selected,
    message = as.character(grid$message),
    stringsAsFactors = FALSE
  )
}

fit_grid_method <- function(method, prepared, lambda, kappa, cutoff) {
  start <- proc.time()[["elapsed"]]
  grid <- if (identical(method, "fixed_fmcp") ||
              identical(method, "smooth")) {
    FMCP:::cpp_fmcp_grid(
      prepared$centered_y, prepared$design, roughness, local_gram,
      lambda, kappa, config$gamma, interval_scale,
      config$max_iter, config$tolerance, config$norm_floor,
      cutoff, config$refit, FALSE
    )
  } else if (identical(method, "fixed_fscad")) {
    FMCP:::cpp_fscad_grid(
      prepared$centered_y, prepared$design, roughness, local_gram,
      lambda, kappa, config$scad_a, interval_scale,
      config$max_iter, config$tolerance, config$norm_floor,
      cutoff, config$refit, FALSE
    )
  } else {
    stop("Unsupported grid method: ", method, call. = FALSE)
  }
  selected <- strict_candidate(grid)
  if (!is.finite(selected)) {
    stop("No strict candidate for ", method, call. = FALSE)
  }
  coefficient <- as.numeric(grid$coefficients[, selected])
  list(
    method = method,
    success = TRUE,
    coefficient = coefficient,
    beta_observed = as.numeric(
      coefficient_basis$basis_observed %*% coefficient
    ),
    beta_dense = as.numeric(dense_basis %*% coefficient),
    validation_prediction = as.numeric(
      prepared$means$y + prepared$validation_design %*% coefficient
    ),
    test_prediction = as.numeric(
      prepared$means$y + prepared$test_design %*% coefficient
    ),
    selected = list(
      M = config$M,
      lambda = as.numeric(grid$lambda[selected]),
      kappa = as.numeric(grid$kappa[selected]),
      bic = as.numeric(grid$bic[selected]),
      df = as.numeric(grid$df[selected]),
      active = as.integer(grid$active[selected]),
      converged = as.logical(grid$converged[selected]),
      iterations = as.integer(grid$iterations[selected])
    ),
    tuning_grid = grid,
    selected_index = selected,
    elapsed = proc.time()[["elapsed"]] - start,
    warnings = "", error = ""
  )
}

fit_ols <- function(prepared) {
  start <- proc.time()[["elapsed"]]
  decomposition <- svd(prepared$design)
  threshold <- max(dim(prepared$design)) * max(decomposition$d) *
    .Machine$double.eps
  keep <- decomposition$d > threshold
  if (!any(keep)) stop("OLS design has zero numerical rank.", call. = FALSE)
  coefficient <- as.numeric(
    decomposition$v[, keep, drop = FALSE] %*%
      ((crossprod(decomposition$u[, keep, drop = FALSE],
                  prepared$centered_y)) / decomposition$d[keep])
  )
  residual <- prepared$centered_y - prepared$design %*% coefficient
  rss <- sum(residual^2)
  df <- sum(keep)
  bic <- length(prepared$centered_y) *
    log(max(rss / length(prepared$centered_y), .Machine$double.xmin)) +
    log(length(prepared$centered_y)) * df
  list(
    method = "ols", success = TRUE, coefficient = coefficient,
    beta_observed = as.numeric(coefficient_basis$basis_observed %*% coefficient),
    beta_dense = as.numeric(dense_basis %*% coefficient),
    validation_prediction = as.numeric(
      prepared$means$y + prepared$validation_design %*% coefficient
    ),
    test_prediction = as.numeric(
      prepared$means$y + prepared$test_design %*% coefficient
    ),
    selected = list(
      M = config$M, lambda = 0, kappa = 0, bic = bic, df = df,
      active = df, converged = TRUE, iterations = 1L
    ),
    tuning = data.frame(
      M = config$M, lambda = 0, kappa = 0, bic = bic, rss = rss,
      df = df, success = TRUE, converged = TRUE,
      refit_success = TRUE, selected = TRUE, message = "",
      stringsAsFactors = FALSE
    ),
    elapsed = proc.time()[["elapsed"]] - start,
    warnings = "", error = ""
  )
}

interpolate_dense <- function(beta_observed) {
  as.numeric(stats::approx(
    argvals, beta_observed, xout = support_grid, rule = 2,
    ties = "ordered"
  )$y)
}

classic_result <- function(selection, method, data, applicable = TRUE) {
  if (!applicable) {
    return(list(
      method = method, applicable = FALSE, success = NA,
      beta_observed = rep(NA_real_, length(argvals)),
      beta_dense = rep(NA_real_, length(support_grid)),
      validation_prediction = rep(NA_real_, data$validation$x |> nrow()),
      test_prediction = rep(NA_real_, data$test$x |> nrow()),
      selected = list(
        M = NA_integer_, lambda = NA_real_, kappa = NA_real_,
        bic = NA_real_, df = NA_real_, active = NA_integer_,
        converged = NA, iterations = NA_integer_
      ),
      tuning = NULL, elapsed = selection$elapsed,
      warnings = selection$warnings,
      error = if (nzchar(selection$error)) selection$error else "Not applicable"
    ))
  }
  if (is.null(selection$fit) || is.null(selection$selected) ||
      any(!is.finite(selection$beta))) {
    return(list(
      method = method, applicable = TRUE, success = FALSE,
      beta_observed = rep(NA_real_, length(argvals)),
      beta_dense = rep(NA_real_, length(support_grid)),
      validation_prediction = rep(NA_real_, nrow(data$validation$x)),
      test_prediction = rep(NA_real_, nrow(data$test$x)),
      selected = list(
        M = NA_integer_, lambda = NA_real_, kappa = NA_real_,
        bic = NA_real_, df = NA_real_, active = NA_integer_,
        converged = FALSE, iterations = NA_integer_
      ),
      tuning = selection$diagnostics, elapsed = selection$elapsed,
      warnings = selection$warnings, error = selection$error
    ))
  }
  validation_prediction <- tryCatch(
    as.numeric(selection$predict(data$validation$x)),
    error = function(condition) rep(NA_real_, nrow(data$validation$x))
  )
  test_prediction <- tryCatch(
    as.numeric(selection$predict(data$test$x)),
    error = function(condition) rep(NA_real_, nrow(data$test$x))
  )
  row <- selection$selected[1L, , drop = FALSE]
  value <- function(name, default = NA_real_) {
    if (name %in% names(row)) row[[name]][1L] else default
  }
  list(
    method = method, applicable = TRUE,
    success = all(is.finite(c(validation_prediction, test_prediction))),
    beta_observed = as.numeric(selection$beta),
    beta_dense = interpolate_dense(selection$beta),
    validation_prediction = validation_prediction,
    test_prediction = test_prediction,
    selected = list(
      M = as.integer(value("M", NA_integer_)),
      lambda = as.numeric(value("lambda")),
      kappa = as.numeric(value("kappa")),
      bic = as.numeric(value("bic")),
      df = as.numeric(value("df")),
      active = as.integer(value("active", NA_integer_)),
      converged = as.logical(value("converged", TRUE)),
      iterations = as.integer(value("iterations", 1L))
    ),
    tuning = selection$diagnostics,
    elapsed = selection$elapsed,
    warnings = selection$warnings, error = selection$error
  )
}

normalize_tuning <- function(result, task) {
  if (!is.null(result$tuning_grid)) {
    return(grid_tuning(
      result$tuning_grid, result$method, result$selected_index, task
    ))
  }
  tuning <- result$tuning
  if (is.null(tuning) || !is.data.frame(tuning) || !nrow(tuning)) {
    tuning <- data.frame(
      M = result$selected$M, lambda = result$selected$lambda,
      kappa = result$selected$kappa, bic = result$selected$bic,
      rss = NA_real_, df = result$selected$df,
      success = isTRUE(result$success),
      converged = result$selected$converged,
      refit_success = isTRUE(result$success), selected = TRUE,
      message = result$error, stringsAsFactors = FALSE
    )
  }
  if (identical(result$method, "pcr") && "components" %in% names(tuning)) {
    tuning$M <- tuning$components
    tuning$lambda <- NA_real_
    tuning$kappa <- NA_real_
    tuning$df <- tuning$components
    tuning$success <- is.finite(tuning$bic)
    tuning$converged <- tuning$success
    tuning$refit_success <- tuning$success
    tuning$selected <- seq_len(nrow(tuning)) == which.min(tuning$bic)
    tuning$message <- ""
  }
  if (identical(result$method, "flirti")) {
    source_tuning <- result$selection_fit_tuning
    if (!is.null(source_tuning)) tuning <- source_tuning
    if ("Mf" %in% names(tuning)) tuning$M <- tuning$Mf
  }
  defaults <- list(
    M = NA_integer_, lambda = NA_real_, kappa = NA_real_, bic = NA_real_,
    rss = NA_real_, df = NA_real_, success = isTRUE(result$success),
    converged = isTRUE(result$success), refit_success = isTRUE(result$success),
    selected = FALSE, message = ""
  )
  for (name in names(defaults)) {
    if (!name %in% names(tuning)) tuning[[name]] <- defaults[[name]]
  }
  if (!any(tuning$selected %in% TRUE) && isTRUE(result$success)) {
    if (any(is.finite(tuning$bic))) {
      tuning$selected[which.min(tuning$bic)] <- TRUE
    } else {
      tuning$selected[1L] <- TRUE
    }
  }
  data.frame(
    task_id = task$task_id, scenario = task$scenario,
    replication = task$replication, method = result$method,
    candidate_index = seq_len(nrow(tuning)),
    M = as.integer(tuning$M), lambda = as.numeric(tuning$lambda),
    kappa = as.numeric(tuning$kappa), bic = as.numeric(tuning$bic),
    rss = as.numeric(tuning$rss), df = as.numeric(tuning$df),
    success = as.logical(tuning$success),
    converged = as.logical(tuning$converged),
    refit_success = as.logical(tuning$refit_success),
    selected = as.logical(tuning$selected),
    message = as.character(tuning$message), stringsAsFactors = FALSE
  )
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

support_metrics <- function(beta, truth, task, structural_sparsity) {
  if (length(beta) != length(support_grid) || any(!is.finite(beta))) {
    return(setNames(rep(NA_real_, length(support_metric_names)),
                    support_metric_names))
  }
  truth_null <- true_null_mask(task$scenario, support_grid)
  truth_nonnull <- !truth_null & abs(truth) > 1e-10
  squared_error <- (beta - truth)^2
  basic <- c(
    ise_total = weighted_integral(squared_error),
    ise_null = weighted_integral(squared_error, truth_null),
    ise_nonnull = weighted_integral(squared_error, truth_nonnull)
  )
  if (!isTRUE(structural_sparsity)) {
    return(c(
      basic,
      setNames(rep(NA_real_, length(support_metric_names) - length(basic)),
               support_metric_names[-seq_along(basic)])
    ))
  }
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
    basic,
    true_zero_rate = if (null_length > 0) intersection / null_length else NA,
    false_zero_rate = if (nonnull_length > 0) false_length / nonnull_length else NA,
    false_active_proportion = if (identical(task$scenario, "case_i_null")) {
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

evaluate_result <- function(result, task, data) {
  metadata <- method_metadata[method_metadata$method == result$method, ]
  structural <- if (isTRUE(result$applicable %||% TRUE) &&
                    isTRUE(result$success)) {
    support_metrics(
      result$beta_dense, data$truth_dense, task,
      metadata$structural_sparsity
    )
  } else {
    setNames(rep(NA_real_, length(support_metric_names)), support_metric_names)
  }
  applicable <- if (is.null(result$applicable)) TRUE else result$applicable
  valid_predictions <- isTRUE(result$success) &&
    all(is.finite(result$validation_prediction)) &&
    all(is.finite(result$test_prediction))
  selected <- result$selected
  metric <- data.frame(
    task_id = task$task_id, scenario = task$scenario,
    case_label = task$case_label, role = task$role,
    replication = task$replication,
    seed_train = task$seed_train,
    seed_validation = task$seed_validation,
    seed_test = task$seed_test,
    method = result$method, label = metadata$label,
    applicable = applicable,
    fit_success = if (applicable) isTRUE(result$success) else NA,
    structural_sparsity = metadata$structural_sparsity,
    truth_aware = metadata$truth_aware,
    M = selected$M, lambda = selected$lambda, kappa = selected$kappa,
    shape = if (identical(result$method, "fixed_fmcp")) {
      config$gamma
    } else if (identical(result$method, "fixed_fscad")) {
      config$scad_a
    } else NA_real_,
    bic = selected$bic, df = selected$df, active = selected$active,
    converged = selected$converged, iterations = selected$iterations,
    elapsed_seconds = result$elapsed,
    warnings = result$warnings, error = result$error,
    stringsAsFactors = FALSE
  )
  metric <- cbind(metric, as.data.frame(as.list(structural)))
  metric$validation_signal_mse <- if (valid_predictions) {
    mean((result$validation_prediction - data$validation$signal)^2)
  } else NA_real_
  metric$validation_response_mse <- if (valid_predictions) {
    mean((result$validation_prediction - data$validation$y)^2)
  } else NA_real_
  metric$signal_prediction_mse <- if (valid_predictions) {
    mean((result$test_prediction - data$test$signal)^2)
  } else NA_real_
  metric$response_prediction_mse <- if (valid_predictions) {
    mean((result$test_prediction - data$test$y)^2)
  } else NA_real_
  metric$test_used_for_tuning <- FALSE
  metric$validation_used_for_tuning <- FALSE
  metric$runtime_is_scientific <- FALSE
  metric
}

`%||%` <- function(left, right) if (is.null(left)) right else left

wide_vector <- function(result, task, values, prefix_name) {
  identity <- data.frame(
    task_id = task$task_id, scenario = task$scenario,
    case_label = task$case_label, replication = task$replication,
    method = result$method, stringsAsFactors = FALSE
  )
  payload <- as.data.frame(t(as.numeric(values)))
  names(payload) <- sprintf("%s_%04d", prefix_name, seq_along(values))
  cbind(identity, payload)
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

  parent_diagnostic <- parent$diagnostics[
    parent$diagnostics$task_id == task$task_id, , drop = FALSE
  ]
  path_kappa_parent <- unique(parent_diagnostic$path_kappa)
  cutoff_parent <- unique(parent_diagnostic$cutoff)
  if (length(path_kappa_parent) != 1L || length(cutoff_parent) != 1L ||
      !is.finite(path_kappa_parent) || !is.finite(cutoff_parent) ||
      cutoff_parent <= 0) {
    stop("Invalid parent tuning for task ", task$task_id, call. = FALSE)
  }

  smooth <- fit_grid_method(
    "smooth", prepared, lambda = 0, kappa = config$smooth_kappa, cutoff = 0
  )
  smooth_l2 <- sqrt(drop(crossprod(
    smooth$coefficient, global_gram %*% smooth$coefficient
  )))
  cutoff_reconstructed <-
    config$relative_cutoff_multiplier * smooth_l2 / domain_length
  provenance_values <- c(
    data$signal_variance - parent_diagnostic$signal_variance[1L],
    data$noise_sd - parent_diagnostic$noise_sd[1L],
    smooth$selected$kappa - path_kappa_parent,
    cutoff_reconstructed - cutoff_parent
  )
  if (is.finite(data$empirical_variance_snr)) {
    provenance_values <- c(
      provenance_values,
      data$empirical_variance_snr -
        parent_diagnostic$empirical_variance_snr[1L]
    )
  }
  provenance_difference <- max(abs(provenance_values))
  if (!is.finite(provenance_difference) || provenance_difference > 1e-10) {
    stop("Task data/tuning did not reproduce parent task ", task$task_id,
         call. = FALSE)
  }

  fixed_fmcp <- fit_grid_method(
    "fixed_fmcp", prepared, config$lambda, path_kappa_parent, cutoff_parent
  )
  fixed_fscad <- fit_grid_method(
    "fixed_fscad", prepared, config$lambda, path_kappa_parent, cutoff_parent
  )
  ols <- fit_ols(prepared)

  pcr_selection <- fit_pcr_bic(
    data$train$y, data$train$x, argvals, integration_weights,
    max_components = config$pcr_max_components
  )
  pcr <- classic_result(pcr_selection, "pcr", data)

  flirti_selection <- wrap_pflr_comparator(
    "FLiRTI", data$train$y, data$train$x, argvals,
    M = config$flirti_M, tuning = config$flirti_tuning
  )
  flirti <- classic_result(flirti_selection, "flirti", data)
  if (!is.null(flirti_selection$fit$extra$tuning)) {
    flirti$tuning <- flirti_selection$fit$extra$tuning
  }

  oracle_applicable <- task$scenario %in%
    c("case_i_null", "case_ii_sparse")
  oracle_selection <- fit_oracle_bic(
    data$train$y, data$train$x, argvals, integration_weights,
    scenario = task$scenario, M_grid = config$M,
    kappa_grid = config$smooth_kappa, degree = config$degree
  )
  oracle <- classic_result(
    oracle_selection, "oracle", data, applicable = oracle_applicable
  )

  results <- list(
    fixed_fmcp = fixed_fmcp, fixed_fscad = fixed_fscad,
    smooth = smooth, ols = ols, pcr = pcr, flirti = flirti,
    oracle = oracle
  )
  results$smooth$applicable <- TRUE
  results$ols$applicable <- TRUE
  results$fixed_fmcp$applicable <- TRUE
  results$fixed_fscad$applicable <- TRUE

  metrics <- do.call(rbind, lapply(
    config$new_method_order,
    function(method) evaluate_result(results[[method]], task, data)
  ))
  tuning <- do.call(rbind, lapply(
    config$new_method_order,
    function(method) normalize_tuning(results[[method]], task)
  ))
  curves <- do.call(rbind, lapply(
    config$new_method_order,
    function(method) wide_vector(
      results[[method]], task, results[[method]]$beta_observed, "beta"
    )
  ))
  validation_predictions <- do.call(rbind, lapply(
    config$new_method_order,
    function(method) wide_vector(
      results[[method]], task, results[[method]]$validation_prediction,
      "validation"
    )
  ))
  test_predictions <- do.call(rbind, lapply(
    config$new_method_order,
    function(method) wide_vector(
      results[[method]], task, results[[method]]$test_prediction, "test"
    )
  ))
  diagnostics <- data.frame(
    task_id = task$task_id, scenario = task$scenario,
    case_label = task$case_label, replication = task$replication,
    seed_train = task$seed_train,
    path_kappa = path_kappa_parent, cutoff = cutoff_parent,
    reconstructed_cutoff = cutoff_reconstructed,
    data_provenance_max_difference = provenance_difference,
    signal_variance = data$signal_variance, noise_sd = data$noise_sd,
    target_variance_snr = data$target_variance_snr,
    empirical_variance_snr = data$empirical_variance_snr,
    applicable_fits = sum(metrics$applicable %in% TRUE),
    successful_applicable_fits = sum(
      metrics$applicable %in% TRUE & metrics$fit_success %in% TRUE
    ),
    stringsAsFactors = FALSE
  )
  list(
    task_id = task$task_id, metrics = metrics, tuning = tuning,
    curves = curves, validation_predictions = validation_predictions,
    test_predictions = test_predictions, diagnostics = diagnostics
  )
}

results_root <- file.path(simulation_dir, "results", prefix)
run_dir <- file.path(results_root, paste0(mode, "_run_v1"))
checkpoint_dir <- file.path(
  results_root, paste0(mode, "_checkpoints_v1"), scientific_signature
)
if (dir.exists(run_dir)) {
  stop("Completed benchmark bundle already exists: ", run_dir,
       call. = FALSE)
}
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
existing_checkpoints <- if (dir.exists(checkpoint_dir)) {
  list.files(checkpoint_dir, pattern = "^task_[0-9]+[.]rds$", full.names = TRUE)
} else character()
if (length(existing_checkpoints) && !resume) {
  stop("Matching checkpoints exist; use --resume.", call. = FALSE)
}
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

atomic_save_rds <- function(object, path) {
  if (file.exists(path)) stop("Refusing to overwrite: ", path, call. = FALSE)
  temporary <- tempfile(basename(path), tmpdir = dirname(path), fileext = ".tmp")
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
  result <- run_task(tasks[index, , drop = FALSE])
  applicable <- result$metrics$applicable %in% TRUE
  if (!all(result$metrics$fit_success[applicable] %in% TRUE)) {
    failed <- result$metrics$method[
      applicable & !(result$metrics$fit_success %in% TRUE)
    ]
    stop("Applicable method failure in task ", task_id, ": ",
         paste(failed, collapse = ", "), call. = FALSE)
  }
  atomic_save_rds(list(
    scientific_signature = scientific_signature,
    task_id = as.integer(task_id), result = result
  ), path)
  if (index %% 10L == 0L || index == nrow(tasks)) {
    cat(sprintf("Completed %d/%d tasks\n", index, nrow(tasks)))
  }
}

checkpoint_files <- vapply(tasks$task_id, checkpoint_path, character(1L))
if (!all(file.exists(checkpoint_files))) {
  stop("Benchmark is incomplete after task loop.", call. = FALSE)
}
checkpoints <- lapply(checkpoint_files, readRDS)
if (!all(vapply(checkpoints, function(value) {
  identical(value$scientific_signature, scientific_signature)
}, logical(1L)))) {
  stop("Mixed scientific signatures in checkpoints.", call. = FALSE)
}
task_results <- lapply(checkpoints, `[[`, "result")
metrics <- do.call(rbind, lapply(task_results, `[[`, "metrics"))
tuning <- do.call(rbind, lapply(task_results, `[[`, "tuning"))
curve_estimates <- do.call(rbind, lapply(task_results, `[[`, "curves"))
validation_predictions <- do.call(rbind, lapply(
  task_results, `[[`, "validation_predictions"
))
test_predictions <- do.call(rbind, lapply(
  task_results, `[[`, "test_predictions"
))
diagnostics <- do.call(rbind, lapply(task_results, `[[`, "diagnostics"))

expected_rows <- nrow(tasks) * length(config$new_method_order)
applicable_expected <- nrow(tasks) * 6L +
  sum(tasks$scenario %in% c("case_i_null", "case_ii_sparse"))
applicable <- metrics$applicable %in% TRUE
if (nrow(metrics) != expected_rows || nrow(curve_estimates) != expected_rows ||
    nrow(validation_predictions) != expected_rows ||
    nrow(test_predictions) != expected_rows || nrow(diagnostics) != nrow(tasks) ||
    sum(applicable) != applicable_expected ||
    !all(metrics$fit_success[applicable] %in% TRUE) ||
    anyDuplicated(metrics[c("task_id", "method")]) ||
    anyDuplicated(curve_estimates[c("task_id", "method")]) ||
    anyDuplicated(validation_predictions[c("task_id", "method")]) ||
    anyDuplicated(test_predictions[c("task_id", "method")]) ||
    max(diagnostics$data_provenance_max_difference) > 1e-10 ||
    any(grepl("overall", c(metrics$scenario, diagnostics$scenario),
              ignore.case = TRUE))) {
  stop("Final benchmark completeness audit failed.", call. = FALSE)
}

run <- list(
  study = config$study, mode = config$mode,
  scientific_signature = scientific_signature,
  config = config, scenarios = scenarios, tasks = tasks,
  source_md5 = source_md5, package_source_md5 = package_source_md5,
  parent_signature = parent$scientific_signature,
  parent_manifest_md5 = parent_manifest_md5,
  method_metadata = method_metadata,
  metrics = metrics, tuning = tuning,
  curve_estimates = curve_estimates,
  validation_predictions = validation_predictions,
  test_predictions = test_predictions,
  diagnostics = diagnostics,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

staging <- file.path(
  results_root, paste0(".", basename(run_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) {
  stop("Staging directory already exists: ", staging, call. = FALSE)
}
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE), add = TRUE)
payloads <- c(
  run = paste0(prefix, "_run_v1.rds"),
  tasks = paste0(prefix, "_tasks_v1.csv"),
  methods = paste0(prefix, "_method_metadata_v1.csv"),
  metrics = paste0(prefix, "_metrics_v1.csv"),
  tuning = paste0(prefix, "_tuning_v1.csv"),
  diagnostics = paste0(prefix, "_diagnostics_v1.csv"),
  session = paste0(prefix, "_sessionInfo_v1.txt")
)
saveRDS(run, file.path(staging, payloads["run"]),
        compress = "xz", version = 3)
utils::write.csv(tasks, file.path(staging, payloads["tasks"]),
                 row.names = FALSE)
utils::write.csv(method_metadata, file.path(staging, payloads["methods"]),
                 row.names = FALSE)
utils::write.csv(metrics, file.path(staging, payloads["metrics"]),
                 row.names = FALSE)
utils::write.csv(tuning, file.path(staging, payloads["tuning"]),
                 row.names = FALSE)
utils::write.csv(diagnostics, file.path(staging, payloads["diagnostics"]),
                 row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(staging, payloads["session"]))
payload_paths <- file.path(staging, unname(payloads))
manifest <- data.frame(
  scientific_signature = scientific_signature,
  file = basename(payload_paths),
  bytes = as.numeric(file.info(payload_paths)$size),
  md5 = unname(as.character(tools::md5sum(payload_paths))),
  stringsAsFactors = FALSE
)
manifest_file <- file.path(
  staging, paste0(prefix, "_manifest_v1.csv")
)
utils::write.csv(manifest, manifest_file, row.names = FALSE)
if (!file.rename(staging, run_dir)) {
  stop("Atomic benchmark bundle install failed: ", run_dir, call. = FALSE)
}
cat("ARFC classical-comparator benchmark complete.\n")
cat("Tasks:", nrow(tasks), "\n")
cat("Applicable strict fits:", sum(applicable), "\n")
cat("Scientific signature:", scientific_signature, "\n")
cat("Bundle:", run_dir, "\n")
