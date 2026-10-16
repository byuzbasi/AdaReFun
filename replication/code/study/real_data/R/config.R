real_data_config <- function(mode = c("full", "smoke"), project_root = NULL) {
  mode <- match.arg(mode)
  if (is.null(project_root)) project_root <- find_fmcp_project_root()
  real_dir <- file.path(project_root, "real_data")
  config <- list(
    version = "real_data_v1",
    mode = mode,
    master_seed = 20260819L,
    datasets = c("Truck", "Tecator"),
    methods = c(
      "FMCP", "fSCAD", "Smooth", "PCR",
      "PenS", "NGR", "FLiRTI", "OLS"
    ),
    primary_sparse_methods = c("FMCP", "fSCAD"),
    M = c(10L, 20L, 30L, 40L, 50L),
    degree = 3L,
    penalty_order = 2L,
    lambda = sort(unique(c(1e-5, 10^seq(-4, 0, length.out = 9L), 10))),
    kappa = sort(unique(c(1e-11, 10^seq(-10, -1, length.out = 7L)))),
    mcp_gamma = 3,
    mcp_sensitivity = c(2.7, 3.7),
    scad_a = 3.7,
    max_iter = 120L,
    tolerance = 1e-8,
    norm_floor = 1e-8,
    relative_cutoff_multiplier = 1e-4,
    zero_tolerance = sqrt(.Machine$double.eps),
    integration = "trapezoid",
    pcr_max_components = 30L,
    pflr_M = 40L,
    pflr_tuning = list(
      PenS = list(alpha = 10^seq(-10, -1, length.out = 7L)),
      NGR = list(
        alphaPS = 10^seq(-9, -3, length.out = 4L),
        kappa = 10^seq(-8, -3, length.out = 4L),
        tau = exp(seq(-36, -24, length.out = 5L)),
        gamma = 0.5,
        niter = 40L
      ),
      FLiRTI = list(
        Mf = c(10L, 20L, 30L),
        lambda = exp(seq(log(5e-4), log(5e-2), length.out = 7L)),
        cons = 4L
      )
    ),
    boundary_rate_threshold = 0.10,
    maximum_grid_expansions = 1L,
    cv_repeats = if (identical(mode, "smoke")) 1L else 20L,
    cv_folds = if (identical(mode, "smoke")) 2L else 5L,
    subject_bootstrap_replicates = if (identical(mode, "smoke")) 20L else 2000L,
    stability_replicates = if (identical(mode, "smoke")) 4L else 500L,
    stability_thresholds = c(0.7, 0.8, 0.9),
    primary_stability_threshold = 0.8,
    permutation_replicates = if (identical(mode, "smoke")) 5L else 999L,
    workers = max(1L, as.integer(Sys.getenv(
      "FMCP_REAL_WORKERS",
      unset = if (.Platform$OS.type == "windows") "1" else "4"
    ))),
    processed_dir = file.path(real_dir, "data", "processed", "v1"),
    results_dir = file.path(real_dir, "results", "v1"),
    figures_dir = file.path(real_dir, "figures", "v1"),
    outputs_dir = file.path(
      project_root, "outputs", "01a01660-203e-7871-8510-13dcd8b1533f",
      "real_data_v1"
    )
  )
  config$package_versions <- vapply(
    c("FMCP", "PFLR", "Rcpp", "RcppArmadillo", "fda", "glmnet", "flare"),
    function(package) {
      if (requireNamespace(package, quietly = TRUE)) {
        as.character(utils::packageVersion(package))
      } else {
        NA_character_
      }
    },
    FUN.VALUE = character(1L)
  )
  config
}

resolved_grid_path <- function(config) {
  file.path(
    config$results_dir,
    versioned_name("resolved_tuning_grid", config$mode, "rds")
  )
}

apply_resolved_grid <- function(config, required = TRUE) {
  path <- resolved_grid_path(config)
  if (!file.exists(path)) {
    if (required) {
      stop(
        "Resolved tuning grid is missing. Run 02_fit_full_models.R first.",
        call. = FALSE
      )
    }
    return(config)
  }
  resolved <- readRDS(path)
  config$M <- resolved$M
  config$lambda <- resolved$lambda
  config$kappa <- resolved$kappa
  config$grid_resolution <- resolved
  config
}

task_seed <- function(config, component, dataset_index = 0L,
                      index = 0L, subindex = 0L) {
  offsets <- c(data = 100000L, cv = 200000L, stability = 300000L,
               permutation = 400000L, analysis = 500000L)
  if (!component %in% names(offsets)) {
    stop("Unknown seed component: ", component, call. = FALSE)
  }
  as.integer(config$master_seed + offsets[[component]] +
               dataset_index * 10000L + index * 10L + subindex)
}

