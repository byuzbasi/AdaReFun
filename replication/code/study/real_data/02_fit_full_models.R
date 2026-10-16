script_dir <- local({
  file_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_argument)) dirname(normalizePath(sub("^--file=", "", file_argument[1L]))) else normalizePath(getwd())
})
project_root <- normalizePath(file.path(script_dir, ".."))
for (file in c("utils.R", "config.R", "data.R", "methods.R")) {
  source(file.path(script_dir, "R", file))
}
activate_fmcp_libraries()
assert_packages(c(
  "FMCP", "PFLR", "Rcpp", "RcppArmadillo", "fda", "glmnet", "flare"
))

mode <- real_mode()
config <- real_data_config(mode, project_root)
datasets <- load_processed_datasets(config)
dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  grid = resolved_grid_path(config),
  pilot = file.path(config$results_dir,
                    versioned_name("boundary_pilot_selection", mode, "csv")),
  boundary = file.path(config$results_dir,
                       versioned_name("boundary_diagnostics", mode, "csv")),
  full = file.path(config$results_dir,
                   versioned_name("full_models", mode, "rds")),
  selection = file.path(config$results_dir,
                        versioned_name("full_model_selection", mode, "csv")),
  coefficients = file.path(config$results_dir,
                           versioned_name("full_coefficient_curves", mode, "csv")),
  gamma = file.path(config$results_dir,
                    versioned_name("gamma_sensitivity", mode, "csv")),
  session = file.path(config$results_dir,
                      versioned_name("full_models_sessionInfo", mode, "txt"))
)
ensure_new_output(unlist(paths), overwrite = real_overwrite())

initialize_real_method_environment(project_root)

cat("Running one boundary pilot on the approved initial grid...\n")
pilot_fits <- lapply(datasets, function(dataset) {
  fit_real_method_set(
    dataset, seq_along(dataset$y), NULL, config,
    methods = config$primary_sparse_methods
  )
})
pilot_selection <- bind_rows_base(Map(
  function(fit, name) selection_rows(fit, name, "boundary_pilot"),
  pilot_fits, names(datasets)
))
pilot_boundary <- grid_boundary_diagnostics(pilot_selection, config)
resolved <- expand_grid_once(config, pilot_boundary)
resolved$initial_M <- config$M
resolved$initial_lambda <- config$lambda
resolved$initial_kappa <- config$kappa
resolved$pilot_boundary_diagnostics <- pilot_boundary
resolved$resolved_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
atomic_save_rds(resolved, paths$grid)
atomic_write_csv(pilot_selection, paths$pilot)

config$M <- resolved$M
config$lambda <- resolved$lambda
config$kappa <- resolved$kappa
cat(
  "Resolved grid: M=", paste(config$M, collapse = ","),
  "; lambda range=", paste(range(config$lambda), collapse = "--"),
  "; kappa range=", paste(range(config$kappa), collapse = "--"), "\n",
  sep = ""
)

cat("Fitting all full-data methods...\n")
full_fits <- lapply(datasets, function(dataset) {
  fit_real_method_set(
    dataset, seq_along(dataset$y), NULL, config,
    methods = config$methods
  )
})
selection <- bind_rows_base(Map(
  function(fit, name) selection_rows(fit, name, "full"),
  full_fits, names(datasets)
))
coefficients <- bind_rows_base(Map(
  function(fit, dataset) coefficient_rows(fit, dataset, "full"),
  full_fits, datasets
))

cat("Fitting FMCP gamma sensitivity models...\n")
gamma_fits <- list()
gamma_rows <- list()
gamma_index <- 1L
for (dataset_name in names(datasets)) {
  dataset <- datasets[[dataset_name]]
  for (gamma in config$mcp_sensitivity) {
    current <- fit_real_method_set(
      dataset, seq_along(dataset$y), NULL, config,
      methods = "FMCP", mcp_gamma = gamma
    )
    key <- paste(dataset_name, gamma, sep = "::")
    gamma_fits[[key]] <- current
    current_curve <- coefficient_rows(
      current, dataset, paste0("gamma_", gamma)
    )
    current_selection <- selection_rows(
      current, dataset_name, paste0("gamma_", gamma)
    )
    gamma_rows[[gamma_index]] <- merge(
      current_curve, current_selection[
        c("dataset", "context", "M", "lambda", "kappa", "shape", "bic")
      ],
      by = c("dataset", "context"), all.x = TRUE, sort = FALSE
    )
    gamma_index <- gamma_index + 1L
  }
}
gamma_curves <- bind_rows_base(gamma_rows)

post_boundary <- grid_boundary_diagnostics(selection, config)
post_boundary$stage <- "resolved_full_fit"
pilot_boundary$stage <- "initial_pilot"
boundary <- bind_rows_base(list(pilot_boundary, post_boundary))

result <- list(
  version = config$version, mode = mode,
  configuration = config, resolved_grid = resolved,
  datasets = lapply(datasets, function(dataset) {
    dataset[c("name", "subject_id", "argvals_u", "argvals_physical",
              "x_axis_label", "y_axis_label", "source", "notes")]
  }),
  fits = full_fits, gamma_fits = gamma_fits,
  selection = selection, coefficients = coefficients,
  gamma_curves = gamma_curves, boundary_diagnostics = boundary
)
atomic_save_rds(result, paths$full)
atomic_write_csv(selection, paths$selection)
atomic_write_csv(coefficients, paths$coefficients)
atomic_write_csv(gamma_curves, paths$gamma)
atomic_write_csv(boundary, paths$boundary)
capture_session_info(paths$session)

cat("Full-data fits completed. Success table:\n")
print(with(selection, table(dataset, method, success)))

