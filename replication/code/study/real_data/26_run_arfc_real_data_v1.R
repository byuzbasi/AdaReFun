#!/usr/bin/env Rscript

# Versioned ARFC-Select real-data runner.
#
# The four-candidate library is selected inside the frozen outer folds.  Long
# production execution is intended for the user's local terminal; The local validation workflow runs
# only --validate-design and --smoke.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("real_data/26_run_arfc_real_data_v1.R", mustWork = TRUE)
}
real_dir <- dirname(script_file)
project_dir <- dirname(real_dir)

arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
resume <- "--resume" %in% arguments
validate_design <- "--validate-design" %in% arguments
unknown <- setdiff(arguments, c("--smoke", "--resume", "--validate-design"))
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}
if (resume && validate_design) {
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
  stop("ARFC real-data v1 requires frozen FMCP 0.2.0.", call. = FALSE)
}

helper_file <- file.path(real_dir, "R", "arfc_real_data_v1.R")
protocol_file <- file.path(
  project_dir, "config", "protocol_arfc_real_data_v1.md"
)
framework_protocol <- file.path(
  project_dir, "config",
  "protocol_adaptive_regional_folded_concave_v1.md"
)
core_file <- file.path(
  project_dir, "simulation", "src", "adaptive_regional_core.cpp"
)
downstream_sources <- file.path(real_dir, c(
  "27_analyze_arfc_real_data_v1.R",
  "28_write_arfc_real_data_v1.R",
  "29_validate_arfc_real_data_v1.R"
))
required_sources <- c(
  helper_file, protocol_file, framework_protocol, core_file, script_file,
  downstream_sources
)
if (!all(file.exists(required_sources))) {
  stop("ARFC real-data source inventory is incomplete.", call. = FALSE)
}
source(helper_file, local = FALSE)
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

config <- arfc_real_config(project_dir, smoke)
datasets <- arfc_real_load_datasets(config)
input_inventory <- arfc_real_input_inventory(config)
folds <- arfc_real_load_folds(config, datasets)
tasks <- arfc_real_build_tasks(config, folds)
secondary_predictions <- arfc_real_load_secondary_predictions(config, tasks)
method_metadata <- arfc_real_method_metadata(config)
contexts <- lapply(datasets, arfc_real_build_context, config = config)
mass_checks <- arfc_real_bind_rows(lapply(contexts, `[[`, "mass_checks"))

expected_tasks <- if (smoke) 4L else 110L
if (nrow(tasks) != expected_tasks || anyDuplicated(tasks$task_id) ||
    !all(mass_checks$passed) ||
    abs(1 / config$gamma - 1 / (config$scad_a - 1)) > 1e-14 ||
    !identical(config$levels, c(8L, 16L, 32L)) ||
    abs(sum(config$level_weights) - 1) > 1e-14 ||
    (!smoke && length(config$lambda) != 25L) ||
    !identical(config$method_order,
               c("asr_mcp", "asr_scad", "amr_mcp", "amr_scad"))) {
  stop("Frozen ARFC real-data design failed its static audit.",
       call. = FALSE)
}

source_md5 <- setNames(
  arfc_real_md5(required_sources), basename(required_sources)
)
package_source_files <- c(
  file.path(project_dir, "FMCP-package", "src", "fmcp_core.cpp"),
  file.path(project_dir, "FMCP-package", "src", "fscad_core.cpp"),
  file.path(project_dir, "FMCP-package", "R", "basis.R")
)
if (!all(file.exists(package_source_files))) {
  stop("FMCP package provenance sources are missing.", call. = FALSE)
}
package_source_md5 <- setNames(
  arfc_real_md5(package_source_files), basename(package_source_files)
)
package_versions <- vapply(
  required_packages,
  function(package) as.character(utils::packageVersion(package)),
  FUN.VALUE = character(1L)
)
dataset_identity <- arfc_real_bind_rows(lapply(datasets, function(value) {
  data.frame(
    dataset = value$dataset, n = length(value$ids),
    id_md5 = arfc_real_object_md5(as.character(value$ids)),
    source_md5 = value$source_md5,
    processed_md5 = value$processed_md5,
    stringsAsFactors = FALSE
  )
}))
signature_config <- config
signature_config$workers <- NULL
scientific_signature <- arfc_real_object_md5(list(
  config = signature_config, tasks = tasks,
  method_metadata = method_metadata,
  dataset_identity = dataset_identity,
  input_inventory = input_inventory[c("role", "bytes", "md5")],
  source_md5 = source_md5,
  package_source_md5 = package_source_md5,
  package_versions = package_versions
))

dataset_metadata <- arfc_real_bind_rows(lapply(datasets, function(value) {
  data.frame(
    dataset = value$dataset, n = length(value$y),
    grid_points = ncol(value$x),
    physical_min = min(value$argvals_physical),
    physical_max = max(value$argvals_physical),
    response = value$response_description,
    sampling_unit = value$sampling_unit,
    id_md5 = arfc_real_object_md5(as.character(value$ids)),
    source_md5 = value$source_md5,
    processed_md5 = value$processed_md5,
    limitation = value$limitation,
    stringsAsFactors = FALSE
  )
}))

if (validate_design) {
  cat("Study:", config$study, "\n")
  cat("Mode:", config$mode, "\n")
  cat("Datasets:", paste(dataset_metadata$dataset, collapse = ", "), "\n")
  cat("Outer tasks:", nrow(tasks), "\n")
  cat("Inner folds:", config$inner_folds, "\n")
  cat("Candidate methods:", paste(config$method_order, collapse = ", "),
      "\n")
  cat("Lambda candidates:", length(config$lambda), "\n")
  cat("Bootstrap replications per dataset:", config$bootstrap_replicates,
      "\n")
  cat("Mass checks:", sum(mass_checks$passed), "/", nrow(mass_checks),
      "PASS\n")
  cat("Scientific signature:", scientific_signature, "\n")
  quit(save = "no", status = 0L)
}

mode_prefix <- if (smoke) "smoke" else "production"
run_dir <- file.path(config$results_root, paste0(mode_prefix, "_run_v1"))
checkpoint_root <- file.path(
  config$results_root, paste0(mode_prefix, "_checkpoints_v1"),
  scientific_signature
)
outer_checkpoint_dir <- file.path(checkpoint_root, "outer")
bootstrap_checkpoint_dir <- file.path(checkpoint_root, "bootstrap")
if (dir.exists(run_dir)) {
  stop("Completed ARFC real-data bundle already exists: ", run_dir,
       call. = FALSE)
}
dir.create(outer_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bootstrap_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
existing <- c(
  list.files(outer_checkpoint_dir, pattern = "[.]rds$", full.names = TRUE),
  list.files(bootstrap_checkpoint_dir, pattern = "[.]rds$", full.names = TRUE)
)
if (length(existing) && !resume) {
  stop("Matching checkpoints exist; use --resume.", call. = FALSE)
}

outer_checkpoint_path <- function(task_id) {
  file.path(outer_checkpoint_dir, paste0(task_id, ".rds"))
}
bootstrap_checkpoint_path <- function(dataset, replication) {
  file.path(
    bootstrap_checkpoint_dir,
    sprintf("%s_bootstrap_%03d.rds", dataset, replication)
  )
}

run_outer_index <- function(index) {
  task <- tasks[index, , drop = FALSE]
  path <- outer_checkpoint_path(task$task_id)
  if (file.exists(path)) {
    checkpoint <- readRDS(path)
    if (!identical(checkpoint$scientific_signature, scientific_signature) ||
        !identical(checkpoint$task_id, task$task_id)) {
      stop("Outer checkpoint signature/task mismatch: ", path,
           call. = FALSE)
    }
    return(path)
  }
  result <- arfc_real_run_outer_task(
    task, datasets[[task$dataset]], folds,
    contexts[[task$dataset]], config
  )
  arfc_real_atomic_save_rds(list(
    scientific_signature = scientific_signature,
    task_id = task$task_id, result = result
  ), path)
  path
}

outer_indices <- seq_len(nrow(tasks))
if (config$workers > 1L && .Platform$OS.type != "windows") {
  outer_paths <- unlist(parallel::mclapply(
    outer_indices, run_outer_index, mc.cores = config$workers,
    mc.preschedule = FALSE
  ), use.names = FALSE)
} else {
  outer_paths <- vapply(outer_indices, function(index) {
    path <- run_outer_index(index)
    cat(sprintf(
      "[%s] outer task %d/%d: %s\n",
      format(Sys.time(), "%H:%M:%S"), index, nrow(tasks),
      tasks$task_id[index]
    ))
    path
  }, character(1L))
}

outer_checkpoints <- lapply(outer_paths, readRDS)
outer_valid <- vapply(seq_along(outer_checkpoints), function(index) {
  value <- outer_checkpoints[[index]]
  identical(value$scientific_signature, scientific_signature) &&
    identical(value$task_id, tasks$task_id[index]) &&
    is.list(value$result)
}, logical(1L))
if (!all(outer_valid)) {
  stop("Mixed or invalid outer checkpoints detected.", call. = FALSE)
}
outer_results <- lapply(outer_checkpoints, `[[`, "result")
candidate_predictions <- arfc_real_bind_rows(lapply(
  outer_results, `[[`, "predictions"
))
all_predictions <- rbind(candidate_predictions, secondary_predictions)

candidate_test_keys <- unique(candidate_predictions[c("task_id", "id")])
for (method in config$secondary_methods) {
  benchmark_keys <- unique(secondary_predictions[
    secondary_predictions$method == method, c("task_id", "id")
  ])
  if (nrow(benchmark_keys) != nrow(candidate_test_keys) ||
      !setequal(paste(benchmark_keys$task_id, benchmark_keys$id),
                paste(candidate_test_keys$task_id, candidate_test_keys$id))) {
    stop("Secondary benchmark outer rows do not match for ", method,
         call. = FALSE)
  }
}

full_fits <- lapply(config$dataset_order, function(dataset_name) {
  arfc_real_run_full_fit(
    datasets[[dataset_name]], contexts[[dataset_name]],
    candidate_predictions, config
  )
})
names(full_fits) <- config$dataset_order

bootstrap_grid <- expand.grid(
  dataset = config$dataset_order,
  replication = seq_len(config$bootstrap_replicates),
  stringsAsFactors = FALSE
)
bootstrap_grid <- bootstrap_grid[
  order(match(bootstrap_grid$dataset, config$dataset_order),
        bootstrap_grid$replication), , drop = FALSE
]
rownames(bootstrap_grid) <- NULL

run_bootstrap_index <- function(index) {
  current <- bootstrap_grid[index, , drop = FALSE]
  path <- bootstrap_checkpoint_path(current$dataset, current$replication)
  key <- sprintf("%s_%03d", current$dataset, current$replication)
  if (file.exists(path)) {
    checkpoint <- readRDS(path)
    if (!identical(checkpoint$scientific_signature, scientific_signature) ||
        !identical(checkpoint$key, key)) {
      stop("Bootstrap checkpoint signature/key mismatch: ", path,
           call. = FALSE)
    }
    return(path)
  }
  result <- arfc_real_run_bootstrap(
    datasets[[current$dataset]], contexts[[current$dataset]],
    current$replication, config
  )
  arfc_real_atomic_save_rds(list(
    scientific_signature = scientific_signature,
    key = key, result = result
  ), path)
  path
}

bootstrap_indices <- seq_len(nrow(bootstrap_grid))
if (config$workers > 1L && .Platform$OS.type != "windows") {
  bootstrap_paths <- unlist(parallel::mclapply(
    bootstrap_indices, run_bootstrap_index,
    mc.cores = config$workers, mc.preschedule = FALSE
  ), use.names = FALSE)
} else {
  bootstrap_paths <- vapply(bootstrap_indices, function(index) {
    path <- run_bootstrap_index(index)
    cat(sprintf(
      "[%s] bootstrap %d/%d: %s %d\n",
      format(Sys.time(), "%H:%M:%S"), index, nrow(bootstrap_grid),
      bootstrap_grid$dataset[index], bootstrap_grid$replication[index]
    ))
    path
  }, character(1L))
}

bootstrap_checkpoints <- lapply(bootstrap_paths, readRDS)
bootstrap_valid <- vapply(seq_along(bootstrap_checkpoints), function(index) {
  expected <- sprintf(
    "%s_%03d", bootstrap_grid$dataset[index],
    bootstrap_grid$replication[index]
  )
  value <- bootstrap_checkpoints[[index]]
  identical(value$scientific_signature, scientific_signature) &&
    identical(value$key, expected) && is.list(value$result)
}, logical(1L))
if (!all(bootstrap_valid)) {
  stop("Mixed or invalid bootstrap checkpoints detected.", call. = FALSE)
}
bootstrap_results <- lapply(bootstrap_checkpoints, `[[`, "result")

inner_validation <- arfc_real_bind_rows(lapply(
  outer_results, `[[`, "inner_validation"
))
outer_path <- arfc_real_bind_rows(lapply(outer_results, `[[`, "path"))
outer_diagnostics <- arfc_real_bind_rows(lapply(
  outer_results, `[[`, "diagnostics"
))
outer_selection <- arfc_real_bind_rows(lapply(
  outer_results, `[[`, "selection"
))
outer_candidate_metrics <- arfc_real_bind_rows(lapply(
  outer_results, `[[`, "outer_candidate_metrics"
))
full_selection <- arfc_real_bind_rows(lapply(full_fits, `[[`, "selection"))
full_curves <- arfc_real_bind_rows(lapply(full_fits, `[[`, "curves"))
full_metrics <- arfc_real_bind_rows(lapply(full_fits, `[[`, "metrics"))
full_path <- arfc_real_bind_rows(lapply(full_fits, `[[`, "path"))
full_diagnostics <- arfc_real_bind_rows(lapply(
  full_fits, `[[`, "diagnostics"
))
bootstrap_selection <- arfc_real_bind_rows(lapply(
  bootstrap_results, `[[`, "selection"
))
bootstrap_curves <- arfc_real_bind_rows(lapply(
  bootstrap_results, `[[`, "curves"
))
bootstrap_metrics <- arfc_real_bind_rows(lapply(
  bootstrap_results, `[[`, "metrics"
))
bootstrap_path <- arfc_real_bind_rows(lapply(
  bootstrap_results, `[[`, "path"
))
bootstrap_diagnostics <- arfc_real_bind_rows(lapply(
  bootstrap_results, `[[`, "diagnostics"
))

expected_outer_test_units <- sum(vapply(seq_len(nrow(tasks)), function(index) {
  task <- tasks[index, , drop = FALSE]
  current <- folds[
    folds$dataset == task$dataset &
      folds$repeat_index == task$repeat_index &
      folds$fold == task$fold, , drop = FALSE
  ]
  nrow(current)
}, integer(1L)))
expected_candidate_predictions <- expected_outer_test_units *
  (length(config$method_order) + 1L)
if (nrow(candidate_predictions) != expected_candidate_predictions ||
    nrow(outer_selection) != nrow(tasks) ||
    nrow(inner_validation) !=
      nrow(tasks) * config$inner_folds * length(config$method_order) ||
    nrow(outer_candidate_metrics) !=
      nrow(tasks) * length(config$method_order) ||
    nrow(full_selection) != length(config$dataset_order) ||
    nrow(full_curves) != length(config$dataset_order) *
      length(config$method_order) * config$support_grid_points ||
    nrow(bootstrap_selection) != nrow(bootstrap_grid) ||
    nrow(bootstrap_curves) != nrow(bootstrap_grid) *
      config$support_grid_points ||
    anyDuplicated(outer_selection$task_id) ||
    any(outer_selection$outer_test_used_for_selection) ||
    !all(bootstrap_selection$full_reselection) ||
    any(bootstrap_selection$n_oob < config$minimum_oob_units)) {
  stop("Final ARFC real-data completeness audit failed.", call. = FALSE)
}

run <- list(
  study = config$study, mode = config$mode,
  scientific_signature = scientific_signature,
  config = config, package_versions = package_versions,
  source_md5 = source_md5,
  package_source_md5 = package_source_md5,
  input_inventory = input_inventory,
  dataset_identity = dataset_identity,
  dataset_metadata = dataset_metadata,
  tasks = tasks, folds = folds,
  method_metadata = method_metadata,
  mass_checks = mass_checks,
  inner_validation = inner_validation,
  outer_path = outer_path,
  outer_diagnostics = outer_diagnostics,
  outer_selection = outer_selection,
  outer_candidate_metrics = outer_candidate_metrics,
  predictions = all_predictions,
  full_selection = full_selection,
  full_curves = full_curves,
  full_metrics = full_metrics,
  full_path = full_path,
  full_diagnostics = full_diagnostics,
  bootstrap_selection = bootstrap_selection,
  bootstrap_curves = bootstrap_curves,
  bootstrap_metrics = bootstrap_metrics,
  bootstrap_path = bootstrap_path,
  bootstrap_diagnostics = bootstrap_diagnostics,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

staging <- file.path(
  config$results_root,
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
prefix <- paste0("arfc_real_data_v1_", mode_prefix)
payloads <- c(
  run = paste0(prefix, "_run.rds"),
  inputs = paste0(prefix, "_input_inventory.csv"),
  datasets = paste0(prefix, "_dataset_metadata.csv"),
  tasks = paste0(prefix, "_tasks.csv"),
  folds = paste0(prefix, "_folds.csv"),
  methods = paste0(prefix, "_method_metadata.csv"),
  mass = paste0(prefix, "_mass_checks.csv"),
  inner = paste0(prefix, "_inner_validation.csv"),
  outer_path = paste0(prefix, "_outer_path.csv"),
  outer_selection = paste0(prefix, "_outer_selection.csv"),
  predictions = paste0(prefix, "_predictions.csv"),
  full_selection = paste0(prefix, "_full_selection.csv"),
  full_curves = paste0(prefix, "_full_curves.csv"),
  bootstrap_selection = paste0(prefix, "_bootstrap_selection.csv"),
  bootstrap_curves = paste0(prefix, "_bootstrap_curves.csv"),
  diagnostics = paste0(prefix, "_diagnostics.csv"),
  session = paste0(prefix, "_sessionInfo.txt")
)
saveRDS(run, file.path(staging, payloads["run"]),
        compress = "xz", version = 3)
write_csv <- function(object, role) utils::write.csv(
  object, file.path(staging, payloads[role]), row.names = FALSE
)
write_csv(input_inventory, "inputs")
write_csv(dataset_metadata, "datasets")
write_csv(tasks, "tasks")
write_csv(folds, "folds")
write_csv(method_metadata, "methods")
write_csv(mass_checks, "mass")
write_csv(inner_validation, "inner")
write_csv(outer_path, "outer_path")
write_csv(outer_selection, "outer_selection")
write_csv(all_predictions, "predictions")
write_csv(full_selection, "full_selection")
write_csv(full_curves, "full_curves")
write_csv(bootstrap_selection, "bootstrap_selection")
write_csv(bootstrap_curves, "bootstrap_curves")
write_csv(arfc_real_bind_rows(list(
  outer_diagnostics, full_diagnostics, bootstrap_diagnostics
)), "diagnostics")
writeLines(capture.output(sessionInfo()),
           file.path(staging, payloads["session"]), useBytes = TRUE)
manifest_name <- paste0(prefix, "_manifest.csv")
invisible(arfc_real_manifest(
  staging, scientific_signature, unname(payloads), manifest_name
))
if (!file.rename(staging, run_dir)) {
  stop("Atomic final bundle installation failed: ", run_dir,
       call. = FALSE)
}
committed <- TRUE

cat("ARFC-Select real-data runner complete.\n")
cat("Mode:", config$mode, "\n")
cat("Outer tasks:", nrow(tasks), "\n")
cat("Bootstrap tasks:", nrow(bootstrap_grid), "\n")
cat("Strict outer candidate fits:", nrow(outer_candidate_metrics), "\n")
cat("Scientific signature:", scientific_signature, "\n")
cat("Bundle:", run_dir, "\n")
