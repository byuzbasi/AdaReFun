#!/usr/bin/env Rscript

# Independent integrity and semantic validator for ARFC real-data v1.

options(stringsAsFactors = FALSE, warn = 1)

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
script_file <- if (length(file_argument)) {
  normalizePath(sub("^--file=", "", file_argument[1L]), mustWork = TRUE)
} else {
  normalizePath("real_data/29_validate_arfc_real_data_v1.R", mustWork = TRUE)
}
real_dir <- dirname(script_file)
project_dir <- dirname(real_dir)
arguments <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% arguments
core_only <- "--core-only" %in% arguments
report_version_argument <- grep("^--report-version=", arguments, value = TRUE)
if (length(report_version_argument) > 1L) {
  stop("Use at most one --report-version argument.", call. = FALSE)
}
report_version <- if (length(report_version_argument)) {
  sub("^--report-version=", "", report_version_argument)
} else {
  "v1"
}
if (!report_version %in% c("v1", "v2", "v3", "v4")) {
  stop("--report-version must be v1, v2, v3, or v4.", call. = FALSE)
}
unknown <- setdiff(
  arguments, c("--smoke", "--core-only", report_version_argument)
)
if (length(unknown)) {
  stop("Unknown argument(s): ", paste(unknown, collapse = ", "),
       call. = FALSE)
}

helper_file <- file.path(real_dir, "R", "arfc_real_data_v1.R")
source(helper_file, local = FALSE)
config_expected <- arfc_real_config(project_dir, smoke)
mode_prefix <- if (smoke) "smoke" else "production"
prefix <- paste0("arfc_real_data_v1_", mode_prefix)
run_dir <- file.path(
  config_expected$results_root, paste0(mode_prefix, "_run_v1")
)
analysis_dir <- file.path(
  config_expected$results_root, paste0(mode_prefix, "_analysis_v1")
)
report_dir <- file.path(
  config_expected$output_root,
  paste0(mode_prefix, "_reporting_", report_version)
)
validation_dir <- file.path(
  config_expected$results_root,
  paste0(mode_prefix, "_validation_", report_version)
)
if (!dir.exists(run_dir) || !dir.exists(analysis_dir) ||
    (!core_only && !dir.exists(report_dir))) {
  stop("Required run/analysis/reporting bundle is missing.", call. = FALSE)
}
if (dir.exists(validation_dir)) {
  stop("Validation bundle already exists: ", validation_dir,
       call. = FALSE)
}

verify_manifest <- function(directory, manifest_name) {
  path <- file.path(directory, manifest_name)
  if (!file.exists(path)) stop("Missing manifest: ", path, call. = FALSE)
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  files <- file.path(directory, manifest$file)
  passed <- !anyDuplicated(manifest$file) && all(file.exists(files)) &&
    all(as.numeric(file.info(files)$size) == manifest$bytes) &&
    all(arfc_real_md5(files) == manifest$md5) &&
    length(unique(manifest$scientific_signature)) == 1L
  list(manifest = manifest, files = files, passed = passed)
}
run_verified <- verify_manifest(
  run_dir, paste0(prefix, "_manifest.csv")
)
analysis_verified <- verify_manifest(
  analysis_dir, paste0(prefix, "_analysis_manifest.csv")
)
report_verified <- if (!core_only) verify_manifest(
  report_dir, paste0(prefix, "_reporting_manifest.csv")
) else NULL
run <- readRDS(file.path(run_dir, paste0(prefix, "_run.rds")))
analysis <- readRDS(file.path(
  analysis_dir, paste0(prefix, "_analysis.rds")
))
config <- run$config
datasets <- arfc_real_load_datasets(config_expected)
folds_rebuilt <- arfc_real_load_folds(config_expected, datasets)
tasks_rebuilt <- arfc_real_build_tasks(config_expected, folds_rebuilt)
inputs_current <- arfc_real_input_inventory(config_expected)

checks <- list()
add <- function(scope, check, passed, detail, applicable = TRUE) {
  checks[[length(checks) + 1L]] <<- data.frame(
    scope = scope, check = check,
    applicable = isTRUE(applicable),
    passed = if (applicable) isTRUE(passed) else NA,
    detail = as.character(detail), stringsAsFactors = FALSE
  )
}

add("integrity", "run_manifest_bytes_md5", run_verified$passed,
    paste(nrow(run_verified$manifest), "payloads"))
add("integrity", "analysis_manifest_bytes_md5", analysis_verified$passed,
    paste(nrow(analysis_verified$manifest), "payloads"))
add("integrity", "report_manifest_bytes_md5",
    if (!core_only) report_verified$passed else FALSE,
    if (core_only) "core-only" else
      paste(nrow(report_verified$manifest), "payloads"),
    applicable = !core_only)
signature_values <- c(
  unique(run_verified$manifest$scientific_signature),
  unique(analysis_verified$manifest$scientific_signature),
  if (!core_only) unique(report_verified$manifest$scientific_signature) else
    run$scientific_signature,
  run$scientific_signature, analysis$scientific_signature
)
add("integrity", "scientific_signature_chain",
    length(unique(signature_values)) == 1L,
    paste(unique(signature_values), collapse = "|"))
add("integrity", "input_md5_unchanged",
    identical(run$input_inventory$role, inputs_current$role) &&
      identical(run$input_inventory$md5, inputs_current$md5) &&
      identical(run$input_inventory$bytes, inputs_current$bytes),
    paste(inputs_current$md5, collapse = "|"))
add("integrity", "source_mode_and_package",
    identical(run$mode, config_expected$mode) &&
      identical(run$package_versions[["FMCP"]], "0.2.0"),
    paste(run$mode, run$package_versions[["FMCP"]]))

add("data", "dataset_dimensions",
    identical(vapply(datasets, function(value) nrow(value$x), integer(1L)),
              c(Truck = 108L, NeoSpectra = 1974L)) &&
      identical(vapply(datasets, function(value) ncol(value$x), integer(1L)),
                c(Truck = 90L, NeoSpectra = 601L)),
    paste(run$dataset_metadata$dataset, run$dataset_metadata$n,
          sep = "=", collapse = "|"))
add("data", "sample_ids_unique",
    all(vapply(datasets, function(value) !anyDuplicated(value$ids),
               logical(1L))),
    "Truck rows and NeoSpectra soil-sample IDs are unique.")
identity_rebuilt <- arfc_real_bind_rows(lapply(datasets, function(value) {
  data.frame(
    dataset = value$dataset, n = length(value$ids),
    id_md5 = arfc_real_object_md5(as.character(value$ids)),
    source_md5 = value$source_md5,
    processed_md5 = value$processed_md5,
    stringsAsFactors = FALSE
  )
}))
add("data", "dataset_identity_hashes",
    identical(run$dataset_identity, identity_rebuilt),
    paste(identity_rebuilt$dataset, identity_rebuilt$id_md5,
          sep = "=", collapse = "|"))
add("data", "frozen_fold_table_reconstructed",
    identical(run$folds, folds_rebuilt),
    paste(nrow(folds_rebuilt), "fold rows"))
add("data", "outer_tasks_reconstructed",
    identical(run$tasks, tasks_rebuilt) &&
      nrow(run$tasks) == if (smoke) 4L else 110L,
    paste(nrow(run$tasks), "outer tasks"))

add("method", "four_candidate_library",
    identical(run$method_metadata$method, config$method_order),
    paste(run$method_metadata$method, collapse = "|"))
add("method", "equal_scalar_maximum_concavity",
    abs(1 / config$gamma - 1 / (config$scad_a - 1)) <= 1e-14,
    sprintf("MCP %.12g; SCAD %.12g",
            1 / config$gamma, 1 / (config$scad_a - 1)))
add("method", "regional_mass_checks",
    nrow(run$mass_checks) == 8L && all(run$mass_checks$passed),
    paste(sum(run$mass_checks$passed), "/8 PASS"))
add("method", "lambda_path_frozen",
    identical(config$lambda,
              if (smoke) c(1e-9, 1e-5, 1e-1, 1e3) else
                10^seq(-9, 3, by = 0.5)),
    paste(length(config$lambda), "ascending candidates"))

expected_path_rows <- nrow(run$tasks) * (config$inner_folds + 1L) *
  length(config$method_order) * length(config$lambda)
add("tuning", "outer_and_inner_path_complete",
    nrow(run$outer_path) == expected_path_rows,
    paste("expected", expected_path_rows, "observed", nrow(run$outer_path)))
path_groups <- split(
  run$outer_path,
  interaction(run$outer_path$fit_id, run$outer_path$method, drop = TRUE)
)
strict_group <- vapply(path_groups, function(value) {
  selected <- value[value$bic_selected, , drop = FALSE]
  eligible <- value[
    value$fit_usable & value$success & value$converged &
      value$refit_success & is.finite(value$bic), , drop = FALSE
  ]
  nrow(value) == length(config$lambda) && nrow(selected) == 1L &&
    nrow(eligible) > 0L &&
    selected$candidate_index == eligible$candidate_index[
      which.min(eligible$bic)
    ]
}, logical(1L))
add("tuning", "strict_bic_reconstructed", all(strict_group),
    paste(sum(strict_group), "/", length(strict_group), "fit-method paths"))
add("tuning", "outer_candidate_fits_strict",
    nrow(run$outer_candidate_metrics) ==
      nrow(run$tasks) * length(config$method_order) &&
      all(with(run$outer_candidate_metrics,
               success & converged & refit_success & fit_usable)),
    paste(nrow(run$outer_candidate_metrics), "strict fits"))

inner_reconstructed <- arfc_real_bind_rows(lapply(
  split(run$inner_validation, run$inner_validation$task_id),
  function(value) {
    sse <- tapply(value$validation_sse, value$method, sum)
    count <- tapply(value$n_validation, value$method, sum)
    selection <- arfc_real_select_method(sse / count, config)
    data.frame(
      task_id = value$task_id[1L], selected_method = selection$method,
      selected_inner_validation_mse = selection$minimum,
      validation_tie_count = selection$tie_count,
      stringsAsFactors = FALSE
    )
  }
))
inner_reconstructed <- inner_reconstructed[
  match(run$outer_selection$task_id, inner_reconstructed$task_id),
  , drop = FALSE
]
add("selection", "inner_selection_reconstructed",
    identical(inner_reconstructed$selected_method,
              run$outer_selection$selected_method) &&
      max(abs(inner_reconstructed$selected_inner_validation_mse -
                run$outer_selection$selected_inner_validation_mse)) < 1e-12,
    paste(nrow(inner_reconstructed), "task selections"))
add("selection", "outer_test_blind",
    !any(run$outer_selection$outer_test_used_for_selection),
    "No outer-test response entered candidate selection.")
add("selection", "practical_refit_explicit",
    all(run$outer_selection$practical_outer_refit),
    "Selected family/resolution refitted on complete outer training data.")

candidate_methods <- c(config$method_order, config$selection_label)
candidate_prediction <- run$predictions[
  run$predictions$method %in% candidate_methods, , drop = FALSE
]
expected_outer_test_units <- sum(vapply(seq_len(nrow(run$tasks)), function(index) {
  task <- run$tasks[index, , drop = FALSE]
  nrow(run$folds[
    run$folds$dataset == task$dataset &
      run$folds$repeat_index == task$repeat_index &
      run$folds$fold == task$fold, , drop = FALSE
  ])
}, integer(1L)))
expected_prediction_rows <- expected_outer_test_units *
  length(candidate_methods)
add("prediction", "candidate_prediction_rows",
    nrow(candidate_prediction) == expected_prediction_rows &&
      !anyDuplicated(candidate_prediction[c("task_id", "id", "method")]),
    paste(nrow(candidate_prediction), "matched predictions"))
selected_match <- arfc_real_bind_rows(lapply(
  seq_len(nrow(run$outer_selection)), function(index) {
    selected <- run$outer_selection[index, ]
    fixed <- run$predictions[
      run$predictions$task_id == selected$task_id &
        run$predictions$method == selected$selected_method, ]
    chosen <- run$predictions[
      run$predictions$task_id == selected$task_id &
        run$predictions$method == config$selection_label, ]
    fixed <- fixed[order(fixed$id), ]
    chosen <- chosen[order(chosen$id), ]
    data.frame(
      task_id = selected$task_id,
      passed = identical(fixed$id, chosen$id) &&
        max(abs(fixed$predicted - chosen$predicted)) <= 1e-14,
      stringsAsFactors = FALSE
    )
  }
))
add("prediction", "selected_predictions_match_fixed_candidate",
    all(selected_match$passed),
    paste(sum(selected_match$passed), "/", nrow(selected_match), "tasks"))
secondary <- run$predictions[
  run$predictions$method %in% config$secondary_methods, ]
add("prediction", "secondary_benchmarks_frozen_and_matched",
    all(secondary$source %in% c("frozen_v1", "frozen_v3")) &&
      all(secondary$success) &&
      !anyDuplicated(secondary[c("task_id", "id", "method")]),
    paste(nrow(secondary), "Smooth/PCR predictions"))

full_reconstructed <- arfc_real_bind_rows(lapply(
  config$dataset_order, function(dataset_name) {
    value <- run$predictions[
      run$predictions$dataset == dataset_name &
        run$predictions$method %in% config$method_order, ]
    mse <- tapply((value$observed - value$predicted)^2,
                  value$method, mean)
    selection <- arfc_real_select_method(mse, config)
    data.frame(dataset = dataset_name,
               selected_method = selection$method,
               pooled_outer_cv_mse = selection$minimum,
               stringsAsFactors = FALSE)
  }
))
full_reconstructed <- full_reconstructed[
  match(run$full_selection$dataset, full_reconstructed$dataset), ]
add("full_fit", "full_candidate_reconstructed",
    identical(full_reconstructed$selected_method,
              run$full_selection$selected_method) &&
      max(abs(full_reconstructed$pooled_outer_cv_mse -
                run$full_selection$pooled_outer_cv_mse)) <= 1e-14,
    paste(run$full_selection$dataset,
          run$full_selection$selected_method, collapse = "|"))
add("full_fit", "full_curve_grid_complete",
    nrow(run$full_curves) == 2L * length(config$method_order) *
      config$support_grid_points &&
      all(vapply(split(run$full_curves,
                       interaction(run$full_curves$dataset,
                                   run$full_curves$method)),
                 nrow, integer(1L)) == config$support_grid_points),
    paste(nrow(run$full_curves), "curve points"))

add("bootstrap", "full_reselection_complete",
    nrow(run$bootstrap_selection) ==
      2L * config$bootstrap_replicates &&
      all(run$bootstrap_selection$full_reselection) &&
      all(run$bootstrap_selection$n_oob >= config$minimum_oob_units),
    paste(nrow(run$bootstrap_selection), "bootstrap selections"))
bootstrap_loss_groups <- split(
  run$bootstrap_metrics, run$bootstrap_metrics$fit_id
)
bootstrap_reconstructed <- vapply(bootstrap_loss_groups, function(value) {
  losses <- setNames(value$validation_mse, value$method)
  arfc_real_select_method(losses, config)$method
}, character(1L))
bootstrap_expected <- setNames(
  run$bootstrap_selection$selected_method,
  sprintf("%s_bootstrap_%03d", run$bootstrap_selection$dataset,
          run$bootstrap_selection$replication)
)
add("bootstrap", "bootstrap_selection_reconstructed",
    identical(bootstrap_reconstructed[names(bootstrap_expected)],
              bootstrap_expected),
    paste(length(bootstrap_expected), "OOB selections"))
add("bootstrap", "bootstrap_curve_grid_complete",
    nrow(run$bootstrap_curves) ==
      2L * config$bootstrap_replicates * config$support_grid_points,
    paste(nrow(run$bootstrap_curves), "selected curve points"))

metric_recomputed <- arfc_real_bind_rows(lapply(
  split(run$predictions,
        interaction(run$predictions$dataset, run$predictions$method,
                    drop = TRUE)),
  function(value) {
    residual <- value$observed - value$predicted
    tss <- sum((value$observed - mean(value$observed))^2)
    data.frame(
      dataset = value$dataset[1L], method = value$method[1L],
      RMSE = sqrt(mean(residual^2)), MAE = mean(abs(residual)),
      R2 = 1 - sum(residual^2) / tss,
      stringsAsFactors = FALSE
    )
  }
))
metric_matched <- merge(
  metric_recomputed, analysis$prediction_metrics,
  by = c("dataset", "method"), suffixes = c("_new", "_saved")
)
metric_deviation <- max(abs(c(
  metric_matched$RMSE_new - metric_matched$RMSE_saved,
  metric_matched$MAE_new - metric_matched$MAE_saved,
  metric_matched$R2_new - metric_matched$R2_saved
)))
add("analysis", "prediction_metrics_reconstructed",
    nrow(metric_matched) == nrow(analysis$prediction_metrics) &&
      metric_deviation <= 1e-14,
    sprintf("maximum deviation %.3g", metric_deviation))
add("analysis", "analysis_checks_pass",
    all(analysis$validation$passed),
    paste(sum(analysis$validation$passed), "/",
          nrow(analysis$validation), "PASS"))
add("analysis", "no_truth_dependent_real_metrics",
    !any(grepl(
      "ise|true_zero|false_zero|boundary_error|oracle",
      names(analysis), ignore.case = TRUE
    )),
    "Real-data analysis reports prediction and selection stability only.")
add("analysis", "dataset_specific_no_overall",
    all(analysis$prediction_metrics$dataset %in% config$dataset_order) &&
      !any(tolower(analysis$prediction_metrics$dataset) == "overall"),
    "Truck and NeoSpectra remain separate.")

if (!core_only) {
  report_files <- basename(report_verified$files)
  expected_semantics <- c(
    "table_prediction_selection[.]tex$",
    "coefficient_stability[.]pdf$",
    "coefficient_stability[.]png$",
    "figure_include[.]tex$", "results_text[.]tex$",
    "limitations[.]tex$"
  )
  semantic_pass <- all(vapply(expected_semantics, function(pattern) {
    sum(grepl(pattern, report_files)) == 1L
  }, logical(1L)))
  add("reporting", "semantic_artifacts_complete", semantic_pass,
      paste(report_files, collapse = "|"))
  pdf_path <- report_verified$files[
    grepl("coefficient_stability[.]pdf$", report_files)
  ]
  png_path <- report_verified$files[
    grepl("coefficient_stability[.]png$", report_files)
  ]
  pdf_pass <- length(pdf_path) == 1L && file.info(pdf_path)$size > 10000L
  if (pdf_pass && nzchar(Sys.which("pdfinfo"))) {
    info <- suppressWarnings(system2("pdfinfo", pdf_path, stdout = TRUE,
                                     stderr = TRUE))
    pdf_pass <- any(grepl("^Pages:[[:space:]]+1", info))
  }
  add("reporting", "figure_pdf_readable", pdf_pass,
      if (length(pdf_path)) basename(pdf_path) else "missing")
  png_pass <- length(png_path) == 1L && file.info(png_path)$size > 10000L
  if (png_pass && requireNamespace("png", quietly = TRUE)) {
    image <- tryCatch(png::readPNG(png_path), error = function(error) NULL)
    png_pass <- !is.null(image) && length(dim(image)) >= 2L &&
      stats::sd(as.numeric(image)) > 0.01
  }
  add("reporting", "figure_png_readable_nonblank", png_pass,
      if (length(png_path)) basename(png_path) else "missing")
  table_path <- report_verified$files[
    grepl("table_prediction_selection[.]csv$", report_files)
  ]
  table_data <- if (length(table_path) == 1L) {
    utils::read.csv(table_path, check.names = FALSE)
  } else {
    data.frame()
  }
  table_rows <- nrow(table_data)
  add("reporting", "compact_table_rows", table_rows == 14L,
      paste(table_rows, "rows; one table only"))
  table_tex_path <- report_verified$files[
    grepl("table_prediction_selection[.]tex$", report_files)
  ]
  table_tex_text <- if (length(table_tex_path) == 1L) {
    paste(readLines(table_tex_path, warn = FALSE), collapse = "\n")
  } else {
    ""
  }
  v4_table_pass <- if (identical(report_version, "v4")) {
    identical(
      names(table_data),
      c(
        "Dataset", "Method", "RMSE", "MAE", "$R^2$",
        "Nested selection", "Bootstrap selection"
      )
    ) && !grepl("Full-data choice", table_tex_text, fixed = TRUE) &&
      grepl("$R^2$", table_tex_text, fixed = TRUE) &&
      grepl("\\multirow", table_tex_text, fixed = TRUE) &&
      length(unique(table_data$Dataset)) == 2L &&
      identical(as.character(unique(table_data$Dataset)),
                config$dataset_order)
  } else {
    TRUE
  }
  add(
    "reporting", "report_version_table_semantics", v4_table_pass,
    if (identical(report_version, "v4")) {
      "R^2 math header; collapsed dataset rows; no full-data-choice column"
    } else {
      paste("legacy", report_version, "table schema")
    }
  )
  text_paths <- report_verified$files[
    grepl("[.]tex$", report_files)
  ]
  report_text <- paste(unlist(lapply(text_paths, readLines, warn = FALSE)),
                       collapse = "\n")
  smoke_label_pass <- if (smoke) {
    grepl("smoke|noninferential", report_text, ignore.case = TRUE)
  } else {
    !grepl("technical smoke|noninferential", report_text,
           ignore.case = TRUE)
  }
  add("reporting", "mode_label_semantics", smoke_label_pass,
      if (smoke) "Smoke is explicitly noninferential." else
        "Production fragments contain no smoke claim.")
} else {
  for (check in c(
    "semantic_artifacts_complete", "figure_pdf_readable",
    "figure_png_readable_nonblank", "compact_table_rows",
    "report_version_table_semantics", "mode_label_semantics"
  )) add("reporting", check, FALSE, "core-only", applicable = FALSE)
}

validation <- arfc_real_bind_rows(checks)
failed <- validation$applicable & !validation$passed
if (any(failed)) {
  stop("ARFC real-data validation failed: ",
       paste(validation$check[failed], collapse = ", "),
       call. = FALSE)
}

staging <- file.path(
  config$results_root,
  paste0(".", basename(validation_dir), "_staging_", Sys.getpid())
)
if (dir.exists(staging)) stop("Validation staging path exists.",
                              call. = FALSE)
dir.create(staging, recursive = TRUE, showWarnings = FALSE)
committed <- FALSE
on.exit(if (!committed && dir.exists(staging)) {
  unlink(staging, recursive = TRUE, force = TRUE)
}, add = TRUE)
validation_rds <- paste0(prefix, "_validation.rds")
validation_csv <- paste0(prefix, "_validation.csv")
validation_tex <- paste0(prefix, "_validation.tex")
session_file <- paste0(prefix, "_validation_sessionInfo.txt")
saveRDS(list(
  study = run$study, mode = run$mode,
  scientific_signature = run$scientific_signature,
  core_only = core_only, report_version = report_version, checks = validation,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
), file.path(staging, validation_rds), compress = "xz", version = 3)
utils::write.csv(validation, file.path(staging, validation_csv),
                 row.names = FALSE)
latex <- c(
  "\\begin{tabular}{lll}", "\\toprule",
  "Scope & Check & Result \\\\", "\\midrule",
  vapply(seq_len(nrow(validation)), function(index) {
    check <- gsub("_", "\\\\_", validation$check[index], fixed = TRUE)
    scope <- gsub("_", "\\\\_", validation$scope[index], fixed = TRUE)
    result <- if (!validation$applicable[index]) "N/A" else if (
      validation$passed[index]
    ) "PASS" else "FAIL"
    paste(scope, "&", check, "&", result, "\\\\")
  }, character(1L)),
  "\\bottomrule", "\\end{tabular}"
)
writeLines(latex, file.path(staging, validation_tex), useBytes = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(staging, session_file), useBytes = TRUE)
payloads <- c(validation_rds, validation_csv, validation_tex, session_file)
manifest_name <- paste0(prefix, "_validation_manifest.csv")
invisible(arfc_real_manifest(
  staging, run$scientific_signature, payloads, manifest_name
))
if (!file.rename(staging, validation_dir)) {
  stop("Atomic validation installation failed: ", validation_dir,
       call. = FALSE)
}
committed <- TRUE

applicable <- validation$applicable
cat("ARFC-Select real-data validation complete.\n")
cat("Mode:", run$mode, "\n")
cat("Applicable checks:", sum(applicable), "\n")
cat("PASS:", sum(validation$passed[applicable]), "\n")
cat("FAIL:", sum(!validation$passed[applicable]), "\n")
cat("N/A:", sum(!applicable), "\n")
cat("Bundle:", validation_dir, "\n")
