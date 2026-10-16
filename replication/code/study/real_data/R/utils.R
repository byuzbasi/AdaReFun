find_fmcp_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "FMCP-package", "DESCRIPTION")) &&
        file.exists(file.path(current, "PFLR-master", "data", "truck.rda"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the FMCP project root.", call. = FALSE)
    }
    current <- parent
  }
}

real_data_script_directory <- function() {
  file_argument <- grep(
    "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
  )
  if (length(file_argument)) {
    return(dirname(normalizePath(sub("^--file=", "", file_argument[1L]))))
  }
  normalizePath(getwd())
}

activate_fmcp_libraries <- function() {
  temporary_library <- Sys.getenv(
    "FMCP_REAL_LIBRARY", unset = ".replication-library"
  )
  if (dir.exists(temporary_library)) {
    .libPaths(unique(c(temporary_library, .libPaths())))
  }
  invisible(.libPaths())
}

assert_packages <- function(packages) {
  missing <- packages[!vapply(
    packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1L)
  )]
  if (length(missing)) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      ". Use the project library documented in real_data/README.md.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

real_mode <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  if ("--smoke" %in% arguments) "smoke" else "full"
}

real_resume <- function(arguments = commandArgs(trailingOnly = TRUE)) {
  "--resume" %in% arguments
}

real_overwrite <- function() {
  identical(Sys.getenv("FMCP_REAL_OVERWRITE"), "1")
}

versioned_name <- function(stem, mode = "full", extension = "rds") {
  suffix <- if (identical(mode, "smoke")) "_smoke" else ""
  paste0(stem, "_v1", suffix, ".", extension)
}

ensure_new_output <- function(paths, resume = FALSE, overwrite = FALSE) {
  existing <- paths[file.exists(paths)]
  if (length(existing) && !resume && !overwrite) {
    stop(
      "Versioned outputs already exist: ",
      paste(basename(existing), collapse = ", "),
      ". Use --resume for checkpoints; replacement requires ",
      "FMCP_REAL_OVERWRITE=1.", call. = FALSE
    )
  }
  invisible(TRUE)
}

atomic_save_rds <- function(object, path, compress = "xz") {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "-"),
                        tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, compress = compress)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "-"),
                        tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = row.names, na = "")
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_lines <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = paste0(basename(path), "-"),
                        tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

capture_session_info <- function(path) {
  atomic_write_lines(capture.output(utils::sessionInfo()), path)
}

trapezoid_weights <- function(argvals) {
  if (!is.numeric(argvals) || length(argvals) < 2L ||
      any(!is.finite(argvals)) || any(diff(argvals) <= 0)) {
    stop("`argvals` must be finite and strictly increasing.", call. = FALSE)
  }
  increments <- diff(argvals)
  c(increments[1L] / 2,
    (increments[-length(increments)] + increments[-1L]) / 2,
    increments[length(increments)] / 2)
}

integral_l2 <- function(values, weights) {
  sqrt(sum(weights * values^2))
}

with_reproducible_seed <- function(seed, expression) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expression)
}

run_task_batch <- function(tasks, task_function, workers = 1L) {
  workers <- max(1L, as.integer(workers))
  if (workers == 1L || .Platform$OS.type == "windows") {
    return(lapply(tasks, task_function))
  }
  parallel::mclapply(
    tasks, task_function, mc.cores = workers,
    mc.preschedule = FALSE, mc.set.seed = FALSE
  )
}

bind_rows_base <- function(objects) {
  objects <- Filter(function(x) !is.null(x) && nrow(x) > 0L, objects)
  if (!length(objects)) return(data.frame())
  columns <- unique(unlist(lapply(objects, names), use.names = FALSE))
  objects <- lapply(objects, function(object) {
    missing <- setdiff(columns, names(object))
    for (column in missing) object[[column]] <- NA
    object[columns]
  })
  rownames_result <- NULL
  result <- do.call(rbind, objects)
  rownames(result) <- rownames_result
  result
}

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

mae <- function(observed, predicted) {
  mean(abs(observed - predicted))
}

r_squared <- function(observed, predicted) {
  denominator <- sum((observed - mean(observed))^2)
  if (denominator <= 0) return(NA_real_)
  1 - sum((observed - predicted)^2) / denominator
}

metric_vector <- function(observed, predicted) {
  complete <- is.finite(observed) & is.finite(predicted)
  if (!all(complete) || !any(complete)) {
    return(c(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_))
  }
  c(
    RMSE = rmse(observed, predicted),
    MAE = mae(observed, predicted),
    R2 = r_squared(observed, predicted)
  )
}

latex_escape <- function(x) {
  x <- as.character(x)
  replacements <- c(
    "\\" = "\\textbackslash{}", "&" = "\\&", "%" = "\\%",
    "$" = "\\$", "#" = "\\#", "_" = "\\_",
    "{" = "\\{", "}" = "\\}", "~" = "\\textasciitilde{}",
    "^" = "\\textasciicircum{}"
  )
  for (pattern in names(replacements)) {
    x <- gsub(pattern, replacements[[pattern]], x, fixed = TRUE)
  }
  x
}

