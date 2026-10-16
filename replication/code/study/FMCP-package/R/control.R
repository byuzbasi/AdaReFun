#' Control numerical details of fMCP and fSCAD fitting
#'
#' @param max_iter Positive integer. Maximum local quadratic approximation
#'   iterations for each tuning-parameter combination.
#' @param tol Positive number. Relative coefficient convergence tolerance.
#' @param norm_floor Positive number used only to stabilize division by an
#'   interval norm near zero.
#' @param cutoff Nonnegative absolute threshold applied to B-spline
#'   coefficients after optimization. Coefficients below this threshold are
#'   set exactly to zero.
#' @param refit Logical. If `TRUE`, refit the roughness-penalized model on the
#'   selected B-spline coefficients after thresholding.
#' @param trace Logical. If `TRUE`, retain per-candidate objective traces.
#' @param keep_matrices Logical. If `TRUE`, retain design and penalty matrices
#'   in the returned fit. This is useful for diagnostics and testing but may
#'   substantially increase object size.
#'
#' @return An object of class `fmcp_control`.
#' @export
fmcp_control <- function(max_iter = 100L,
                         tol = 1e-8,
                         norm_floor = 1e-8,
                         cutoff = 1e-6,
                         refit = TRUE,
                         trace = FALSE,
                         keep_matrices = FALSE) {
  max_iter <- .fmcp_scalar_integer(max_iter, "max_iter", lower = 1L)
  tol <- .fmcp_scalar_numeric(tol, "tol", lower = 0, strict = TRUE)
  norm_floor <- .fmcp_scalar_numeric(
    norm_floor, "norm_floor", lower = 0, strict = TRUE
  )
  cutoff <- .fmcp_scalar_numeric(cutoff, "cutoff", lower = 0)
  refit <- .fmcp_scalar_logical(refit, "refit")
  trace <- .fmcp_scalar_logical(trace, "trace")
  keep_matrices <- .fmcp_scalar_logical(keep_matrices, "keep_matrices")

  structure(
    list(
      max_iter = max_iter,
      tol = tol,
      norm_floor = norm_floor,
      cutoff = cutoff,
      refit = refit,
      trace = trace,
      keep_matrices = keep_matrices
    ),
    class = "fmcp_control"
  )
}

.fmcp_validate_control <- function(control) {
  if (missing(control) || is.null(control)) {
    return(fmcp_control())
  }
  if (!is.list(control)) {
    stop("`control` must be created by `fmcp_control()` or be a named list.",
         call. = FALSE)
  }
  if (is.null(names(control)) || any(names(control) == "")) {
    stop("Every supplied `control` element must be named.", call. = FALSE)
  }
  defaults <- unclass(fmcp_control())
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown `control` element(s): ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  do.call(fmcp_control, utils::modifyList(defaults, unclass(control)))
}

.fmcp_scalar_numeric <- function(x, name, lower = -Inf, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("`", name, "` must be one finite numeric value.", call. = FALSE)
  }
  invalid <- if (strict) x <= lower else x < lower
  if (invalid) {
    relation <- if (strict) "greater than" else "at least"
    stop("`", name, "` must be ", relation, " ", lower, ".", call. = FALSE)
  }
  as.numeric(x)
}

.fmcp_scalar_integer <- function(x, name, lower = -Inf) {
  value <- .fmcp_scalar_numeric(x, name, lower = lower)
  if (value != round(value)) {
    stop("`", name, "` must be an integer.", call. = FALSE)
  }
  as.integer(value)
}

.fmcp_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  x
}

.fmcp_numeric_grid <- function(x, name, lower = 0, strict = FALSE) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x))) {
    stop("`", name, "` must be a nonempty finite numeric vector.", call. = FALSE)
  }
  invalid <- if (strict) x <= lower else x < lower
  if (any(invalid)) {
    relation <- if (strict) "greater than" else "at least"
    stop("Every `", name, "` value must be ", relation, " ", lower, ".",
         call. = FALSE)
  }
  unique(as.numeric(x))
}
