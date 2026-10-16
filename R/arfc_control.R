#' Control adaptive regional folded-concave fitting
#'
#' `arfc_control()` collects numerical and adaptive-weight controls used by
#' [arfc()] and [arfc_select()]. The defaults reproduce the computational
#' convention used in the accompanying ARFC study. No training/validation
#' split is created by these controls.
#'
#' @param max_iter Positive integer giving the maximum LQA iterations for each
#'   tuning candidate.
#' @param tol Positive relative coefficient convergence tolerance.
#' @param norm_floor Positive stabilizer for regional norms near zero.
#' @param relative_cutoff_multiplier Nonnegative multiplier applied to the
#'   largest absolute smooth-pilot coefficient. The resulting threshold is
#'   used for exact coefficient thresholding.
#' @param cutoff_floor Positive lower bound for the absolute threshold.
#' @param refit Logical; refit the roughness-penalized model on retained basis
#'   coefficients after thresholding.
#' @param trace Logical; retain exact-objective traces for tuning candidates.
#' @param keep_matrices Logical; retain design and regional matrices in fitted
#'   objects. This can substantially increase object size.
#' @param adaptive_alpha Positive exponent in the pilot-adaptive regional
#'   weights.
#' @param adaptive_epsilon_fraction Positive fraction of the median positive
#'   pilot norm used in the stabilizing offset.
#' @param adaptive_clip Two positive increasing numbers bounding normalized
#'   adaptive weights.
#' @param zero_tolerance Nonnegative tolerance used only to display estimated
#'   zero regions in plots.
#' @param tie_relative_tolerance Nonnegative relative tolerance for validation
#'   loss ties in [arfc_select()].
#' @param strict Logical. If `TRUE`, only successful, converged, successfully
#'   refitted candidates with finite coefficients and BIC are eligible.
#'
#' @return An object of class `arfc_control`.
#' @export
arfc_control <- function(max_iter = 120L,
                         tol = 1e-8,
                         norm_floor = 1e-8,
                         relative_cutoff_multiplier = 1e-4,
                         cutoff_floor = 1e-8,
                         refit = TRUE,
                         trace = FALSE,
                         keep_matrices = FALSE,
                         adaptive_alpha = 1,
                         adaptive_epsilon_fraction = 0.05,
                         adaptive_clip = c(0.25, 4),
                         zero_tolerance = 1e-8,
                         tie_relative_tolerance = 1e-12,
                         strict = TRUE) {
  max_iter <- .fmcp_scalar_integer(max_iter, "max_iter", lower = 1L)
  tol <- .fmcp_scalar_numeric(tol, "tol", lower = 0, strict = TRUE)
  norm_floor <- .fmcp_scalar_numeric(
    norm_floor, "norm_floor", lower = 0, strict = TRUE
  )
  relative_cutoff_multiplier <- .fmcp_scalar_numeric(
    relative_cutoff_multiplier, "relative_cutoff_multiplier", lower = 0
  )
  cutoff_floor <- .fmcp_scalar_numeric(
    cutoff_floor, "cutoff_floor", lower = 0, strict = TRUE
  )
  refit <- .fmcp_scalar_logical(refit, "refit")
  trace <- .fmcp_scalar_logical(trace, "trace")
  keep_matrices <- .fmcp_scalar_logical(keep_matrices, "keep_matrices")
  adaptive_alpha <- .fmcp_scalar_numeric(
    adaptive_alpha, "adaptive_alpha", lower = 0, strict = TRUE
  )
  adaptive_epsilon_fraction <- .fmcp_scalar_numeric(
    adaptive_epsilon_fraction,
    "adaptive_epsilon_fraction",
    lower = 0,
    strict = TRUE
  )
  if (!is.numeric(adaptive_clip) || length(adaptive_clip) != 2L ||
      anyNA(adaptive_clip) || any(!is.finite(adaptive_clip)) ||
      adaptive_clip[1L] <= 0 || adaptive_clip[2L] <= adaptive_clip[1L]) {
    stop(
      "`adaptive_clip` must contain two finite positive increasing values.",
      call. = FALSE
    )
  }
  zero_tolerance <- .fmcp_scalar_numeric(
    zero_tolerance, "zero_tolerance", lower = 0
  )
  tie_relative_tolerance <- .fmcp_scalar_numeric(
    tie_relative_tolerance, "tie_relative_tolerance", lower = 0
  )
  strict <- .fmcp_scalar_logical(strict, "strict")

  structure(
    list(
      max_iter = max_iter,
      tol = tol,
      norm_floor = norm_floor,
      relative_cutoff_multiplier = relative_cutoff_multiplier,
      cutoff_floor = cutoff_floor,
      refit = refit,
      trace = trace,
      keep_matrices = keep_matrices,
      adaptive_alpha = adaptive_alpha,
      adaptive_epsilon_fraction = adaptive_epsilon_fraction,
      adaptive_clip = as.numeric(adaptive_clip),
      zero_tolerance = zero_tolerance,
      tie_relative_tolerance = tie_relative_tolerance,
      strict = strict
    ),
    class = "arfc_control"
  )
}

.arfc_validate_control <- function(control) {
  if (missing(control) || is.null(control)) {
    return(arfc_control())
  }
  if (!is.list(control) || is.null(names(control)) ||
      any(names(control) == "")) {
    stop(
      "`control` must be created by `arfc_control()` or be a named list.",
      call. = FALSE
    )
  }
  defaults <- unclass(arfc_control())
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop(
      "Unknown `control` element(s): ", paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  do.call(arfc_control, utils::modifyList(defaults, unclass(control)))
}
