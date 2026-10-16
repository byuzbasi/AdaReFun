#' Fit adaptive regional folded-concave functional regression
#'
#' `arfc()` fits a scalar-on-function linear model with pilot-adaptive regional
#' MCP or SCAD penalties. A single-resolution fit uses the finest partition;
#' a multi-resolution fit uses a nested collection of partitions in the same
#' B-spline coefficient space. Regional multipliers are normalized by domain
#' mass, and level weights must sum to one.
#'
#' The roughness-only pilot chooses one `kappa` value by strict BIC. Its
#' regional norms define level-normalized and clipped adaptive weights. The
#' folded-concave path is then evaluated over `lambda` at that common
#' smoothness value. No response or predictor standardization is performed;
#' centering estimates the intercept.
#'
#' @param y Finite numeric response vector.
#' @param x Finite numeric matrix with subjects in rows and functional grid
#'   points in columns.
#' @param argvals Strictly increasing functional argument values.
#' @param lambda Nonempty nonnegative folded-concave tuning grid.
#' @param kappa Nonempty nonnegative roughness grid for the common smooth
#'   pilot. The BIC-selected value is used for the ARFC path.
#' @param resolution Either `"single"` (ASR) or `"multi"` (AMR).
#' @param penalty Either `"MCP"` or `"SCAD"`.
#' @param levels Strictly increasing nested region counts. Every value must
#'   divide the largest value. The default reproduces the study specification.
#' @param level_weights Positive weights, one per `levels`, summing to one.
#'   `NULL` gives equal weights. They are not silently normalized.
#' @param adaptive Logical. Use smooth-pilot adaptive regional weights. Setting
#'   this to `FALSE` gives the corresponding nonadaptive regional fit.
#' @param degree B-spline polynomial degree.
#' @param penalty_order Derivative order for roughness penalization.
#' @param mcp_gamma MCP shape parameter greater than one.
#' @param scad_a SCAD shape parameter greater than two. The defaults
#'   `mcp_gamma = 2.7` and `scad_a = 3.7` have equal scalar maximum concavity.
#' @param integration Numerical integration rule for the observed predictor.
#' @param control Numerical and adaptive-weight settings from [arfc_control()].
#'
#' @return An object of class `arfc` containing the training-BIC-tuned fit for
#'   the one requested penalty-resolution specification: its coefficient
#'   function, regional norms and weights, complete tuning diagnostics, smooth
#'   pilot, fitted values, and prediction metadata. This function does not
#'   choose between ASR and AMR or between MCP and SCAD; use [arfc_select()] for
#'   that four-candidate comparison.
#'
#' @references
#' Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized
#' likelihood and its oracle properties. *Journal of the American Statistical
#' Association*, 96(456), 1348-1360.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *The Annals of Statistics*, 38(2), 894-942.
#'
#' Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally sparse estimator
#' for functional linear regression models. *Journal of Computational and
#' Graphical Statistics*, 26(2), 306-318.
#'
#' @examples
#' set.seed(2026)
#' n <- 36
#' argvals <- seq(0, 1, length.out = 21)
#' x <- matrix(rnorm(n * length(argvals)), n)
#' beta <- ifelse(argvals < 0.35, 1.5 * sin(2 * pi * argvals), 0)
#' w <- c(0.025, rep(0.05, 19), 0.025)
#' y <- as.numeric(x %*% (w * beta)) + rnorm(n, sd = 0.2)
#' fit <- arfc(
#'   y, x, argvals,
#'   lambda = c(0.002, 0.02),
#'   kappa = c(1e-5, 1e-3),
#'   resolution = "single",
#'   penalty = "MCP",
#'   levels = c(4L, 8L),
#'   control = arfc_control(max_iter = 40L)
#' )
#' coef(fit, type = "basis")
#' predict(fit, x[1:2, , drop = FALSE])
#'
#' @export
arfc <- function(y,
                 x,
                 argvals,
                 lambda,
                 kappa,
                 resolution = c("single", "multi"),
                 penalty = c("MCP", "SCAD"),
                 levels = c(8L, 16L, 32L),
                 level_weights = NULL,
                 adaptive = TRUE,
                 degree = 3L,
                 penalty_order = 2L,
                 mcp_gamma = 2.7,
                 scad_a = 3.7,
                 integration = c("simpson", "trapezoid"),
                 control = arfc_control()) {
  call <- match.call()
  resolution <- match.arg(tolower(resolution), c("single", "multi"))
  penalty <- match.arg(toupper(penalty), c("MCP", "SCAD"))
  adaptive <- .fmcp_scalar_logical(adaptive, "adaptive")
  degree <- .fmcp_scalar_integer(degree, "degree", lower = 1L)
  penalty_order <- .fmcp_scalar_integer(
    penalty_order, "penalty_order", lower = 1L
  )
  if (penalty_order > degree) {
    stop("`penalty_order` cannot exceed `degree`.", call. = FALSE)
  }
  lambda <- .fmcp_numeric_grid(lambda, "lambda", lower = 0)
  kappa <- .fmcp_numeric_grid(kappa, "kappa", lower = 0)
  mcp_gamma <- .fmcp_scalar_numeric(
    mcp_gamma, "mcp_gamma", lower = 1, strict = TRUE
  )
  scad_a <- .fmcp_scalar_numeric(
    scad_a, "scad_a", lower = 2, strict = TRUE
  )
  integration <- match.arg(integration)
  control <- .arfc_validate_control(control)
  specification <- .arfc_validate_levels(
    levels, level_weights, resolution
  )
  context <- .arfc_prepare_context(
    y, x, argvals,
    specification$levels,
    specification$level_weights,
    degree,
    penalty_order,
    integration
  )
  pilot <- .arfc_smooth_pilot(context, kappa, mcp_gamma, control)
  weights <- .arfc_adaptive_weights(
    pilot$coefficients, context, adaptive, control
  )
  .arfc_fit_prepared(
    context = context,
    resolution = resolution,
    penalty = penalty,
    lambda = lambda,
    mcp_gamma = mcp_gamma,
    scad_a = scad_a,
    adaptive = adaptive,
    pilot = pilot,
    adaptive_weights = weights,
    control = control,
    call = call
  )
}

#' Select the ARFC penalty family and resolution by validation loss
#'
#' `arfc_select()` fits the four prespecified candidates ASR-MCP, ASR-SCAD,
#' AMR-MCP, and AMR-SCAD to the same training data, pilot, tuning grids, and
#' preprocessing. It selects the smallest response MSE on an explicitly
#' supplied validation sample. The validation responses are never used to fit
#' a candidate. Numerical ties follow the documented order above.
#'
#' @inheritParams arfc
#' @param validation Named list with finite `y` and numeric matrix `x` on the
#'   same functional grid as the training data. A split is never created
#'   silently.
#'
#' @return An object of class `arfc_select`. It contains one `selected_fit`, all
#'   four candidate fits, validation losses, tie diagnostics, and an explicit
#'   no-leakage record. The candidates are not averaged or ensembled. A final
#'   combined-data refit is not performed automatically.
#'
#' @examples
#' set.seed(2027)
#' argvals <- seq(0, 1, length.out = 21)
#' x <- matrix(rnorm(48 * length(argvals)), 48)
#' beta <- ifelse(argvals >= 0.35 & argvals <= 0.65, 0,
#'                sin(2 * pi * argvals))
#' w <- c(0.025, rep(0.05, 19), 0.025)
#' y <- as.numeric(x %*% (w * beta)) + rnorm(48, sd = 0.2)
#' selected <- arfc_select(
#'   y[1:32], x[1:32, , drop = FALSE], argvals,
#'   validation = list(y = y[33:48], x = x[33:48, , drop = FALSE]),
#'   lambda = c(0.002, 0.02),
#'   kappa = c(1e-5, 1e-3),
#'   levels = c(4L, 8L),
#'   control = arfc_control(max_iter = 40L)
#' )
#' selected$selection
#'
#' @export
arfc_select <- function(y,
                        x,
                        argvals,
                        validation,
                        lambda,
                        kappa,
                        levels = c(8L, 16L, 32L),
                        level_weights = NULL,
                        adaptive = TRUE,
                        degree = 3L,
                        penalty_order = 2L,
                        mcp_gamma = 2.7,
                        scad_a = 3.7,
                        integration = c("simpson", "trapezoid"),
                        control = arfc_control()) {
  call <- match.call()
  if (!is.list(validation) ||
      !all(c("y", "x") %in% names(validation))) {
    stop(
      "`validation` must be a named list containing `y` and `x`.",
      call. = FALSE
    )
  }
  adaptive <- .fmcp_scalar_logical(adaptive, "adaptive")
  degree <- .fmcp_scalar_integer(degree, "degree", lower = 1L)
  penalty_order <- .fmcp_scalar_integer(
    penalty_order, "penalty_order", lower = 1L
  )
  if (penalty_order > degree) {
    stop("`penalty_order` cannot exceed `degree`.", call. = FALSE)
  }
  lambda <- .fmcp_numeric_grid(lambda, "lambda", lower = 0)
  kappa <- .fmcp_numeric_grid(kappa, "kappa", lower = 0)
  mcp_gamma <- .fmcp_scalar_numeric(
    mcp_gamma, "mcp_gamma", lower = 1, strict = TRUE
  )
  scad_a <- .fmcp_scalar_numeric(
    scad_a, "scad_a", lower = 2, strict = TRUE
  )
  integration <- match.arg(integration)
  control <- .arfc_validate_control(control)
  specification <- .arfc_validate_levels(levels, level_weights, "multi")
  context <- .arfc_prepare_context(
    y, x, argvals,
    specification$levels,
    specification$level_weights,
    degree,
    penalty_order,
    integration
  )
  validation_data <- .fmcp_validate_data(
    validation$y, validation$x, context$argvals
  )
  pilot <- .arfc_smooth_pilot(context, kappa, mcp_gamma, control)
  weights <- .arfc_adaptive_weights(
    pilot$coefficients, context, adaptive, control
  )
  candidate_specifications <- data.frame(
    method = c("asr_mcp", "asr_scad", "amr_mcp", "amr_scad"),
    label = c("ASR--MCP", "ASR--SCAD", "AMR--MCP", "AMR--SCAD"),
    resolution = c("single", "single", "multi", "multi"),
    penalty = c("MCP", "SCAD", "MCP", "SCAD"),
    tie_priority = 1:4,
    stringsAsFactors = FALSE
  )
  candidates <- stats::setNames(
    vector("list", nrow(candidate_specifications)),
    candidate_specifications$method
  )
  for (index in seq_len(nrow(candidate_specifications))) {
    specification_row <- candidate_specifications[index, ]
    candidates[[specification_row$method]] <- .arfc_fit_prepared(
      context = context,
      resolution = specification_row$resolution,
      penalty = specification_row$penalty,
      lambda = lambda,
      mcp_gamma = mcp_gamma,
      scad_a = scad_a,
      adaptive = adaptive,
      pilot = pilot,
      adaptive_weights = weights,
      control = control,
      call = call
    )
  }
  centered_validation <- sweep(
    validation_data$x, 2L, context$means$x, FUN = "-"
  )
  validation_design <- cpp_design_matrix(
    centered_validation,
    context$basis$basis_observed,
    context$integration_weights
  )
  losses <- vapply(candidates, function(fit) {
    prediction <- as.numeric(
      context$means$y + validation_design %*% fit$coefficients
    )
    mean((validation_data$y - prediction)^2)
  }, numeric(1L))
  if (any(!is.finite(losses))) {
    stop("At least one ARFC validation loss is nonfinite.", call. = FALSE)
  }
  minimum <- min(losses)
  tolerance <- control$tie_relative_tolerance * (1 + abs(minimum))
  tied <- names(losses)[losses <= minimum + tolerance]
  selected_method <- candidate_specifications$method[
    candidate_specifications$method %in% tied
  ][1L]
  second <- sort(losses, partial = 2L)[2L]
  selection <- candidate_specifications
  selection$validation_response_mse <- unname(losses[selection$method])
  selection$difference_from_minimum <-
    selection$validation_response_mse - minimum
  selection$tied <- selection$method %in% tied
  selection$selected <- selection$method == selected_method
  structure(
    list(
      call = call,
      selected_method = selected_method,
      selected_label = selection$label[selection$selected],
      selected_fit = candidates[[selected_method]],
      candidates = candidates,
      selection = selection,
      diagnostics = list(
        minimum_validation_mse = minimum,
        margin_to_second = second - minimum,
        tie_tolerance = tolerance,
        tie_count = length(tied),
        tie_order = candidate_specifications$method,
        validation_used_for_fitting = FALSE,
        validation_sample_size = length(validation_data$y)
      ),
      training_dimensions = c(
        n = nrow(context$x), p = ncol(context$x),
        basis = ncol(context$basis$basis_observed)
      ),
      validation_dimensions = c(
        n = nrow(validation_data$x), p = ncol(validation_data$x)
      ),
      control = control
    ),
    class = "arfc_select"
  )
}
