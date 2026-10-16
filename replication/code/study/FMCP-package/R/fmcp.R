#' Fit functional MCP (fMCP) regression
#'
#' `fmcp()` fits the scalar-on-function linear model
#' \deqn{Y_i = \mu + \int X_i(t)\beta(t)\,dt + \varepsilon_i}
#' with a B-spline representation of \eqn{\beta(t)}, a derivative roughness
#' penalty, and an interval-wise functional minimax concave penalty (fMCP).
#'
#' For interval \eqn{j}, let \eqn{W_j} be its local B-spline Gram matrix and
#' \deqn{q_j = \{(M/T)b^T W_j b\}^{1/2}.}
#' The fitted objective is
#' \deqn{n^{-1}\|y-Ub\|_2^2 + \kappa b^T Vb +
#'       \sum_j p_{\lambda,a}(q_j),}
#' where
#' \deqn{p_{\lambda,a}(u) = \lambda u-u^2/(2a)}
#' for \eqn{0 \le u \le a\lambda}, and
#' \eqn{p_{\lambda,a}(u)=a\lambda^2/2} otherwise. The nonconvex objective is
#' optimized with a local quadratic approximation initialized from the
#' roughness-only fit. No predictor or response standardization is performed;
#' the response and functional predictor are centered internally to estimate
#' an intercept.
#'
#' @param y Finite numeric response vector of length `n`.
#' @param x Finite numeric `n` by `p` matrix. Rows are subjects and columns are
#'   evaluations of a functional predictor on the common `argvals` grid.
#' @param argvals Strictly increasing finite numeric vector of length `p`.
#' @param M Positive integer giving the number of equal B-spline subintervals.
#' @param lambda Nonempty nonnegative grid for the fMCP shrinkage parameter.
#' @param kappa Nonempty nonnegative grid for the derivative roughness penalty.
#' @param degree B-spline polynomial degree. Cubic splines (`degree = 3`) are
#'   used by default.
#' @param penalty_order Positive derivative order used in the roughness
#'   penalty. It cannot exceed `degree`.
#' @param mcp_gamma One or more MCP concavity parameters, all strictly greater
#'   than one. If several values are supplied, they are included in BIC model
#'   selection.
#' @param integration Numerical integration rule for the observed functional
#'   predictor. Simpson integration requires an equally spaced grid with an
#'   odd number of points; trapezoid integration accepts irregular grids.
#' @param control Numerical controls from [fmcp_control()].
#' @param ... Arguments forwarded by the `fMCP()` and `FMCP()` aliases to
#'   `fmcp()`.
#'
#' @return An object of class `fmcp`. It contains the selected basis
#'   coefficients, coefficient function on `argvals`, fitted values,
#'   residuals, interval norms, complete tuning table, and convergence
#'   diagnostics.
#'
#' @references
#' Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally Sparse Estimator
#' for Functional Linear Regression Models. *Journal of Computational and
#' Graphical Statistics*, 26(2), 306-318.
#'
#' Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
#' concave penalty. *The Annals of Statistics*, 38(2), 894-942.
#'
#' @examples
#' set.seed(42)
#' n <- 40
#' argvals <- seq(0, 1, length.out = 21)
#' scores <- matrix(rnorm(n * 6), n, 6)
#' x <- matrix(0, n, length(argvals))
#' for (j in seq_len(6)) {
#'   x <- x + scores[, j] %o% sin(j * pi * argvals)
#' }
#' x <- x + matrix(rnorm(n * length(argvals), sd = 0.05), n)
#' beta <- ifelse(argvals <= 0.55, 2 * sin(2 * pi * argvals), 0)
#' y <- rowMeans(x * rep(beta, each = n)) + rnorm(n, sd = 0.1)
#' fit <- fmcp(
#'   y, x, argvals,
#'   M = 8,
#'   lambda = c(0.005, 0.02),
#'   kappa = c(1e-5, 1e-4),
#'   control = fmcp_control(max_iter = 50)
#' )
#' coef(fit, type = "basis")
#' predict(fit, x[1:3, , drop = FALSE])
#'
#' @export
fmcp <- function(y,
                 x,
                 argvals,
                 M,
                 lambda,
                 kappa,
                 degree = 3L,
                 penalty_order = 2L,
                 mcp_gamma = 3,
                 integration = c("simpson", "trapezoid"),
                 control = fmcp_control()) {
  call <- match.call()
  validated <- .fmcp_validate_data(y, x, argvals)
  y <- validated$y
  x <- validated$x
  argvals <- validated$argvals

  M <- .fmcp_scalar_integer(M, "M", lower = 2L)
  degree <- .fmcp_scalar_integer(degree, "degree", lower = 1L)
  penalty_order <- .fmcp_scalar_integer(
    penalty_order, "penalty_order", lower = 1L
  )
  if (penalty_order > degree) {
    stop("`penalty_order` cannot exceed `degree`.", call. = FALSE)
  }
  lambda <- .fmcp_numeric_grid(lambda, "lambda", lower = 0)
  kappa <- .fmcp_numeric_grid(kappa, "kappa", lower = 0)
  mcp_gamma <- .fmcp_numeric_grid(
    mcp_gamma, "mcp_gamma", lower = 1, strict = TRUE
  )
  integration <- match.arg(integration)
  control <- .fmcp_validate_control(control)

  integration_weights <- .fmcp_integration_weights(argvals, integration)
  basis <- .fmcp_basis_setup(argvals, M, degree, penalty_order)

  y_mean <- mean(y)
  x_mean <- colMeans(x)
  centered_y <- y - y_mean
  centered_x <- sweep(x, 2L, x_mean, FUN = "-")

  design <- cpp_design_matrix(
    centered_x, basis$basis_observed, integration_weights
  )
  local_gram <- cpp_local_gram(
    basis$basis_quad,
    basis$quad_weights,
    basis$quad_intervals,
    M,
    degree + 1L
  )
  roughness <- cpp_weighted_crossprod(
    basis$derivative_quad, basis$quad_weights
  )
  interval_scale <- sqrt(M / diff(basis$domain))

  grid <- cpp_fmcp_grid(
    centered_y,
    design,
    roughness,
    local_gram,
    lambda,
    kappa,
    mcp_gamma,
    interval_scale,
    control$max_iter,
    control$tol,
    control$norm_floor,
    control$cutoff,
    control$refit,
    control$trace
  )

  tuning <- data.frame(
    lambda = grid$lambda,
    kappa = grid$kappa,
    mcp_gamma = grid$mcp_gamma,
    bic = grid$bic,
    rss = grid$rss,
    df = grid$df,
    objective = grid$objective,
    iterations = grid$iterations,
    active = grid$active,
    success = grid$success,
    converged = grid$converged,
    refit_success = grid$refit_success,
    objective_increased = grid$objective_increased,
    message = as.character(grid$message),
    stringsAsFactors = FALSE
  )

  valid <- which(tuning$success %in% TRUE & is.finite(tuning$bic))
  if (!length(valid)) {
    details <- unique(tuning$message[nzchar(tuning$message)])
    stop(
      "No tuning candidate produced a valid fit.",
      if (length(details)) paste0(" ", paste(details, collapse = " ")) else "",
      call. = FALSE
    )
  }
  converged_valid <- valid[tuning$converged[valid] %in% TRUE]
  eligible <- if (length(converged_valid)) converged_valid else valid
  if (!length(converged_valid)) {
    warning(
      "No tuning candidate converged before `max_iter`; selecting the valid ",
      "candidate with minimum BIC.",
      call. = FALSE
    )
  }
  selected <- eligible[which.min(tuning$bic[eligible])]

  coefficients <- as.numeric(grid$coefficients[, selected])
  names(coefficients) <- paste0("B", seq_along(coefficients))
  beta <- as.numeric(basis$basis_observed %*% coefficients)
  fitted_values <- as.numeric(y_mean + design %*% coefficients)
  residual_values <- y - fitted_values
  interval_norms <- as.numeric(
    cpp_interval_norms(coefficients, local_gram, interval_scale)
  )
  selected_row <- tuning[selected, , drop = FALSE]

  if (!isTRUE(selected_row$converged)) {
    warning("The selected fMCP fit did not satisfy the convergence tolerance.",
            call. = FALSE)
  }
  if (isTRUE(selected_row$objective_increased)) {
    warning(
      "The exact objective increased during at least one selected-fit LQA step; ",
      "inspect `diagnostics` and the tuning table.",
      call. = FALSE
    )
  }

  internal <- if (control$keep_matrices) {
    list(
      centered_y = centered_y,
      design = design,
      roughness = roughness,
      local_gram = local_gram
    )
  } else {
    NULL
  }

  structure(
    list(
      call = call,
      coefficients = coefficients,
      beta = beta,
      argvals = argvals,
      fitted.values = fitted_values,
      residuals = residual_values,
      interval_norms = interval_norms,
      selected = selected_row,
      tuning = tuning,
      diagnostics = list(
        selected_index = selected,
        failed_candidates = sum(!tuning$success),
        nonconverged_candidates = sum(tuning$success & !tuning$converged),
        objective_trace = if (control$trace) grid$traces[[selected]] else NULL
      ),
      means = list(y = y_mean, x = x_mean),
      dimensions = c(n = nrow(x), p = ncol(x), basis = length(coefficients)),
      basis = list(
        M = M,
        degree = degree,
        penalty_order = penalty_order,
        domain = basis$domain,
        breaks = basis$breaks,
        full_knots = basis$full_knots,
        observed = basis$basis_observed,
        integration = integration,
        integration_weights = integration_weights,
        interval_scale = interval_scale
      ),
      control = control,
      internal = internal
    ),
    class = "fmcp"
  )
}

#' @rdname fmcp
#' @export
fMCP <- function(...) {
  fmcp(...)
}

#' @rdname fmcp
#' @export
FMCP <- function(...) {
  fmcp(...)
}
