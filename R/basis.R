.fmcp_integration_weights <- function(argvals, method) {
  method <- match.arg(method, c("simpson", "trapezoid"))
  p <- length(argvals)

  if (method == "simpson") {
    if (p < 3L || p %% 2L == 0L) {
      stop(
        "Simpson integration requires an odd number of at least three grid points.",
        call. = FALSE
      )
    }
    gaps <- diff(argvals)
    scale <- max(1, abs(argvals[length(argvals)] - argvals[1L]))
    if (max(abs(gaps - gaps[1L])) > 100 * .Machine$double.eps * scale) {
      stop("Simpson integration requires an equally spaced `argvals` grid.",
           call. = FALSE)
    }
    pattern <- rep(c(4, 2), length.out = p - 2L)
    return(gaps[1L] / 3 * c(1, pattern, 1))
  }

  gaps <- diff(argvals)
  c(gaps[1L] / 2, (gaps[-length(gaps)] + gaps[-1L]) / 2,
    gaps[length(gaps)] / 2)
}

.fmcp_gauss_legendre <- function(order) {
  order <- .fmcp_scalar_integer(order, "quadrature order", lower = 1L)
  if (order == 1L) {
    return(list(nodes = 0, weights = 2))
  }
  index <- seq_len(order - 1L)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi <- matrix(0, order, order)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  ordering <- order(decomposition$values)
  list(
    nodes = decomposition$values[ordering],
    weights = 2 * decomposition$vectors[1L, ordering]^2
  )
}

.fmcp_basis_setup <- function(argvals, M, degree, penalty_order) {
  domain <- range(argvals)
  breaks <- seq(domain[1L], domain[2L], length.out = M + 1L)
  spline_order <- degree + 1L
  internal <- if (M > 1L) breaks[2L:M] else numeric()
  full_knots <- c(
    rep(domain[1L], spline_order),
    internal,
    rep(domain[2L], spline_order)
  )

  basis_observed <- splines::splineDesign(
    knots = full_knots,
    x = argvals,
    ord = spline_order,
    derivs = rep.int(0L, length(argvals)),
    outer.ok = TRUE
  )

  quadrature_order <- max(degree + 1L, 2L)
  rule <- .fmcp_gauss_legendre(quadrature_order)
  quad_values <- numeric(M * quadrature_order)
  quad_weights <- numeric(M * quadrature_order)
  quad_intervals <- integer(M * quadrature_order)

  offset <- 0L
  for (j in seq_len(M)) {
    half_width <- (breaks[j + 1L] - breaks[j]) / 2
    midpoint <- (breaks[j + 1L] + breaks[j]) / 2
    indices <- offset + seq_len(quadrature_order)
    quad_values[indices] <- midpoint + half_width * rule$nodes
    quad_weights[indices] <- half_width * rule$weights
    quad_intervals[indices] <- j
    offset <- offset + quadrature_order
  }

  basis_quad <- splines::splineDesign(
    knots = full_knots,
    x = quad_values,
    ord = spline_order,
    derivs = rep.int(0L, length(quad_values)),
    outer.ok = TRUE
  )
  derivative_quad <- splines::splineDesign(
    knots = full_knots,
    x = quad_values,
    ord = spline_order,
    derivs = rep.int(penalty_order, length(quad_values)),
    outer.ok = TRUE
  )

  list(
    domain = domain,
    breaks = breaks,
    full_knots = full_knots,
    degree = degree,
    order = spline_order,
    M = M,
    basis_observed = unclass(basis_observed),
    basis_quad = unclass(basis_quad),
    derivative_quad = unclass(derivative_quad),
    quad_weights = quad_weights,
    quad_intervals = quad_intervals
  )
}

.fmcp_validate_data <- function(y, x, argvals) {
  if (!is.numeric(y) || !length(y) || anyNA(y) || any(!is.finite(y))) {
    stop("`y` must be a nonempty finite numeric vector.", call. = FALSE)
  }
  if (!is.matrix(x) || !is.numeric(x) || anyNA(x) || any(!is.finite(x))) {
    stop("`x` must be a finite numeric matrix with subjects in rows.",
         call. = FALSE)
  }
  if (nrow(x) != length(y)) {
    stop("`nrow(x)` must equal `length(y)`.", call. = FALSE)
  }
  if (nrow(x) < 2L || ncol(x) < 3L) {
    stop("`x` must contain at least two subjects and three grid points.",
         call. = FALSE)
  }
  if (!is.numeric(argvals) || length(argvals) != ncol(x) || anyNA(argvals) ||
      any(!is.finite(argvals))) {
    stop("`argvals` must be a finite numeric vector of length `ncol(x)`.",
         call. = FALSE)
  }
  if (any(diff(argvals) <= 0)) {
    stop("`argvals` must be strictly increasing.", call. = FALSE)
  }
  list(y = as.numeric(y), x = unclass(x), argvals = as.numeric(argvals))
}
