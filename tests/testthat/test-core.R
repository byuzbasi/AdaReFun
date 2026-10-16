test_that("MCP penalty and derivative match their analytic definitions", {
  values <- c(0, 0.5, 1.5, 3, 4)
  lambda <- 1
  gamma <- 3
  expected_penalty <- ifelse(
    values <= gamma * lambda,
    lambda * values - values^2 / (2 * gamma),
    gamma * lambda^2 / 2
  )
  expected_derivative <- pmax(lambda - values / gamma, 0)

  expect_equal(
    as.numeric(AdaReFun:::cpp_mcp_penalty(values, lambda, gamma)),
    expected_penalty,
    tolerance = 1e-14
  )
  expect_equal(
    as.numeric(AdaReFun:::cpp_mcp_derivative(values, lambda, gamma)),
    expected_derivative,
    tolerance = 1e-14
  )

  point <- 1.2
  step <- 1e-6
  numerical <- (
    AdaReFun:::cpp_mcp_penalty(point + step, lambda, gamma) -
      AdaReFun:::cpp_mcp_penalty(point - step, lambda, gamma)
  ) / (2 * step)
  expect_equal(
    as.numeric(numerical),
    as.numeric(AdaReFun:::cpp_mcp_derivative(point, lambda, gamma)),
    tolerance = 1e-8
  )
})

test_that("Rcpp matrix construction matches independent R calculations", {
  argvals <- seq(0, 1, length.out = 11)
  setup <- AdaReFun:::.fmcp_basis_setup(argvals, M = 4L, degree = 3L,
                                    penalty_order = 2L)
  weights <- AdaReFun:::.fmcp_integration_weights(argvals, "simpson")
  set.seed(12)
  x <- matrix(rnorm(7 * length(argvals)), 7, length(argvals))

  design_cpp <- AdaReFun:::cpp_design_matrix(x, setup$basis_observed, weights)
  design_r <- x %*% (weights * setup$basis_observed)
  expect_equal(design_cpp, design_r, tolerance = 1e-13)

  roughness_cpp <- AdaReFun:::cpp_weighted_crossprod(
    setup$derivative_quad, setup$quad_weights
  )
  roughness_r <- crossprod(
    setup$derivative_quad,
    setup$quad_weights * setup$derivative_quad
  )
  expect_equal(roughness_cpp, roughness_r, tolerance = 1e-11)

  local_cpp <- AdaReFun:::cpp_local_gram(
    setup$basis_quad,
    setup$quad_weights,
    setup$quad_intervals,
    4L,
    4L
  )
  local_r <- array(0, dim = c(4L, 4L, 4L))
  for (j in seq_len(4L)) {
    rows <- setup$quad_intervals == j
    columns <- j:(j + 3L)
    local_basis <- setup$basis_quad[rows, columns, drop = FALSE]
    local_r[, , j] <- crossprod(
      local_basis,
      setup$quad_weights[rows] * local_basis
    )
  }
  expect_equal(local_cpp, local_r, tolerance = 1e-13)
})

test_that("B-spline and Gram constructions have the expected structure", {
  argvals <- seq(-2, 3, length.out = 21)
  setup <- AdaReFun:::.fmcp_basis_setup(argvals, M = 6L, degree = 3L,
                                    penalty_order = 2L)
  expect_equal(rowSums(setup$basis_observed), rep(1, length(argvals)),
               tolerance = 1e-13)
  expect_equal(ncol(setup$basis_observed), 9L)

  roughness <- AdaReFun:::cpp_weighted_crossprod(
    setup$derivative_quad, setup$quad_weights
  )
  expect_equal(roughness, t(roughness), tolerance = 1e-13)
  expect_gte(min(eigen(roughness, symmetric = TRUE, only.values = TRUE)$values),
             -1e-9)
})
