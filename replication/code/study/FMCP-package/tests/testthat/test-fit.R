test_that("lambda zero equals the roughness-only B-spline solution", {
  dat <- make_fmcp_test_data()
  fit <- fmcp(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = 0,
    kappa = 1e-3,
    control = fmcp_control(keep_matrices = TRUE)
  )
  internal <- fit$internal
  n <- length(dat$y)
  expected <- solve(
    crossprod(internal$design) + n * 1e-3 * internal$roughness,
    crossprod(internal$design, internal$centered_y)
  )
  expect_equal(unname(coef(fit)), as.numeric(expected), tolerance = 1e-9)
})

test_that("Rcpp LQA matches an independent R implementation", {
  dat <- make_fmcp_test_data(n = 60L, seed = 99L)
  controls <- fmcp_control(
    max_iter = 150L,
    tol = 1e-10,
    norm_floor = 1e-8,
    cutoff = 0,
    refit = FALSE,
    keep_matrices = TRUE
  )
  fit <- fmcp(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = 0.01,
    kappa = 1e-3,
    mcp_gamma = 3,
    control = controls
  )
  expected <- reference_fmcp_lqa(
    fit$internal$centered_y,
    fit$internal$design,
    fit$internal$roughness,
    fit$internal$local_gram,
    lambda = 0.01,
    kappa = 1e-3,
    mcp_gamma = 3,
    scale = fit$basis$interval_scale,
    max_iter = controls$max_iter,
    tolerance = controls$tol,
    norm_floor = controls$norm_floor
  )
  expect_equal(unname(coef(fit)), expected, tolerance = 1e-8)
})

test_that("a centered zero signal produces an exact zero coefficient function", {
  set.seed(8)
  n <- 35L
  argvals <- seq(0, 1, length.out = 21)
  x <- matrix(rnorm(n * length(argvals)), n, length(argvals))
  fit <- fmcp(
    rep(2.5, n), x, argvals,
    M = 6,
    lambda = 1,
    kappa = 1e-3
  )
  expect_true(all(coef(fit) == 0))
  expect_true(all(fit$beta == 0))
  expect_true(all(fit$interval_norms == 0))
  expect_equal(fitted(fit), rep(2.5, n), tolerance = 1e-14)
})

test_that("a mixed signal can produce exact local null intervals", {
  set.seed(314)
  n <- 100L
  argvals <- seq(0, 1, length.out = 51)
  scores <- matrix(rnorm(n * 12L), n, 12L)
  x <- matrix(0, n, length(argvals))
  for (j in seq_len(12L)) {
    x <- x + scores[, j] %o% (sqrt(2) * sin(j * pi * argvals) / j)
  }
  x <- x + matrix(rnorm(n * length(argvals), sd = 0.01), n)
  beta <- ifelse(argvals <= 0.5, 3 * sin(2 * pi * argvals), 0)
  weights <- FMCP:::.fmcp_integration_weights(argvals, "simpson")
  y <- as.numeric(x %*% (weights * beta)) + rnorm(n, sd = 0.08)

  fit <- fmcp(
    y, x, argvals,
    M = 14,
    lambda = 0.03,
    kappa = 1e-5,
    mcp_gamma = 3
  )
  tolerance <- max(fit$control$norm_floor, sqrt(.Machine$double.eps))
  zero_intervals <- fit$interval_norms <= tolerance
  expect_true(any(zero_intervals))
  expect_false(all(zero_intervals))
  expect_true(all(tail(zero_intervals, 3L)))
})

test_that("public methods and tuning-grid bookkeeping are coherent", {
  dat <- make_fmcp_test_data(n = 55L, seed = 101L)
  fit <- fmcp(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = c(0.005, 0.02),
    kappa = c(1e-4, 1e-3),
    mcp_gamma = c(2, 3),
    control = fmcp_control(max_iter = 100L)
  )

  expect_s3_class(fit, "fmcp")
  expect_equal(nrow(fit$tuning), 8L)
  expect_equal(nobs(fit), 55L)
  expect_length(coef(fit, type = "basis"), 11L)
  expect_named(coef(fit, type = "function"), c("argvals", "beta"))
  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, dat$x), fitted(fit), tolerance = 1e-12)
  expect_equal(residuals(fit), dat$y - fitted(fit), tolerance = 1e-14)
  expect_s3_class(summary(fit), "summary.fmcp")
})

test_that("invalid preprocessing choices fail explicitly", {
  dat <- make_fmcp_test_data(n = 20L)
  expect_error(
    fmcp(dat$y, dat$x[, -1L], dat$argvals[-1L], M = 6,
         lambda = 0.1, kappa = 0.01),
    "odd number"
  )
  expect_error(
    fmcp(dat$y, dat$x, dat$argvals, M = 6,
         lambda = 0.1, kappa = 0.01, mcp_gamma = 1),
    "greater than"
  )
  bad_x <- dat$x
  bad_x[1, 1] <- NA_real_
  expect_error(
    fmcp(dat$y, bad_x, dat$argvals, M = 6,
         lambda = 0.1, kappa = 0.01),
    "finite numeric matrix"
  )
  expect_error(
    predict(
      fmcp(dat$y, dat$x, dat$argvals, M = 6,
           lambda = 0.1, kappa = 0.01),
      dat$x[, -1L, drop = FALSE]
    ),
    "same grid columns"
  )
})
