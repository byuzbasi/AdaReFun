test_that("SCAD penalty, derivative, knots, and concavity are analytic", {
  lambda <- 1.25
  scad_a <- 3.7
  values <- c(-5, -scad_a * lambda, -lambda, 0, lambda,
              2 * lambda, scad_a * lambda, 5)
  magnitudes <- abs(values)
  expected_penalty <- ifelse(
    magnitudes <= lambda,
    lambda * magnitudes,
    ifelse(
      magnitudes < scad_a * lambda,
      (-magnitudes^2 + 2 * scad_a * lambda * magnitudes - lambda^2) /
        (2 * (scad_a - 1)),
      0.5 * (scad_a + 1) * lambda^2
    )
  )
  expected_derivative <- ifelse(
    magnitudes <= lambda,
    lambda,
    ifelse(
      magnitudes < scad_a * lambda,
      (scad_a * lambda - magnitudes) / (scad_a - 1),
      0
    )
  )

  expect_equal(
    as.numeric(AdaReFun:::cpp_scad_penalty(values, lambda, scad_a)),
    expected_penalty,
    tolerance = 1e-14
  )
  expect_equal(
    as.numeric(AdaReFun:::cpp_scad_derivative(values, lambda, scad_a)),
    expected_derivative,
    tolerance = 1e-14
  )

  point <- 2 * lambda
  step <- 1e-6
  numerical <- (
    AdaReFun:::cpp_scad_penalty(point + step, lambda, scad_a) -
      AdaReFun:::cpp_scad_penalty(point - step, lambda, scad_a)
  ) / (2 * step)
  expect_equal(
    as.numeric(numerical),
    as.numeric(AdaReFun:::cpp_scad_derivative(point, lambda, scad_a)),
    tolerance = 1e-8
  )

  middle_left <- AdaReFun:::cpp_scad_derivative(
    c(1.5, 2.5) * lambda, lambda, scad_a
  )
  observed_slope <- diff(middle_left) / lambda
  expect_equal(as.numeric(observed_slope), -1 / (scad_a - 1),
               tolerance = 1e-14)
})

test_that("fSCAD lambda zero equals the shared roughness-only solution", {
  dat <- make_fmcp_test_data()
  control <- fmcp_control(keep_matrices = TRUE)
  fit_scad <- fscad(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = 0,
    kappa = 1e-3,
    control = control
  )
  fit_mcp <- fmcp(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = 0,
    kappa = 1e-3,
    control = control
  )
  internal <- fit_scad$internal
  n <- length(dat$y)
  expected <- solve(
    crossprod(internal$design) + n * 1e-3 * internal$roughness,
    crossprod(internal$design, internal$centered_y)
  )

  expect_equal(unname(coef(fit_scad)), as.numeric(expected), tolerance = 1e-9)
  expect_equal(unname(coef(fit_scad)), unname(coef(fit_mcp)), tolerance = 1e-13)
})

test_that("fSCAD Rcpp LQA matches an independent R implementation", {
  dat <- make_fmcp_test_data(n = 60L, seed = 207L)
  controls <- fmcp_control(
    max_iter = 150L,
    tol = 1e-10,
    norm_floor = 1e-8,
    cutoff = 0,
    refit = FALSE,
    keep_matrices = TRUE
  )
  fit <- fscad(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = 0.01,
    kappa = 1e-3,
    scad_a = 3.7,
    control = controls
  )
  expected <- reference_fscad_lqa(
    fit$internal$centered_y,
    fit$internal$design,
    fit$internal$roughness,
    fit$internal$local_gram,
    lambda = 0.01,
    kappa = 1e-3,
    scad_a = 3.7,
    scale = fit$basis$interval_scale,
    max_iter = controls$max_iter,
    tolerance = controls$tol,
    norm_floor = controls$norm_floor
  )
  expect_equal(unname(coef(fit)), expected, tolerance = 1e-8)
})

test_that("fSCAD returns an exact zero function for a centered zero signal", {
  set.seed(210)
  n <- 35L
  argvals <- seq(0, 1, length.out = 21)
  x <- matrix(rnorm(n * length(argvals)), n, length(argvals))
  fit <- fscad(
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

test_that("fSCAD full grid and public methods are coherent", {
  dat <- make_fmcp_test_data(n = 55L, seed = 208L)
  fit <- fSCAD(
    dat$y, dat$x, dat$argvals,
    M = 8,
    lambda = c(0.005, 0.02),
    kappa = c(1e-4, 1e-3),
    scad_a = 3.7,
    control = fmcp_control(max_iter = 100L, keep_matrices = TRUE)
  )

  expect_s3_class(fit, "fscad")
  expect_false(inherits(fit, "fmcp"))
  expect_equal(nrow(fit$tuning), 4L)
  expect_equal(unique(fit$tuning$scad_a), 3.7)
  expect_equal(nobs(fit), 55L)
  expect_length(coef(fit, type = "basis"), 11L)
  expect_named(coef(fit, type = "function"), c("argvals", "beta"))
  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, dat$x), fitted(fit), tolerance = 1e-12)
  expect_equal(residuals(fit), dat$y - fitted(fit), tolerance = 1e-14)
  expect_s3_class(summary(fit), "summary.fscad")
  expect_match(paste(capture.output(print(fit)), collapse = "\n"), "fSCAD")

  internal <- fit$internal
  grid <- AdaReFun:::cpp_fscad_grid(
    internal$centered_y,
    internal$design,
    internal$roughness,
    internal$local_gram,
    c(0.005, 0.02),
    c(1e-4, 1e-3),
    3.7,
    fit$basis$interval_scale,
    fit$control$max_iter,
    fit$control$tol,
    fit$control$norm_floor,
    fit$control$cutoff,
    fit$control$refit,
    fit$control$trace
  )
  expect_equal(dim(grid$coefficients), c(11L, 4L))
  expect_equal(as.numeric(grid$lambda), fit$tuning$lambda)
  expect_equal(as.numeric(grid$kappa), fit$tuning$kappa)
})

test_that("new scientific aliases preserve legacy fMCP behavior", {
  dat <- make_fmcp_test_data(n = 30L, seed = 209L)
  arguments <- list(
    y = dat$y,
    x = dat$x,
    argvals = dat$argvals,
    M = 6,
    lambda = 0.01,
    kappa = 1e-3,
    mcp_gamma = 3
  )
  legacy_lower <- do.call(fmcp, arguments)
  legacy_upper <- do.call(FMCP, arguments)
  scientific <- do.call(fMCP, arguments)

  expect_equal(coef(scientific), coef(legacy_lower), tolerance = 0)
  expect_equal(scientific$tuning, legacy_lower$tuning, tolerance = 0)
  expect_equal(coef(legacy_upper), coef(legacy_lower), tolerance = 0)

  scad_lower <- fscad(
    dat$y, dat$x, dat$argvals,
    M = 6, lambda = 0.01, kappa = 1e-3
  )
  scad_upper <- fSCAD(
    dat$y, dat$x, dat$argvals,
    M = 6, lambda = 0.01, kappa = 1e-3
  )
  expect_equal(coef(scad_upper), coef(scad_lower), tolerance = 0)
  expect_equal(scad_upper$tuning, scad_lower$tuning, tolerance = 0)
})

test_that("fSCAD rejects invalid shape and preprocessing inputs", {
  dat <- make_fmcp_test_data(n = 20L)
  expect_error(
    fscad(dat$y, dat$x, dat$argvals, M = 6,
          lambda = 0.1, kappa = 0.01, scad_a = 2),
    "greater than 2"
  )
  expect_error(
    fscad(dat$y, dat$x[, -1L], dat$argvals[-1L], M = 6,
          lambda = 0.1, kappa = 0.01),
    "odd number"
  )
})
