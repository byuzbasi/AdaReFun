make_arfc_test_data <- function(n = 48L, p = 21L, seed = 771L) {
  set.seed(seed)
  argvals <- seq(0, 1, length.out = p)
  scores <- matrix(rnorm(n * 7L), n, 7L)
  x <- matrix(0, n, p)
  for (index in seq_len(7L)) {
    x <- x + scores[, index] %o%
      (sqrt(2) * sin(index * pi * argvals) / index)
  }
  beta <- ifelse(
    argvals >= 0.35 & argvals <= 0.65,
    0,
    1.5 * sin(2 * pi * argvals)
  )
  weights <- AdaReFun:::.fmcp_integration_weights(argvals, "simpson")
  y <- 0.4 + as.numeric(x %*% (weights * beta)) + rnorm(n, sd = 0.12)
  list(y = y, x = x, argvals = argvals, beta = beta)
}

test_that("regional MCP and SCAD scalar penalties are analytic", {
  values <- c(0, 0.5, 1.5, 3, 4)
  mcp <- AdaReFun:::cpp_arfc_penalty_scalar_v1(values, 1, "mcp", 3)
  expect_equal(
    mcp$penalty,
    ifelse(values <= 3, values - values^2 / 6, 1.5),
    tolerance = 1e-14
  )
  expect_equal(mcp$derivative, pmax(1 - values / 3, 0), tolerance = 1e-14)

  scad <- AdaReFun:::cpp_arfc_penalty_scalar_v1(values, 1, "scad", 3.7)
  expected_derivative <- ifelse(
    values <= 1, 1,
    ifelse(values <= 3.7, (3.7 - values) / 2.7, 0)
  )
  expect_equal(scad$derivative, expected_derivative, tolerance = 1e-14)
  expect_equal(1 / 2.7, 1 / (3.7 - 1), tolerance = 1e-14)
})

test_that("canonical single and multi regional mass matrices agree", {
  dat <- make_arfc_test_data(n = 24L)
  context <- AdaReFun:::.arfc_prepare_context(
    dat$y, dat$x, dat$argvals,
    levels = c(4L, 8L), level_weights = c(0.5, 0.5),
    degree = 3L, penalty_order = 2L, integration = "simpson"
  )
  single <- AdaReFun:::cpp_arfc_curvature_envelope_v1(
    context$fine_gram, context$fine_length, context$fine_multiplier,
    "mcp", 2.7
  )
  multi <- AdaReFun:::cpp_arfc_curvature_envelope_v1(
    context$region_tree$gram,
    context$region_tree$metadata$region_length,
    context$multi_multiplier,
    "mcp", 2.7
  )
  scad <- AdaReFun:::cpp_arfc_curvature_envelope_v1(
    context$region_tree$gram,
    context$region_tree$metadata$region_length,
    context$multi_multiplier,
    "scad", 3.7
  )
  expect_equal(single$mass_matrix, multi$mass_matrix, tolerance = 1e-11)
  expect_equal(
    multi$curvature_envelope,
    scad$curvature_envelope,
    tolerance = 1e-11
  )
  expect_equal(sum(context$fine_multiplier), 1, tolerance = 1e-14)
  expect_equal(sum(context$multi_multiplier), 1, tolerance = 1e-14)
})

test_that("all four ARFC specifications fit through one public API", {
  dat <- make_arfc_test_data()
  specifications <- expand.grid(
    resolution = c("single", "multi"),
    penalty = c("MCP", "SCAD"),
    stringsAsFactors = FALSE
  )
  fits <- lapply(seq_len(nrow(specifications)), function(index) {
    arfc(
      dat$y, dat$x, dat$argvals,
      lambda = c(0.002, 0.01),
      kappa = c(1e-5, 1e-3),
      resolution = specifications$resolution[index],
      penalty = specifications$penalty[index],
      levels = c(4L, 8L),
      control = arfc_control(max_iter = 80L)
    )
  })
  expect_true(all(vapply(fits, inherits, logical(1L), "arfc")))
  expect_equal(
    vapply(fits, function(value) nrow(value$tuning), integer(1L)),
    rep.int(2L, 4L)
  )
  expect_true(all(vapply(
    fits, function(value) isTRUE(value$selected$converged), logical(1L)
  )))
  expect_equal(predict(fits[[1L]], dat$x), fitted(fits[[1L]]),
               tolerance = 1e-12)
  expect_s3_class(summary(fits[[1L]]), "summary.arfc")
})

test_that("ARFC-Select uses explicit validation without fitting on it", {
  dat <- make_arfc_test_data(n = 56L, seed = 812L)
  selected <- arfc_select(
    dat$y[1:40], dat$x[1:40, , drop = FALSE], dat$argvals,
    validation = list(
      y = dat$y[41:56],
      x = dat$x[41:56, , drop = FALSE]
    ),
    lambda = c(0.002, 0.01),
    kappa = c(1e-5, 1e-3),
    levels = c(4L, 8L),
    control = arfc_control(max_iter = 80L)
  )
  expect_s3_class(selected, "arfc_select")
  expect_equal(nrow(selected$selection), 4L)
  expect_equal(sum(selected$selection$selected), 1L)
  expect_false(selected$diagnostics$validation_used_for_fitting)
  expect_equal(selected$diagnostics$validation_sample_size, 16L)
  expect_equal(nobs(selected), 40L)
  expect_equal(coef(selected), coef(selected$selected_fit))
  expect_equal(predict(selected), fitted(selected$selected_fit))
  expect_s3_class(summary(selected), "summary.arfc_select")
})

test_that("ARFC preprocessing and selection choices fail explicitly", {
  dat <- make_arfc_test_data(n = 24L)
  expect_error(
    arfc(
      dat$y, dat$x, dat$argvals,
      lambda = 0.01, kappa = 1e-3,
      levels = c(5L, 8L)
    ),
    "divide"
  )
  expect_error(
    arfc(
      dat$y, dat$x, dat$argvals,
      lambda = 0.01, kappa = 1e-3,
      levels = c(4L, 8L), level_weights = c(1, 1)
    ),
    "sum to one"
  )
  expect_error(
    arfc_select(
      dat$y, dat$x, dat$argvals,
      validation = NULL,
      lambda = 0.01, kappa = 1e-3,
      levels = c(4L, 8L)
    ),
    "named list"
  )
})

test_that("ARFC plotting methods render all documented views", {
  dat <- make_arfc_test_data(n = 40L, seed = 945L)
  selected <- arfc_select(
    dat$y[1:28], dat$x[1:28, , drop = FALSE], dat$argvals,
    validation = list(
      y = dat$y[29:40],
      x = dat$x[29:40, , drop = FALSE]
    ),
    lambda = c(0.002, 0.01),
    kappa = c(1e-5, 1e-3),
    levels = c(4L, 8L),
    control = arfc_control(max_iter = 80L)
  )
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(output)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(output)
  }, add = TRUE)
  expect_silent(plot(selected$selected_fit, type = "coefficient"))
  expect_silent(plot(selected$selected_fit, type = "regions"))
  expect_silent(plot(selected$selected_fit, type = "weights"))
  expect_silent(plot(selected$selected_fit, type = "tuning"))
  expect_silent(plot(selected, type = "validation"))
  expect_silent(plot(selected, type = "coefficient"))
  expect_silent(plot(selected, type = "candidates"))
  grDevices::dev.off()
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
})
