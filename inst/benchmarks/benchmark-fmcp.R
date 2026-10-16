# Run after installing AdaReFun. This script writes no files.
library(AdaReFun)

set.seed(2026)
n <- 250L
p <- 101L
argvals <- seq(0, 1, length.out = p)
x <- matrix(rnorm(n * p), n, p)
basis <- AdaReFun:::.fmcp_basis_setup(
  argvals, M = 30L, degree = 3L, penalty_order = 2L
)
weights <- AdaReFun:::.fmcp_integration_weights(argvals, "simpson")

r_design_time <- system.time({
  for (iteration in seq_len(100L)) {
    design_r <- x %*% (weights * basis$basis_observed)
  }
})

cpp_design_time <- system.time({
  for (iteration in seq_len(100L)) {
    design_cpp <- AdaReFun:::cpp_design_matrix(
      x, basis$basis_observed, weights
    )
  }
})

stopifnot(isTRUE(all.equal(design_r, design_cpp, tolerance = 1e-12)))
print(rbind(R = r_design_time, RcppArmadillo = cpp_design_time))

beta <- ifelse(argvals <= 0.55, 2 * sin(2 * pi * argvals), 0)
y <- as.numeric(x %*% (weights * beta)) + rnorm(n, sd = 0.2)
fit_time <- system.time({
  fit <- fmcp(
    y, x, argvals,
    M = 30L,
    lambda = c(0.001, 0.003, 0.01, 0.03),
    kappa = 10^seq(-6, -3, length.out = 4L),
    mcp_gamma = c(2, 3),
    control = fmcp_control(max_iter = 100L)
  )
})
print(fit_time)
print(fit)
