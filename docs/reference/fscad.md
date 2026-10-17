# Fit functional SCAD regression

`fscad()` fits the scalar-on-function linear model \$\$Y_i = \mu + \int
X_i(t)\beta(t)\\dt + \varepsilon_i\$\$ with a B-spline representation of
\\\beta(t)\\, a derivative roughness penalty, and an interval-wise
functional smoothly clipped absolute deviation penalty (fSCAD).

## Usage

``` r
fscad(
  y,
  x,
  argvals,
  M,
  lambda,
  kappa,
  degree = 3L,
  penalty_order = 2L,
  scad_a = 3.7,
  integration = c("simpson", "trapezoid"),
  control = fmcp_control()
)

fSCAD(...)
```

## Arguments

- y:

  Finite numeric response vector of length `n`.

- x:

  Finite numeric `n` by `p` matrix. Rows are subjects and columns are
  evaluations of a functional predictor on the common `argvals` grid.

- argvals:

  Strictly increasing finite numeric vector of length `p`.

- M:

  Positive integer giving the number of equal B-spline subintervals.

- lambda:

  Nonempty nonnegative grid for the fSCAD shrinkage parameter.

- kappa:

  Nonempty nonnegative grid for the derivative roughness penalty.

- degree:

  B-spline polynomial degree. Cubic splines (`degree = 3`) are used by
  default.

- penalty_order:

  Positive derivative order used in the roughness penalty. It cannot
  exceed `degree`.

- scad_a:

  SCAD shape parameter strictly greater than two.

- integration:

  Numerical integration rule for the observed functional predictor.
  Simpson integration requires an equally spaced grid with an odd number
  of points; trapezoid integration accepts irregular grids.

- control:

  Numerical controls from
  [`fmcp_control()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp_control.md).

- ...:

  Arguments forwarded by the `fSCAD()` alias to `fscad()`.

## Value

An object of class `fscad`. It contains the selected basis coefficients,
coefficient function on `argvals`, fitted values, residuals, interval
norms, complete tuning table, and convergence diagnostics.

## Details

For interval \\j\\, let \\W_j\\ be its local B-spline Gram matrix and
\$\$q_j = \\(M/T)b^T W_j b\\^{1/2}.\$\$ The fitted objective is
\$\$n^{-1}\\y-Ub\\\_2^2 + \kappa b^T Vb + \sum_j
p\_{\lambda,a}(q_j),\$\$ where the SCAD derivative for \\u \ge 0\\ is
\$\$p'\_{\lambda,a}(u) = \lambda I(u \le \lambda) +
\frac{(a\lambda-u)\_+}{a-1}I(u \> \lambda), \qquad a\>2.\$\$ The
nonconvex objective is optimized with the same local quadratic
approximation, roughness-only initialization, active-set refit,
effective degrees of freedom, and BIC convention used by
[`fmcp()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp.md).

## References

Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association*, 96(456), 1348-1360.

Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally Sparse
Estimator for Functional Linear Regression Models. *Journal of
Computational and Graphical Statistics*, 26(2), 306-318.

## Examples

``` r
set.seed(42)
n <- 40
argvals <- seq(0, 1, length.out = 21)
x <- matrix(rnorm(n * length(argvals)), n, length(argvals))
beta <- ifelse(argvals <= 0.55, 2 * sin(2 * pi * argvals), 0)
weights <- c(1, rep(c(4, 2), length.out = 19), 1) / (3 * 20)
y <- as.numeric(x %*% (weights * beta)) + rnorm(n, sd = 0.1)
fit <- fscad(
  y, x, argvals,
  M = 8,
  lambda = c(0.005, 0.02),
  kappa = c(1e-5, 1e-4),
  control = fmcp_control(max_iter = 50)
)
coef(fit, type = "basis")
predict(fit, x[1:3, , drop = FALSE])
```
