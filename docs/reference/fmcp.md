# Fit functional MCP (fMCP) regression

`fmcp()` fits the scalar-on-function linear model \$\$Y_i = \mu + \int
X_i(t)\beta(t)\\dt + \varepsilon_i\$\$ with a B-spline representation of
\\\beta(t)\\, a derivative roughness penalty, and an interval-wise
functional minimax concave penalty (fMCP).

## Usage

``` r
fmcp(
  y,
  x,
  argvals,
  M,
  lambda,
  kappa,
  degree = 3L,
  penalty_order = 2L,
  mcp_gamma = 3,
  integration = c("simpson", "trapezoid"),
  control = fmcp_control()
)

fMCP(...)

FMCP(...)
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

  Nonempty nonnegative grid for the fMCP shrinkage parameter.

- kappa:

  Nonempty nonnegative grid for the derivative roughness penalty.

- degree:

  B-spline polynomial degree. Cubic splines (`degree = 3`) are used by
  default.

- penalty_order:

  Positive derivative order used in the roughness penalty. It cannot
  exceed `degree`.

- mcp_gamma:

  One or more MCP concavity parameters, all strictly greater than one.
  If several values are supplied, they are included in BIC model
  selection.

- integration:

  Numerical integration rule for the observed functional predictor.
  Simpson integration requires an equally spaced grid with an odd number
  of points; trapezoid integration accepts irregular grids.

- control:

  Numerical controls from
  [`fmcp_control()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp_control.md).

- ...:

  Arguments forwarded by the `fMCP()` and `FMCP()` aliases to `fmcp()`.

## Value

An object of class `fmcp`. It contains the selected basis coefficients,
coefficient function on `argvals`, fitted values, residuals, interval
norms, complete tuning table, and convergence diagnostics.

## Details

For interval \\j\\, let \\W_j\\ be its local B-spline Gram matrix and
\$\$q_j = \\(M/T)b^T W_j b\\^{1/2}.\$\$ The fitted objective is
\$\$n^{-1}\\y-Ub\\\_2^2 + \kappa b^T Vb + \sum_j
p\_{\lambda,a}(q_j),\$\$ where \$\$p\_{\lambda,a}(u) = \lambda
u-u^2/(2a)\$\$ for \\0 \le u \le a\lambda\\, and
\\p\_{\lambda,a}(u)=a\lambda^2/2\\ otherwise. The nonconvex objective is
optimized with a local quadratic approximation initialized from the
roughness-only fit. No predictor or response standardization is
performed; the response and functional predictor are centered internally
to estimate an intercept.

## References

Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally Sparse
Estimator for Functional Linear Regression Models. *Journal of
Computational and Graphical Statistics*, 26(2), 306-318.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *The Annals of Statistics*, 38(2), 894-942.

## Examples

``` r
set.seed(42)
n <- 40
argvals <- seq(0, 1, length.out = 21)
scores <- matrix(rnorm(n * 6), n, 6)
x <- matrix(0, n, length(argvals))
for (j in seq_len(6)) {
  x <- x + scores[, j] %o% sin(j * pi * argvals)
}
x <- x + matrix(rnorm(n * length(argvals), sd = 0.05), n)
beta <- ifelse(argvals <= 0.55, 2 * sin(2 * pi * argvals), 0)
y <- rowMeans(x * rep(beta, each = n)) + rnorm(n, sd = 0.1)
fit <- fmcp(
  y, x, argvals,
  M = 8,
  lambda = c(0.005, 0.02),
  kappa = c(1e-5, 1e-4),
  control = fmcp_control(max_iter = 50)
)
coef(fit, type = "basis")
predict(fit, x[1:3, , drop = FALSE])
```
