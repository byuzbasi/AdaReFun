# FMCP

`FMCP` fits smooth and locally sparse scalar-on-function linear regression
models with functional minimax concave penalty (fMCP) and functional smoothly
clipped absolute deviation (fSCAD) penalties. The package name remains `FMCP`
for backward compatibility. The implementation is independent: it does not
modify or import code from `PFLR`.

The model is

\[
Y_i = \mu + \int X_i(t)\beta(t)\,dt + \varepsilon_i,
\]

with a B-spline representation \(\beta(t)=B(t)^\top b\). For each of \(M\)
B-spline intervals, the functional magnitude is

\[
q_j = \sqrt{\frac{M}{T}b^\top W_jb}.
\]

The fitted objective combines mean squared error, a derivative roughness
penalty, and either MCP or SCAD applied to the interval norms:

\[
\frac{1}{n}\lVert y-Ub\rVert_2^2 + \kappa b^\top Vb +
\sum_{j=1}^M p_\lambda(q_j).
\]

The matrix-intensive parts of integration, Gram-matrix construction, local
quadratic approximation, linear solving, active-set refitting, effective
degrees of freedom, and BIC evaluation are implemented in C++17 with
RcppArmadillo.

## Installation

From the package source directory:

```r
install.packages(".", repos = NULL, type = "source")
```

## Example

```r
library(FMCP)

set.seed(42)
n <- 60
argvals <- seq(0, 1, length.out = 51)
x <- matrix(rnorm(n * length(argvals)), n, length(argvals))
beta <- ifelse(argvals <= 0.6, 1.5 * sin(2 * pi * argvals), 0)
weights <- c(1, rep(c(4, 2), length.out = 49), 1) / (3 * 50)
y <- 0.5 + as.numeric(x %*% (weights * beta)) + rnorm(n, sd = 0.1)

fit_mcp <- fMCP(
  y = y,
  x = x,
  argvals = argvals,
  M = 20,
  lambda = c(0.002, 0.01, 0.05),
  kappa = 10^seq(-5, -3),
  mcp_gamma = 3
)

fit_scad <- fSCAD(
  y = y,
  x = x,
  argvals = argvals,
  M = 20,
  lambda = c(0.002, 0.01, 0.05),
  kappa = 10^seq(-5, -3),
  scad_a = 3.7
)

fit_mcp
summary(fit_scad)
plot(fit_mcp)
predict(fit_scad, x[1:5, , drop = FALSE])
```

## Explicit analysis choices

- `lambda`, `kappa`, and `M` are required. The package does not silently pick
  scale-dependent tuning grids.
- The response and functional predictor are centered to fit an intercept, but
  they are not standardized.
- Missing values are rejected; no observations are removed or imputed.
- Simpson integration is the default and requires an equally spaced grid with
  an odd number of points. Set `integration = "trapezoid"` explicitly for an
  irregular common grid.
- Each fMCP and fSCAD tuning candidate starts from the same roughness-only
  solution for its value of `kappa`. This avoids path-order dependence in the
  nonconvex fit.
- Failed and nonconverged candidates remain visible in the returned tuning
  table. No hidden diagonal jitter is added to singular systems.

MCP and SCAD are nonconvex, so reported solutions are local solutions
associated with the documented initialization. `fmcp()` and `FMCP()` remain
available as backward-compatible aliases; new analyses may use the scientific
label `fMCP()` alongside `fSCAD()`.

## References

- Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally Sparse Estimator for
  Functional Linear Regression Models. *Journal of Computational and Graphical
  Statistics*, 26(2), 306-318.
- Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized
  likelihood and its oracle properties. *Journal of the American Statistical
  Association*, 96(456), 1348-1360.
- Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax concave
  penalty. *The Annals of Statistics*, 38(2), 894-942.
