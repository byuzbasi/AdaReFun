# Fit adaptive regional folded-concave functional regression

`arfc()` fits a scalar-on-function linear model with pilot-adaptive
regional MCP or SCAD penalties. A single-resolution fit uses the finest
partition; a multi-resolution fit uses a nested collection of partitions
in the same B-spline coefficient space. Regional multipliers are
normalized by domain mass, and level weights must sum to one.

## Usage

``` r
arfc(
  y,
  x,
  argvals,
  lambda,
  kappa,
  resolution = c("single", "multi"),
  penalty = c("MCP", "SCAD"),
  levels = c(8L, 16L, 32L),
  level_weights = NULL,
  adaptive = TRUE,
  degree = 3L,
  penalty_order = 2L,
  mcp_gamma = 2.7,
  scad_a = 3.7,
  integration = c("simpson", "trapezoid"),
  control = arfc_control()
)
```

## Arguments

- y:

  Finite numeric response vector.

- x:

  Finite numeric matrix with subjects in rows and functional grid points
  in columns.

- argvals:

  Strictly increasing functional argument values.

- lambda:

  Nonempty nonnegative folded-concave tuning grid.

- kappa:

  Nonempty nonnegative roughness grid for the common smooth pilot. The
  BIC-selected value is used for the ARFC path.

- resolution:

  Either `"single"` (ASR) or `"multi"` (AMR).

- penalty:

  Either `"MCP"` or `"SCAD"`.

- levels:

  Strictly increasing nested region counts. Every value must divide the
  largest value. The default reproduces the study specification.

- level_weights:

  Positive weights, one per `levels`, summing to one. `NULL` gives equal
  weights. They are not silently normalized.

- adaptive:

  Logical. Use smooth-pilot adaptive regional weights. Setting this to
  `FALSE` gives the corresponding nonadaptive regional fit.

- degree:

  B-spline polynomial degree.

- penalty_order:

  Derivative order for roughness penalization.

- mcp_gamma:

  MCP shape parameter greater than one.

- scad_a:

  SCAD shape parameter greater than two. The defaults `mcp_gamma = 2.7`
  and `scad_a = 3.7` have equal scalar maximum concavity.

- integration:

  Numerical integration rule for the observed predictor.

- control:

  Numerical and adaptive-weight settings from
  [`arfc_control()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_control.md).

## Value

An object of class `arfc` containing the training-BIC-tuned fit for the
one requested penalty-resolution specification: its coefficient
function, regional norms and weights, complete tuning diagnostics,
smooth pilot, fitted values, and prediction metadata. This function does
not choose between ASR and AMR or between MCP and SCAD; use
[`arfc_select()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_select.md)
for that four-candidate comparison.

## Details

The roughness-only pilot chooses one `kappa` value by strict BIC. Its
regional norms define level-normalized and clipped adaptive weights. The
folded-concave path is then evaluated over `lambda` at that common
smoothness value. No response or predictor standardization is performed;
centering estimates the intercept.

## References

Fan, J., and Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. *Journal of the American
Statistical Association*, 96(456), 1348-1360.

Zhang, C.-H. (2010). Nearly unbiased variable selection under minimax
concave penalty. *The Annals of Statistics*, 38(2), 894-942.

Lin, Z., Cao, J., Wang, L., and Wang, H. (2017). Locally sparse
estimator for functional linear regression models. *Journal of
Computational and Graphical Statistics*, 26(2), 306-318.

## Examples

``` r
set.seed(2026)
n <- 36
argvals <- seq(0, 1, length.out = 21)
x <- matrix(rnorm(n * length(argvals)), n)
beta <- ifelse(argvals < 0.35, 1.5 * sin(2 * pi * argvals), 0)
w <- c(0.025, rep(0.05, 19), 0.025)
y <- as.numeric(x %*% (w * beta)) + rnorm(n, sd = 0.2)
fit <- arfc(
  y, x, argvals,
  lambda = c(0.002, 0.02),
  kappa = c(1e-5, 1e-3),
  resolution = "single",
  penalty = "MCP",
  levels = c(4L, 8L),
  control = arfc_control(max_iter = 40L)
)
coef(fit, type = "basis")
predict(fit, x[1:2, , drop = FALSE])
```
