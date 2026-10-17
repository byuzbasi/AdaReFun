# Select the ARFC penalty family and resolution by validation loss

`arfc_select()` fits the four prespecified candidates ASR-MCP, ASR-SCAD,
AMR-MCP, and AMR-SCAD to the same training data, pilot, tuning grids,
and preprocessing. It selects the smallest response MSE on an explicitly
supplied validation sample. The validation responses are never used to
fit a candidate. Numerical ties follow the documented order above.

## Usage

``` r
arfc_select(
  y,
  x,
  argvals,
  validation,
  lambda,
  kappa,
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

- validation:

  Named list with finite `y` and numeric matrix `x` on the same
  functional grid as the training data. A split is never created
  silently.

- lambda:

  Nonempty nonnegative folded-concave tuning grid.

- kappa:

  Nonempty nonnegative roughness grid for the common smooth pilot. The
  BIC-selected value is used for the ARFC path.

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

An object of class `arfc_select`. It contains one `selected_fit`, all
four candidate fits, validation losses, tie diagnostics, and an explicit
no-leakage record. The candidates are not averaged or ensembled. A final
combined-data refit is not performed automatically.

## Examples

``` r
set.seed(2027)
argvals <- seq(0, 1, length.out = 21)
x <- matrix(rnorm(48 * length(argvals)), 48)
beta <- ifelse(argvals >= 0.35 & argvals <= 0.65, 0,
               sin(2 * pi * argvals))
w <- c(0.025, rep(0.05, 19), 0.025)
y <- as.numeric(x %*% (w * beta)) + rnorm(48, sd = 0.2)
selected <- arfc_select(
  y[1:32], x[1:32, , drop = FALSE], argvals,
  validation = list(y = y[33:48], x = x[33:48, , drop = FALSE]),
  lambda = c(0.002, 0.02),
  kappa = c(1e-5, 1e-3),
  levels = c(4L, 8L),
  control = arfc_control(max_iter = 40L)
)
selected$selection
```
