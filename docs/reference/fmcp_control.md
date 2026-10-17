# Control numerical details of fMCP and fSCAD fitting

Control numerical details of fMCP and fSCAD fitting

## Usage

``` r
fmcp_control(
  max_iter = 100L,
  tol = 1e-08,
  norm_floor = 1e-08,
  cutoff = 1e-06,
  refit = TRUE,
  trace = FALSE,
  keep_matrices = FALSE
)
```

## Arguments

- max_iter:

  Positive integer. Maximum local quadratic approximation iterations for
  each tuning-parameter combination.

- tol:

  Positive number. Relative coefficient convergence tolerance.

- norm_floor:

  Positive number used only to stabilize division by an interval norm
  near zero.

- cutoff:

  Nonnegative absolute threshold applied to B-spline coefficients after
  optimization. Coefficients below this threshold are set exactly to
  zero.

- refit:

  Logical. If `TRUE`, refit the roughness-penalized model on the
  selected B-spline coefficients after thresholding.

- trace:

  Logical. If `TRUE`, retain per-candidate objective traces.

- keep_matrices:

  Logical. If `TRUE`, retain design and penalty matrices in the returned
  fit. This is useful for diagnostics and testing but may substantially
  increase object size.

## Value

An object of class `fmcp_control`.
