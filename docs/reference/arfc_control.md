# Control adaptive regional folded-concave fitting

`arfc_control()` collects numerical and adaptive-weight controls used by
[`arfc()`](https://byuzbasi.github.io/AdaReFun/reference/arfc.md) and
[`arfc_select()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_select.md).
The defaults reproduce the computational convention used in the
accompanying ARFC study. No training/validation split is created by
these controls.

## Usage

``` r
arfc_control(
  max_iter = 120L,
  tol = 1e-08,
  norm_floor = 1e-08,
  relative_cutoff_multiplier = 1e-04,
  cutoff_floor = 1e-08,
  refit = TRUE,
  trace = FALSE,
  keep_matrices = FALSE,
  adaptive_alpha = 1,
  adaptive_epsilon_fraction = 0.05,
  adaptive_clip = c(0.25, 4),
  zero_tolerance = 1e-08,
  tie_relative_tolerance = 1e-12,
  strict = TRUE
)
```

## Arguments

- max_iter:

  Positive integer giving the maximum LQA iterations for each tuning
  candidate.

- tol:

  Positive relative coefficient convergence tolerance.

- norm_floor:

  Positive stabilizer for regional norms near zero.

- relative_cutoff_multiplier:

  Nonnegative multiplier applied to the largest absolute smooth-pilot
  coefficient. The resulting threshold is used for exact coefficient
  thresholding.

- cutoff_floor:

  Positive lower bound for the absolute threshold.

- refit:

  Logical; refit the roughness-penalized model on retained basis
  coefficients after thresholding.

- trace:

  Logical; retain exact-objective traces for tuning candidates.

- keep_matrices:

  Logical; retain design and regional matrices in fitted objects. This
  can substantially increase object size.

- adaptive_alpha:

  Positive exponent in the pilot-adaptive regional weights.

- adaptive_epsilon_fraction:

  Positive fraction of the median positive pilot norm used in the
  stabilizing offset.

- adaptive_clip:

  Two positive increasing numbers bounding normalized adaptive weights.

- zero_tolerance:

  Nonnegative tolerance used only to display estimated zero regions in
  plots.

- tie_relative_tolerance:

  Nonnegative relative tolerance for validation loss ties in
  [`arfc_select()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_select.md).

- strict:

  Logical. If `TRUE`, only successful, converged, successfully refitted
  candidates with finite coefficients and BIC are eligible.

## Value

An object of class `arfc_control`.
