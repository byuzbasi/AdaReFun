# FMCP 0.2.0

- Added `fscad()` and `fSCAD()` with a separate C++17/RcppArmadillo fSCAD
  fitting, active-set refitting, effective-degrees-of-freedom, and BIC core.
- Added `fMCP()` as the preferred scientific-label alias while retaining
  `fmcp()` and `FMCP()` for backward compatibility.
- Added independent `fscad` fit and summary classes with coefficient,
  prediction, residual, observation-count, print, summary, and plot methods.
- Added analytic SCAD, independent LQA, API, and backward-compatibility tests;
  verified numerical parity with the controlled simulation comparator.

# FMCP 0.1.0

- Added scalar-on-function regression with interval-wise functional MCP.
- Added a C++17/RcppArmadillo fitting and BIC-selection core.
- Added Simpson and trapezoid integration on a common observation grid.
- Added active-set refitting, convergence diagnostics, and complete candidate
  bookkeeping.
- Added `coef()`, `fitted()`, `residuals()`, `predict()`, `print()`,
  `summary()`, `nobs()`, and `plot()` methods.
- Added analytic, reference-implementation, zero-signal, prediction, and input
  validation tests.
