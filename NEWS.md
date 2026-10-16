# AdaReFun 0.1.2

* Adds an English package showcase, original vector identity, local pkgdown
  documentation, worked examples, result gallery and reproduction guidance.
* Records the GitHub address and adds citation metadata and contribution guidance.
* Explains the distinction between independent validation, post-selection refits,
  structural recovery, and conditional prediction-record uncertainty.
* Preserves all 0.1.1 statistical interfaces, defaults, R code and C++ kernels.
  Manuscript v18 corrects regional ISE reporting separately from estimation.

# AdaReFun 0.1.1

* Records Bahadır Yüzbaşı and Jiguo Cao as package authors, with verified
  institutional email addresses and ORCID identifiers.
* Aligns the package citation, overview, README, and vignette with the ARFC
  manuscript terminology and the non-ensemble role of ARFC-Select.
* Documents the review-stage availability policy without claiming a public
  repository or archive DOI before those resources exist.
* Leaves the validated statistical interfaces and RcppArmadillo numerical
  implementation unchanged from 0.1.0.

# AdaReFun 0.1.0

* Defines **AdaReFun** as **Adaptive Regional Functional Regression** and
  distinguishes the ARFC estimator family, a single `arfc()` specification,
  and the non-ensemble four-candidate `arfc_select()` procedure.
* Introduces adaptive regional folded-concave regression through `arfc()` and
  four-candidate validation selection through `arfc_select()`.
* Supports adaptive single- and multi-resolution MCP and SCAD penalties with a
  shared RcppArmadillo numerical core.
* Adds coefficient, regional-norm, adaptive-weight, tuning-path, and candidate
  comparison plots.
* Retains the documented `fmcp()`/`fMCP()` and `fscad()`/`fSCAD()` interfaces
  from the earlier research package.

# FMCP 0.2.0 (legacy source)

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
