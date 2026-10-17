# Changelog

## AdaReFun 0.1.2

- Adds an English package showcase, original vector identity, local
  pkgdown documentation, worked examples, result gallery and
  reproduction guidance.
- Records the GitHub address and adds citation metadata and contribution
  guidance.
- Explains the distinction between independent validation,
  post-selection refits, structural recovery, and conditional
  prediction-record uncertainty.
- Preserves all 0.1.1 statistical interfaces, defaults, R code and C++
  kernels. Manuscript v18 corrects regional ISE reporting separately
  from estimation.

## AdaReFun 0.1.1

- Records Bahadır Yüzbaşı and Jiguo Cao as package authors, with
  verified institutional email addresses and ORCID identifiers.
- Aligns the package citation, overview, README, and vignette with the
  ARFC manuscript terminology and the non-ensemble role of ARFC-Select.
- Documents the review-stage availability policy without claiming a
  public repository or archive DOI before those resources exist.
- Leaves the validated statistical interfaces and RcppArmadillo
  numerical implementation unchanged from 0.1.0.

## AdaReFun 0.1.0

- Defines **AdaReFun** as **Adaptive Regional Functional Regression**
  and distinguishes the ARFC estimator family, a single
  [`arfc()`](https://byuzbasi.github.io/AdaReFun/reference/arfc.md)
  specification, and the non-ensemble four-candidate
  [`arfc_select()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_select.md)
  procedure.
- Introduces adaptive regional folded-concave regression through
  [`arfc()`](https://byuzbasi.github.io/AdaReFun/reference/arfc.md) and
  four-candidate validation selection through
  [`arfc_select()`](https://byuzbasi.github.io/AdaReFun/reference/arfc_select.md).
- Supports adaptive single- and multi-resolution MCP and SCAD penalties
  with a shared RcppArmadillo numerical core.
- Adds coefficient, regional-norm, adaptive-weight, tuning-path, and
  candidate comparison plots.
- Retains the documented
  [`fmcp()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp.md)/[`fMCP()`](https://byuzbasi.github.io/AdaReFun/reference/fmcp.md)
  and
  [`fscad()`](https://byuzbasi.github.io/AdaReFun/reference/fscad.md)/[`fSCAD()`](https://byuzbasi.github.io/AdaReFun/reference/fscad.md)
  interfaces from the earlier research package.
