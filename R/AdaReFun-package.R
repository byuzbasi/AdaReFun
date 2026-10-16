#' AdaReFun: Adaptive Regional Functional Regression
#'
#' The package name **AdaReFun** is a mnemonic for **Adaptive Regional
#' Functional Regression**. It fits smooth and locally sparse
#' scalar-on-function linear regression models using functional MCP (fMCP),
#' functional SCAD (fSCAD), and the adaptive regional folded-concave (ARFC)
#' family. One ARFC specification fixes a resolution structure (adaptive
#' single resolution, ASR, or adaptive multi resolution, AMR) and a penalty
#' family (MCP or SCAD). [arfc()] fits one such specification, whereas
#' [arfc_select()] fits the four ASR-MCP, ASR-SCAD, AMR-MCP, and AMR-SCAD
#' candidates and selects exactly one with an explicitly supplied validation
#' sample; it is not an ensemble. B-spline basis construction and user-facing
#' methods are implemented in R; matrix-intensive fitting is implemented in
#' C++ with RcppArmadillo.
#'
#' The package is developed by Bahadır Yüzbaşı and Jiguo Cao. Bahadır Yüzbaşı
#' is the maintainer and corresponding author for the accompanying work.
#' During peer review, the versioned package source and replication archive
#' are available from the corresponding author. A public repository URL and
#' archival DOI will be added only after those resources are released.
#'
#' @keywords internal
#' @useDynLib AdaReFun, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics plot
#' @importFrom stats coef fitted nobs predict residuals
"_PACKAGE"
