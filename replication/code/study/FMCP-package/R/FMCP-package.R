#' FMCP: Functional MCP and SCAD Regression
#'
#' FMCP fits smooth and locally sparse scalar-on-function linear regression
#' models using functional MCP (fMCP) or functional SCAD (fSCAD). B-spline
#' basis construction and user-facing model methods are implemented in R;
#' matrix-intensive fitting is implemented in C++ with RcppArmadillo.
#'
#' @keywords internal
#' @useDynLib FMCP, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics plot
#' @importFrom stats coef fitted nobs predict residuals
"_PACKAGE"
