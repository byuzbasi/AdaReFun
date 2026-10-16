#' Model coefficients for an fMCP fit
#'
#' @param object An object of class `fmcp`.
#' @param type Return B-spline coefficients (`"basis"`) or coefficient-function
#'   values on the training grid (`"function"`).
#' @param ... Unused.
#'
#' @return A named numeric vector for `type = "basis"`, or a data frame with
#'   `argvals` and `beta` for `type = "function"`.
#' @export
coef.fmcp <- function(object, type = c("basis", "function"), ...) {
  type <- match.arg(type)
  if (type == "basis") {
    return(object$coefficients)
  }
  data.frame(argvals = object$argvals, beta = object$beta)
}

#' Predict from an fMCP fit
#'
#' @param object An object of class `fmcp`.
#' @param newdata Optional finite numeric matrix whose columns correspond
#'   exactly to the training `argvals`. If omitted, training fitted values are
#'   returned.
#' @param ... Unused.
#'
#' @return Numeric vector of predicted responses.
#' @export
predict.fmcp <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    return(object$fitted.values)
  }
  if (!is.matrix(newdata) || !is.numeric(newdata) || anyNA(newdata) ||
      any(!is.finite(newdata))) {
    stop("`newdata` must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(newdata) != object$dimensions[["p"]]) {
    stop("`newdata` must have the same grid columns as the training data.",
         call. = FALSE)
  }
  centered <- sweep(newdata, 2L, object$means$x, FUN = "-")
  design <- cpp_design_matrix(
    centered,
    object$basis$observed,
    object$basis$integration_weights
  )
  as.numeric(object$means$y + design %*% object$coefficients)
}

#' @rdname predict.fmcp
#' @export
fitted.fmcp <- function(object, ...) {
  object$fitted.values
}

#' @rdname predict.fmcp
#' @export
residuals.fmcp <- function(object, ...) {
  object$residuals
}

#' Number of observations in an fMCP fit
#'
#' @param object An object of class `fmcp`.
#' @param ... Unused.
#'
#' @return Integer number of training subjects.
#' @export
nobs.fmcp <- function(object, ...) {
  unname(object$dimensions[["n"]])
}

#' Print an fMCP fit
#'
#' @param x An object of class `fmcp`.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.fmcp <- function(x, ...) {
  selected <- x$selected
  cat("Functional MCP (fMCP) regression fit\n")
  cat("  Subjects:", x$dimensions[["n"]],
      " Grid points:", x$dimensions[["p"]],
      " Basis functions:", x$dimensions[["basis"]], "\n")
  cat("  lambda:", format(selected$lambda),
      " kappa:", format(selected$kappa),
      " MCP gamma:", format(selected$mcp_gamma), "\n")
  cat("  BIC:", format(selected$bic),
      " Active coefficients:", selected$active,
      " Converged:", if (selected$converged) "yes" else "no", "\n")
  invisible(x)
}

#' Summarize an fMCP fit
#'
#' @param object An object of class `fmcp`.
#' @param x An object of class `summary.fmcp` for printing.
#' @param ... Unused.
#'
#' @return An object of class `summary.fmcp`.
#' @export
summary.fmcp <- function(object, ...) {
  selected <- object$selected
  residual_df <- object$dimensions[["n"]] - selected$df
  residual_scale <- if (is.finite(residual_df) && residual_df > 0) {
    sqrt(selected$rss / residual_df)
  } else {
    NA_real_
  }
  null_tolerance <- max(object$control$norm_floor, .Machine$double.eps^0.5)
  structure(
    list(
      call = object$call,
      dimensions = object$dimensions,
      selected = selected,
      residual_scale = residual_scale,
      coefficient_summary = summary(object$beta),
      null_intervals = sum(object$interval_norms <= null_tolerance),
      total_intervals = length(object$interval_norms),
      diagnostics = object$diagnostics
    ),
    class = "summary.fmcp"
  )
}

#' @rdname summary.fmcp
#' @export
print.summary.fmcp <- function(x, ...) {
  cat("Functional MCP (fMCP) regression summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDimensions:\n")
  cat("  Subjects:", x$dimensions[["n"]],
      " Grid points:", x$dimensions[["p"]],
      " Basis functions:", x$dimensions[["basis"]], "\n")
  cat("\nSelected tuning parameters:\n")
  cat("  lambda =", format(x$selected$lambda),
      " kappa =", format(x$selected$kappa),
      " MCP gamma =", format(x$selected$mcp_gamma), "\n")
  cat("  BIC =", format(x$selected$bic),
      " df =", format(x$selected$df),
      " residual scale =", format(x$residual_scale), "\n")
  cat("  Active coefficients:", x$selected$active,
      " LQA iterations:", x$selected$iterations,
      " Converged:", if (x$selected$converged) "yes" else "no", "\n")
  cat("  Null intervals:", x$null_intervals, "of", x$total_intervals, "\n")
  cat("\nCoefficient-function values:\n")
  print(x$coefficient_summary)
  invisible(x)
}

#' Plot an fMCP coefficient function
#'
#' @param x An object of class `fmcp`.
#' @param shade_null Logical. Shade intervals whose fitted functional norm is
#'   numerically zero.
#' @param null_col Color for shaded null intervals.
#' @param xlab,ylab,main,lwd Standard plotting parameters.
#' @param ... Additional arguments passed to [graphics::plot.default()].
#'
#' @return `x`, invisibly.
#' @export
plot.fmcp <- function(x,
                      shade_null = TRUE,
                      null_col = "grey90",
                      xlab = "Argument",
                      ylab = expression(beta(t)),
                      main = "fMCP coefficient function",
                      lwd = 2,
                      ...) {
  limits <- range(c(0, x$beta), finite = TRUE)
  if (diff(limits) == 0) {
    limits <- limits + c(-1, 1)
  }
  graphics::plot.default(
    x$argvals, x$beta,
    type = "n", xlab = xlab, ylab = ylab, main = main,
    ylim = limits, ...
  )
  if (isTRUE(shade_null)) {
    tolerance <- max(x$control$norm_floor, .Machine$double.eps^0.5)
    null <- which(x$interval_norms <= tolerance)
    for (j in null) {
      graphics::rect(
        x$basis$breaks[j], limits[1L],
        x$basis$breaks[j + 1L], limits[2L],
        col = null_col, border = NA
      )
    }
  }
  graphics::abline(h = 0, col = "grey60", lty = 2)
  graphics::lines(x$argvals, x$beta, lwd = lwd)
  invisible(x)
}
