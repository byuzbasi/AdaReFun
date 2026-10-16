#' Model coefficients for an fSCAD fit
#'
#' @param object An object of class `fscad`.
#' @param type Return B-spline coefficients (`"basis"`) or coefficient-function
#'   values on the training grid (`"function"`).
#' @param ... Unused.
#'
#' @return A named numeric vector for `type = "basis"`, or a data frame with
#'   `argvals` and `beta` for `type = "function"`.
#' @export
coef.fscad <- function(object, type = c("basis", "function"), ...) {
  coef.fmcp(object, type = type, ...)
}

#' Predict from an fSCAD fit
#'
#' @param object An object of class `fscad`.
#' @param newdata Optional finite numeric matrix whose columns correspond
#'   exactly to the training `argvals`. If omitted, training fitted values are
#'   returned.
#' @param ... Unused.
#'
#' @return Numeric vector of predicted responses.
#' @export
predict.fscad <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(predict.fmcp(object, ...))
  }
  predict.fmcp(object, newdata = newdata, ...)
}

#' @rdname predict.fscad
#' @export
fitted.fscad <- function(object, ...) {
  object$fitted.values
}

#' @rdname predict.fscad
#' @export
residuals.fscad <- function(object, ...) {
  object$residuals
}

#' Number of observations in an fSCAD fit
#'
#' @param object An object of class `fscad`.
#' @param ... Unused.
#'
#' @return Integer number of training subjects.
#' @export
nobs.fscad <- function(object, ...) {
  unname(object$dimensions[["n"]])
}

#' Print an fSCAD fit
#'
#' @param x An object of class `fscad`.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#' @export
print.fscad <- function(x, ...) {
  selected <- x$selected
  cat("Functional SCAD (fSCAD) regression fit\n")
  cat("  Subjects:", x$dimensions[["n"]],
      " Grid points:", x$dimensions[["p"]],
      " Basis functions:", x$dimensions[["basis"]], "\n")
  cat("  lambda:", format(selected$lambda),
      " kappa:", format(selected$kappa),
      " SCAD a:", format(selected$scad_a), "\n")
  cat("  BIC:", format(selected$bic),
      " Active coefficients:", selected$active,
      " Converged:", if (selected$converged) "yes" else "no", "\n")
  invisible(x)
}

#' Summarize an fSCAD fit
#'
#' @param object An object of class `fscad`.
#' @param x An object of class `summary.fscad` for printing.
#' @param ... Unused.
#'
#' @return An object of class `summary.fscad`.
#' @export
summary.fscad <- function(object, ...) {
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
    class = "summary.fscad"
  )
}

#' @rdname summary.fscad
#' @export
print.summary.fscad <- function(x, ...) {
  cat("Functional SCAD (fSCAD) regression summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nDimensions:\n")
  cat("  Subjects:", x$dimensions[["n"]],
      " Grid points:", x$dimensions[["p"]],
      " Basis functions:", x$dimensions[["basis"]], "\n")
  cat("\nSelected tuning parameters:\n")
  cat("  lambda =", format(x$selected$lambda),
      " kappa =", format(x$selected$kappa),
      " SCAD a =", format(x$selected$scad_a), "\n")
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

#' Plot an fSCAD coefficient function
#'
#' @param x An object of class `fscad`.
#' @param shade_null Logical. Shade intervals whose fitted functional norm is
#'   numerically zero.
#' @param null_col Color for shaded null intervals.
#' @param xlab,ylab,main,lwd Standard plotting parameters.
#' @param ... Additional arguments passed to [graphics::plot.default()].
#'
#' @return `x`, invisibly.
#' @export
plot.fscad <- function(x,
                       shade_null = TRUE,
                       null_col = "grey90",
                       xlab = "Argument",
                       ylab = expression(beta(t)),
                       main = "fSCAD coefficient function",
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
