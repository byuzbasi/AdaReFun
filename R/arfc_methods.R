#' Methods for an ARFC fit
#'
#' @param object,x An object of class `arfc`.
#' @param type For `coef()`, return basis coefficients or coefficient-function
#'   values. For `plot()`, choose the coefficient curve, regional norms,
#'   adaptive weights, or tuning path.
#' @param newdata Optional finite predictor matrix on the training grid.
#' @param shade_null Shade finest-level regions with numerically zero fitted
#'   regional norm in coefficient plots.
#' @param null_col Color used for zero-region shading.
#' @param ... Additional graphical arguments where applicable.
#'
#' @return `coef()` returns a vector or data frame; `predict()`, `fitted()`, and
#'   `residuals()` return numeric vectors; `nobs()` returns an integer;
#'   printing and plotting methods return their object invisibly.
#' @name arfc-methods
NULL

#' @rdname arfc-methods
#' @export
coef.arfc <- function(object, type = c("basis", "function"), ...) {
  type <- match.arg(type)
  if (type == "basis") {
    return(object$coefficients)
  }
  data.frame(argvals = object$argvals, beta = object$beta)
}

.arfc_predict_matrix <- function(object, newdata) {
  if (!is.matrix(newdata) || !is.numeric(newdata) || anyNA(newdata) ||
      any(!is.finite(newdata))) {
    stop("`newdata` must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(newdata) != object$dimensions[["p"]]) {
    stop("`newdata` must have the same functional grid as the training data.",
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

#' @rdname arfc-methods
#' @export
predict.arfc <- function(object, newdata, ...) {
  if (missing(newdata) || is.null(newdata)) {
    return(object$fitted.values)
  }
  .arfc_predict_matrix(object, newdata)
}

#' @rdname arfc-methods
#' @export
fitted.arfc <- function(object, ...) object$fitted.values

#' @rdname arfc-methods
#' @export
residuals.arfc <- function(object, ...) object$residuals

#' @rdname arfc-methods
#' @export
nobs.arfc <- function(object, ...) unname(object$dimensions[["n"]])

#' @rdname arfc-methods
#' @export
print.arfc <- function(x, ...) {
  cat("Adaptive regional folded-concave functional regression\n")
  cat("  Method:", x$method$label,
      " Adaptive:", if (x$method$adaptive) "yes" else "no", "\n")
  cat("  Subjects:", x$dimensions[["n"]],
      " Grid points:", x$dimensions[["p"]],
      " Regions:", x$dimensions[["regions"]], "\n")
  cat("  lambda:", format(x$selected$lambda),
      " kappa:", format(x$selected$kappa),
      " shape:", format(x$method$shape), "\n")
  cat("  BIC:", format(x$selected$bic),
      " Active coefficients:", x$selected$active,
      " Converged:", if (x$selected$converged) "yes" else "no", "\n")
  invisible(x)
}

#' @rdname arfc-methods
#' @export
summary.arfc <- function(object, ...) {
  residual_df <- object$dimensions[["n"]] - object$selected$df
  residual_scale <- if (is.finite(residual_df) && residual_df > 0) {
    sqrt(object$selected$rss / residual_df)
  } else {
    NA_real_
  }
  fine <- object$regions$level_M == max(object$regions$level_M)
  tolerance <- max(object$control$norm_floor, .Machine$double.eps^0.5)
  structure(
    list(
      call = object$call,
      method = object$method,
      dimensions = object$dimensions,
      selected = object$selected,
      residual_scale = residual_scale,
      coefficient_summary = summary(object$beta),
      zero_finest_regions = sum(
        object$regions$selected_region_norm[fine] <= tolerance
      ),
      total_finest_regions = sum(fine),
      diagnostics = object$diagnostics
    ),
    class = "summary.arfc"
  )
}

#' @rdname arfc-methods
#' @export
print.summary.arfc <- function(x, ...) {
  cat("Adaptive regional folded-concave regression summary\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nMethod:", x$method$label,
      " adaptive =", x$method$adaptive, "\n")
  cat("Dimensions: n =", x$dimensions[["n"]],
      ", p =", x$dimensions[["p"]],
      ", basis =", x$dimensions[["basis"]],
      ", regions =", x$dimensions[["regions"]], "\n", sep = "")
  cat("Selected: lambda =", format(x$selected$lambda),
      ", kappa =", format(x$selected$kappa),
      ", BIC =", format(x$selected$bic), "\n", sep = "")
  cat("Residual scale:", format(x$residual_scale), "\n")
  cat("Finest-level zero regions:", x$zero_finest_regions, "of",
      x$total_finest_regions, "\n")
  cat("Coefficient-function values:\n")
  print(x$coefficient_summary)
  invisible(x)
}

.arfc_plot_coefficient <- function(x, shade_null, null_col, ...) {
  limits <- range(c(0, x$beta), finite = TRUE)
  if (diff(limits) == 0) limits <- limits + c(-1, 1)
  graphics::plot.default(
    x$argvals, x$beta, type = "n", xlab = "Functional argument",
    ylab = expression(beta(t)), main = paste(x$method$label, "coefficient"),
    ylim = limits, ...
  )
  if (isTRUE(shade_null)) {
    fine <- x$regions$level_M == max(x$regions$level_M)
    tolerance <- max(x$control$norm_floor, .Machine$double.eps^0.5)
    zero <- which(fine & x$regions$selected_region_norm <= tolerance)
    for (index in zero) {
      graphics::rect(
        x$regions$left[index], limits[1L],
        x$regions$right[index], limits[2L],
        col = null_col, border = NA
      )
    }
  }
  graphics::abline(h = 0, col = "grey60", lty = 2)
  graphics::lines(x$argvals, x$beta, lwd = 2, col = "#0072B2")
}

.arfc_plot_regional <- function(x, variable, ylab, main, ...) {
  level <- factor(x$regions$level_M, levels = unique(x$regions$level_M))
  colors <- grDevices::hcl.colors(nlevels(level), "Dark 3")
  graphics::plot.default(
    x$regions$midpoint, x$regions[[variable]],
    pch = 19, col = colors[level], xlab = "Region midpoint", ylab = ylab,
    main = main, ...
  )
  graphics::legend(
    "topright", legend = paste("M =", levels(level)), col = colors,
    pch = 19, bty = "n"
  )
}

#' @rdname arfc-methods
#' @export
plot.arfc <- function(x,
                      type = c("coefficient", "regions", "weights", "tuning"),
                      shade_null = TRUE,
                      null_col = "grey90",
                      ...) {
  type <- match.arg(type)
  if (type == "coefficient") {
    .arfc_plot_coefficient(x, shade_null, null_col, ...)
  } else if (type == "regions") {
    .arfc_plot_regional(
      x, "selected_region_norm", "Selected regional norm",
      paste(x$method$label, "regional norms"), ...
    )
  } else if (type == "weights") {
    .arfc_plot_regional(
      x, "adaptive_weight", "Adaptive weight",
      paste(x$method$label, "adaptive weights"), ...
    )
  } else {
    graphics::plot.default(
      x$tuning$lambda, x$tuning$bic, type = "b", pch = 19,
      xlab = expression(lambda), ylab = "BIC",
      main = paste(x$method$label, "tuning path"), ...
    )
    graphics::points(
      x$selected$lambda, x$selected$bic, pch = 21,
      bg = "#D55E00", cex = 1.4
    )
  }
  invisible(x)
}

#' Methods for an ARFC-Select fit
#'
#' @param object,x An object of class `arfc_select`.
#' @param type For `plot()`, show validation losses, the selected coefficient,
#'   or all candidate coefficient functions. For `coef()`, passed to the
#'   selected [arfc()] fit.
#' @param newdata Optional predictor matrix on the training grid.
#' @param ... Additional arguments.
#'
#' @return Methods delegate numerical results to the selected `arfc` fit;
#'   plotting and printing return `x` invisibly.
#' @name arfc-select-methods
NULL

#' @rdname arfc-select-methods
#' @export
coef.arfc_select <- function(object, type = c("basis", "function"), ...) {
  coef(object$selected_fit, type = match.arg(type), ...)
}

#' @rdname arfc-select-methods
#' @export
predict.arfc_select <- function(object, newdata, ...) {
  if (missing(newdata)) return(predict(object$selected_fit, ...))
  predict(object$selected_fit, newdata = newdata, ...)
}

#' @rdname arfc-select-methods
#' @export
fitted.arfc_select <- function(object, ...) fitted(object$selected_fit, ...)

#' @rdname arfc-select-methods
#' @export
residuals.arfc_select <- function(object, ...) residuals(object$selected_fit, ...)

#' @rdname arfc-select-methods
#' @export
nobs.arfc_select <- function(object, ...) nobs(object$selected_fit, ...)

#' @rdname arfc-select-methods
#' @export
print.arfc_select <- function(x, ...) {
  cat("ARFC-Select validation selection\n")
  cat("  Selected:", x$selected_label, "(", x$selected_method, ")\n", sep = "")
  cat("  Validation MSE:",
      format(x$diagnostics$minimum_validation_mse),
      " Margin to second:", format(x$diagnostics$margin_to_second), "\n")
  cat("  Validation observations:",
      x$diagnostics$validation_sample_size,
      " Used for fitting: no\n")
  invisible(x)
}

#' @rdname arfc-select-methods
#' @export
summary.arfc_select <- function(object, ...) {
  structure(
    list(
      call = object$call,
      selected_method = object$selected_method,
      selected_label = object$selected_label,
      selection = object$selection,
      diagnostics = object$diagnostics,
      selected_fit = summary(object$selected_fit)
    ),
    class = "summary.arfc_select"
  )
}

#' @rdname arfc-select-methods
#' @export
print.summary.arfc_select <- function(x, ...) {
  cat("ARFC-Select summary\n\n")
  cat("Selected candidate:", x$selected_label, "\n")
  print(x$selection[, c(
    "label", "validation_response_mse", "difference_from_minimum",
    "tied", "selected"
  )], row.names = FALSE)
  cat("\nValidation used for fitting:",
      x$diagnostics$validation_used_for_fitting, "\n")
  invisible(x)
}

#' @rdname arfc-select-methods
#' @export
plot.arfc_select <- function(x,
                             type = c("validation", "coefficient", "candidates"),
                             ...) {
  type <- match.arg(type)
  if (type == "validation") {
    colors <- ifelse(x$selection$selected, "#0072B2", "grey70")
    graphics::barplot(
      x$selection$validation_response_mse,
      names.arg = x$selection$label,
      col = colors,
      ylab = "Validation response MSE",
      main = "ARFC-Select candidate comparison",
      ...
    )
  } else if (type == "coefficient") {
    plot(x$selected_fit, type = "coefficient", ...)
  } else {
    fits <- x$candidates
    values <- do.call(cbind, lapply(fits, `[[`, "beta"))
    colors <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7")
    graphics::matplot(
      x$selected_fit$argvals, values, type = "l", lty = 1,
      lwd = 2, col = colors, xlab = "Functional argument",
      ylab = expression(beta(t)), main = "ARFC candidate coefficients", ...
    )
    graphics::abline(h = 0, col = "grey60", lty = 2)
    graphics::legend(
      "topright", legend = x$selection$label, col = colors,
      lty = 1, lwd = 2, bty = "n"
    )
  }
  invisible(x)
}
