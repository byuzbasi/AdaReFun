# Method adapters used by the controlled version-3 and confirmatory version-4
# FSCAD comparison studies.
# They deliberately do not modify the FMCP or PFLR package sources.

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "FMCP-package", "DESCRIPTION")) &&
        file.exists(file.path(current, "simulation", "src",
                              "fscad_comparator.cpp"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the FMCP project root.", call. = FALSE)
    }
    current <- parent
  }
}

compile_fscad_comparator <- function(project_root = find_project_root(),
                                     rebuild = FALSE) {
  source_file <- file.path(
    project_root, "simulation", "src", "fscad_comparator.cpp"
  )
  cache_dir <- file.path(tempdir(), "fmcp-v3-rcpp-cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  Rcpp::sourceCpp(source_file, cacheDir = cache_dir, rebuild = rebuild)
  invisible(TRUE)
}

capture_method <- function(expression) {
  warnings <- character()
  start <- proc.time()[["elapsed"]]
  value <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(condition) {
        warnings <<- c(warnings, conditionMessage(condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) condition
  )
  elapsed <- proc.time()[["elapsed"]] - start
  list(
    value = if (inherits(value, "error")) NULL else value,
    success = !inherits(value, "error"),
    elapsed = unname(elapsed),
    warnings = paste(unique(warnings), collapse = " | "),
    error = if (inherits(value, "error")) conditionMessage(value) else ""
  )
}

prepare_common_design <- function(y, x, argvals, M, degree = 3L,
                                  penalty_order = 2L,
                                  integration = "simpson") {
  validated <- FMCP:::.fmcp_validate_data(y, x, argvals)
  y <- validated$y
  x <- validated$x
  argvals <- validated$argvals
  integration_weights <- FMCP:::.fmcp_integration_weights(
    argvals, integration
  )
  basis <- FMCP:::.fmcp_basis_setup(
    argvals, as.integer(M), as.integer(degree), as.integer(penalty_order)
  )
  y_mean <- mean(y)
  x_mean <- colMeans(x)
  centered_y <- y - y_mean
  centered_x <- sweep(x, 2L, x_mean, FUN = "-")
  design <- FMCP:::cpp_design_matrix(
    centered_x, basis$basis_observed, integration_weights
  )
  local_gram <- FMCP:::cpp_local_gram(
    basis$basis_quad,
    basis$quad_weights,
    basis$quad_intervals,
    as.integer(M),
    as.integer(degree + 1L)
  )
  roughness <- FMCP:::cpp_weighted_crossprod(
    basis$derivative_quad, basis$quad_weights
  )
  list(
    y = y,
    x = x,
    argvals = argvals,
    centered_y = centered_y,
    design = design,
    roughness = roughness,
    local_gram = local_gram,
    interval_scale = sqrt(M / diff(basis$domain)),
    integration_weights = integration_weights,
    basis = basis,
    means = list(y = y_mean, x = x_mean)
  )
}

fit_controlled_fscad <- function(y, x, argvals, M, lambda, kappa,
                                 scad_a = 3.7, degree = 3L,
                                 penalty_order = 2L,
                                 integration = "simpson",
                                 max_iter = 150L, tolerance = 1e-8,
                                 norm_floor = 1e-8, cutoff = 0,
                                 refit = TRUE, keep_trace = FALSE) {
  if (!exists("cpp_fscad_grid", mode = "function", inherits = TRUE)) {
    compile_fscad_comparator()
  }
  prepared <- prepare_common_design(
    y, x, argvals, M, degree, penalty_order, integration
  )
  grid <- cpp_fscad_grid(
    prepared$centered_y,
    prepared$design,
    prepared$roughness,
    prepared$local_gram,
    as.numeric(lambda),
    as.numeric(kappa),
    as.numeric(scad_a),
    prepared$interval_scale,
    as.integer(max_iter),
    as.numeric(tolerance),
    as.numeric(norm_floor),
    as.numeric(cutoff),
    isTRUE(refit),
    isTRUE(keep_trace)
  )
  tuning <- data.frame(
    lambda = grid$lambda,
    kappa = grid$kappa,
    scad_a = rep(grid$scad_a, length(grid$lambda)),
    bic = grid$bic,
    rss = grid$rss,
    df = grid$df,
    objective = grid$objective,
    iterations = grid$iterations,
    active = grid$active,
    success = grid$success,
    converged = grid$converged,
    refit_success = grid$refit_success,
    objective_increased = grid$objective_increased,
    message = as.character(grid$message),
    stringsAsFactors = FALSE
  )
  valid <- which(tuning$success %in% TRUE & is.finite(tuning$bic))
  if (!length(valid)) {
    details <- unique(tuning$message[nzchar(tuning$message)])
    stop(
      "No controlled fSCAD tuning candidate produced a valid fit. ",
      paste(details, collapse = " "), call. = FALSE
    )
  }
  converged <- valid[tuning$converged[valid] %in% TRUE]
  eligible <- if (length(converged)) converged else valid
  selected_index <- eligible[which.min(tuning$bic[eligible])]
  coefficients <- as.numeric(grid$coefficients[, selected_index])
  beta <- as.numeric(prepared$basis$basis_observed %*% coefficients)
  fitted_values <- as.numeric(
    prepared$means$y + prepared$design %*% coefficients
  )
  selected <- tuning[selected_index, , drop = FALSE]

  structure(
    list(
      method = "fSCAD",
      coefficients = coefficients,
      beta = beta,
      argvals = prepared$argvals,
      fitted.values = fitted_values,
      residuals = prepared$y - fitted_values,
      selected = selected,
      tuning = tuning,
      means = prepared$means,
      basis = list(
        M = as.integer(M),
        degree = as.integer(degree),
        penalty_order = as.integer(penalty_order),
        observed = prepared$basis$basis_observed,
        integration = integration,
        integration_weights = prepared$integration_weights
      ),
      diagnostics = list(
        selected_index = selected_index,
        objective_trace = if (keep_trace) {
          grid$traces[[selected_index]]
        } else {
          NULL
        }
      )
    ),
    class = "controlled_fscad"
  )
}

predict.controlled_fscad <- function(object, newx, ...) {
  if (!is.matrix(newx) || ncol(newx) != length(object$argvals)) {
    stop("`newx` must have the training grid in columns.", call. = FALSE)
  }
  centered <- sweep(newx, 2L, object$means$x, FUN = "-")
  design <- FMCP:::cpp_design_matrix(
    centered,
    object$basis$observed,
    object$basis$integration_weights
  )
  as.numeric(object$means$y + design %*% object$coefficients)
}

resample_curves <- function(x, from, to) {
  if (length(from) == length(to) && isTRUE(all.equal(from, to))) {
    return(x)
  }
  result <- matrix(NA_real_, nrow(x), length(to))
  for (row in seq_len(nrow(x))) {
    result[row, ] <- stats::approx(
      from, x[row, ], xout = to, rule = 2, ties = "ordered"
    )$y
  }
  result
}

evaluate_pflr_basis <- function(coefficients, M, degree, argvals,
                                domain = range(argvals)) {
  knots <- seq(domain[1L], domain[2L], length.out = M + 1L)
  basis <- fda::create.bspline.basis(
    rangeval = domain,
    breaks = knots,
    nbasis = M + degree,
    norder = degree + 1L
  )
  as.numeric(fda::eval.basis(argvals, basis) %*% coefficients)
}

fit_reference_slos <- function(y, x, argvals, M, lambda, kappa,
                               cutoff, degree = 3L, max_iter = 150L,
                               tolerance = 1e-8) {
  reference_grid <- seq(min(argvals), max(argvals), length.out = M + 1L)
  if (length(reference_grid) %% 2L == 0L) {
    stop("Reference SLoS needs even `M` for Simpson integration.",
         call. = FALSE)
  }
  reference_x <- resample_curves(x, argvals, reference_grid)
  knots <- seq(min(argvals), max(argvals), length.out = M + 1L)
  basis <- fda::create.bspline.basis(
    rangeval = range(argvals),
    breaks = knots,
    nbasis = M + degree,
    norder = degree + 1L
  )
  centered_y <- y - mean(y)
  centered_x <- sweep(reference_x, 2L, colMeans(reference_x), FUN = "-")
  core <- PFLR:::slos_temp(
    Y = centered_y,
    X = centered_x,
    Maxiter = as.integer(max_iter),
    lambda = as.numeric(lambda),
    gamma = as.numeric(kappa),
    beta.basis = basis,
    absTol = as.numeric(tolerance),
    Cutoff = as.numeric(cutoff)
  )
  fit <- list(
    beta = as.numeric(fda::eval.fd(knots, core$beta)),
    extra = list(
      X = t(reference_x), Y = y, M = M, d = degree,
      domain = range(argvals), b = core$b,
      U = core$U,
      V = fda::eval.penalty(basis, fda::int2Lfd(2)),
      delta = PFLR:::T1.hat(as.numeric(fda::eval.fd(knots, core$beta))),
      Optgamma = kappa, Optlambda = lambda
    )
  )
  class(fit) <- "slos"
  fit$beta_common <- evaluate_pflr_basis(
    fit$extra$b, M, degree, argvals, range(argvals)
  )
  fit$reference_grid <- reference_grid
  fit$argvals_common <- argvals
  fit
}

fit_reference_slos_grid <- function(y, x, argvals, M, lambda, kappa,
                                     cutoff, degree = 3L,
                                     max_iter = 150L,
                                     tolerance = 1e-10) {
  candidates <- list()
  tuning_rows <- list()
  candidate <- 1L
  for (current_kappa in kappa) {
    for (current_lambda in lambda) {
      captured <- capture_method(fit_reference_slos(
        y, x, argvals, M,
        lambda = current_lambda,
        kappa = current_kappa,
        cutoff = cutoff,
        degree = degree,
        max_iter = max_iter,
        tolerance = tolerance
      ))
      if (captured$success) {
        fit <- captured$value
        b <- as.numeric(fit$extra$b)
        fitted <- mean(y) + as.numeric(fit$extra$U %*% b)
        rss <- sum((y - fitted)^2)
        active <- which(is.finite(b) & b != 0)
        if (!length(active)) {
          df <- 0
        } else {
          active_U <- fit$extra$U[, active, drop = FALSE]
          active_V <- fit$extra$V[active, active, drop = FALSE]
          system <- crossprod(active_U) +
            length(y) * current_kappa * active_V
          influence <- tryCatch(
            solve(system, crossprod(active_U)),
            error = function(condition) NULL
          )
          df <- if (is.null(influence)) NA_real_ else sum(diag(influence))
        }
        bic <- if (is.finite(df)) {
          length(y) * log(max(rss / length(y), .Machine$double.xmin)) +
            log(length(y)) * df
        } else {
          Inf
        }
      } else {
        fit <- NULL
        rss <- df <- NA_real_
        bic <- Inf
      }
      candidates[[candidate]] <- fit
      tuning_rows[[candidate]] <- data.frame(
        lambda = current_lambda,
        kappa = current_kappa,
        bic = bic,
        rss = rss,
        df = df,
        success = captured$success && is.finite(bic),
        warnings = captured$warnings,
        error = captured$error,
        stringsAsFactors = FALSE
      )
      candidate <- candidate + 1L
    }
  }
  tuning <- do.call(rbind, tuning_rows)
  valid <- which(tuning$success)
  if (!length(valid)) {
    stop(
      "Every reference SLoS tuning candidate failed: ",
      paste(unique(tuning$error[nzchar(tuning$error)]), collapse = " | "),
      call. = FALSE
    )
  }
  selected_index <- valid[which.min(tuning$bic[valid])]
  fit <- candidates[[selected_index]]
  fit$extra$Optlambda <- tuning$lambda[selected_index]
  fit$extra$Optgamma <- tuning$kappa[selected_index]
  fit$extra$tuning <- tuning
  fit
}

predict_reference_slos <- function(object, newx) {
  newx_reference <- resample_curves(
    newx, object$argvals_common, object$reference_grid
  )
  as.numeric(stats::predict(object, newx_reference))
}

safe_effective_df <- function(U, V, n, penalty, coefficients) {
  active <- which(is.finite(coefficients) & coefficients != 0)
  if (!length(active)) {
    return(0)
  }
  active_U <- U[, active, drop = FALSE]
  active_V <- V[active, active, drop = FALSE]
  system <- crossprod(active_U) + n * penalty * active_V
  influence_core <- tryCatch(
    solve(system, crossprod(active_U)),
    error = function(condition) NULL
  )
  if (is.null(influence_core)) NA_real_ else sum(diag(influence_core))
}

fit_ngr_safe <- function(y, x, M, degree, domain, tuning) {
  knots <- seq(domain[1L], domain[2L], length.out = M + 1L)
  basis <- fda::create.bspline.basis(
    rangeval = domain, breaks = knots,
    nbasis = M + degree, norder = degree + 1L
  )
  V <- fda::eval.penalty(basis, fda::int2Lfd(2))
  argvals <- seq(domain[1L], domain[2L], length.out = ncol(x))
  basismat <- fda::eval.basis(argvals, basis)
  candidates <- list()
  tuning_rows <- list()
  candidate <- 1L
  for (kappa in tuning$kappa) {
    for (tau in tuning$tau) {
      captured <- capture_method(PFLR:::ngr.algorithm(
        Y = y, X = x, V = V, domain = domain, basismat = basismat,
        alphaPS = tuning$alphaPS, kappa = kappa, tau = tau,
        gamma = if (is.null(tuning$gamma)) 0.5 else tuning$gamma,
        M = M, d = degree,
        niter = if (is.null(tuning$niter)) 100L else tuning$niter
      ))
      if (captured$success) {
        fit <- captured$value
        residual <- y - as.numeric(fit$fitted.values)
        df <- safe_effective_df(fit$U, V, length(y), kappa, fit$b)
        bic <- if (is.finite(df)) {
          length(y) * log(max(mean(residual^2), .Machine$double.xmin)) +
            log(length(y)) * df
        } else {
          Inf
        }
      } else {
        fit <- NULL
        df <- NA_real_
        bic <- Inf
      }
      candidates[[candidate]] <- fit
      tuning_rows[[candidate]] <- data.frame(
        kappa = kappa, tau = tau, bic = bic, df = df,
        success = captured$success && is.finite(bic),
        warnings = captured$warnings, error = captured$error,
        stringsAsFactors = FALSE
      )
      candidate <- candidate + 1L
    }
  }
  tuning_table <- do.call(rbind, tuning_rows)
  valid <- which(tuning_table$success)
  if (!length(valid)) {
    stop(
      "Every safe NGR tuning candidate failed: ",
      paste(unique(tuning_table$error[nzchar(tuning_table$error)]),
            collapse = " | "),
      call. = FALSE
    )
  }
  selected_index <- valid[which.min(tuning_table$bic[valid])]
  selected_fit <- candidates[[selected_index]]
  selected_row <- tuning_table[selected_index, , drop = FALSE]
  result <- list(
    beta = selected_fit$beta,
    extra = list(
      b = selected_fit$b,
      delta = selected_fit$delta,
      Ymean = selected_fit$Ymean,
      Xmean = selected_fit$Xmean,
      Optkappa = selected_row$kappa,
      Opttau = selected_row$tau,
      M = M, d = degree, domain = domain,
      tuning = tuning_table
    )
  )
  class(result) <- "ngr"
  result
}

fit_flirti_safe <- function(y, x, degree, domain, tuning) {
  n <- length(y)
  p <- ncol(x)
  centered_y <- y - mean(y)
  centered_x <- sweep(x, 2L, colMeans(x), FUN = "-")
  argvals <- seq(domain[1L], domain[2L], length.out = p)
  h <- diff(domain) / (p - 1L)
  simpson <- c(1, rep(c(4, 2), length.out = p - 2L), 1)
  candidates <- list()
  tuning_rows <- list()
  candidate <- 1L
  cons <- if (is.null(tuning$cons)) 4L else tuning$cons

  for (Mf in tuning$Mf) {
    knots <- seq(domain[1L], domain[2L], length.out = Mf + 1L)
    basis <- fda::create.bspline.basis(
      rangeval = domain, breaks = knots,
      nbasis = Mf + degree, norder = degree + 1L
    )
    basismat <- fda::eval.basis(argvals, basis)
    U <- h / 3 * centered_x %*% (basismat * simpson)
    for (lambda in tuning$lambda) {
      captured <- capture_method(PFLR:::flirty1_dz2(
        Y = centered_y, U = U, basis = basis,
        lambda = lambda, cons = cons
      ))
      if (captured$success) {
        current <- captured$value
        fitted <- as.numeric(U %*% solve(current$A, current$gamma))
        rss <- sum((centered_y - fitted)^2)
        bic <- n * log(max(rss / n, .Machine$double.xmin)) +
          log(n) * current$df
      } else {
        current <- NULL
        bic <- Inf
      }
      candidates[[candidate]] <- current
      tuning_rows[[candidate]] <- data.frame(
        Mf = Mf, lambda = lambda, bic = bic,
        success = captured$success && is.finite(bic),
        warnings = captured$warnings, error = captured$error,
        stringsAsFactors = FALSE
      )
      candidate <- candidate + 1L
    }
  }
  tuning_table <- do.call(rbind, tuning_rows)
  valid <- which(tuning_table$success)
  if (!length(valid)) {
    stop(
      "Every safe FLiRTI tuning candidate failed: ",
      paste(unique(tuning_table$error[nzchar(tuning_table$error)]),
            collapse = " | "),
      call. = FALSE
    )
  }
  selected_index <- valid[which.min(tuning_table$bic[valid])]
  selected <- candidates[[selected_index]]
  selected_row <- tuning_table[selected_index, , drop = FALSE]
  result <- list(
    beta = selected$beta,
    extra = list(
      X = t(x), Y = y, domain = domain,
      delta = PFLR:::T1.hat(selected$beta),
      OptM = selected_row$Mf,
      Optlambda = selected_row$lambda,
      tuning = tuning_table
    )
  )
  class(result) <- "flirti"
  result
}

fit_pflr_method <- function(method, y, x, argvals, M = 40L,
                            degree = 3L, tuning = list()) {
  domain <- range(argvals)
  method <- match.arg(method, c("NGR", "PenS", "FLiRTI"))
  if (identical(method, "NGR")) {
    fit <- fit_ngr_safe(
      y, x, as.integer(M), as.integer(degree), domain, tuning
    )
    beta <- as.numeric(fit$beta)
  } else if (identical(method, "PenS")) {
    fit <- PFLR::PenS(
      Y = y, X = x, alpha = tuning$alpha,
      M = as.integer(M), d = as.integer(degree), domain = domain
    )
    beta <- as.numeric(fit$beta)
  } else {
    fit <- fit_flirti_safe(
      y, x, as.integer(degree), domain, tuning
    )
    beta <- as.numeric(stats::approx(
      seq(domain[1L], domain[2L], length.out = length(fit$beta)),
      fit$beta,
      xout = argvals,
      rule = 2,
      ties = "ordered"
    )$y)
  }
  fit$beta_common <- beta
  fit$argvals_common <- argvals
  fit$method_label <- method
  fit
}

predict_pflr_method <- function(object, newx) {
  as.numeric(stats::predict(object, newx))
}
