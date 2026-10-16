initialize_real_method_environment <- function(project_root, rebuild = FALSE) {
  source(file.path(
    project_root, "simulation", "comparators",
    "fscad_and_pflr_adapters.R"
  ))
  source(file.path(
    project_root, "simulation", "comparators", "classic_methods.R"
  ))
  compile_fscad_comparator(project_root, rebuild = rebuild)
  invisible(TRUE)
}

fit_training_scaler <- function(y, x, argvals_u) {
  weights <- trapezoid_weights(argvals_u)
  centered_x <- sweep(x, 2L, colMeans(x), FUN = "-")
  y_scale <- stats::sd(y)
  x_scale <- sqrt(mean(rowSums(sweep(centered_x^2, 2L, weights, `*`))))
  if (!is.finite(y_scale) || y_scale <= sqrt(.Machine$double.eps)) {
    stop("Training-response scale is numerically zero.", call. = FALSE)
  }
  if (!is.finite(x_scale) || x_scale <= sqrt(.Machine$double.eps)) {
    stop("Training functional-RMS scale is numerically zero.", call. = FALSE)
  }
  list(y_scale = y_scale, x_scale = x_scale, weights = weights)
}

apply_training_scaler <- function(y, x, scaler) {
  list(y = y / scaler$y_scale, x = x / scaler$x_scale)
}

back_transform_beta <- function(beta_scaled, scaler, argvals_physical) {
  physical_length <- diff(range(argvals_physical))
  if (!is.finite(physical_length) || physical_length <= 0) {
    stop("Physical argument domain must have positive length.", call. = FALSE)
  }
  beta_scaled * scaler$y_scale / (scaler$x_scale * physical_length)
}

common_engine_fit_over_M <- function(y, x, argvals, config, method,
                                     cutoff = 0, mcp_gamma = NULL) {
  if (is.null(mcp_gamma)) mcp_gamma <- config$mcp_gamma
  candidates <- vector("list", length(config$M))
  rows <- vector("list", length(config$M))
  start <- proc.time()[["elapsed"]]

  for (index in seq_along(config$M)) {
    current_M <- as.integer(config$M[index])
    captured <- if (identical(method, "FMCP")) {
      capture_method(FMCP::fmcp(
        y, x, argvals, M = current_M,
        lambda = config$lambda, kappa = config$kappa,
        degree = config$degree, penalty_order = config$penalty_order,
        mcp_gamma = mcp_gamma, integration = config$integration,
        control = FMCP::fmcp_control(
          max_iter = config$max_iter, tol = config$tolerance,
          norm_floor = config$norm_floor, cutoff = cutoff,
          refit = TRUE, trace = FALSE, keep_matrices = FALSE
        )
      ))
    } else if (identical(method, "fSCAD")) {
      capture_method(fit_controlled_fscad(
        y, x, argvals, M = current_M,
        lambda = config$lambda, kappa = config$kappa,
        scad_a = config$scad_a, degree = config$degree,
        penalty_order = config$penalty_order,
        integration = config$integration,
        max_iter = config$max_iter, tolerance = config$tolerance,
        norm_floor = config$norm_floor, cutoff = cutoff,
        refit = TRUE, keep_trace = FALSE
      ))
    } else if (identical(method, "Smooth")) {
      capture_method(FMCP::fmcp(
        y, x, argvals, M = current_M,
        lambda = 0, kappa = config$kappa,
        degree = config$degree, penalty_order = config$penalty_order,
        mcp_gamma = config$mcp_gamma, integration = config$integration,
        control = FMCP::fmcp_control(
          max_iter = 1L, tol = config$tolerance,
          norm_floor = config$norm_floor, cutoff = 0,
          refit = FALSE, trace = FALSE, keep_matrices = FALSE
        )
      ))
    } else if (identical(method, "OLS")) {
      capture_method(FMCP::fmcp(
        y, x, argvals, M = current_M,
        lambda = 0, kappa = 0,
        degree = config$degree, penalty_order = config$penalty_order,
        mcp_gamma = config$mcp_gamma, integration = config$integration,
        control = FMCP::fmcp_control(
          max_iter = 1L, tol = config$tolerance,
          norm_floor = config$norm_floor, cutoff = 0,
          refit = FALSE, trace = FALSE, keep_matrices = FALSE
        )
      ))
    } else {
      stop("Unknown common-engine method: ", method, call. = FALSE)
    }

    candidates[[index]] <- captured
    if (captured$success) {
      selected <- captured$value$selected[1L, , drop = FALSE]
      rows[[index]] <- data.frame(
        method = method, M = current_M, fit_success = TRUE,
        bic = selected$bic, rss = selected$rss,
        lambda = selected$lambda, kappa = selected$kappa,
        shape = if (identical(method, "FMCP")) selected$mcp_gamma else
          if (identical(method, "fSCAD")) selected$scad_a else NA_real_,
        df = selected$df, active = selected$active,
        converged = selected$converged,
        iterations = selected$iterations,
        warnings = captured$warnings, error = captured$error,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[index]] <- data.frame(
        method = method, M = current_M, fit_success = FALSE,
        bic = NA_real_, rss = NA_real_, lambda = NA_real_,
        kappa = NA_real_, shape = NA_real_, df = NA_real_,
        active = NA_real_, converged = FALSE, iterations = NA_real_,
        warnings = captured$warnings, error = captured$error,
        stringsAsFactors = FALSE
      )
    }
  }

  diagnostics <- do.call(rbind, rows)
  valid <- which(diagnostics$fit_success & is.finite(diagnostics$bic))
  elapsed <- proc.time()[["elapsed"]] - start
  if (!length(valid)) {
    return(list(
      method = method, fit = NULL,
      beta = rep(NA_real_, length(argvals)),
      predict = function(newx) rep(NA_real_, nrow(newx)),
      selected = NULL, diagnostics = diagnostics,
      elapsed = unname(elapsed),
      warnings = paste(unique(diagnostics$warnings), collapse = " | "),
      error = paste(unique(diagnostics$error), collapse = " | ")
    ))
  }
  chosen <- valid[which.min(diagnostics$bic[valid])]
  fit <- candidates[[chosen]]$value
  list(
    method = method, fit = fit, beta = as.numeric(fit$beta),
    predict = function(newx) as.numeric(stats::predict(fit, newx)),
    selected = diagnostics[chosen, , drop = FALSE],
    diagnostics = diagnostics, elapsed = unname(elapsed),
    warnings = candidates[[chosen]]$warnings,
    error = candidates[[chosen]]$error
  )
}

make_functional_basis <- function(argvals, M, degree) {
  knots <- seq(min(argvals), max(argvals), length.out = M + 1L)
  basis <- fda::create.bspline.basis(
    rangeval = range(argvals), breaks = knots,
    nbasis = M + degree, norder = degree + 1L
  )
  list(
    object = basis,
    observed = unclass(fda::eval.basis(argvals, basis)),
    roughness = unclass(fda::eval.penalty(basis, fda::int2Lfd(2)))
  )
}

effective_df <- function(design, roughness, n, penalty) {
  system <- crossprod(design) + n * penalty * roughness
  influence_core <- tryCatch(
    solve(system, crossprod(design)), error = function(condition) NULL
  )
  if (is.null(influence_core)) NA_real_ else sum(diag(influence_core))
}

fit_pens_trapezoid <- function(y, x, argvals, weights, config) {
  captured <- capture_method({
    basis <- make_functional_basis(argvals, config$pflr_M, config$degree)
    y_mean <- mean(y)
    x_mean <- colMeans(x)
    yc <- y - y_mean
    xc <- sweep(x, 2L, x_mean, FUN = "-")
    design <- xc %*% sweep(basis$observed, 1L, weights, `*`)
    candidates <- vector("list", length(config$pflr_tuning$PenS$alpha))
    rows <- vector("list", length(candidates))
    for (index in seq_along(candidates)) {
      alpha <- config$pflr_tuning$PenS$alpha[index]
      system <- crossprod(design) + length(y) * alpha * basis$roughness
      coefficient <- tryCatch(
        as.numeric(solve(system, crossprod(design, yc))),
        error = function(condition) rep(NA_real_, ncol(design))
      )
      if (all(is.finite(coefficient))) {
        residual <- yc - as.numeric(design %*% coefficient)
        rss <- sum(residual^2)
        df <- effective_df(design, basis$roughness, length(y), alpha)
        bic <- length(y) * log(max(rss / length(y), .Machine$double.xmin)) +
          log(length(y)) * df
      } else {
        rss <- df <- NA_real_
        bic <- Inf
      }
      candidates[[index]] <- coefficient
      rows[[index]] <- data.frame(
        alpha = alpha, bic = bic, rss = rss, df = df,
        success = is.finite(bic), stringsAsFactors = FALSE
      )
    }
    tuning <- do.call(rbind, rows)
    valid <- which(tuning$success)
    if (!length(valid)) stop("Every trapezoid PenS candidate failed.")
    selected_index <- valid[which.min(tuning$bic[valid])]
    coefficient <- candidates[[selected_index]]
    list(
      coefficient = coefficient,
      beta = as.numeric(basis$observed %*% coefficient),
      basis = basis$observed, weights = weights,
      x_mean = x_mean, y_mean = y_mean,
      tuning = tuning, selected = tuning[selected_index, , drop = FALSE]
    )
  })
  if (!captured$success) return(failed_method("PenS", length(argvals), captured))
  fit <- captured$value
  list(
    method = "PenS", fit = fit, beta = fit$beta,
    predict = function(newx) {
      centered <- sweep(newx, 2L, fit$x_mean, FUN = "-")
      design <- centered %*% sweep(fit$basis, 1L, fit$weights, `*`)
      as.numeric(fit$y_mean + design %*% fit$coefficient)
    },
    selected = data.frame(
      method = "PenS", M = config$pflr_M, fit_success = TRUE,
      bic = fit$selected$bic, rss = fit$selected$rss,
      lambda = NA_real_, kappa = fit$selected$alpha,
      shape = NA_real_, df = fit$selected$df,
      active = sum(abs(fit$coefficient) > config$zero_tolerance),
      converged = TRUE, iterations = 1L,
      warnings = captured$warnings, error = "", stringsAsFactors = FALSE
    ),
    diagnostics = fit$tuning, elapsed = captured$elapsed,
    warnings = captured$warnings, error = ""
  )
}

run_ngr_core_trapezoid <- function(y, x, design, basis_observed,
                                   roughness, tuning, kappa, tau,
                                   tolerance, norm_floor) {
  n <- length(y)
  yc <- y - mean(y)
  ps_tune <- PFLR:::PS.tune.BIC(
    Y = yc, U = design, V = roughness, n = n,
    alpha = tuning$alphaPS, M = ncol(design) - 3L, d = 3L,
    plot.it = FALSE
  )
  ps_fit <- PFLR:::PS_ngr(
    Y = yc, U = design, V = roughness, n = n,
    alpha = ps_tune$Optalpha,
    M = ncol(design) - 3L, d = 3L
  )
  b0 <- as.numeric(ps_fit$bsmoothhat)
  coefficient <- b0
  iterations <- 1L
  M <- length(coefficient) - 3L
  gamma <- tuning$gamma
  for (iteration in seq_len(tuning$niter)) {
    if (iteration == 1L) next
    theta <- numeric(M)
    hn <- numeric(M)
    for (j in seq_len(M)) {
      current_cj <- PFLR:::cj(
        M = M, d = 3L, gamma = gamma, j = j, bsmooth = b0
      )
      theta[j] <- PFLR:::theta_js(
        b_sminus1 = coefficient, cj = current_cj, n = n,
        gamma = gamma, tau = tau, M = M, d = 3L, j = j
      )
      hn[j] <- PFLR:::h_js(theta[j], current_cj, gamma)
    }
    g <- vapply(
      seq_along(coefficient),
      function(k) PFLR:::g_ks(h_js = hn, k = k, M = M, d = 3L),
      FUN.VALUE = numeric(1L)
    )
    g[!is.finite(g) | g <= norm_floor] <- norm_floor
    augmented_x <- PFLR:::Ustar(
      U = design, V = roughness, n = n,
      alpha = kappa, M = M, d = 3L
    ) %*% diag(1 / (n * g))
    augmented_y <- PFLR:::Ystar(
      Y = yc, alpha = kappa, M = M, d = 3L
    )
    lasso <- glmnet::glmnet(
      x = augmented_x, y = augmented_y, standardize = FALSE,
      alpha = 1, lambda = 0.5 / n, family = "gaussian",
      intercept = FALSE
    )
    updated <- as.numeric(stats::coef(lasso))[-1L] / (n * g)
    if (any(!is.finite(updated))) stop("NGR produced non-finite coefficients.")
    denominator <- pmax(abs(coefficient), norm_floor)
    relative_change <- max(abs(updated - coefficient) / denominator)
    coefficient <- updated
    iterations <- iteration
    if (relative_change <= tolerance) break
  }
  list(
    coefficient = coefficient,
    beta = as.numeric(basis_observed %*% coefficient),
    fitted = as.numeric(mean(y) + design %*% coefficient),
    iterations = iterations,
    ps_alpha = ps_tune$Optalpha
  )
}

fit_ngr_trapezoid <- function(y, x, argvals, weights, config) {
  start <- proc.time()[["elapsed"]]
  basis <- make_functional_basis(argvals, config$pflr_M, config$degree)
  x_mean <- colMeans(x)
  y_mean <- mean(y)
  centered <- sweep(x, 2L, x_mean, FUN = "-")
  design <- centered %*% sweep(basis$observed, 1L, weights, `*`)
  candidates <- list()
  rows <- list()
  index <- 1L
  captured_warnings <- character()
  for (kappa in config$pflr_tuning$NGR$kappa) {
    for (tau in config$pflr_tuning$NGR$tau) {
      captured <- capture_method(run_ngr_core_trapezoid(
        y, x, design, basis$observed, basis$roughness,
        config$pflr_tuning$NGR, kappa, tau, config$tolerance,
        config$norm_floor
      ))
      captured_warnings <- c(captured_warnings, captured$warnings)
      if (captured$success) {
        current <- captured$value
        residual <- y - current$fitted
        rss <- sum(residual^2)
        df <- safe_effective_df(
          design, basis$roughness, length(y), kappa,
          current$coefficient
        )
        bic <- if (is.finite(df)) {
          length(y) * log(max(rss / length(y), .Machine$double.xmin)) +
            log(length(y)) * df
        } else Inf
      } else {
        current <- NULL
        rss <- df <- NA_real_
        bic <- Inf
      }
      candidates[[index]] <- current
      rows[[index]] <- data.frame(
        kappa = kappa, tau = tau, bic = bic, rss = rss, df = df,
        success = captured$success && is.finite(bic),
        iterations = if (captured$success) current$iterations else NA_integer_,
        warnings = captured$warnings, error = captured$error,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  tuning <- do.call(rbind, rows)
  valid <- which(tuning$success)
  if (!length(valid)) {
    return(list(
      method = "NGR", fit = NULL, beta = rep(NA_real_, length(argvals)),
      predict = function(newx) rep(NA_real_, nrow(newx)),
      selected = NULL, diagnostics = tuning,
      elapsed = proc.time()[["elapsed"]] - start,
      warnings = paste(unique(captured_warnings), collapse = " | "),
      error = paste(unique(tuning$error[nzchar(tuning$error)]), collapse = " | ")
    ))
  }
  selected_index <- valid[which.min(tuning$bic[valid])]
  core <- candidates[[selected_index]]
  selected_tuning <- tuning[selected_index, , drop = FALSE]
  fit <- list(
    coefficient = core$coefficient, beta = core$beta,
    basis = basis$observed, weights = weights,
    x_mean = x_mean, y_mean = y_mean,
    tuning = tuning, selected = selected_tuning
  )
  list(
    method = "NGR", fit = fit, beta = fit$beta,
    predict = function(newx) {
      centered_new <- sweep(newx, 2L, fit$x_mean, FUN = "-")
      design_new <- centered_new %*% sweep(fit$basis, 1L, fit$weights, `*`)
      as.numeric(fit$y_mean + design_new %*% fit$coefficient)
    },
    selected = data.frame(
      method = "NGR", M = config$pflr_M, fit_success = TRUE,
      bic = selected_tuning$bic, rss = selected_tuning$rss,
      lambda = selected_tuning$tau, kappa = selected_tuning$kappa,
      shape = config$pflr_tuning$NGR$gamma, df = selected_tuning$df,
      active = sum(abs(core$coefficient) > config$zero_tolerance),
      converged = selected_tuning$iterations < config$pflr_tuning$NGR$niter,
      iterations = selected_tuning$iterations,
      warnings = paste(unique(captured_warnings), collapse = " | "),
      error = "", stringsAsFactors = FALSE
    ),
    diagnostics = tuning,
    elapsed = unname(proc.time()[["elapsed"]] - start),
    warnings = paste(unique(captured_warnings), collapse = " | "),
    error = ""
  )
}

fit_flirti_trapezoid <- function(y, x, argvals, weights, config) {
  start <- proc.time()[["elapsed"]]
  y_mean <- mean(y)
  x_mean <- colMeans(x)
  yc <- y - y_mean
  xc <- sweep(x, 2L, x_mean, FUN = "-")
  candidates <- list()
  rows <- list()
  index <- 1L
  for (Mf in config$pflr_tuning$FLiRTI$Mf) {
    basis <- make_functional_basis(argvals, Mf, config$degree)
    design <- xc %*% sweep(basis$observed, 1L, weights, `*`)
    for (lambda in config$pflr_tuning$FLiRTI$lambda) {
      captured <- capture_method(PFLR:::flirty1_dz2(
        Y = yc, U = design, basis = basis$object,
        lambda = lambda, cons = config$pflr_tuning$FLiRTI$cons
      ))
      if (captured$success) {
        current <- captured$value
        coefficient <- tryCatch(
          as.numeric(solve(current$A, current$gamma)),
          error = function(condition) rep(NA_real_, ncol(design))
        )
        if (all(is.finite(coefficient))) {
          fitted <- as.numeric(design %*% coefficient)
          rss <- sum((yc - fitted)^2)
          bic <- length(y) * log(max(rss / length(y), .Machine$double.xmin)) +
            log(length(y)) * current$df
        } else {
          rss <- NA_real_
          bic <- Inf
        }
      } else {
        current <- NULL
        coefficient <- NULL
        rss <- NA_real_
        bic <- Inf
      }
      candidates[[index]] <- list(
        core = current, coefficient = coefficient,
        basis = basis$observed, design = design
      )
      rows[[index]] <- data.frame(
        Mf = Mf, lambda = lambda, bic = bic, rss = rss,
        df = if (captured$success) current$df else NA_real_,
        success = captured$success && is.finite(bic),
        warnings = captured$warnings, error = captured$error,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  tuning <- do.call(rbind, rows)
  valid <- which(tuning$success)
  if (!length(valid)) {
    captured <- list(elapsed = proc.time()[["elapsed"]] - start,
                     warnings = paste(unique(tuning$warnings), collapse = " | "),
                     error = paste(unique(tuning$error), collapse = " | "))
    return(failed_method("FLiRTI", length(argvals), captured, tuning))
  }
  selected_index <- valid[which.min(tuning$bic[valid])]
  selected_tuning <- tuning[selected_index, , drop = FALSE]
  selected <- candidates[[selected_index]]
  fit <- list(
    coefficient = selected$coefficient,
    beta = as.numeric(selected$basis %*% selected$coefficient),
    basis = selected$basis, weights = weights,
    x_mean = x_mean, y_mean = y_mean,
    tuning = tuning, selected = selected_tuning
  )
  list(
    method = "FLiRTI", fit = fit, beta = fit$beta,
    predict = function(newx) {
      centered_new <- sweep(newx, 2L, fit$x_mean, FUN = "-")
      design_new <- centered_new %*% sweep(fit$basis, 1L, fit$weights, `*`)
      as.numeric(fit$y_mean + design_new %*% fit$coefficient)
    },
    selected = data.frame(
      method = "FLiRTI", M = selected_tuning$Mf,
      fit_success = TRUE, bic = selected_tuning$bic,
      rss = selected_tuning$rss, lambda = selected_tuning$lambda,
      kappa = NA_real_, shape = NA_real_, df = selected_tuning$df,
      active = sum(abs(fit$coefficient) > config$zero_tolerance),
      converged = TRUE, iterations = 1L,
      warnings = paste(unique(tuning$warnings), collapse = " | "),
      error = "", stringsAsFactors = FALSE
    ),
    diagnostics = tuning,
    elapsed = unname(proc.time()[["elapsed"]] - start),
    warnings = paste(unique(tuning$warnings), collapse = " | "), error = ""
  )
}

failed_method <- function(method, grid_length, captured,
                          diagnostics = NULL) {
  list(
    method = method, fit = NULL, beta = rep(NA_real_, grid_length),
    predict = function(newx) rep(NA_real_, nrow(newx)),
    selected = NULL, diagnostics = diagnostics,
    elapsed = unname(captured$elapsed), warnings = captured$warnings,
    error = captured$error
  )
}

fit_single_scaled_method <- function(method, y, x, argvals, weights,
                                     config, cutoff = 0,
                                     mcp_gamma = NULL) {
  if (method %in% c("FMCP", "fSCAD", "Smooth", "OLS")) {
    return(common_engine_fit_over_M(
      y, x, argvals, config, method, cutoff, mcp_gamma
    ))
  }
  if (identical(method, "PCR")) {
    return(fit_pcr_bic(
      y, x, argvals, weights,
      max_components = min(config$pcr_max_components, length(y) - 2L)
    ))
  }
  if (identical(method, "PenS")) {
    return(fit_pens_trapezoid(y, x, argvals, weights, config))
  }
  if (identical(method, "NGR")) {
    return(fit_ngr_trapezoid(y, x, argvals, weights, config))
  }
  if (identical(method, "FLiRTI")) {
    return(fit_flirti_trapezoid(y, x, argvals, weights, config))
  }
  stop("Unknown method: ", method, call. = FALSE)
}

compact_real_fit <- function(fit, method, dataset, scaler, cutoff,
                             prediction_scaled = NULL) {
  success <- !is.null(fit$fit) && all(is.finite(fit$beta))
  prediction <- if (is.null(prediction_scaled)) NULL else
    prediction_scaled * scaler$y_scale
  selected <- fit$selected
  if (is.null(selected)) {
    selected <- data.frame(
      method = method, M = NA_integer_, fit_success = FALSE,
      bic = NA_real_, rss = NA_real_, lambda = NA_real_, kappa = NA_real_,
      shape = NA_real_, df = NA_real_, active = NA_real_,
      converged = FALSE, iterations = NA_real_,
      warnings = fit$warnings, error = fit$error,
      stringsAsFactors = FALSE
    )
  }
  selected$cutoff_scaled <- cutoff
  selected$y_scale <- scaler$y_scale
  selected$x_functional_rms_scale <- scaler$x_scale
  selected$integration <- "trapezoid"
  list(
    method = method, success = success,
    beta_scaled = fit$beta,
    beta_physical = if (success) {
      back_transform_beta(fit$beta, scaler, dataset$argvals_physical)
    } else rep(NA_real_, length(dataset$argvals_u)),
    prediction = prediction,
    selected = selected,
    diagnostics = fit$diagnostics,
    elapsed = fit$elapsed, warnings = fit$warnings, error = fit$error
  )
}

fit_real_method_set <- function(dataset, train_index, test_index = NULL,
                                config, methods = config$methods,
                                mcp_gamma = NULL) {
  train_y <- dataset$y[train_index]
  train_x <- dataset$x[train_index, , drop = FALSE]
  scaler <- fit_training_scaler(train_y, train_x, dataset$argvals_u)
  scaled <- apply_training_scaler(train_y, train_x, scaler)
  test_scaled_x <- if (is.null(test_index)) NULL else
    dataset$x[test_index, , drop = FALSE] / scaler$x_scale

  smooth_fit <- fit_single_scaled_method(
    "Smooth", scaled$y, scaled$x, dataset$argvals_u,
    scaler$weights, config, cutoff = 0
  )
  if (is.null(smooth_fit$fit) || any(!is.finite(smooth_fit$beta))) {
    cutoff <- 0
  } else {
    cutoff <- config$relative_cutoff_multiplier *
      integral_l2(smooth_fit$beta, scaler$weights)
  }

  results <- vector("list", length(methods))
  names(results) <- methods
  for (method in methods) {
    fit <- if (identical(method, "Smooth")) smooth_fit else
      fit_single_scaled_method(
        method, scaled$y, scaled$x, dataset$argvals_u,
        scaler$weights, config,
        cutoff = if (method %in% config$primary_sparse_methods) cutoff else 0,
        mcp_gamma = mcp_gamma
      )
    prediction_scaled <- if (is.null(test_scaled_x)) NULL else
      tryCatch(fit$predict(test_scaled_x), error = function(condition) {
        fit$error <- paste(fit$error, conditionMessage(condition), sep = " | ")
        rep(NA_real_, nrow(test_scaled_x))
      })
    results[[method]] <- compact_real_fit(
      fit, method, dataset, scaler, cutoff, prediction_scaled
    )
  }
  list(
    fits = results, scaler = scaler, cutoff_scaled = cutoff,
    train_index = train_index, test_index = test_index
  )
}

selection_rows <- function(fit_set, dataset_name, context = "full") {
  rows <- lapply(fit_set$fits, function(fit) {
    row <- fit$selected
    row$dataset <- dataset_name
    row$context <- context
    row$elapsed_seconds <- fit$elapsed
    row$success <- fit$success
    row
  })
  result <- bind_rows_base(rows)
  result[c("dataset", "context", setdiff(names(result), c("dataset", "context")))]
}

coefficient_rows <- function(fit_set, dataset, context = "full") {
  bind_rows_base(lapply(fit_set$fits, function(fit) {
    data.frame(
      dataset = dataset$name, context = context, method = fit$method,
      grid_index = seq_along(dataset$argvals_u),
      argval_u = dataset$argvals_u,
      argval_physical = dataset$argvals_physical,
      beta_scaled = fit$beta_scaled,
      beta_physical = fit$beta_physical,
      estimated_null = abs(fit$beta_scaled) <= sqrt(.Machine$double.eps),
      success = fit$success,
      stringsAsFactors = FALSE
    )
  }))
}

grid_boundary_diagnostics <- function(selection, config) {
  selection <- selection[
    selection$method %in% config$primary_sparse_methods & selection$success,
    , drop = FALSE
  ]
  parameter_grid <- list(M = config$M, lambda = config$lambda,
                         kappa = config$kappa)
  rows <- list()
  index <- 1L
  for (parameter in names(parameter_grid)) {
    values <- selection[[parameter]]
    values <- values[is.finite(values)]
    grid <- parameter_grid[[parameter]]
    for (side in c("lower", "upper")) {
      boundary <- if (identical(side, "lower")) min(grid) else max(grid)
      hits <- sum(abs(values - boundary) <=
                    10 * .Machine$double.eps * max(1, abs(boundary)))
      rows[[index]] <- data.frame(
        parameter = parameter, side = side, boundary = boundary,
        selections = length(values), hits = hits,
        rate = if (length(values)) hits / length(values) else NA_real_,
        exceeds_threshold = if (length(values)) {
          hits / length(values) > config$boundary_rate_threshold
        } else FALSE,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }
  do.call(rbind, rows)
}

expand_grid_once <- function(config, diagnostics) {
  resolved <- list(M = config$M, lambda = config$lambda,
                   kappa = config$kappa)
  applied <- character()
  for (parameter in c("M", "lambda", "kappa")) {
    current <- resolved[[parameter]]
    lower <- diagnostics$exceeds_threshold[
      diagnostics$parameter == parameter & diagnostics$side == "lower"
    ]
    upper <- diagnostics$exceeds_threshold[
      diagnostics$parameter == parameter & diagnostics$side == "upper"
    ]
    additions <- numeric()
    if (any(lower)) {
      additions <- c(additions, if (identical(parameter, "M")) {
        max(2L, min(current) - 5L)
      } else min(current) / 10)
      applied <- c(applied, paste0(parameter, "_lower"))
    }
    if (any(upper)) {
      additions <- c(additions, if (identical(parameter, "M")) {
        max(current) + 10L
      } else max(current) * 10)
      applied <- c(applied, paste0(parameter, "_upper"))
    }
    resolved[[parameter]] <- sort(unique(c(current, additions)))
    if (identical(parameter, "M")) {
      resolved[[parameter]] <- as.integer(resolved[[parameter]])
    }
  }
  resolved$expansion_applied <- length(applied) > 0L
  resolved$expanded_boundaries <- applied
  resolved$boundary_rate_threshold <- config$boundary_rate_threshold
  resolved$maximum_grid_expansions <- 1L
  resolved
}
