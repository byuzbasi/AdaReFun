# Efficient comparison methods for the version-3 FSCAD simulation.

select_fit_over_M <- function(y, x, argvals, M_grid, method,
                              lambda, kappa, cutoff = 0,
                              mcp_gamma = 3, scad_a = 3.7,
                              degree = 3L, penalty_order = 2L,
                              max_iter = 150L, tolerance = 1e-8,
                              norm_floor = 1e-8,
                              package_fscad = FALSE) {
  start <- proc.time()[["elapsed"]]
  candidates <- vector("list", length(M_grid))
  rows <- vector("list", length(M_grid))

  for (index in seq_along(M_grid)) {
    current_M <- as.integer(M_grid[index])
    captured <- if (identical(method, "FMCP")) {
      capture_method(FMCP::fmcp(
        y, x, argvals,
        M = current_M,
        lambda = lambda,
        kappa = kappa,
        degree = degree,
        penalty_order = penalty_order,
        mcp_gamma = mcp_gamma,
        integration = "simpson",
        control = FMCP::fmcp_control(
          max_iter = max_iter,
          tol = tolerance,
          norm_floor = norm_floor,
          cutoff = cutoff,
          refit = TRUE,
          trace = FALSE,
          keep_matrices = FALSE
        )
      ))
    } else if (identical(method, "fSCAD")) {
      if (isTRUE(package_fscad)) {
        capture_method(FMCP::fSCAD(
          y, x, argvals,
          M = current_M,
          lambda = lambda,
          kappa = kappa,
          scad_a = scad_a,
          degree = degree,
          penalty_order = penalty_order,
          integration = "simpson",
          control = FMCP::fmcp_control(
            max_iter = max_iter,
            tol = tolerance,
            norm_floor = norm_floor,
            cutoff = cutoff,
            refit = TRUE,
            trace = FALSE,
            keep_matrices = FALSE
          )
        ))
      } else {
        capture_method(fit_controlled_fscad(
          y, x, argvals,
          M = current_M,
          lambda = lambda,
          kappa = kappa,
          scad_a = scad_a,
          degree = degree,
          penalty_order = penalty_order,
          integration = "simpson",
          max_iter = max_iter,
          tolerance = tolerance,
          norm_floor = norm_floor,
          cutoff = cutoff,
          refit = TRUE,
          keep_trace = FALSE
        ))
      }
    } else if (identical(method, "Smooth")) {
      capture_method(FMCP::fmcp(
        y, x, argvals,
        M = current_M,
        lambda = 0,
        kappa = kappa,
        degree = degree,
        penalty_order = penalty_order,
        mcp_gamma = mcp_gamma,
        integration = "simpson",
        control = FMCP::fmcp_control(
          max_iter = 1L,
          tol = tolerance,
          norm_floor = norm_floor,
          cutoff = 0,
          refit = FALSE,
          keep_matrices = FALSE
        )
      ))
    } else if (identical(method, "OLS")) {
      capture_method(FMCP::fmcp(
        y, x, argvals,
        M = current_M,
        lambda = 0,
        kappa = 0,
        degree = degree,
        penalty_order = penalty_order,
        mcp_gamma = mcp_gamma,
        integration = "simpson",
        control = FMCP::fmcp_control(
          max_iter = 1L,
          tol = tolerance,
          norm_floor = norm_floor,
          cutoff = 0,
          refit = FALSE,
          keep_matrices = FALSE
        )
      ))
    } else {
      stop("Unknown common-engine method: ", method, call. = FALSE)
    }

    candidates[[index]] <- captured
    if (captured$success) {
      selected <- captured$value$selected[1L, , drop = FALSE]
      rows[[index]] <- data.frame(
        method = method,
        M = current_M,
        fit_success = TRUE,
        bic = selected$bic,
        lambda = selected$lambda,
        kappa = selected$kappa,
        shape = if (identical(method, "FMCP")) {
          selected$mcp_gamma
        } else if (identical(method, "fSCAD")) {
          selected$scad_a
        } else {
          NA_real_
        },
        df = selected$df,
        active = selected$active,
        converged = selected$converged,
        iterations = selected$iterations,
        warnings = captured$warnings,
        error = captured$error,
        stringsAsFactors = FALSE
      )
    } else {
      rows[[index]] <- data.frame(
        method = method,
        M = current_M,
        fit_success = FALSE,
        bic = NA_real_,
        lambda = NA_real_,
        kappa = NA_real_,
        shape = NA_real_,
        df = NA_real_,
        active = NA_integer_,
        converged = FALSE,
        iterations = NA_integer_,
        warnings = captured$warnings,
        error = captured$error,
        stringsAsFactors = FALSE
      )
    }
  }

  diagnostics <- do.call(rbind, rows)
  valid <- which(diagnostics$fit_success & is.finite(diagnostics$bic))
  elapsed <- proc.time()[["elapsed"]] - start
  if (!length(valid)) {
    return(list(
      method = method,
      fit = NULL,
      beta = rep(NA_real_, length(argvals)),
      predict = function(newx) rep(NA_real_, nrow(newx)),
      selected = NULL,
      diagnostics = diagnostics,
      elapsed = unname(elapsed),
      warnings = paste(unique(diagnostics$warnings), collapse = " | "),
      error = paste(unique(diagnostics$error), collapse = " | ")
    ))
  }
  chosen <- valid[which.min(diagnostics$bic[valid])]
  fit <- candidates[[chosen]]$value
  list(
    method = method,
    fit = fit,
    beta = as.numeric(fit$beta),
    predict = function(newx) as.numeric(stats::predict(fit, newx)),
    selected = diagnostics[chosen, , drop = FALSE],
    diagnostics = diagnostics,
    elapsed = unname(elapsed),
    warnings = candidates[[chosen]]$warnings,
    error = candidates[[chosen]]$error
  )
}

fit_pcr_bic <- function(y, x, argvals, integration_weights,
                        max_components = 30L) {
  captured <- capture_method({
    y_mean <- mean(y)
    x_mean <- colMeans(x)
    yc <- y - y_mean
    xc <- sweep(x, 2L, x_mean, FUN = "-")
    sqrt_weights <- sqrt(integration_weights)
    weighted_x <- sweep(xc, 2L, sqrt_weights, FUN = "*")
    decomposition <- svd(weighted_x, nu = 0L, nv = min(
      as.integer(max_components), nrow(x) - 1L, ncol(x)
    ))
    maximum <- ncol(decomposition$v)
    scores <- weighted_x %*% decomposition$v
    tuning <- vector("list", maximum)
    coefficients <- vector("list", maximum)
    for (components in seq_len(maximum)) {
      current_scores <- scores[, seq_len(components), drop = FALSE]
      coefficient <- qr.solve(current_scores, yc, tol = 1e-10)
      residual <- yc - as.numeric(current_scores %*% coefficient)
      rss <- sum(residual^2)
      bic <- length(y) * log(max(rss / length(y), .Machine$double.xmin)) +
        log(length(y)) * components
      tuning[[components]] <- data.frame(
        components = components, bic = bic, rss = rss,
        stringsAsFactors = FALSE
      )
      coefficients[[components]] <- coefficient
    }
    tuning <- do.call(rbind, tuning)
    selected_index <- which.min(tuning$bic)
    components <- tuning$components[selected_index]
    regression_coefficient <- coefficients[[selected_index]]
    weighted_beta <- as.numeric(
      decomposition$v[, seq_len(components), drop = FALSE] %*%
        regression_coefficient
    )
    beta <- weighted_beta / sqrt_weights
    list(
      beta = beta,
      y_mean = y_mean,
      x_mean = x_mean,
      sqrt_weights = sqrt_weights,
      rotation = decomposition$v[, seq_len(components), drop = FALSE],
      regression_coefficient = regression_coefficient,
      tuning = tuning,
      selected = tuning[selected_index, , drop = FALSE]
    )
  })

  if (!captured$success) {
    return(list(
      method = "PCR", fit = NULL,
      beta = rep(NA_real_, length(argvals)),
      predict = function(newx) rep(NA_real_, nrow(newx)),
      selected = NULL, diagnostics = NULL,
      elapsed = captured$elapsed, warnings = captured$warnings,
      error = captured$error
    ))
  }
  fit <- captured$value
  list(
    method = "PCR",
    fit = fit,
    beta = fit$beta,
    predict = function(newx) {
      centered <- sweep(newx, 2L, fit$x_mean, FUN = "-")
      weighted <- sweep(centered, 2L, fit$sqrt_weights, FUN = "*")
      scores <- weighted %*% fit$rotation
      as.numeric(fit$y_mean + scores %*% fit$regression_coefficient)
    },
    selected = data.frame(
      method = "PCR", M = NA_integer_, fit_success = TRUE,
      bic = fit$selected$bic, lambda = NA_real_, kappa = NA_real_,
      shape = NA_real_, df = fit$selected$components,
      active = fit$selected$components, converged = TRUE,
      iterations = 1L, warnings = captured$warnings, error = "",
      stringsAsFactors = FALSE
    ),
    diagnostics = fit$tuning,
    elapsed = captured$elapsed,
    warnings = captured$warnings,
    error = captured$error
  )
}

make_block_basis <- function(argvals, blocks, M, degree = 3L) {
  total_length <- sum(vapply(blocks, function(block) diff(block), numeric(1L)))
  basis_parts <- list()
  roughness_parts <- list()
  interval_counts <- pmax(
    2L,
    as.integer(round(M * vapply(
      blocks, function(block) diff(block) / total_length, numeric(1L)
    )))
  )
  rule <- FMCP:::.fmcp_gauss_legendre(max(degree + 1L, 4L))

  for (index in seq_along(blocks)) {
    block <- blocks[[index]]
    interval_count <- interval_counts[index]
    breaks <- seq(block[1L], block[2L], length.out = interval_count + 1L)
    order <- degree + 1L
    full_knots <- c(
      rep(block[1L], order),
      breaks[-c(1L, length(breaks))],
      rep(block[2L], order)
    )
    observed <- splines::splineDesign(
      full_knots,
      pmin(pmax(argvals, block[1L]), block[2L]),
      ord = order,
      derivs = rep.int(0L, length(argvals)),
      outer.ok = TRUE
    )
    outside <- argvals < block[1L] | argvals > block[2L]
    observed[outside, ] <- 0

    quad_values <- numeric(interval_count * length(rule$nodes))
    quad_weights <- numeric(length(quad_values))
    offset <- 0L
    for (j in seq_len(interval_count)) {
      indices <- offset + seq_along(rule$nodes)
      half_width <- diff(breaks[j + 0:1]) / 2
      midpoint <- mean(breaks[j + 0:1])
      quad_values[indices] <- midpoint + half_width * rule$nodes
      quad_weights[indices] <- half_width * rule$weights
      offset <- offset + length(rule$nodes)
    }
    derivative <- splines::splineDesign(
      full_knots, quad_values, ord = order,
      derivs = rep.int(2L, length(quad_values)), outer.ok = TRUE
    )
    basis_parts[[index]] <- unclass(observed)
    roughness_parts[[index]] <- FMCP:::cpp_weighted_crossprod(
      unclass(derivative), quad_weights
    )
  }
  basis <- do.call(cbind, basis_parts)
  roughness <- matrix(0, ncol(basis), ncol(basis))
  offset <- 0L
  for (part in roughness_parts) {
    indices <- offset + seq_len(nrow(part))
    roughness[indices, indices] <- part
    offset <- offset + nrow(part)
  }
  list(basis = basis, roughness = roughness, intervals = interval_counts)
}

fit_oracle_bic <- function(y, x, argvals, integration_weights,
                           scenario, M_grid, kappa_grid,
                           degree = 3L) {
  captured <- capture_method({
    y_mean <- mean(y)
    x_mean <- colMeans(x)
    yc <- y - y_mean
    xc <- sweep(x, 2L, x_mean, FUN = "-")

    if (identical(scenario, "case_i_null")) {
      fit <- list(
        beta = numeric(length(argvals)), y_mean = y_mean,
        x_mean = x_mean, coefficients = numeric(), basis = NULL
      )
      selected <- data.frame(
        method = "Oracle", M = 0L, fit_success = TRUE,
        bic = length(y) * log(max(mean(yc^2), .Machine$double.xmin)),
        lambda = NA_real_, kappa = NA_real_, shape = NA_real_,
        df = 0, active = 0L, converged = TRUE, iterations = 1L,
        warnings = "", error = "", stringsAsFactors = FALSE
      )
      list(fit = fit, selected = selected, tuning = selected)
    } else if (identical(scenario, "case_ii_sparse")) {
      tuning_rows <- list()
      candidate_fits <- list()
      candidate <- 1L
      for (M in M_grid) {
        block <- make_block_basis(
          argvals, list(c(0, 0.3), c(0.7, 1)), M, degree
        )
        design <- FMCP:::cpp_design_matrix(
          xc, block$basis, integration_weights
        )
        crossprod <- crossprod(design)
        response <- crossprod(design, yc)
        for (kappa in kappa_grid) {
          system <- crossprod + length(y) * kappa * block$roughness
          coefficient <- tryCatch(
            as.numeric(solve(system, response)),
            error = function(condition) rep(NA_real_, ncol(design))
          )
          if (all(is.finite(coefficient))) {
            residual <- yc - as.numeric(design %*% coefficient)
            rss <- sum(residual^2)
            influence <- tryCatch(
              solve(system, crossprod),
              error = function(condition) NULL
            )
            df <- if (is.null(influence)) NA_real_ else sum(diag(influence))
            bic <- if (is.finite(df)) {
              length(y) * log(max(rss / length(y), .Machine$double.xmin)) +
                log(length(y)) * df
            } else {
              Inf
            }
          } else {
            rss <- df <- NA_real_
            bic <- Inf
          }
          tuning_rows[[candidate]] <- data.frame(
            method = "Oracle", M = M, fit_success = is.finite(bic),
            bic = bic, lambda = NA_real_, kappa = kappa,
            shape = NA_real_, df = df,
            active = if (all(is.finite(coefficient))) {
              sum(coefficient != 0)
            } else {
              NA_integer_
            },
            converged = is.finite(bic), iterations = 1L,
            warnings = "", error = "", stringsAsFactors = FALSE
          )
          candidate_fits[[candidate]] <- list(
            coefficient = coefficient, basis = block$basis,
            y_mean = y_mean, x_mean = x_mean
          )
          candidate <- candidate + 1L
        }
      }
      tuning <- do.call(rbind, tuning_rows)
      selected_index <- which.min(tuning$bic)
      fit <- candidate_fits[[selected_index]]
      fit$beta <- as.numeric(fit$basis %*% fit$coefficient)
      list(
        fit = fit,
        selected = tuning[selected_index, , drop = FALSE],
        tuning = tuning
      )
    } else {
      list(fit = NULL, selected = NULL, tuning = NULL)
    }
  })

  if (!captured$success || is.null(captured$value$fit)) {
    return(list(
      method = "Oracle", fit = NULL,
      beta = rep(NA_real_, length(argvals)),
      predict = function(newx) rep(NA_real_, nrow(newx)),
      selected = NULL,
      diagnostics = if (captured$success) captured$value$tuning else NULL,
      elapsed = captured$elapsed, warnings = captured$warnings,
      error = if (captured$success) "Not applicable" else captured$error
    ))
  }
  fit <- captured$value$fit
  list(
    method = "Oracle", fit = fit, beta = fit$beta,
    predict = function(newx) {
      if (is.null(fit$basis)) {
        return(rep(fit$y_mean, nrow(newx)))
      }
      centered <- sweep(newx, 2L, fit$x_mean, FUN = "-")
      design <- FMCP:::cpp_design_matrix(
        centered, fit$basis, integration_weights
      )
      as.numeric(fit$y_mean + design %*% fit$coefficient)
    },
    selected = captured$value$selected,
    diagnostics = captured$value$tuning,
    elapsed = captured$elapsed,
    warnings = captured$warnings,
    error = captured$error
  )
}

wrap_pflr_comparator <- function(method, y, x, argvals, M, tuning) {
  captured <- capture_method(fit_pflr_method(
    method, y, x, argvals, M = M, degree = 3L, tuning = tuning
  ))
  if (!captured$success) {
    return(list(
      method = method, fit = NULL,
      beta = rep(NA_real_, length(argvals)),
      predict = function(newx) rep(NA_real_, nrow(newx)),
      selected = NULL, diagnostics = NULL,
      elapsed = captured$elapsed, warnings = captured$warnings,
      error = captured$error
    ))
  }
  fit <- captured$value
  selected <- if (identical(method, "NGR")) {
    data.frame(
      method = method, M = M, fit_success = TRUE, bic = NA_real_,
      lambda = fit$extra$Opttau, kappa = fit$extra$Optkappa,
      shape = 0.5, df = sum(fit$extra$b != 0),
      active = sum(fit$extra$b != 0), converged = TRUE,
      iterations = if (is.null(tuning$niter)) 100L else tuning$niter,
      warnings = captured$warnings, error = "", stringsAsFactors = FALSE
    )
  } else if (identical(method, "PenS")) {
    data.frame(
      method = method, M = M, fit_success = TRUE, bic = NA_real_,
      lambda = NA_real_, kappa = fit$extra$Optalpha,
      shape = NA_real_, df = sum(fit$extra$b != 0),
      active = sum(fit$extra$b != 0), converged = TRUE,
      iterations = 1L, warnings = captured$warnings, error = "",
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      method = method, M = fit$extra$OptM, fit_success = TRUE,
      bic = NA_real_, lambda = fit$extra$Optlambda, kappa = NA_real_,
      shape = NA_real_, df = NA_real_, active = sum(fit$beta == 0),
      converged = TRUE, iterations = 1L, warnings = captured$warnings,
      error = "", stringsAsFactors = FALSE
    )
  }
  list(
    method = method,
    fit = fit,
    beta = fit$beta_common,
    predict = function(newx) predict_pflr_method(fit, newx),
    selected = selected,
    diagnostics = NULL,
    elapsed = captured$elapsed,
    warnings = captured$warnings,
    error = captured$error
  )
}
