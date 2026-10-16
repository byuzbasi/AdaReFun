# Common-design adapters for the version-5 MCP-versus-SCAD penalty-geometry
# study.  This file is intentionally separate from both package cores and the
# frozen version-4 comparison workflow.

penalty_path_require_context <- function() {
  required <- c("prepare_common_design", "capture_method")
  missing <- required[!vapply(
    required, exists, logical(1L), mode = "function", inherits = TRUE
  )]
  if (length(missing)) {
    stop(
      "Source fscad_and_pflr_adapters.R before penalty_path_adapters.R. ",
      "Missing: ", paste(missing, collapse = ", "), call. = FALSE
    )
  }
  if (!requireNamespace("FMCP", quietly = TRUE)) {
    stop("The installed FMCP package is required.", call. = FALSE)
  }
  invisible(TRUE)
}

penalty_grid_selected_index <- function(tuning) {
  valid <- which(tuning$success %in% TRUE & is.finite(tuning$bic))
  if (!length(valid)) return(NA_integer_)
  converged <- valid[tuning$converged[valid] %in% TRUE]
  eligible <- if (length(converged)) converged else valid
  eligible[which.min(tuning$bic[eligible])]
}

penalty_grid_tuning <- function(grid, method, M, shape, elapsed_seconds,
                                warnings = "", error = "") {
  if (is.null(grid)) {
    return(data.frame(
      method = method, M = as.integer(M), candidate = NA_integer_,
      lambda = NA_real_, kappa = NA_real_, shape = as.numeric(shape),
      bic = NA_real_, rss = NA_real_, df = NA_real_, objective = NA_real_,
      iterations = NA_integer_, active = NA_integer_, success = FALSE,
      converged = FALSE, refit_success = FALSE,
      objective_increased = NA, message = error,
      grid_elapsed_seconds = elapsed_seconds, warnings = warnings,
      error = error, stringsAsFactors = FALSE
    ))
  }
  n_candidates <- length(grid$lambda)
  data.frame(
    method = method,
    M = rep(as.integer(M), n_candidates),
    candidate = seq_len(n_candidates),
    lambda = as.numeric(grid$lambda),
    kappa = as.numeric(grid$kappa),
    shape = rep(as.numeric(shape), n_candidates),
    bic = as.numeric(grid$bic),
    rss = as.numeric(grid$rss),
    df = as.numeric(grid$df),
    objective = as.numeric(grid$objective),
    iterations = as.integer(grid$iterations),
    active = as.integer(grid$active),
    success = as.logical(grid$success),
    converged = as.logical(grid$converged),
    refit_success = as.logical(grid$refit_success),
    objective_increased = as.logical(grid$objective_increased),
    message = as.character(grid$message),
    grid_elapsed_seconds = rep(elapsed_seconds, n_candidates),
    warnings = rep(warnings, n_candidates),
    error = rep(error, n_candidates),
    stringsAsFactors = FALSE
  )
}

fit_common_smooth_reference <- function(
    y, x, argvals, M, kappa, degree = 3L, penalty_order = 2L,
    integration = "simpson") {
  penalty_path_require_context()
  M <- unique(as.integer(M))
  kappa <- unique(as.numeric(kappa))
  if (!length(M) || any(!is.finite(M)) || any(M < 1L)) {
    stop("`M` must contain positive finite integers.", call. = FALSE)
  }
  if (!length(kappa) || any(!is.finite(kappa)) || any(kappa < 0)) {
    stop("`kappa` must contain nonnegative finite values.", call. = FALSE)
  }

  prepared_by_M <- setNames(vector("list", length(M)), as.character(M))
  fit_store <- list()
  tuning_rows <- list()
  row_index <- 1L
  fit_index <- 1L

  for (current_M in M) {
    prepared <- prepare_common_design(
      y, x, argvals, current_M, degree, penalty_order, integration
    )
    prepared_by_M[[as.character(current_M)]] <- prepared
    xtx <- crossprod(prepared$design)
    xty <- crossprod(prepared$design, prepared$centered_y)
    n <- length(prepared$centered_y)

    for (current_kappa in kappa) {
      captured <- capture_method({
        system <- xtx + n * current_kappa * prepared$roughness
        coefficient <- as.numeric(solve(system, xty))
        fitted <- as.numeric(prepared$design %*% coefficient)
        residual <- prepared$centered_y - fitted
        rss <- sum(residual^2)
        effective_df <- sum(diag(solve(system, xtx)))
        bic <- if (is.finite(rss) && rss > 0 && is.finite(effective_df)) {
          n * log(rss / n) + log(n) * effective_df
        } else {
          Inf
        }
        list(
          coefficient = coefficient, rss = rss, df = effective_df, bic = bic
        )
      })
      fit <- captured$value
      success <- captured$success && !is.null(fit) && is.finite(fit$bic)
      if (success) fit_store[[fit_index]] <- fit
      tuning_rows[[row_index]] <- data.frame(
        M = current_M,
        kappa = current_kappa,
        bic = if (success) fit$bic else Inf,
        rss = if (success) fit$rss else NA_real_,
        df = if (success) fit$df else NA_real_,
        success = success,
        converged = success,
        elapsed_seconds = captured$elapsed,
        warnings = captured$warnings,
        error = captured$error,
        fit_index = if (success) fit_index else NA_integer_,
        stringsAsFactors = FALSE
      )
      if (success) fit_index <- fit_index + 1L
      row_index <- row_index + 1L
    }
  }

  tuning <- do.call(rbind, tuning_rows)
  selected_index <- penalty_grid_selected_index(tuning)
  if (!is.finite(selected_index)) {
    stop("No common smooth-reference fit had a finite BIC.", call. = FALSE)
  }
  selected_row <- tuning[selected_index, , drop = FALSE]
  selected_fit <- fit_store[[selected_row$fit_index]]
  selected_prepared <- prepared_by_M[[as.character(selected_row$M)]]
  selected_beta <- as.numeric(
    selected_prepared$basis$basis_observed %*% selected_fit$coefficient
  )

  selected_by_M <- setNames(vector("list", length(M)), as.character(M))
  for (current_M in M) {
    candidates <- which(
      tuning$M == current_M & tuning$success %in% TRUE & is.finite(tuning$bic)
    )
    if (!length(candidates)) next
    chosen <- candidates[which.min(tuning$bic[candidates])]
    row <- tuning[chosen, , drop = FALSE]
    selected_by_M[[as.character(current_M)]] <- list(
      tuning = row,
      coefficient = fit_store[[row$fit_index]]$coefficient
    )
  }

  list(
    tuning = tuning,
    selected = selected_row,
    coefficient = selected_fit$coefficient,
    beta = selected_beta,
    prepared_by_M = prepared_by_M,
    selected_by_M = selected_by_M
  )
}

smooth_reference_cutoff <- function(reference, multiplier = 1e-4) {
  if (!is.list(reference) || is.null(reference$beta) ||
      is.null(reference$prepared_by_M)) {
    stop("`reference` is not a common smooth-reference result.", call. = FALSE)
  }
  selected_M <- as.character(reference$selected$M[1L])
  prepared <- reference$prepared_by_M[[selected_M]]
  domain_length <- diff(range(prepared$argvals))
  beta_norm <- sqrt(sum(
    prepared$integration_weights * as.numeric(reference$beta)^2
  ))
  as.numeric(multiplier) * beta_norm / domain_length
}

run_matched_penalty_grid <- function(
    prepared, lambda, kappa, mcp_gamma = 3.7, scad_a = 3.7,
    max_iter = 120L, tolerance = 1e-8, norm_floor = 1e-8,
    cutoff = 0, refit = TRUE, keep_trace = FALSE) {
  penalty_path_require_context()
  if (!exists("cpp_fscad_grid", mode = "function", inherits = TRUE)) {
    stop("Compile the controlled fSCAD comparator before grid fitting.",
         call. = FALSE)
  }
  lambda <- as.numeric(lambda)
  kappa <- as.numeric(kappa)
  if (!length(lambda) || any(!is.finite(lambda)) || any(lambda < 0)) {
    stop("`lambda` must contain nonnegative finite values.", call. = FALSE)
  }
  if (!length(kappa) || any(!is.finite(kappa)) || any(kappa < 0)) {
    stop("`kappa` must contain nonnegative finite values.", call. = FALSE)
  }
  M <- as.integer(prepared$basis$M)
  if (!length(M) || !is.finite(M) || M < 1L) {
    M <- as.integer(ncol(prepared$basis$basis_observed) - 3L)
  }

  mcp_gamma <- unique(as.numeric(mcp_gamma))
  if (!length(mcp_gamma) || any(!is.finite(mcp_gamma)) ||
      any(mcp_gamma <= 1)) {
    stop("`mcp_gamma` must contain finite values greater than one.",
         call. = FALSE)
  }
  fmcp <- lapply(mcp_gamma, function(current_gamma) {
    capture_method(FMCP:::cpp_fmcp_grid(
      prepared$centered_y, prepared$design, prepared$roughness,
      prepared$local_gram, lambda, kappa, current_gamma,
      prepared$interval_scale, as.integer(max_iter), as.numeric(tolerance),
      as.numeric(norm_floor), as.numeric(cutoff), isTRUE(refit),
      isTRUE(keep_trace)
    ))
  })
  fscad <- capture_method(cpp_fscad_grid(
    prepared$centered_y, prepared$design, prepared$roughness,
    prepared$local_gram, lambda, kappa, as.numeric(scad_a),
    prepared$interval_scale, as.integer(max_iter), as.numeric(tolerance),
    as.numeric(norm_floor), as.numeric(cutoff), isTRUE(refit),
    isTRUE(keep_trace)
  ))

  fscad_tuning <- penalty_grid_tuning(
    fscad$value, "fSCAD", M, scad_a, fscad$elapsed,
    fscad$warnings, fscad$error
  )
  components <- list()
  for (index in seq_along(mcp_gamma)) {
    current_gamma <- mcp_gamma[index]
    method <- sprintf("FMCP_gamma_%.1f", current_gamma)
    captured <- fmcp[[index]]
    tuning <- penalty_grid_tuning(
      captured$value, method, M, current_gamma, captured$elapsed,
      captured$warnings, captured$error
    )
    if (captured$success && fscad$success) {
      same_lambda <- isTRUE(all.equal(
        as.numeric(captured$value$lambda), as.numeric(fscad$value$lambda),
        tolerance = 0
      ))
      same_kappa <- isTRUE(all.equal(
        as.numeric(captured$value$kappa), as.numeric(fscad$value$kappa),
        tolerance = 0
      ))
      if (!same_lambda || !same_kappa) {
        stop("FMCP and fSCAD returned different candidate-grid ordering.",
             call. = FALSE)
      }
    }
    components[[method]] <- list(
      grid = captured$value, tuning = tuning, capture = captured,
      shape = current_gamma
    )
  }
  components$fSCAD <- list(
      grid = fscad$value, tuning = fscad_tuning, capture = fscad,
      shape = scad_a
  )
  components
}
