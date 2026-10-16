.arfc_validate_levels <- function(levels, level_weights, resolution) {
  if (!is.numeric(levels) || !length(levels) || anyNA(levels) ||
      any(!is.finite(levels)) || any(levels < 2) ||
      any(levels != round(levels))) {
    stop("`levels` must contain finite integers of at least two.", call. = FALSE)
  }
  levels <- as.integer(levels)
  if (is.unsorted(levels, strictly = TRUE) || anyDuplicated(levels)) {
    stop("`levels` must be strictly increasing and unique.", call. = FALSE)
  }
  finest <- max(levels)
  if (any(finest %% levels != 0L)) {
    stop("Every resolution level must divide the finest level.", call. = FALSE)
  }
  if (identical(resolution, "multi") && length(levels) < 2L) {
    stop("Multi-resolution ARFC requires at least two `levels`.", call. = FALSE)
  }
  if (is.null(level_weights)) {
    level_weights <- rep.int(1 / length(levels), length(levels))
  }
  if (!is.numeric(level_weights) || length(level_weights) != length(levels) ||
      anyNA(level_weights) || any(!is.finite(level_weights)) ||
      any(level_weights <= 0)) {
    stop(
      "`level_weights` must contain one finite positive value per level.",
      call. = FALSE
    )
  }
  if (abs(sum(level_weights) - 1) > 1e-12) {
    stop("`level_weights` must sum to one; they are not normalized silently.",
         call. = FALSE)
  }
  list(
    levels = levels,
    level_weights = as.numeric(level_weights),
    finest = finest
  )
}

.arfc_build_region_tree <- function(basis, levels, level_weights) {
  basis_count <- ncol(basis$basis_quad)
  region_count <- sum(levels)
  gram <- array(0, dim = c(basis_count, basis_count, region_count))
  metadata <- vector("list", region_count)
  domain_length <- diff(basis$domain)
  cursor <- 0L
  for (level_index in seq_along(levels)) {
    level_M <- levels[level_index]
    descendants <- basis$M %/% level_M
    for (region_index in seq_len(level_M)) {
      cursor <- cursor + 1L
      first_fine <- (region_index - 1L) * descendants + 1L
      last_fine <- region_index * descendants
      rows <- basis$quad_intervals >= first_fine &
        basis$quad_intervals <= last_fine
      regional_basis <- basis$basis_quad[rows, , drop = FALSE]
      regional_weights <- basis$quad_weights[rows]
      gram[, , cursor] <- crossprod(
        regional_basis,
        regional_basis * regional_weights
      )
      left <- basis$breaks[first_fine]
      right <- basis$breaks[last_fine + 1L]
      length <- right - left
      metadata[[cursor]] <- data.frame(
        region_id = cursor,
        level_index = level_index,
        level_M = level_M,
        region_index = region_index,
        first_fine_interval = first_fine,
        last_fine_interval = last_fine,
        left = left,
        right = right,
        midpoint = (left + right) / 2,
        region_length = length,
        level_weight = level_weights[level_index],
        canonical_multiplier = level_weights[level_index] *
          length / domain_length,
        stringsAsFactors = FALSE
      )
    }
  }
  list(gram = gram, metadata = do.call(rbind, metadata))
}

.arfc_prepare_context <- function(y,
                                  x,
                                  argvals,
                                  levels,
                                  level_weights,
                                  degree,
                                  penalty_order,
                                  integration) {
  validated <- .fmcp_validate_data(y, x, argvals)
  y <- validated$y
  x <- validated$x
  argvals <- validated$argvals
  basis <- .fmcp_basis_setup(
    argvals, max(levels), degree, penalty_order
  )
  integration_weights <- .fmcp_integration_weights(argvals, integration)
  y_mean <- mean(y)
  x_mean <- colMeans(x)
  centered_y <- y - y_mean
  centered_x <- sweep(x, 2L, x_mean, FUN = "-")
  design <- cpp_design_matrix(
    centered_x, basis$basis_observed, integration_weights
  )
  roughness <- cpp_weighted_crossprod(
    basis$derivative_quad, basis$quad_weights
  )
  local_gram <- cpp_local_gram(
    basis$basis_quad,
    basis$quad_weights,
    basis$quad_intervals,
    basis$M,
    basis$order
  )
  region_tree <- .arfc_build_region_tree(basis, levels, level_weights)
  fine <- region_tree$metadata$level_M == basis$M
  domain_length <- diff(basis$domain)
  fine_length <- region_tree$metadata$region_length[fine]
  list(
    y = y,
    x = x,
    argvals = argvals,
    centered_y = centered_y,
    design = design,
    roughness = roughness,
    local_gram = local_gram,
    basis = basis,
    penalty_order = penalty_order,
    integration = integration,
    integration_weights = integration_weights,
    means = list(y = y_mean, x = x_mean),
    region_tree = region_tree,
    fine = fine,
    fine_gram = region_tree$gram[, , fine, drop = FALSE],
    fine_length = fine_length,
    fine_multiplier = fine_length / domain_length,
    multi_multiplier = region_tree$metadata$canonical_multiplier,
    interval_scale = sqrt(basis$M / domain_length)
  )
}

.arfc_eligible_candidates <- function(grid, strict = TRUE) {
  finite_coefficients <- apply(
    grid$coefficients, 2L, function(value) all(is.finite(value))
  )
  usable <- grid$success %in% TRUE & is.finite(grid$bic) &
    finite_coefficients
  if (isTRUE(strict)) {
    usable <- usable & grid$converged %in% TRUE &
      grid$refit_success %in% TRUE
  }
  which(usable)
}

.arfc_select_bic <- function(grid, strict = TRUE, label = "ARFC") {
  eligible <- .arfc_eligible_candidates(grid, strict)
  if (!length(eligible)) {
    details <- unique(as.character(grid$message[nzchar(grid$message)]))
    stop(
      "No eligible ", label, " tuning candidate was available.",
      if (length(details)) paste0(" ", paste(details, collapse = " ")) else "",
      call. = FALSE
    )
  }
  eligible[which.min(grid$bic[eligible])]
}

.arfc_smooth_pilot <- function(context, kappa, mcp_gamma, control) {
  grid <- cpp_fmcp_grid(
    context$centered_y,
    context$design,
    context$roughness,
    context$local_gram,
    0,
    kappa,
    mcp_gamma,
    context$interval_scale,
    control$max_iter,
    control$tol,
    control$norm_floor,
    0,
    control$refit,
    control$trace
  )
  selected <- .arfc_select_bic(grid, control$strict, "smooth-pilot")
  tuning <- data.frame(
    candidate_index = seq_along(grid$kappa),
    kappa = as.numeric(grid$kappa),
    bic = as.numeric(grid$bic),
    rss = as.numeric(grid$rss),
    df = as.numeric(grid$df),
    success = as.logical(grid$success),
    converged = as.logical(grid$converged),
    refit_success = as.logical(grid$refit_success),
    selected = seq_along(grid$kappa) == selected,
    message = as.character(grid$message),
    stringsAsFactors = FALSE
  )
  list(
    coefficients = as.numeric(grid$coefficients[, selected]),
    selected_kappa = as.numeric(grid$kappa[selected]),
    selected_index = selected,
    tuning = tuning,
    grid = if (control$keep_matrices) grid else NULL
  )
}

.arfc_adaptive_weights <- function(pilot, context, adaptive, control) {
  metadata <- context$region_tree$metadata
  norms <- as.numeric(cpp_arfc_region_norms_v1(
    pilot, context$region_tree$gram, metadata$region_length
  ))
  if (!isTRUE(adaptive)) {
    return(list(
      pilot_norm = norms,
      epsilon = NA_real_,
      raw = rep.int(1, length(norms)),
      normalized = rep.int(1, length(norms)),
      weights = rep.int(1, length(norms))
    ))
  }
  positive <- norms[is.finite(norms) & norms > 0]
  if (!length(positive)) {
    stop(
      "The smooth pilot has no positive regional norm; adaptive weights are undefined.",
      call. = FALSE
    )
  }
  epsilon <- max(
    control$norm_floor,
    control$adaptive_epsilon_fraction * stats::median(positive)
  )
  raw <- (norms + epsilon)^(-control$adaptive_alpha)
  normalized <- raw
  for (level_M in unique(metadata$level_M)) {
    index <- which(metadata$level_M == level_M)
    normalized[index] <- raw[index] / stats::median(raw[index])
  }
  weights <- pmin(
    control$adaptive_clip[2L],
    pmax(control$adaptive_clip[1L], normalized)
  )
  list(
    pilot_norm = norms,
    epsilon = epsilon,
    raw = raw,
    normalized = normalized,
    weights = weights
  )
}

.arfc_method_code <- function(resolution, penalty) {
  paste0(if (resolution == "single") "asr" else "amr", "_",
         tolower(penalty))
}

.arfc_method_label <- function(resolution, penalty) {
  paste0(if (resolution == "single") "ASR" else "AMR", "--", penalty)
}

.arfc_fit_prepared <- function(context,
                               resolution,
                               penalty,
                               lambda,
                               mcp_gamma,
                               scad_a,
                               adaptive,
                               pilot,
                               adaptive_weights,
                               control,
                               call = NULL) {
  multi <- identical(resolution, "multi")
  family <- tolower(penalty)
  shape <- if (identical(penalty, "MCP")) mcp_gamma else scad_a
  region_gram <- if (multi) context$region_tree$gram else context$fine_gram
  region_metadata <- if (multi) {
    context$region_tree$metadata
  } else {
    context$region_tree$metadata[context$fine, , drop = FALSE]
  }
  region_length <- region_metadata$region_length
  region_multiplier <- if (multi) {
    context$multi_multiplier
  } else {
    context$fine_multiplier
  }
  weights <- if (multi) {
    adaptive_weights$weights
  } else {
    adaptive_weights$weights[context$fine]
  }
  cutoff <- max(
    control$cutoff_floor,
    control$relative_cutoff_multiplier * max(abs(pilot$coefficients))
  )
  grid <- cpp_arfc_grid_v1(
    context$centered_y,
    context$design,
    context$roughness,
    region_gram,
    region_length,
    region_multiplier,
    weights,
    lambda,
    pilot$selected_kappa,
    family,
    shape,
    control$max_iter,
    control$tol,
    control$norm_floor,
    cutoff,
    control$refit,
    control$trace
  )
  selected <- .arfc_select_bic(
    grid, control$strict, .arfc_method_label(resolution, penalty)
  )
  tuning <- data.frame(
    candidate_index = seq_along(grid$lambda),
    lambda = as.numeric(grid$lambda),
    kappa = as.numeric(grid$kappa),
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
    eligible = seq_along(grid$lambda) %in%
      .arfc_eligible_candidates(grid, control$strict),
    selected = seq_along(grid$lambda) == selected,
    message = as.character(grid$message),
    stringsAsFactors = FALSE
  )
  coefficients <- as.numeric(grid$coefficients[, selected])
  names(coefficients) <- paste0("B", seq_along(coefficients))
  beta <- as.numeric(context$basis$basis_observed %*% coefficients)
  fitted_values <- as.numeric(
    context$means$y + context$design %*% coefficients
  )
  selected_region_norms <- as.numeric(grid$region_norms[, selected])
  weight_table <- cbind(
    region_metadata,
    pilot_norm = if (multi) {
      adaptive_weights$pilot_norm
    } else {
      adaptive_weights$pilot_norm[context$fine]
    },
    raw_adaptive_weight = if (multi) {
      adaptive_weights$raw
    } else {
      adaptive_weights$raw[context$fine]
    },
    normalized_adaptive_weight = if (multi) {
      adaptive_weights$normalized
    } else {
      adaptive_weights$normalized[context$fine]
    },
    adaptive_weight = weights,
    selected_region_norm = selected_region_norms
  )
  curvature <- cpp_arfc_curvature_envelope_v1(
    region_gram, region_length, region_multiplier, family, shape
  )
  internal <- if (control$keep_matrices) {
    list(
      centered_y = context$centered_y,
      design = context$design,
      roughness = context$roughness,
      region_gram = region_gram,
      region_multiplier = region_multiplier,
      grid = grid
    )
  } else {
    NULL
  }
  structure(
    list(
      call = call,
      coefficients = coefficients,
      beta = beta,
      argvals = context$argvals,
      fitted.values = fitted_values,
      residuals = context$y - fitted_values,
      interval_norms = selected_region_norms,
      region_norms = selected_region_norms,
      regions = weight_table,
      adaptive_weights = weight_table,
      selected = tuning[selected, , drop = FALSE],
      tuning = tuning,
      smooth_pilot = list(
        coefficients = pilot$coefficients,
        selected_kappa = pilot$selected_kappa,
        tuning = pilot$tuning,
        adaptive_epsilon = adaptive_weights$epsilon
      ),
      diagnostics = list(
        selected_index = selected,
        cutoff = cutoff,
        eligible_candidates = sum(tuning$eligible),
        failed_candidates = sum(!tuning$success),
        nonconverged_candidates = sum(tuning$success & !tuning$converged),
        objective_increase_candidates = sum(
          tuning$objective_increased %in% TRUE
        ),
        objective_trace = if (control$trace) grid$traces[[selected]] else NULL,
        maximum_scalar_concavity = curvature$maximum_scalar_concavity,
        maximum_envelope_eigenvalue =
          curvature$maximum_envelope_eigenvalue,
        canonical_multiplier_total = sum(region_multiplier)
      ),
      method = list(
        code = .arfc_method_code(resolution, penalty),
        label = .arfc_method_label(resolution, penalty),
        resolution = resolution,
        penalty = penalty,
        shape = shape,
        adaptive = adaptive,
        levels = unique(region_metadata$level_M),
        region_count = nrow(region_metadata)
      ),
      means = context$means,
      dimensions = c(
        n = nrow(context$x), p = ncol(context$x),
        basis = length(coefficients), regions = nrow(region_metadata)
      ),
      basis = list(
        M = context$basis$M,
        degree = context$basis$degree,
        penalty_order = context$penalty_order,
        domain = context$basis$domain,
        breaks = context$basis$breaks,
        full_knots = context$basis$full_knots,
        observed = context$basis$basis_observed,
        integration = context$integration,
        integration_weights = context$integration_weights
      ),
      control = control,
      internal = internal
    ),
    class = "arfc"
  )
}
