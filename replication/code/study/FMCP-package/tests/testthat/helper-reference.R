make_fmcp_test_data <- function(n = 50L, p = 21L, seed = 2026L) {
  set.seed(seed)
  argvals <- seq(0, 1, length.out = p)
  component_count <- 8L
  scores <- matrix(rnorm(n * component_count), n, component_count)
  x <- matrix(0, n, p)
  for (j in seq_len(component_count)) {
    x <- x + scores[, j] %o% sin(j * pi * argvals)
  }
  x <- x + matrix(rnorm(n * p, sd = 0.05), n, p)
  beta <- ifelse(argvals <= 0.6, 1.5 * sin(2 * pi * argvals), 0)
  integration_weights <- FMCP:::.fmcp_integration_weights(argvals, "simpson")
  signal <- as.numeric(x %*% (integration_weights * beta))
  y <- 0.75 + signal + rnorm(n, sd = 0.08)
  list(y = y, x = x, argvals = argvals, beta = beta)
}

reference_fmcp_lqa <- function(y,
                               design,
                               roughness,
                               local_gram,
                               lambda,
                               kappa,
                               mcp_gamma,
                               scale,
                               max_iter,
                               tolerance,
                               norm_floor) {
  n <- length(y)
  crossprod_design <- crossprod(design)
  design_response <- crossprod(design, y)
  coefficients <- as.numeric(solve(
    crossprod_design + n * kappa * roughness,
    design_response
  ))
  if (lambda == 0) {
    return(coefficients)
  }

  interval_count <- dim(local_gram)[3L]
  local_order <- dim(local_gram)[1L]
  for (iteration in seq_len(max_iter)) {
    lqa <- matrix(0, ncol(design), ncol(design))
    for (j in seq_len(interval_count)) {
      index <- j:(j + local_order - 1L)
      local <- coefficients[index]
      squared_norm <- drop(
        crossprod(local, local_gram[, , j] %*% local)
      )
      interval_norm <- scale * sqrt(max(squared_norm, 0))
      derivative <- max(lambda - interval_norm / mcp_gamma, 0)
      if (derivative > 0) {
        multiplier <- 0.5 * derivative * scale^2 /
          max(interval_norm, norm_floor)
        lqa[index, index] <- lqa[index, index] +
          multiplier * local_gram[, , j]
      }
    }
    updated <- as.numeric(solve(
      crossprod_design + n * kappa * roughness + n * lqa,
      design_response
    ))
    change <- max(abs(updated - coefficients))
    reference <- 1 + max(abs(coefficients))
    coefficients <- updated
    if (change <= tolerance * reference) {
      break
    }
  }
  coefficients
}

reference_fscad_lqa <- function(y,
                                design,
                                roughness,
                                local_gram,
                                lambda,
                                kappa,
                                scad_a,
                                scale,
                                max_iter,
                                tolerance,
                                norm_floor) {
  n <- length(y)
  crossprod_design <- crossprod(design)
  design_response <- crossprod(design, y)
  coefficients <- as.numeric(solve(
    crossprod_design + n * kappa * roughness,
    design_response
  ))
  if (lambda == 0) {
    return(coefficients)
  }

  interval_count <- dim(local_gram)[3L]
  local_order <- dim(local_gram)[1L]
  for (iteration in seq_len(max_iter)) {
    lqa <- matrix(0, ncol(design), ncol(design))
    for (j in seq_len(interval_count)) {
      index <- j:(j + local_order - 1L)
      local <- coefficients[index]
      squared_norm <- drop(
        crossprod(local, local_gram[, , j] %*% local)
      )
      interval_norm <- scale * sqrt(max(squared_norm, 0))
      derivative <- if (interval_norm <= lambda) {
        lambda
      } else if (interval_norm < scad_a * lambda) {
        (scad_a * lambda - interval_norm) / (scad_a - 1)
      } else {
        0
      }
      if (derivative > 0) {
        multiplier <- 0.5 * derivative * scale^2 /
          max(interval_norm, norm_floor)
        lqa[index, index] <- lqa[index, index] +
          multiplier * local_gram[, , j]
      }
    }
    updated <- as.numeric(solve(
      crossprod_design + n * kappa * roughness + n * lqa,
      design_response
    ))
    change <- max(abs(updated - coefficients))
    reference <- 1 + max(abs(coefficients))
    coefficients <- updated
    if (change <= tolerance * reference) {
      break
    }
  }
  coefficients
}
