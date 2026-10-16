// fSCAD numerical core for locally sparse functional linear regression.
//
// This implementation mirrors the package's fMCP numerical scheme while
// replacing only the MCP penalty and derivative with SCAD. The separate
// translation unit deliberately leaves the established fMCP core unchanged.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct SingleFit {
  arma::vec coefficients;
  bool success = false;
  bool converged = false;
  bool refit_success = true;
  bool objective_increased = false;
  int iterations = 0;
  int active = 0;
  double objective = NA_REAL;
  double rss = NA_REAL;
  double df = NA_REAL;
  double bic = R_PosInf;
  std::string message;
  std::vector<double> objective_trace;
};

bool solve_vector(const arma::mat& system,
                  const arma::vec& rhs,
                  arma::vec& solution) {
  const arma::mat symmetric_system = arma::symmatu(system);
  bool solved = arma::solve(
    solution,
    symmetric_system,
    rhs,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (!solved) {
    solved = arma::solve(
      solution,
      symmetric_system,
      rhs,
      arma::solve_opts::no_approx
    );
  }
  return solved && solution.is_finite();
}

bool solve_matrix(const arma::mat& system,
                  const arma::mat& rhs,
                  arma::mat& solution) {
  const arma::mat symmetric_system = arma::symmatu(system);
  bool solved = arma::solve(
    solution,
    symmetric_system,
    rhs,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (!solved) {
    solved = arma::solve(
      solution,
      symmetric_system,
      rhs,
      arma::solve_opts::no_approx
    );
  }
  return solved && solution.is_finite();
}

double scad_penalty_scalar(const double value,
                           const double lambda,
                           const double shape) {
  const double magnitude = std::abs(value);
  if (magnitude <= lambda) {
    return lambda * magnitude;
  }
  if (magnitude < shape * lambda) {
    return (-magnitude * magnitude + 2.0 * shape * lambda * magnitude -
            lambda * lambda) / (2.0 * (shape - 1.0));
  }
  return 0.5 * (shape + 1.0) * lambda * lambda;
}

double scad_derivative_scalar(const double value,
                              const double lambda,
                              const double shape) {
  const double magnitude = std::abs(value);
  if (magnitude <= lambda) {
    return lambda;
  }
  if (magnitude < shape * lambda) {
    return (shape * lambda - magnitude) / (shape - 1.0);
  }
  return 0.0;
}

arma::vec interval_norms_internal(const arma::vec& coefficients,
                                  const arma::cube& local_gram,
                                  const double scale) {
  const arma::uword local_order = local_gram.n_rows;
  const arma::uword intervals = local_gram.n_slices;
  arma::vec norms(intervals, arma::fill::zeros);

  for (arma::uword j = 0; j < intervals; ++j) {
    const arma::vec local = coefficients.subvec(j, j + local_order - 1);
    double squared_norm = arma::dot(local, local_gram.slice(j) * local);
    if (squared_norm < 0.0 && squared_norm > -1e-12) {
      squared_norm = 0.0;
    }
    if (squared_norm < 0.0) {
      Rcpp::stop("An interval Gram matrix produced a negative squared norm.");
    }
    norms[j] = scale * std::sqrt(squared_norm);
  }
  return norms;
}

double exact_objective(const arma::vec& y,
                       const arma::mat& design,
                       const arma::mat& roughness,
                       const arma::cube& local_gram,
                       const arma::vec& coefficients,
                       const double kappa,
                       const double lambda,
                       const double shape,
                       const double scale) {
  const arma::vec residual = y - design * coefficients;
  const double loss = arma::dot(residual, residual) /
    static_cast<double>(y.n_elem);
  const double smoothness = kappa *
    arma::dot(coefficients, roughness * coefficients);
  const arma::vec norms = interval_norms_internal(
    coefficients, local_gram, scale
  );
  double penalty = 0.0;
  for (arma::uword j = 0; j < norms.n_elem; ++j) {
    penalty += scad_penalty_scalar(norms[j], lambda, shape);
  }
  return loss + smoothness + penalty;
}

arma::mat lqa_penalty_matrix(const arma::vec& coefficients,
                             const arma::cube& local_gram,
                             const double lambda,
                             const double shape,
                             const double scale,
                             const double norm_floor) {
  const arma::uword local_order = local_gram.n_rows;
  const arma::uword intervals = local_gram.n_slices;
  const arma::uword basis_count = intervals + local_order - 1;
  arma::mat penalty(basis_count, basis_count, arma::fill::zeros);
  const double scale_squared = scale * scale;

  for (arma::uword j = 0; j < intervals; ++j) {
    const arma::vec local = coefficients.subvec(j, j + local_order - 1);
    double squared_norm = arma::dot(local, local_gram.slice(j) * local);
    if (squared_norm < 0.0 && squared_norm > -1e-12) {
      squared_norm = 0.0;
    }
    if (squared_norm < 0.0) {
      Rcpp::stop("An interval Gram matrix produced a negative squared norm.");
    }

    const double interval_norm = scale * std::sqrt(squared_norm);
    const double derivative = scad_derivative_scalar(
      interval_norm, lambda, shape
    );
    if (derivative == 0.0) {
      continue;
    }

    const double denominator = std::max(interval_norm, norm_floor);
    const double multiplier = 0.5 * derivative * scale_squared / denominator;
    penalty.submat(j, j, j + local_order - 1, j + local_order - 1) +=
      multiplier * local_gram.slice(j);
  }

  return arma::symmatu(penalty);
}

SingleFit fit_single(const arma::vec& y,
                     const arma::mat& design,
                     const arma::mat& roughness,
                     const arma::cube& local_gram,
                     const arma::mat& design_crossprod,
                     const arma::vec& design_response,
                     const arma::vec& initial,
                     const double kappa,
                     const double lambda,
                     const double shape,
                     const double scale,
                     const int max_iter,
                     const double tolerance,
                     const double norm_floor,
                     const double cutoff,
                     const bool refit,
                     const bool keep_trace) {
  SingleFit result;
  const arma::uword basis_count = design.n_cols;
  const double sample_size = static_cast<double>(y.n_elem);
  arma::vec coefficients = initial;

  if (lambda == 0.0) {
    result.converged = true;
    if (keep_trace) {
      result.objective_trace.push_back(exact_objective(
        y, design, roughness, local_gram, coefficients,
        kappa, lambda, shape, scale
      ));
    }
  } else {
    double previous_objective = exact_objective(
      y, design, roughness, local_gram, coefficients,
      kappa, lambda, shape, scale
    );
    if (keep_trace) {
      result.objective_trace.push_back(previous_objective);
    }

    for (int iteration = 1; iteration <= max_iter; ++iteration) {
      const arma::mat lqa = lqa_penalty_matrix(
        coefficients, local_gram, lambda, shape, scale, norm_floor
      );
      const arma::mat system = design_crossprod +
        sample_size * kappa * roughness + sample_size * lqa;
      arma::vec updated;
      if (!solve_vector(system, design_response, updated)) {
        result.message = "The fSCAD LQA system could not be solved.";
        result.coefficients = coefficients;
        return result;
      }

      const double current_objective = exact_objective(
        y, design, roughness, local_gram, updated,
        kappa, lambda, shape, scale
      );
      if (keep_trace) {
        result.objective_trace.push_back(current_objective);
      }
      if (current_objective > previous_objective +
          tolerance * (1.0 + std::abs(previous_objective))) {
        result.objective_increased = true;
      }

      const double change = arma::abs(updated - coefficients).max();
      const double reference = 1.0 + arma::abs(coefficients).max();
      coefficients = updated;
      result.iterations = iteration;
      previous_objective = current_objective;
      if (change <= tolerance * reference) {
        result.converged = true;
        break;
      }
    }
  }

  if (lambda > 0.0 && cutoff > 0.0) {
    coefficients.elem(arma::find(arma::abs(coefficients) < cutoff)).zeros();
  }

  arma::uvec active;
  if (lambda > 0.0 && cutoff == 0.0) {
    active = arma::regspace<arma::uvec>(0, basis_count - 1);
  } else {
    active = arma::find(arma::abs(coefficients) > 0.0);
  }

  if (lambda > 0.0 && refit) {
    if (active.n_elem == 0) {
      coefficients.zeros();
    } else {
      const arma::mat active_crossprod = design_crossprod.submat(active, active);
      const arma::mat active_roughness = roughness.submat(active, active);
      const arma::mat active_system = active_crossprod +
        sample_size * kappa * active_roughness;
      arma::vec refitted;
      if (!solve_vector(active_system, design_response.elem(active), refitted)) {
        result.refit_success = false;
        result.message = "The fSCAD active-set refit could not be solved.";
        result.coefficients = coefficients;
        return result;
      }
      coefficients.zeros();
      coefficients.elem(active) = refitted;
    }
  }

  active = arma::find(arma::abs(coefficients) > 0.0);
  const arma::vec residual = y - design * coefficients;
  const double rss = arma::dot(residual, residual);
  double degrees_freedom = 0.0;

  if (active.n_elem > 0) {
    const arma::mat active_crossprod = design_crossprod.submat(active, active);
    const arma::mat active_roughness = roughness.submat(active, active);
    const arma::mat active_system = active_crossprod +
      sample_size * kappa * active_roughness;
    arma::mat influence_core;
    if (!solve_matrix(active_system, active_crossprod, influence_core)) {
      result.message = "fSCAD degrees of freedom could not be computed.";
      result.coefficients = coefficients;
      return result;
    }
    degrees_freedom = arma::trace(influence_core);
  }

  const double mean_rss = std::max(
    rss / sample_size, std::numeric_limits<double>::min()
  );
  result.coefficients = coefficients;
  result.active = static_cast<int>(active.n_elem);
  result.rss = rss;
  result.df = degrees_freedom;
  result.bic = sample_size * std::log(mean_rss) +
    std::log(sample_size) * degrees_freedom;
  result.objective = exact_objective(
    y, design, roughness, local_gram, coefficients,
    kappa, lambda, shape, scale
  );
  result.success = true;
  if (!result.converged) {
    result.message = "Maximum fSCAD LQA iterations reached.";
  }
  return result;
}

}  // namespace


// [[Rcpp::export]]
arma::vec cpp_scad_penalty(const arma::vec& values,
                           const double lambda,
                           const double shape = 3.7) {
  if (lambda < 0.0 || shape <= 2.0 || !std::isfinite(lambda) ||
      !std::isfinite(shape) || !values.is_finite()) {
    Rcpp::stop("Invalid SCAD penalty inputs.");
  }
  arma::vec result(values.n_elem);
  for (arma::uword i = 0; i < values.n_elem; ++i) {
    result[i] = scad_penalty_scalar(values[i], lambda, shape);
  }
  return result;
}


// [[Rcpp::export]]
arma::vec cpp_scad_derivative(const arma::vec& values,
                              const double lambda,
                              const double shape = 3.7) {
  if (lambda < 0.0 || shape <= 2.0 || !std::isfinite(lambda) ||
      !std::isfinite(shape) || !values.is_finite()) {
    Rcpp::stop("Invalid SCAD derivative inputs.");
  }
  arma::vec result(values.n_elem);
  for (arma::uword i = 0; i < values.n_elem; ++i) {
    result[i] = scad_derivative_scalar(values[i], lambda, shape);
  }
  return result;
}


// [[Rcpp::export]]
Rcpp::List cpp_fscad_grid(const arma::vec& y,
                          const arma::mat& design,
                          const arma::mat& roughness,
                          const arma::cube& local_gram,
                          const arma::vec& lambda,
                          const arma::vec& kappa,
                          const double shape,
                          const double scale,
                          const int max_iter,
                          const double tolerance,
                          const double norm_floor,
                          const double cutoff,
                          const bool refit,
                          const bool keep_trace) {
  const arma::uword sample_size = y.n_elem;
  const arma::uword basis_count = design.n_cols;
  const arma::uword expected_basis = local_gram.n_slices +
    local_gram.n_rows - 1;

  if (sample_size < 2 || design.n_rows != sample_size ||
      roughness.n_rows != basis_count || roughness.n_cols != basis_count ||
      local_gram.n_rows != local_gram.n_cols || expected_basis != basis_count) {
    Rcpp::stop("fSCAD core matrices have incompatible dimensions.");
  }
  if (!y.is_finite() || !design.is_finite() || !roughness.is_finite() ||
      !local_gram.is_finite()) {
    Rcpp::stop("fSCAD core matrices must be finite.");
  }
  if (lambda.n_elem == 0 || kappa.n_elem == 0 ||
      arma::any(lambda < 0.0) || arma::any(kappa < 0.0) ||
      !lambda.is_finite() || !kappa.is_finite() || shape <= 2.0 ||
      !std::isfinite(shape)) {
    Rcpp::stop("Invalid fSCAD tuning grid.");
  }
  if (scale <= 0.0 || tolerance <= 0.0 || norm_floor <= 0.0 ||
      cutoff < 0.0 || max_iter < 1) {
    Rcpp::stop("Invalid fSCAD numerical controls.");
  }

  const arma::mat design_crossprod = arma::symmatu(design.t() * design);
  const arma::vec design_response = design.t() * y;
  const int candidate_count = static_cast<int>(lambda.n_elem * kappa.n_elem);

  arma::mat coefficients(basis_count, candidate_count);
  coefficients.fill(NA_REAL);
  Rcpp::NumericVector lambda_out(candidate_count);
  Rcpp::NumericVector kappa_out(candidate_count);
  Rcpp::NumericVector bic(candidate_count, R_PosInf);
  Rcpp::NumericVector rss(candidate_count, NA_REAL);
  Rcpp::NumericVector df(candidate_count, NA_REAL);
  Rcpp::NumericVector objective(candidate_count, NA_REAL);
  Rcpp::IntegerVector iterations(candidate_count);
  Rcpp::IntegerVector active(candidate_count);
  Rcpp::LogicalVector success(candidate_count, false);
  Rcpp::LogicalVector converged(candidate_count, false);
  Rcpp::LogicalVector refit_success(candidate_count, false);
  Rcpp::LogicalVector objective_increased(candidate_count, false);
  Rcpp::CharacterVector message(candidate_count);
  Rcpp::List traces(candidate_count);

  int candidate = 0;
  for (arma::uword k = 0; k < kappa.n_elem; ++k) {
    const arma::mat initial_system = design_crossprod +
      static_cast<double>(sample_size) * kappa[k] * roughness;
    arma::vec initial;
    const bool initial_success = solve_vector(
      initial_system, design_response, initial
    );

    for (arma::uword l = 0; l < lambda.n_elem; ++l) {
      Rcpp::checkUserInterrupt();
      lambda_out[candidate] = lambda[l];
      kappa_out[candidate] = kappa[k];
      if (!initial_success) {
        message[candidate] =
          "The roughness-only initialization could not be solved.";
        traces[candidate] = R_NilValue;
        ++candidate;
        continue;
      }

      const SingleFit fit = fit_single(
        y, design, roughness, local_gram,
        design_crossprod, design_response, initial,
        kappa[k], lambda[l], shape, scale,
        max_iter, tolerance, norm_floor, cutoff, refit, keep_trace
      );

      coefficients.col(candidate) = fit.coefficients;
      bic[candidate] = fit.bic;
      rss[candidate] = fit.rss;
      df[candidate] = fit.df;
      objective[candidate] = fit.objective;
      iterations[candidate] = fit.iterations;
      active[candidate] = fit.active;
      success[candidate] = fit.success;
      converged[candidate] = fit.converged;
      refit_success[candidate] = fit.refit_success;
      objective_increased[candidate] = fit.objective_increased;
      message[candidate] = fit.message;
      traces[candidate] = keep_trace ?
        Rcpp::wrap(fit.objective_trace) : R_NilValue;
      ++candidate;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("lambda") = lambda_out,
    Rcpp::Named("kappa") = kappa_out,
    Rcpp::Named("scad_a") = shape,
    Rcpp::Named("bic") = bic,
    Rcpp::Named("rss") = rss,
    Rcpp::Named("df") = df,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("active") = active,
    Rcpp::Named("success") = success,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("refit_success") = refit_success,
    Rcpp::Named("objective_increased") = objective_increased,
    Rcpp::Named("message") = message,
    Rcpp::Named("traces") = traces
  );
}

