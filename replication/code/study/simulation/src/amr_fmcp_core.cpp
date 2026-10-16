// Exploratory AMR-fMCP numerical core (pilot v0).
//
// This file is intentionally outside FMCP-package.  It implements the frozen
// feasibility protocol in manuscript/protocol_amr_fmcp_pilot_v0.md without
// changing the released fMCP or fSCAD APIs.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct AmrFit {
  arma::vec coefficients;
  arma::vec region_norms;
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

double mcp_penalty_scalar(const double value,
                          const double lambda,
                          const double gamma) {
  const double magnitude = std::abs(value);
  if (magnitude <= gamma * lambda) {
    return lambda * magnitude - magnitude * magnitude / (2.0 * gamma);
  }
  return 0.5 * gamma * lambda * lambda;
}

double mcp_derivative_scalar(const double value,
                             const double lambda,
                             const double gamma) {
  const double magnitude = std::abs(value);
  return std::max(lambda - magnitude / gamma, 0.0);
}

arma::vec region_norms_internal(const arma::vec& coefficients,
                                const arma::cube& region_gram,
                                const arma::vec& region_length) {
  const arma::uword region_count = region_gram.n_slices;
  arma::vec norms(region_count, arma::fill::zeros);

  for (arma::uword region = 0; region < region_count; ++region) {
    double squared_integral = arma::dot(
      coefficients,
      region_gram.slice(region) * coefficients
    );
    if (squared_integral < 0.0 && squared_integral > -1e-12) {
      squared_integral = 0.0;
    }
    if (squared_integral < 0.0) {
      Rcpp::stop("A regional Gram matrix produced a negative squared norm.");
    }
    norms[region] = std::sqrt(squared_integral / region_length[region]);
  }
  return norms;
}

double exact_objective(const arma::vec& y,
                       const arma::mat& design,
                       const arma::mat& roughness,
                       const arma::cube& region_gram,
                       const arma::vec& region_length,
                       const arma::vec& level_weight,
                       const arma::vec& adaptive_weight,
                       const arma::vec& coefficients,
                       const double kappa,
                       const double lambda,
                       const double gamma) {
  const arma::vec residual = y - design * coefficients;
  const double loss = arma::dot(residual, residual) /
    static_cast<double>(y.n_elem);
  const double smoothness = kappa *
    arma::dot(coefficients, roughness * coefficients);
  const arma::vec norms = region_norms_internal(
    coefficients, region_gram, region_length
  );
  double penalty = 0.0;
  for (arma::uword region = 0; region < norms.n_elem; ++region) {
    const double regional_lambda = lambda * adaptive_weight[region];
    penalty += level_weight[region] * mcp_penalty_scalar(
      norms[region], regional_lambda, gamma
    );
  }
  return loss + smoothness + penalty;
}

arma::mat lqa_penalty_matrix(const arma::vec& coefficients,
                             const arma::cube& region_gram,
                             const arma::vec& region_length,
                             const arma::vec& level_weight,
                             const arma::vec& adaptive_weight,
                             const double lambda,
                             const double gamma,
                             const double norm_floor) {
  const arma::uword basis_count = coefficients.n_elem;
  const arma::uword region_count = region_gram.n_slices;
  arma::mat penalty(basis_count, basis_count, arma::fill::zeros);

  for (arma::uword region = 0; region < region_count; ++region) {
    double squared_integral = arma::dot(
      coefficients,
      region_gram.slice(region) * coefficients
    );
    if (squared_integral < 0.0 && squared_integral > -1e-12) {
      squared_integral = 0.0;
    }
    if (squared_integral < 0.0) {
      Rcpp::stop("A regional Gram matrix produced a negative squared norm.");
    }

    const double region_norm = std::sqrt(
      squared_integral / region_length[region]
    );
    const double regional_lambda = lambda * adaptive_weight[region];
    const double derivative = level_weight[region] *
      mcp_derivative_scalar(region_norm, regional_lambda, gamma);
    if (derivative == 0.0) {
      continue;
    }

    const double denominator = std::max(region_norm, norm_floor);
    const double multiplier = 0.5 * derivative /
      (denominator * region_length[region]);
    penalty += multiplier * region_gram.slice(region);
  }

  return arma::symmatu(penalty);
}

AmrFit fit_single(const arma::vec& y,
                  const arma::mat& design,
                  const arma::mat& roughness,
                  const arma::cube& region_gram,
                  const arma::vec& region_length,
                  const arma::vec& level_weight,
                  const arma::vec& adaptive_weight,
                  const arma::mat& design_crossprod,
                  const arma::vec& design_response,
                  const arma::vec& initial,
                  const double kappa,
                  const double lambda,
                  const double gamma,
                  const int max_iter,
                  const double tolerance,
                  const double norm_floor,
                  const double cutoff,
                  const bool refit,
                  const bool keep_trace) {
  AmrFit result;
  const double sample_size = static_cast<double>(y.n_elem);
  arma::vec coefficients = initial;

  if (lambda == 0.0) {
    result.converged = true;
    if (keep_trace) {
      result.objective_trace.push_back(exact_objective(
        y, design, roughness, region_gram, region_length,
        level_weight, adaptive_weight, coefficients,
        kappa, lambda, gamma
      ));
    }
  } else {
    double previous_objective = exact_objective(
      y, design, roughness, region_gram, region_length,
      level_weight, adaptive_weight, coefficients,
      kappa, lambda, gamma
    );
    if (keep_trace) {
      result.objective_trace.push_back(previous_objective);
    }

    for (int iteration = 1; iteration <= max_iter; ++iteration) {
      const arma::mat lqa = lqa_penalty_matrix(
        coefficients, region_gram, region_length,
        level_weight, adaptive_weight, lambda, gamma, norm_floor
      );
      const arma::mat system = design_crossprod +
        sample_size * kappa * roughness + sample_size * lqa;
      arma::vec updated;
      if (!solve_vector(system, design_response, updated)) {
        result.message = "The AMR-fMCP LQA system could not be solved.";
        result.coefficients = coefficients;
        return result;
      }

      const double current_objective = exact_objective(
        y, design, roughness, region_gram, region_length,
        level_weight, adaptive_weight, updated,
        kappa, lambda, gamma
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

  if (lambda > 0.0) {
    const arma::uvec below_cutoff = arma::find(
      arma::abs(coefficients) < cutoff
    );
    coefficients.elem(below_cutoff).zeros();
  }

  arma::uvec active = arma::find(arma::abs(coefficients) > 0.0);
  if (lambda > 0.0 && refit) {
    if (active.n_elem == 0) {
      coefficients.zeros();
    } else {
      const arma::mat active_crossprod = design_crossprod.submat(
        active, active
      );
      const arma::mat active_roughness = roughness.submat(active, active);
      const arma::mat active_system = active_crossprod +
        sample_size * kappa * active_roughness;
      arma::vec refitted;
      if (!solve_vector(
          active_system, design_response.elem(active), refitted)) {
        result.refit_success = false;
        result.message =
          "The AMR-fMCP active-set refit system could not be solved.";
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
      result.message = "AMR-fMCP degrees of freedom could not be computed.";
      result.coefficients = coefficients;
      return result;
    }
    degrees_freedom = arma::trace(influence_core);
  }

  const double mean_rss = std::max(
    rss / sample_size, std::numeric_limits<double>::min()
  );
  result.coefficients = coefficients;
  result.region_norms = region_norms_internal(
    coefficients, region_gram, region_length
  );
  result.active = static_cast<int>(active.n_elem);
  result.rss = rss;
  result.df = degrees_freedom;
  result.bic = sample_size * std::log(mean_rss) +
    std::log(sample_size) * degrees_freedom;
  result.objective = exact_objective(
    y, design, roughness, region_gram, region_length,
    level_weight, adaptive_weight, coefficients,
    kappa, lambda, gamma
  );
  result.success = true;
  if (!result.converged) {
    result.message = "Maximum AMR-fMCP LQA iterations reached.";
  }
  return result;
}

void validate_common_inputs(const arma::vec& y,
                            const arma::mat& design,
                            const arma::mat& roughness,
                            const arma::cube& region_gram,
                            const arma::vec& region_length,
                            const arma::vec& level_weight,
                            const arma::vec& adaptive_weight) {
  const arma::uword sample_size = y.n_elem;
  const arma::uword basis_count = design.n_cols;
  const arma::uword region_count = region_gram.n_slices;

  if (sample_size < 2 || design.n_rows != sample_size || basis_count < 1 ||
      roughness.n_rows != basis_count || roughness.n_cols != basis_count ||
      region_gram.n_rows != basis_count ||
      region_gram.n_cols != basis_count || region_count < 1 ||
      region_length.n_elem != region_count ||
      level_weight.n_elem != region_count ||
      adaptive_weight.n_elem != region_count) {
    Rcpp::stop("AMR-fMCP core matrices have incompatible dimensions.");
  }
  if (!y.is_finite() || !design.is_finite() || !roughness.is_finite() ||
      !region_gram.is_finite() || !region_length.is_finite() ||
      !level_weight.is_finite() || !adaptive_weight.is_finite()) {
    Rcpp::stop("AMR-fMCP core inputs must be finite.");
  }
  if (arma::any(region_length <= 0.0) ||
      arma::any(level_weight <= 0.0) ||
      arma::any(adaptive_weight <= 0.0)) {
    Rcpp::stop("AMR-fMCP regional lengths and weights must be positive.");
  }
}

}  // namespace


// [[Rcpp::export]]
arma::vec cpp_amr_region_norms_v0(const arma::vec& coefficients,
                                  const arma::cube& region_gram,
                                  const arma::vec& region_length) {
  if (region_gram.n_rows != coefficients.n_elem ||
      region_gram.n_cols != coefficients.n_elem ||
      region_gram.n_slices != region_length.n_elem ||
      !coefficients.is_finite() || !region_gram.is_finite() ||
      !region_length.is_finite() || arma::any(region_length <= 0.0)) {
    Rcpp::stop("Invalid AMR regional-norm inputs.");
  }
  return region_norms_internal(coefficients, region_gram, region_length);
}


// [[Rcpp::export]]
Rcpp::List cpp_amr_fmcp_grid_v0(const arma::vec& y,
                                const arma::mat& design,
                                const arma::mat& roughness,
                                const arma::cube& region_gram,
                                const arma::vec& region_length,
                                const arma::vec& level_weight,
                                const arma::vec& adaptive_weight,
                                const arma::vec& lambda,
                                const arma::vec& kappa,
                                const double gamma,
                                const int max_iter,
                                const double tolerance,
                                const double norm_floor,
                                const double cutoff,
                                const bool refit,
                                const bool keep_trace) {
  validate_common_inputs(
    y, design, roughness, region_gram, region_length,
    level_weight, adaptive_weight
  );
  if (lambda.n_elem == 0 || kappa.n_elem == 0 ||
      arma::any(lambda < 0.0) || arma::any(kappa < 0.0) ||
      !lambda.is_finite() || !kappa.is_finite() ||
      !std::isfinite(gamma) || gamma <= 1.0) {
    Rcpp::stop("Invalid AMR-fMCP tuning inputs.");
  }
  if (max_iter < 1 || tolerance <= 0.0 || norm_floor <= 0.0 ||
      cutoff < 0.0 || !std::isfinite(tolerance) ||
      !std::isfinite(norm_floor) || !std::isfinite(cutoff)) {
    Rcpp::stop("Invalid AMR-fMCP numerical controls.");
  }
  if (refit && arma::any(lambda > 0.0) && cutoff <= 0.0) {
    Rcpp::stop(
      "AMR-fMCP forbids `refit = TRUE` with a zero coefficient cutoff."
    );
  }

  const arma::uword sample_size = y.n_elem;
  const arma::uword basis_count = design.n_cols;
  const arma::uword region_count = region_gram.n_slices;
  const arma::mat design_crossprod = arma::symmatu(design.t() * design);
  const arma::vec design_response = design.t() * y;
  const int candidate_count = static_cast<int>(lambda.n_elem * kappa.n_elem);

  arma::mat coefficients(basis_count, candidate_count);
  coefficients.fill(NA_REAL);
  arma::mat region_norms(region_count, candidate_count);
  region_norms.fill(NA_REAL);
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
          "The AMR-fMCP roughness initialization could not be solved.";
        traces[candidate] = R_NilValue;
        ++candidate;
        continue;
      }

      const AmrFit fit = fit_single(
        y, design, roughness, region_gram, region_length,
        level_weight, adaptive_weight, design_crossprod, design_response,
        initial, kappa[k], lambda[l], gamma, max_iter, tolerance,
        norm_floor, cutoff, refit, keep_trace
      );

      coefficients.col(candidate) = fit.coefficients;
      if (fit.region_norms.n_elem == region_count) {
        region_norms.col(candidate) = fit.region_norms;
      }
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
    Rcpp::Named("region_norms") = region_norms,
    Rcpp::Named("lambda") = lambda_out,
    Rcpp::Named("kappa") = kappa_out,
    Rcpp::Named("mcp_gamma") = gamma,
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
