// Adaptive regional folded-concave numerical core (ARFC-v1).
//
// This prototype is intentionally separate from FMCP-package.  It implements
// manuscript/protocol_adaptive_regional_folded_concave_v1.md without changing
// any legacy fMCP/fSCAD compiled symbol or high-level API.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

enum class PenaltyFamily {
  MCP,
  SCAD
};

struct ArfcFit {
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

std::string lower_string(std::string value) {
  std::transform(
    value.begin(), value.end(), value.begin(),
    [](unsigned char character) {
      return static_cast<char>(std::tolower(character));
    }
  );
  return value;
}

PenaltyFamily parse_penalty_family(const std::string& family) {
  const std::string normalized = lower_string(family);
  if (normalized == "mcp") {
    return PenaltyFamily::MCP;
  }
  if (normalized == "scad") {
    return PenaltyFamily::SCAD;
  }
  Rcpp::stop("`penalty_family` must be either 'mcp' or 'scad'.");
  return PenaltyFamily::MCP;
}

std::string penalty_family_name(const PenaltyFamily family) {
  return family == PenaltyFamily::MCP ? "mcp" : "scad";
}

void validate_penalty_shape(const PenaltyFamily family,
                            const double shape) {
  if (!std::isfinite(shape)) {
    Rcpp::stop("The folded-concave shape parameter must be finite.");
  }
  if (family == PenaltyFamily::MCP && shape <= 1.0) {
    Rcpp::stop("MCP requires `shape > 1`.");
  }
  if (family == PenaltyFamily::SCAD && shape <= 2.0) {
    Rcpp::stop("SCAD requires `shape > 2`.");
  }
}

double scalar_maximum_concavity(const PenaltyFamily family,
                                const double shape) {
  return family == PenaltyFamily::MCP ?
    1.0 / shape : 1.0 / (shape - 1.0);
}

double penalty_value_scalar(const double value,
                            const double lambda,
                            const PenaltyFamily family,
                            const double shape) {
  const double magnitude = std::abs(value);
  if (lambda <= 0.0) {
    return 0.0;
  }

  if (family == PenaltyFamily::MCP) {
    if (magnitude <= shape * lambda) {
      return lambda * magnitude -
        magnitude * magnitude / (2.0 * shape);
    }
    return 0.5 * shape * lambda * lambda;
  }

  if (magnitude <= lambda) {
    return lambda * magnitude;
  }
  if (magnitude <= shape * lambda) {
    return (
      -magnitude * magnitude +
      2.0 * shape * lambda * magnitude -
      lambda * lambda
    ) / (2.0 * (shape - 1.0));
  }
  return 0.5 * (shape + 1.0) * lambda * lambda;
}

double penalty_derivative_scalar(const double value,
                                 const double lambda,
                                 const PenaltyFamily family,
                                 const double shape) {
  const double magnitude = std::abs(value);
  if (lambda <= 0.0) {
    return 0.0;
  }

  if (family == PenaltyFamily::MCP) {
    return std::max(lambda - magnitude / shape, 0.0);
  }

  if (magnitude <= lambda) {
    return lambda;
  }
  if (magnitude <= shape * lambda) {
    return (shape * lambda - magnitude) / (shape - 1.0);
  }
  return 0.0;
}

double penalty_second_derivative_scalar(const double value,
                                        const double lambda,
                                        const PenaltyFamily family,
                                        const double shape,
                                        bool& at_knot) {
  const double magnitude = std::abs(value);
  at_knot = false;
  if (lambda <= 0.0) {
    return 0.0;
  }

  const double knot_tolerance = 1e-12 * (1.0 + magnitude + lambda);
  if (family == PenaltyFamily::MCP) {
    const double upper = shape * lambda;
    if (std::abs(magnitude - upper) <= knot_tolerance) {
      at_knot = true;
      return NA_REAL;
    }
    return magnitude < upper ? -1.0 / shape : 0.0;
  }

  const double lower = lambda;
  const double upper = shape * lambda;
  if (std::abs(magnitude - lower) <= knot_tolerance ||
      std::abs(magnitude - upper) <= knot_tolerance) {
    at_knot = true;
    return NA_REAL;
  }
  if (magnitude > lower && magnitude < upper) {
    return -1.0 / (shape - 1.0);
  }
  return 0.0;
}

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

void validate_region_inputs(const arma::vec& coefficients,
                            const arma::cube& region_gram,
                            const arma::vec& region_length) {
  if (coefficients.n_elem < 1 ||
      region_gram.n_rows != coefficients.n_elem ||
      region_gram.n_cols != coefficients.n_elem ||
      region_gram.n_slices < 1 ||
      region_gram.n_slices != region_length.n_elem ||
      !coefficients.is_finite() || !region_gram.is_finite() ||
      !region_length.is_finite() || arma::any(region_length <= 0.0)) {
    Rcpp::stop("Invalid ARFC regional-norm inputs.");
  }
}

void validate_common_inputs(const arma::vec& y,
                            const arma::mat& design,
                            const arma::mat& roughness,
                            const arma::cube& region_gram,
                            const arma::vec& region_length,
                            const arma::vec& region_multiplier,
                            const arma::vec& adaptive_weight) {
  const arma::uword sample_size = y.n_elem;
  const arma::uword basis_count = design.n_cols;
  const arma::uword region_count = region_gram.n_slices;

  if (sample_size < 2 || design.n_rows != sample_size || basis_count < 1 ||
      roughness.n_rows != basis_count || roughness.n_cols != basis_count ||
      region_gram.n_rows != basis_count ||
      region_gram.n_cols != basis_count || region_count < 1 ||
      region_length.n_elem != region_count ||
      region_multiplier.n_elem != region_count ||
      adaptive_weight.n_elem != region_count) {
    Rcpp::stop("ARFC core matrices have incompatible dimensions.");
  }
  if (!y.is_finite() || !design.is_finite() || !roughness.is_finite() ||
      !region_gram.is_finite() || !region_length.is_finite() ||
      !region_multiplier.is_finite() || !adaptive_weight.is_finite()) {
    Rcpp::stop("ARFC core inputs must be finite.");
  }
  if (arma::any(region_length <= 0.0) ||
      arma::any(region_multiplier < 0.0) ||
      arma::all(region_multiplier == 0.0) ||
      arma::any(adaptive_weight <= 0.0)) {
    Rcpp::stop(
      "ARFC lengths/weights are invalid: lengths and adaptive weights must "
      "be positive and regional multipliers must be nonnegative with a "
      "positive total."
    );
  }
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

arma::mat normalized_region_matrix(const arma::cube& region_gram,
                                   const arma::vec& region_length,
                                   const arma::vec& region_multiplier) {
  const arma::uword basis_count = region_gram.n_rows;
  arma::mat result(basis_count, basis_count, arma::fill::zeros);
  for (arma::uword region = 0; region < region_gram.n_slices; ++region) {
    if (region_multiplier[region] == 0.0) {
      continue;
    }
    result += (region_multiplier[region] / region_length[region]) *
      region_gram.slice(region);
  }
  return arma::symmatu(result);
}

double penalty_sum_internal(const arma::vec& norms,
                            const arma::vec& region_multiplier,
                            const arma::vec& adaptive_weight,
                            const double lambda,
                            const PenaltyFamily family,
                            const double shape) {
  double penalty = 0.0;
  for (arma::uword region = 0; region < norms.n_elem; ++region) {
    if (region_multiplier[region] == 0.0) {
      continue;
    }
    const double regional_lambda = lambda * adaptive_weight[region];
    penalty += region_multiplier[region] * penalty_value_scalar(
      norms[region], regional_lambda, family, shape
    );
  }
  return penalty;
}

double exact_objective(const arma::vec& y,
                       const arma::mat& design,
                       const arma::mat& roughness,
                       const arma::cube& region_gram,
                       const arma::vec& region_length,
                       const arma::vec& region_multiplier,
                       const arma::vec& adaptive_weight,
                       const arma::vec& coefficients,
                       const double kappa,
                       const double lambda,
                       const PenaltyFamily family,
                       const double shape) {
  const arma::vec residual = y - design * coefficients;
  const double loss = arma::dot(residual, residual) /
    static_cast<double>(y.n_elem);
  const double smoothness = kappa *
    arma::dot(coefficients, roughness * coefficients);
  const arma::vec norms = region_norms_internal(
    coefficients, region_gram, region_length
  );
  return loss + smoothness + penalty_sum_internal(
    norms, region_multiplier, adaptive_weight, lambda, family, shape
  );
}

arma::mat lqa_penalty_matrix(const arma::vec& coefficients,
                             const arma::cube& region_gram,
                             const arma::vec& region_length,
                             const arma::vec& region_multiplier,
                             const arma::vec& adaptive_weight,
                             const double lambda,
                             const PenaltyFamily family,
                             const double shape,
                             const double norm_floor,
                             const bool exact_denominator) {
  const arma::uword basis_count = coefficients.n_elem;
  const arma::uword region_count = region_gram.n_slices;
  arma::mat penalty(basis_count, basis_count, arma::fill::zeros);

  for (arma::uword region = 0; region < region_count; ++region) {
    if (region_multiplier[region] == 0.0) {
      continue;
    }

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
    if (exact_denominator && region_norm <= 0.0) {
      Rcpp::stop(
        "The ideal ARFC MM step requires every active regional norm to be "
        "strictly positive."
      );
    }

    const double regional_lambda = lambda * adaptive_weight[region];
    const double derivative = region_multiplier[region] *
      penalty_derivative_scalar(
        region_norm, regional_lambda, family, shape
      );
    if (derivative == 0.0) {
      continue;
    }

    const double denominator = exact_denominator ?
      region_norm : std::max(region_norm, norm_floor);
    const double multiplier = 0.5 * derivative /
      (denominator * region_length[region]);
    penalty += multiplier * region_gram.slice(region);
  }

  return arma::symmatu(penalty);
}

ArfcFit fit_single(const arma::vec& y,
                   const arma::mat& design,
                   const arma::mat& roughness,
                   const arma::cube& region_gram,
                   const arma::vec& region_length,
                   const arma::vec& region_multiplier,
                   const arma::vec& adaptive_weight,
                   const arma::mat& design_crossprod,
                   const arma::vec& design_response,
                   const arma::vec& initial,
                   const double kappa,
                   const double lambda,
                   const PenaltyFamily family,
                   const double shape,
                   const int max_iter,
                   const double tolerance,
                   const double norm_floor,
                   const double cutoff,
                   const bool refit,
                   const bool keep_trace) {
  ArfcFit result;
  const double sample_size = static_cast<double>(y.n_elem);
  arma::vec coefficients = initial;

  if (lambda == 0.0) {
    result.converged = true;
    if (keep_trace) {
      result.objective_trace.push_back(exact_objective(
        y, design, roughness, region_gram, region_length,
        region_multiplier, adaptive_weight, coefficients,
        kappa, lambda, family, shape
      ));
    }
  } else {
    double previous_objective = exact_objective(
      y, design, roughness, region_gram, region_length,
      region_multiplier, adaptive_weight, coefficients,
      kappa, lambda, family, shape
    );
    if (keep_trace) {
      result.objective_trace.push_back(previous_objective);
    }

    for (int iteration = 1; iteration <= max_iter; ++iteration) {
      const arma::mat lqa = lqa_penalty_matrix(
        coefficients, region_gram, region_length,
        region_multiplier, adaptive_weight, lambda, family, shape,
        norm_floor, false
      );
      const arma::mat system = design_crossprod +
        sample_size * kappa * roughness + sample_size * lqa;
      arma::vec updated;
      if (!solve_vector(system, design_response, updated)) {
        result.message = "The ARFC LQA system could not be solved.";
        result.coefficients = coefficients;
        return result;
      }

      const double current_objective = exact_objective(
        y, design, roughness, region_gram, region_length,
        region_multiplier, adaptive_weight, updated,
        kappa, lambda, family, shape
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
          "The ARFC active-set refit system could not be solved.";
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
      result.message = "ARFC degrees of freedom could not be computed.";
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
    region_multiplier, adaptive_weight, coefficients,
    kappa, lambda, family, shape
  );
  result.success = true;
  if (!result.converged) {
    result.message = "Maximum ARFC LQA iterations reached.";
  }
  return result;
}

}  // namespace


// [[Rcpp::export]]
Rcpp::DataFrame cpp_arfc_penalty_scalar_v1(const arma::vec& value,
                                           const double lambda,
                                           const std::string penalty_family,
                                           const double shape) {
  if (!value.is_finite() || !std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("ARFC scalar penalty inputs must be finite and lambda >= 0.");
  }
  const PenaltyFamily family = parse_penalty_family(penalty_family);
  validate_penalty_shape(family, shape);

  Rcpp::NumericVector penalty(value.n_elem);
  Rcpp::NumericVector derivative(value.n_elem);
  Rcpp::NumericVector second_derivative(value.n_elem);
  Rcpp::LogicalVector at_knot(value.n_elem);
  for (arma::uword index = 0; index < value.n_elem; ++index) {
    bool knot = false;
    penalty[index] = penalty_value_scalar(
      value[index], lambda, family, shape
    );
    derivative[index] = penalty_derivative_scalar(
      value[index], lambda, family, shape
    );
    second_derivative[index] = penalty_second_derivative_scalar(
      value[index], lambda, family, shape, knot
    );
    at_knot[index] = knot;
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("value") = value,
    Rcpp::Named("penalty") = penalty,
    Rcpp::Named("derivative") = derivative,
    Rcpp::Named("second_derivative") = second_derivative,
    Rcpp::Named("at_knot") = at_knot,
    Rcpp::Named("penalty_family") = penalty_family_name(family),
    Rcpp::Named("shape") = shape
  );
}


// [[Rcpp::export]]
arma::vec cpp_arfc_region_norms_v1(const arma::vec& coefficients,
                                   const arma::cube& region_gram,
                                   const arma::vec& region_length) {
  validate_region_inputs(coefficients, region_gram, region_length);
  return region_norms_internal(coefficients, region_gram, region_length);
}


// [[Rcpp::export]]
Rcpp::List cpp_arfc_curvature_envelope_v1(
    const arma::cube& region_gram,
    const arma::vec& region_length,
    const arma::vec& region_multiplier,
    const std::string penalty_family,
    const double shape) {
  arma::vec dummy(region_gram.n_rows, arma::fill::zeros);
  validate_region_inputs(dummy, region_gram, region_length);
  if (region_multiplier.n_elem != region_gram.n_slices ||
      !region_multiplier.is_finite() ||
      arma::any(region_multiplier < 0.0) ||
      arma::all(region_multiplier == 0.0)) {
    Rcpp::stop("Invalid ARFC regional multipliers.");
  }
  const PenaltyFamily family = parse_penalty_family(penalty_family);
  validate_penalty_shape(family, shape);

  const arma::mat mass_matrix = normalized_region_matrix(
    region_gram, region_length, region_multiplier
  );
  const double maximum_concavity = scalar_maximum_concavity(family, shape);
  const arma::mat envelope = maximum_concavity * mass_matrix;
  const arma::vec eigenvalues = arma::eig_sym(envelope);
  return Rcpp::List::create(
    Rcpp::Named("mass_matrix") = mass_matrix,
    Rcpp::Named("curvature_envelope") = envelope,
    Rcpp::Named("maximum_scalar_concavity") = maximum_concavity,
    Rcpp::Named("maximum_envelope_eigenvalue") = eigenvalues.max(),
    Rcpp::Named("penalty_family") = penalty_family_name(family),
    Rcpp::Named("shape") = shape
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_arfc_penalty_diagnostics_v1(
    const arma::vec& coefficients,
    const arma::cube& region_gram,
    const arma::vec& region_length,
    const arma::vec& region_multiplier,
    const arma::vec& adaptive_weight,
    const double lambda,
    const std::string penalty_family,
    const double shape) {
  validate_region_inputs(coefficients, region_gram, region_length);
  if (region_multiplier.n_elem != region_gram.n_slices ||
      adaptive_weight.n_elem != region_gram.n_slices ||
      !region_multiplier.is_finite() || !adaptive_weight.is_finite() ||
      arma::any(region_multiplier < 0.0) ||
      arma::all(region_multiplier == 0.0) ||
      arma::any(adaptive_weight <= 0.0) ||
      !std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("Invalid ARFC penalty diagnostic inputs.");
  }
  const PenaltyFamily family = parse_penalty_family(penalty_family);
  validate_penalty_shape(family, shape);

  const arma::uword basis_count = coefficients.n_elem;
  const arma::vec norms = region_norms_internal(
    coefficients, region_gram, region_length
  );
  arma::vec gradient(basis_count, arma::fill::zeros);
  arma::mat hessian(basis_count, basis_count, arma::fill::zeros);
  bool any_knot = false;

  for (arma::uword region = 0; region < norms.n_elem; ++region) {
    if (region_multiplier[region] == 0.0) {
      continue;
    }
    const double q = norms[region];
    if (q <= 0.0) {
      Rcpp::stop(
        "Penalty derivatives require strictly positive regional norms."
      );
    }
    const arma::mat A = region_gram.slice(region) / region_length[region];
    const arma::vec Aq = A * coefficients;
    const arma::vec grad_q = Aq / q;
    const arma::mat hess_q = A / q - (grad_q * grad_q.t()) / q;
    const double regional_lambda = lambda * adaptive_weight[region];
    const double first = penalty_derivative_scalar(
      q, regional_lambda, family, shape
    );
    bool at_knot = false;
    const double second = penalty_second_derivative_scalar(
      q, regional_lambda, family, shape, at_knot
    );
    if (at_knot) {
      any_knot = true;
      continue;
    }
    gradient += region_multiplier[region] * first * grad_q;
    hessian += region_multiplier[region] * (
      second * (grad_q * grad_q.t()) + first * hess_q
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("penalty") = penalty_sum_internal(
      norms, region_multiplier, adaptive_weight, lambda, family, shape
    ),
    Rcpp::Named("gradient") = gradient,
    Rcpp::Named("hessian") = arma::symmatu(hessian),
    Rcpp::Named("region_norms") = norms,
    Rcpp::Named("at_penalty_knot") = any_knot,
    Rcpp::Named("penalty_family") = penalty_family_name(family),
    Rcpp::Named("shape") = shape
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_arfc_mm_step_v1(
    const arma::vec& y,
    const arma::mat& design,
    const arma::mat& roughness,
    const arma::cube& region_gram,
    const arma::vec& region_length,
    const arma::vec& region_multiplier,
    const arma::vec& adaptive_weight,
    const arma::vec& current,
    const double lambda,
    const double kappa,
    const std::string penalty_family,
    const double shape) {
  validate_common_inputs(
    y, design, roughness, region_gram, region_length,
    region_multiplier, adaptive_weight
  );
  if (current.n_elem != design.n_cols || !current.is_finite() ||
      !std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(kappa) || kappa < 0.0) {
    Rcpp::stop("Invalid ARFC ideal-MM inputs.");
  }
  const PenaltyFamily family = parse_penalty_family(penalty_family);
  validate_penalty_shape(family, shape);

  const arma::mat design_crossprod = arma::symmatu(design.t() * design);
  const arma::vec design_response = design.t() * y;
  const arma::mat lqa = lqa_penalty_matrix(
    current, region_gram, region_length, region_multiplier,
    adaptive_weight, lambda, family, shape, 0.0, true
  );
  const arma::mat system = design_crossprod +
    static_cast<double>(y.n_elem) * kappa * roughness +
    static_cast<double>(y.n_elem) * lqa;
  arma::vec updated;
  if (!solve_vector(system, design_response, updated)) {
    Rcpp::stop("The ideal ARFC MM system could not be solved.");
  }
  const double before = exact_objective(
    y, design, roughness, region_gram, region_length,
    region_multiplier, adaptive_weight, current,
    kappa, lambda, family, shape
  );
  const double after = exact_objective(
    y, design, roughness, region_gram, region_length,
    region_multiplier, adaptive_weight, updated,
    kappa, lambda, family, shape
  );
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = updated,
    Rcpp::Named("objective_before") = before,
    Rcpp::Named("objective_after") = after,
    Rcpp::Named("nonincrease") = after <= before +
      1e-12 * (1.0 + std::abs(before)),
    Rcpp::Named("lqa_matrix") = lqa,
    Rcpp::Named("penalty_family") = penalty_family_name(family),
    Rcpp::Named("shape") = shape
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_arfc_grid_v1(const arma::vec& y,
                            const arma::mat& design,
                            const arma::mat& roughness,
                            const arma::cube& region_gram,
                            const arma::vec& region_length,
                            const arma::vec& region_multiplier,
                            const arma::vec& adaptive_weight,
                            const arma::vec& lambda,
                            const arma::vec& kappa,
                            const std::string penalty_family,
                            const double shape,
                            const int max_iter,
                            const double tolerance,
                            const double norm_floor,
                            const double cutoff,
                            const bool refit,
                            const bool keep_trace) {
  validate_common_inputs(
    y, design, roughness, region_gram, region_length,
    region_multiplier, adaptive_weight
  );
  const PenaltyFamily family = parse_penalty_family(penalty_family);
  validate_penalty_shape(family, shape);
  if (lambda.n_elem == 0 || kappa.n_elem == 0 ||
      arma::any(lambda < 0.0) || arma::any(kappa < 0.0) ||
      !lambda.is_finite() || !kappa.is_finite()) {
    Rcpp::stop("Invalid ARFC tuning grids.");
  }
  if (max_iter < 1 || tolerance <= 0.0 || norm_floor <= 0.0 ||
      cutoff < 0.0 || !std::isfinite(tolerance) ||
      !std::isfinite(norm_floor) || !std::isfinite(cutoff)) {
    Rcpp::stop("Invalid ARFC numerical controls.");
  }
  if (refit && arma::any(lambda > 0.0) && cutoff <= 0.0) {
    Rcpp::stop(
      "ARFC forbids `refit = TRUE` with a zero coefficient cutoff."
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
          "The ARFC roughness initialization could not be solved.";
        traces[candidate] = R_NilValue;
        ++candidate;
        continue;
      }

      const ArfcFit fit = fit_single(
        y, design, roughness, region_gram, region_length,
        region_multiplier, adaptive_weight, design_crossprod,
        design_response, initial, kappa[k], lambda[l], family, shape,
        max_iter, tolerance, norm_floor, cutoff, refit, keep_trace
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

  const arma::mat mass_matrix = normalized_region_matrix(
    region_gram, region_length, region_multiplier
  );
  const double maximum_concavity = scalar_maximum_concavity(family, shape);
  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("region_norms") = region_norms,
    Rcpp::Named("lambda") = lambda_out,
    Rcpp::Named("kappa") = kappa_out,
    Rcpp::Named("penalty_family") = penalty_family_name(family),
    Rcpp::Named("shape") = shape,
    Rcpp::Named("maximum_scalar_concavity") = maximum_concavity,
    Rcpp::Named("regional_mass_matrix") = mass_matrix,
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
