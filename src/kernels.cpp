#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <vector>

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {

constexpr double kInvSqrtTwoPi = 0.39894228040143267794;
constexpr double kPi = 3.14159265358979323846;
constexpr int kMaximumCompactMoment = 9;
constexpr int kCompactSummaryColumns = 14;

double binomial_coefficient(const int n, const int k) {
  static constexpr double coefficients[10][10] = {
    {1},
    {1, 1},
    {1, 2, 1},
    {1, 3, 3, 1},
    {1, 4, 6, 4, 1},
    {1, 5, 10, 10, 5, 1},
    {1, 6, 15, 20, 15, 6, 1},
    {1, 7, 21, 35, 35, 21, 7, 1},
    {1, 8, 28, 56, 70, 56, 28, 8, 1},
    {1, 9, 36, 84, 126, 126, 84, 36, 9, 1}
  };
  return coefficients[n][k];
}

double translated_moment(const NumericMatrix& summaries, const int node,
                         const int degree, const double shift) {
  double out = 0.0;
  double shift_power = 1.0;
  for (int index = degree; index >= 0; --index) {
    out += binomial_coefficient(degree, index) * shift_power *
      summaries(node, index);
    shift_power *= shift;
  }
  return out;
}

void add_translated_summary(NumericMatrix& summaries, const int parent,
                            const int child, const double shift) {
  for (int degree = 0; degree <= kMaximumCompactMoment; ++degree) {
    double translated = 0.0;
    double shift_power = 1.0;
    for (int index = degree; index >= 0; --index) {
      translated += binomial_coefficient(degree, index) * shift_power *
        summaries(child, index);
      shift_power *= shift;
    }
    summaries(parent, degree) += translated;
  }

  const double angle1 = kPi * shift;
  const double cosine1 = std::cos(angle1);
  const double sine1 = std::sin(angle1);
  summaries(parent, 10) += summaries(child, 10) * cosine1 -
    summaries(child, 11) * sine1;
  summaries(parent, 11) += summaries(child, 11) * cosine1 +
    summaries(child, 10) * sine1;

  const double angle2 = 0.5 * kPi * shift;
  const double cosine2 = std::cos(angle2);
  const double sine2 = std::sin(angle2);
  summaries(parent, 12) += summaries(child, 12) * cosine2 -
    summaries(child, 13) * sine2;
  summaries(parent, 13) += summaries(child, 13) * cosine2 +
    summaries(child, 12) * sine2;
}

double compact_node_sum(const int node, const int kernel_id,
                        const int side, const double y, const double h,
                        const int base, const NumericVector& centers,
                        const NumericMatrix& summaries) {
  if (!std::isfinite(summaries(node, 0))) {
    if (node >= base) {
      return 0.0;
    }
    return compact_node_sum(
      2 * node, kernel_id, side, y, h, base, centers, summaries
    ) + compact_node_sum(
      2 * node + 1, kernel_id, side, y, h, base, centers, summaries
    );
  }

  const double shift = (centers[node] - y) / h;
  const double moment0 = summaries(node, 0);

  if (kernel_id == 3) {
    return 0.5 * moment0;
  }
  if (kernel_id == 4) {
    const double moment1 = translated_moment(
      summaries, node, 1, shift
    );
    return side < 0 ? moment0 + moment1 : moment0 - moment1;
  }
  if (kernel_id == 2) {
    const double moment2 = translated_moment(
      summaries, node, 2, shift
    );
    return 0.75 * (moment0 - moment2);
  }
  if (kernel_id == 5) {
    const double moment2 = translated_moment(
      summaries, node, 2, shift
    );
    const double moment4 = translated_moment(
      summaries, node, 4, shift
    );
    return 15.0 / 16.0 * (moment0 - 2.0 * moment2 + moment4);
  }
  if (kernel_id == 8) {
    const double moment2 = translated_moment(
      summaries, node, 2, shift
    );
    const double moment4 = translated_moment(
      summaries, node, 4, shift
    );
    const double moment6 = translated_moment(
      summaries, node, 6, shift
    );
    return 35.0 / 32.0 *
      (moment0 - 3.0 * moment2 + 3.0 * moment4 - moment6);
  }
  if (kernel_id == 9) {
    const double moment3 = translated_moment(
      summaries, node, 3, shift
    );
    const double moment6 = translated_moment(
      summaries, node, 6, shift
    );
    const double moment9 = translated_moment(
      summaries, node, 9, shift
    );
    return side < 0
      ? 70.0 / 81.0 *
        (moment0 + 3.0 * moment3 + 3.0 * moment6 + moment9)
      : 70.0 / 81.0 *
        (moment0 - 3.0 * moment3 + 3.0 * moment6 - moment9);
  }
  if (kernel_id == 6) {
    const double angle = kPi * shift;
    const double cosine_sum = summaries(node, 10) * std::cos(angle) -
      summaries(node, 11) * std::sin(angle);
    return 0.5 * (moment0 + cosine_sum);
  }
  if (kernel_id == 7) {
    const double angle = 0.5 * kPi * shift;
    const double cosine_sum = summaries(node, 12) * std::cos(angle) -
      summaries(node, 13) * std::sin(angle);
    return kPi / 4.0 * cosine_sum;
  }

  // All kernel ids accepted by the sorted evaluator are handled above.
  // # nocov start
  return NA_REAL;
  // # nocov end
}

double compact_range_sum(int left, int right, const int kernel_id,
                         const int side, const double y, const double h,
                         const int base, const NumericVector& centers,
                         const NumericMatrix& summaries) {
  double out = 0.0;
  left += base;
  right += base;

  while (left < right) {
    if (left & 1) {
      out += compact_node_sum(
        left, kernel_id, side, y, h, base, centers, summaries
      );
      ++left;
    }
    if (right & 1) {
      --right;
      out += compact_node_sum(
        right, kernel_id, side, y, h, base, centers, summaries
      );
    }
    left /= 2;
    right /= 2;
  }
  return out;
}

// Gamma-kernel shape at evaluation point y; the scale is always h.
// kernel_id 13 is Chen's (2000) first kernel, shape y / h + 1 (mean y + h).
// kernel_id 12 is his modified kernel, shape rho_b(y) = y / h for y >= 2h and
// (y / (2h))^2 + 1 below, whose interior mean is exactly y.
double gamma_kernel_shape(const double y, const double h,
                          const int kernel_id) {
  if (kernel_id == 13) {
    return y / h + 1.0;
  }
  const double ratio = y / h;
  return y >= 2.0 * h ? ratio : 0.25 * ratio * ratio + 1.0;
}

double beta_boundary_parameter(const double y, const double h,
                               const double rho_base) {
  return 2.0 * h * h + 2.5 -
    std::sqrt(rho_base - y * y - y / h);
}

int start_parameter_count(const int start_id) {
  switch (start_id) {
    case 1: return 2;  // Normal: mean, sd.
    case 2: return 1;  // Exponential: rate.
    case 3: return 2;  // Lognormal: meanlog, sdlog.
    case 4: return 2;  // Gamma: shape, rate.
    case 5: return 2;  // Weibull: shape, scale.
    case 6: return 2;  // Beta: shape1, shape2.
    // Callers validate start ids before reaching this helper.
    // # nocov start
    default: return 0;
    // # nocov end
  }
}

double start_log_constant(const int start_id,
                          const NumericMatrix& parameters,
                          const R_xlen_t row) {
  const double p0 = parameters(row, 0);

  switch (start_id) {
    case 1:
      return std::log(kInvSqrtTwoPi) - std::log(parameters(row, 1));
    case 2:
      return std::log(p0);
    case 3:
      return std::log(kInvSqrtTwoPi) - std::log(parameters(row, 1));
    case 4:
      return p0 * std::log(parameters(row, 1)) - std::lgamma(p0);
    case 5:
      return std::log(p0) - std::log(parameters(row, 1));
    case 6:
      return std::lgamma(p0 + parameters(row, 1)) -
        std::lgamma(p0) - std::lgamma(parameters(row, 1));
    // Callers validate start ids before reaching this helper.
    // # nocov start
    default:
      return NA_REAL;
    // # nocov end
  }
}

double start_log_density(const double x, const double log_x,
                         const double log1m_x, const int start_id,
                         const NumericMatrix& parameters,
                         const R_xlen_t row, const double constant) {
  const double p0 = parameters(row, 0);

  switch (start_id) {
    case 1: {
      const double z = (x - p0) / parameters(row, 1);
      return constant - 0.5 * z * z;
    }
    case 2:
      return x >= 0.0 ? constant - p0 * x : R_NegInf;
    case 3: {
      if (!(x > 0.0)) {
        return R::dlnorm(x, p0, parameters(row, 1), true);
      }
      const double z = (log_x - p0) / parameters(row, 1);
      return constant - 0.5 * z * z - log_x;
    }
    case 4:
      return x > 0.0
        ? constant + (p0 - 1.0) * log_x - parameters(row, 1) * x
        : R::dgamma(x, p0, 1.0 / parameters(row, 1), true);
    case 5: {
      if (!(x > 0.0)) {
        return R::dweibull(x, p0, parameters(row, 1), true);
      }
      const double z = log_x - std::log(parameters(row, 1));
      return constant + (p0 - 1.0) * z - std::exp(p0 * z);
    }
    case 6:
      return x > 0.0 && x < 1.0
        ? constant + (p0 - 1.0) * log_x +
          (parameters(row, 1) - 1.0) * log1m_x
        : R::dbeta(x, p0, parameters(row, 1), true);
    // Callers validate start ids before reaching this helper.
    // # nocov start
    default:
      return NA_REAL;
    // # nocov end
  }
}

// Shape and normalization terms depend only on the evaluation point.
struct PairKernelConstants {
  double shape1 = 0.0;
  double shape2 = 0.0;
  double log_h = 0.0;
  double log_gamma1 = 0.0;
  double log_gamma2 = 0.0;
  double log_gamma_sum = 0.0;
};

PairKernelConstants pair_kernel_constants(double y, double h, int kernel_id) {
  PairKernelConstants c;
  if (kernel_id < 12) return c;
  c.log_h = std::log(h);
  if (kernel_id <= 13) {
    c.shape1 = gamma_kernel_shape(y, h, kernel_id);
    c.log_gamma1 = std::lgamma(c.shape1);
  } else {
    c.shape1 = y / h + 1.0;
    c.shape2 = (1.0 - y) / h + 1.0;
    if (kernel_id == 14) {
      const double h2 = h * h;
      const double rho_base = 4.0 * h2 * h2 + 6.0 * h2 + 2.25;
      if (y < 2.0 * h) {
        c.shape1 = beta_boundary_parameter(y, h, rho_base);
        c.shape2 = (1.0 - y) / h;
      } else if (y < 1.0 - 2.0 * h) {
        c.shape1 = y / h;
        c.shape2 = (1.0 - y) / h;
      } else {
        c.shape1 = y / h;
        c.shape2 = beta_boundary_parameter(1.0 - y, h, rho_base);
      }
    }
    c.log_gamma1 = std::lgamma(c.shape1);
    c.log_gamma2 = std::lgamma(c.shape2);
    c.log_gamma_sum = std::lgamma(c.shape1 + c.shape2);
  }
  return c;
}

double pair_log_kernel(const double y, const double x, const double h,
                       const int kernel_id, const double qy,
                       const double qx, const double log_x,
                       const double log1m_x, const PairKernelConstants& c) {
  const double u = (x - y) / h;

  if (kernel_id == 1) {
    return std::log(kInvSqrtTwoPi) - 0.5 * u * u;
  }
  if (kernel_id == 2) {
    const double value = 0.75 * (1.0 - u * u);
    return std::abs(u) < 1.0 ? std::log(value) : R_NegInf;
  }
  if (kernel_id == 3) {
    return std::abs(u) <= 1.0 ? std::log(0.5) : R_NegInf;
  }
  if (kernel_id == 4) {
    const double abs_u = std::abs(u);
    return abs_u < 1.0 ? std::log1p(-abs_u) : R_NegInf;
  }
  if (kernel_id == 5) {
    const double one_minus_u2 = 1.0 - u * u;
    return std::abs(u) < 1.0
      ? std::log(15.0 / 16.0) + 2.0 * std::log(one_minus_u2)
      : R_NegInf;
  }
  if (kernel_id == 6) {
    const double value = (1.0 + std::cos(kPi * u)) / 2.0;
    return std::abs(u) < 1.0 ? std::log(value) : R_NegInf;
  }
  if (kernel_id == 7) {
    return std::abs(u) < 1.0
      ? std::log(kPi / 4.0) + std::log(std::cos(kPi / 2.0 * u))
      : R_NegInf;
  }
  if (kernel_id == 8) {
    const double one_minus_u2 = 1.0 - u * u;
    return std::abs(u) < 1.0
      ? std::log(35.0 / 32.0) + 3.0 * std::log(one_minus_u2)
      : R_NegInf;
  }
  if (kernel_id == 9) {
    const double abs_u = std::abs(u);
    const double one_minus_u3 = 1.0 - abs_u * abs_u * abs_u;
    return abs_u < 1.0
      ? std::log(70.0 / 81.0) + 3.0 * std::log(one_minus_u3)
      : R_NegInf;
  }
  if (kernel_id == 10) {
    return std::log(0.5) - std::abs(u);
  }
  if (kernel_id == 11) {
    const double rho = 1.0 - h * h;
    const double inside =
      rho * rho * (qy * qy + qx * qx) - 2.0 * rho * qy * qx;
    return std::log(h) - 0.5 * std::log(1.0 - rho * rho) -
      inside / (2.0 * (1.0 - rho * rho));
  }
  if (kernel_id == 12 || kernel_id == 13) {
    return x > 0.0
      ? (c.shape1 - 1.0) * (log_x - c.log_h) - x / h - c.log_gamma1
      : c.log_h + R::dgamma(x, c.shape1, h, true);
  }
  if (kernel_id == 14 || kernel_id == 15) {
    if (x > 0.0 && x < 1.0) {
      return c.log_h + (c.shape1 - 1.0) * log_x +
        (c.shape2 - 1.0) * log1m_x - c.log_gamma1 -
        c.log_gamma2 + c.log_gamma_sum;
    }
    return c.log_h + R::dbeta(x, c.shape1, c.shape2, true);
  }

  // Callers validate kernel ids before reaching this helper.
  // # nocov start
  return NA_REAL;
  // # nocov end
}

}  // namespace

// Build a segment-tree cache of local weighted moments for compact kernels.
// A node is cached only when its data span is no wider than 2h. Such nodes
// have stable local coordinates whenever they can occur in a kernel window.
//
// [[Rcpp::export]]
List compact_kernel_cache_cpp(const NumericVector& sorted_x,
                              const NumericVector& sorted_weights,
                              const double h) {
  const R_xlen_t n = sorted_x.size();
  if (n == 0 || sorted_weights.size() != n) {
    Rcpp::stop("'sorted_x' and 'sorted_weights' must have equal positive length.");
  }
  if (!std::isfinite(h) || h <= 0.0) {
    Rcpp::stop("'h' must be finite and positive.");
  }
  for (R_xlen_t index = 0; index < n; ++index) {
    if (!std::isfinite(sorted_x[index]) ||
        !std::isfinite(sorted_weights[index])) {
      Rcpp::stop("'sorted_x' and 'sorted_weights' must be finite.");
    }
  }
  if (!std::is_sorted(sorted_x.begin(), sorted_x.end())) {
    Rcpp::stop("'sorted_x' must be sorted.");
  }

  int base = 1;
  while (base < n) {
    if (base > std::numeric_limits<int>::max() / 2) {
      // R vectors cannot reach this size on supported platforms.
      // # nocov start
      Rcpp::stop("The compact-kernel cache is too large.");
      // # nocov end
    }
    base *= 2;
  }
  const int node_count = 2 * base;
  NumericVector centers(node_count, NA_REAL);
  NumericMatrix summaries(node_count, kCompactSummaryColumns);
  std::vector<double> minimum(node_count, NA_REAL);
  std::vector<double> maximum(node_count, NA_REAL);
  std::vector<int> counts(node_count, 0);

  for (R_xlen_t index = 0; index < n; ++index) {
    const int node = base + static_cast<int>(index);
    const double weight = sorted_weights[index];
    centers[node] = sorted_x[index];
    minimum[node] = sorted_x[index];
    maximum[node] = sorted_x[index];
    counts[node] = 1;
    summaries(node, 0) = weight;
    summaries(node, 10) = weight;
    summaries(node, 12) = weight;
  }

  for (int node = base - 1; node >= 1; --node) {
    if ((node & 16383) == 0) {
      // Polling is operational rather than data-dependent behavior.
      // # nocov start
      Rcpp::checkUserInterrupt();
      // # nocov end
    }
    const int left = 2 * node;
    const int right = left + 1;
    counts[node] = counts[left] + counts[right];
    if (counts[node] == 0) {
      continue;
    }
    if (counts[left] == 0 || counts[right] == 0) {
      const int child = counts[left] > 0 ? left : right;
      centers[node] = centers[child];
      minimum[node] = minimum[child];
      maximum[node] = maximum[child];
      for (int column = 0; column < kCompactSummaryColumns; ++column) {
        summaries(node, column) = summaries(child, column);
      }
      continue;
    }

    minimum[node] = minimum[left];
    maximum[node] = maximum[right];
    centers[node] = minimum[node] +
      0.5 * (maximum[node] - minimum[node]);
    const bool stable =
      (maximum[node] - minimum[node]) / h <=
      2.0 + 32.0 * std::numeric_limits<double>::epsilon();
    if (!stable || !std::isfinite(summaries(left, 0)) ||
        !std::isfinite(summaries(right, 0))) {
      summaries(node, 0) = NA_REAL;
      continue;
    }

    add_translated_summary(
      summaries,
      node,
      left,
      (centers[left] - centers[node]) / h
    );
    add_translated_summary(
      summaries,
      node,
      right,
      (centers[right] - centers[node]) / h
    );
  }

  return List::create(
    Rcpp::Named("base") = base,
    Rcpp::Named("centers") = centers,
    Rcpp::Named("summaries") = summaries
  );
}

// Exact sorted-data evaluation for compact polynomial/trigonometric kernels
// and for the Laplace kernel. The weights already include kdensity's
// 1 / (n * h) factor.
//
// [[Rcpp::export]]
NumericVector kernel_weighted_sum_sorted_cpp(
    const NumericVector& y,
    const NumericVector& sorted_x,
    const NumericVector& sorted_weights,
    const double h,
    const int kernel_id,
    const List& cache) {
  const R_xlen_t nx = sorted_x.size();
  const R_xlen_t ny = y.size();
  if (nx == 0 || sorted_weights.size() != nx) {
    Rcpp::stop("'sorted_x' and 'sorted_weights' must have equal positive length.");
  }
  if (!std::isfinite(h) || h <= 0.0) {
    Rcpp::stop("'h' must be finite and positive.");
  }
  for (R_xlen_t index = 0; index < nx; ++index) {
    if (!std::isfinite(sorted_x[index]) ||
        !std::isfinite(sorted_weights[index])) {
      Rcpp::stop("'sorted_x' and 'sorted_weights' must be finite.");
    }
  }
  if (kernel_id < 2 || kernel_id > 10) {
    Rcpp::stop("The sorted evaluator supports kernel ids 2 through 10.");
  }
  if (!std::is_sorted(sorted_x.begin(), sorted_x.end())) {
    Rcpp::stop("'sorted_x' must be sorted.");
  }

  NumericVector out(ny);
  if (kernel_id == 10) {
    NumericVector left_values(ny);
    std::vector<R_xlen_t> order;
    order.reserve(ny);
    for (R_xlen_t index = 0; index < ny; ++index) {
      if (std::isnan(y[index])) {
        out[index] = y[index];
      } else if (!std::isfinite(y[index])) {
        out[index] = 0.0;
      } else {
        order.push_back(index);
      }
    }
    std::stable_sort(
      order.begin(),
      order.end(),
      [&y](const R_xlen_t left, const R_xlen_t right) {
        return y[left] < y[right];
      }
    );

    R_xlen_t data_index = 0;
    double left_sum = 0.0;
    double previous = 0.0;
    bool first = true;
    for (const R_xlen_t index : order) {
      const double point = y[index];
      if (!first) {
        left_sum *= std::exp(-(point - previous) / h);
      }
      while (data_index < nx && sorted_x[data_index] <= point) {
        left_sum += sorted_weights[data_index] *
          std::exp(-(point - sorted_x[data_index]) / h);
        ++data_index;
      }
      left_values[index] = left_sum;
      previous = point;
      first = false;
    }

    data_index = nx;
    double right_sum = 0.0;
    first = true;
    for (auto iterator = order.rbegin(); iterator != order.rend(); ++iterator) {
      const R_xlen_t index = *iterator;
      const double point = y[index];
      if (!first) {
        right_sum *= std::exp(-(previous - point) / h);
      }
      while (data_index > 0 && sorted_x[data_index - 1] > point) {
        --data_index;
        right_sum += sorted_weights[data_index] *
          std::exp(-(sorted_x[data_index] - point) / h);
      }
      out[index] = 0.5 * (left_values[index] + right_sum);
      previous = point;
      first = false;
    }
  } else {
    if (!cache.containsElementNamed("base") ||
        !cache.containsElementNamed("centers") ||
        !cache.containsElementNamed("summaries")) {
      Rcpp::stop("The compact-kernel cache is incomplete.");
    }
    const int base = Rcpp::as<int>(cache["base"]);
    const NumericVector centers = cache["centers"];
    const NumericMatrix summaries = cache["summaries"];
    if (centers.size() != 2 * base || summaries.nrow() != 2 * base ||
        summaries.ncol() != kCompactSummaryColumns || base < nx) {
      Rcpp::stop("The compact-kernel cache has incompatible dimensions.");
    }

    for (R_xlen_t index = 0; index < ny; ++index) {
      if ((index & 1023) == 0) {
        Rcpp::checkUserInterrupt();
      }
      const double point = y[index];
      if (std::isnan(point)) {
        out[index] = point;
        continue;
      }
      if (!std::isfinite(point)) {
        out[index] = 0.0;
        continue;
      }

      // Select the window with the SAME divided predicate the direct path and
      // stats::density use: include x iff abs((x - point) / h) <= 1. The
      // subtraction bounds point +/- h are not float-equivalent at the
      // boundary, which matters for the uniform kernel (id 3), whose value
      // jumps to 0.5 there rather than vanishing like the polynomial/trig
      // kernels. (x - point) / h is monotone increasing in x, so binary search
      // over the sorted data stays valid.
      const int left = static_cast<int>(
        std::lower_bound(
          sorted_x.begin(), sorted_x.end(), point,
          [h](double xv, double pt) { return (xv - pt) / h < -1.0; }
        ) - sorted_x.begin()
      );
      const int right = static_cast<int>(
        std::lower_bound(
          sorted_x.begin(), sorted_x.end(), point,
          [h](double xv, double pt) { return (xv - pt) / h <= 1.0; }
        ) - sorted_x.begin()
      );
      double value = 0.0;
      if (kernel_id == 4 || kernel_id == 9) {
        const int middle = static_cast<int>(
          std::upper_bound(
            sorted_x.begin() + left,
            sorted_x.begin() + right,
            point
          ) - sorted_x.begin()
        );
        value += compact_range_sum(
          left, middle, kernel_id, -1, point, h,
          base, centers, summaries
        );
        value += compact_range_sum(
          middle, right, kernel_id, 1, point, h,
          base, centers, summaries
        );
      } else {
        value = compact_range_sum(
          left, right, kernel_id, 0, point, h,
          base, centers, summaries
        );
      }
      out[index] = std::max(0.0, value);
    }
  }

  if (y.hasAttribute("names")) {
    out.attr("names") = y.attr("names");
  }
  return out;
}

// Exact weighted sums for the built-in kernels. The weights include the
// 1 / (n * h) factor used by kdensity.
//
// [[Rcpp::export]]
NumericVector kernel_weighted_sum_cpp(const NumericVector& y,
                                      const NumericVector& x,
                                      const NumericVector& weights,
                                      const double h,
                                      const int kernel_id,
                                      const NumericVector& cache1,
                                      const NumericVector& cache2) {
  const R_xlen_t nx = x.size();
  const R_xlen_t ny = y.size();
  NumericVector out = Rcpp::no_init(ny);

  if (weights.size() != nx) {
    Rcpp::stop("'weights' and 'x' must have the same length.");
  }
  if (kernel_id < 1 || kernel_id > 15) {
    Rcpp::stop("Unknown compiled kernel id.");
  }
  if (kernel_id == 11 && cache1.size() != nx) {
    Rcpp::stop("The Gaussian copula cache has the wrong length.");
  }
  if (kernel_id >= 12 && kernel_id <= 15 && cache1.size() != nx) {
    Rcpp::stop("The log-data cache has the wrong length.");
  }
  if (kernel_id >= 14 && cache2.size() != nx) {
    Rcpp::stop("The beta-kernel cache has the wrong length.");
  }

  const double h2 = h * h;
  const double log_h = std::log(h);
  const double rho = 1.0 - h2;
  const double copula_denominator = 2.0 * (1.0 - rho * rho);
  const double lower_cut = 2.0 * h;
  const double upper_cut = 1.0 - lower_cut;
  const double rho_base = 4.0 * h2 * h2 + 6.0 * h2 + 2.25;

  for (R_xlen_t i = 0; i < ny; ++i) {
    if ((i & 63) == 0) {
      Rcpp::checkUserInterrupt();
    }

    const double yi = y[i];
    if (std::isnan(yi)) {
      out[i] = yi;
      continue;
    }

    double shape1 = 0.0;
    double shape2 = 0.0;
    double qy = 0.0;

    if (kernel_id == 11) {
      qy = R::qnorm5(yi, 0.0, 1.0, true, false);
    } else if (kernel_id == 12 || kernel_id == 13) {
      shape1 = gamma_kernel_shape(yi, h, kernel_id);
    } else if (kernel_id == 14) {
      if (yi < lower_cut) {
        shape1 = beta_boundary_parameter(yi, h, rho_base);
        shape2 = (1.0 - yi) / h;
      } else if (yi < upper_cut) {
        shape1 = yi / h;
        shape2 = (1.0 - yi) / h;
      } else {
        shape1 = yi / h;
        shape2 = beta_boundary_parameter(1.0 - yi, h, rho_base);
      }
    } else if (kernel_id == 15) {
      shape1 = yi / h + 1.0;
      shape2 = (1.0 - yi) / h + 1.0;
    }

    const double gamma_constant = (kernel_id == 12 || kernel_id == 13)
      ? -std::lgamma(shape1) - (shape1 - 1.0) * log_h
      : 0.0;
    const double beta_constant = (kernel_id == 14 || kernel_id == 15)
      ? log_h - std::lgamma(shape1) - std::lgamma(shape2) +
        std::lgamma(shape1 + shape2)
      : 0.0;

    double sum = 0.0;
    for (R_xlen_t j = 0; j < nx; ++j) {
      if (i == 0 && j != 0 && (j & 1048575) == 0) {
        // Polling is operational rather than data-dependent behavior.
        // # nocov start
        Rcpp::checkUserInterrupt();
        // # nocov end
      }

      const double xj = x[j];
      double value = 0.0;

      if (kernel_id == 1) {
        const double u = (xj - yi) / h;
        value = kInvSqrtTwoPi * std::exp(-0.5 * u * u);
      } else if (kernel_id == 2) {
        const double u = (xj - yi) / h;
        value = std::abs(u) <= 1.0 ? 0.75 * (1.0 - u * u) : 0.0;
      } else if (kernel_id == 3) {
        const double u = (xj - yi) / h;
        value = std::abs(u) <= 1.0 ? 0.5 : 0.0;
      } else if (kernel_id == 4) {
        const double abs_u = std::abs((xj - yi) / h);
        value = abs_u <= 1.0 ? 1.0 - abs_u : 0.0;
      } else if (kernel_id == 5) {
        const double u = (xj - yi) / h;
        const double one_minus_u2 = 1.0 - u * u;
        value = std::abs(u) <= 1.0
          ? 15.0 / 16.0 * one_minus_u2 * one_minus_u2
          : 0.0;
      } else if (kernel_id == 6) {
        const double u = (xj - yi) / h;
        value = std::abs(u) <= 1.0
          ? (1.0 + std::cos(kPi * u)) / 2.0
          : 0.0;
      } else if (kernel_id == 7) {
        const double u = (xj - yi) / h;
        value = std::abs(u) <= 1.0
          ? kPi / 4.0 * std::cos(kPi / 2.0 * u)
          : 0.0;
      } else if (kernel_id == 8) {
        const double u = (xj - yi) / h;
        const double one_minus_u2 = 1.0 - u * u;
        value = std::abs(u) <= 1.0
          ? 35.0 / 32.0 * one_minus_u2 * one_minus_u2 * one_minus_u2
          : 0.0;
      } else if (kernel_id == 9) {
        const double abs_u = std::abs((xj - yi) / h);
        const double one_minus_u3 = 1.0 - abs_u * abs_u * abs_u;
        value = abs_u <= 1.0
          ? 70.0 / 81.0 * one_minus_u3 * one_minus_u3 * one_minus_u3
          : 0.0;
      } else if (kernel_id == 10) {
        value = 0.5 * std::exp(-std::abs((xj - yi) / h));
      } else if (kernel_id == 11) {
        const double qx = cache1[j];
        const double inside =
          rho * rho * (qy * qy + qx * qx) - 2.0 * rho * qy * qx;
        value = h / std::sqrt(1.0 - rho * rho) *
          std::exp(-inside / copula_denominator);
      } else if (kernel_id == 12 || kernel_id == 13) {
        value = xj > 0.0
          ? std::exp((shape1 - 1.0) * cache1[j] - xj / h +
                     gamma_constant)
          : h * R::dgamma(xj, shape1, h, false);
      } else if (kernel_id == 14 || kernel_id == 15) {
        value = xj > 0.0 && xj < 1.0
          ? std::exp((shape1 - 1.0) * cache1[j] +
                     (shape2 - 1.0) * cache2[j] + beta_constant)
          : h * R::dbeta(xj, shape1, shape2, false);
      }

      sum += weights[j] * value;
    }
    out[i] = sum;
  }

  if (y.hasAttribute("names")) {
    out.attr("names") = y.attr("names");
  }

  return out;
}

namespace {

// The embedded Newton fits below report failure the "Result<>" way -- via an
// empty std::optional -- instead of throwing. A leave-one-out fit that hits a
// NaN or a degenerate configuration returns std::nullopt; the caller writes NaN
// into that row and the R side (make_compiled_loo_parameters / ucv) falls back
// to the univariateML path. std::optional needs no error payload here because
// the only failure information the caller acts on is "did it converge".
struct BetaShapes {
  double alpha;
  double exp_b;
};

std::optional<double> fit_gamma_shape_loo(const double s, const double tolerance,
                                          const int iteration_limit) {
  double shape = 1.0 / (12.0 * s) *
    (3.0 - s + std::sqrt((s - 3.0) * (s - 3.0) + 24.0 * s)) *
    0.9989 + 0.0102;

  for (int iteration = 0; iteration < iteration_limit; ++iteration) {
    const double step =
      (std::log(shape) - R::digamma(shape) - s) /
      (1.0 / shape - R::trigamma(shape));
    const double next = shape - step;
    if (std::isnan(next)) {
      return std::nullopt;
    }
    if (std::abs((shape - next) / shape) < tolerance) {
      shape = next;
      break;
    }
    shape = next;
  }

  return shape;
}

std::optional<BetaShapes> fit_beta_shapes_loo(const double s, const double r,
                                              const double tolerance,
                                              const int iteration_limit) {
  const double g1 = std::exp(s);
  const double g2 = std::exp(r);
  const double b_init = std::log(0.5 + 0.5 * g2 / (1.0 - (g1 + g2)));
  if (!std::isfinite(b_init)) {
    return std::nullopt;
  }

  double b = b_init;
  double alpha = NA_REAL;
  double exp_b = NA_REAL;

  for (int iteration = 0; iteration < iteration_limit; ++iteration) {
    exp_b = std::exp(b);
    const double digamma_b_r = R::digamma(exp_b) - r;
    const double trigamma_b = R::trigamma(exp_b);
    const double y = s + digamma_b_r;
    alpha = 1.0 / std::log1p(std::exp(-y));
    for (int inner = 0; inner < 3; ++inner) {
      alpha += (y - R::digamma(alpha)) / R::trigamma(alpha);
    }
    const double f =
      exp_b * (digamma_b_r - R::digamma(alpha + exp_b));
    const double df =
      f + exp_b * exp_b *
        (trigamma_b - R::trigamma(alpha + exp_b) *
          (trigamma_b / R::trigamma(alpha) + 1.0));
    const double next = b - f / df;
    if (std::isnan(next)) {
      return std::nullopt;
    }
    if (std::abs((next - b) / next) < tolerance) {
      break;
    }
    b = next;
  }

  return BetaShapes{alpha, exp_b};
}

std::optional<double> fit_weibull_sigma_loo(const NumericVector& negative_log_x,
                                            const double mean_z,
                                            const R_xlen_t i, const R_xlen_t n,
                                            const double n_loo,
                                            const double sigma_init,
                                            const double tolerance,
                                            const int iteration_limit) {
  double sigma = sigma_init;

  for (int iteration = 0; iteration < iteration_limit; ++iteration) {
    const double negative_sigma_inverse = -1.0 / sigma;
    double psi0 = 0.0;
    double psi1 = 0.0;
    double psi2 = 0.0;

    for (R_xlen_t j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      const double centered = negative_log_x[j] - mean_z;
      const double exponential =
        std::exp(centered * negative_sigma_inverse);
      psi0 += exponential;
      psi1 += centered * exponential;
      psi2 += centered * centered * exponential;
    }
    psi1 /= psi0;
    psi2 = psi2 / psi0 - psi1 * psi1;
    const double step =
      (sigma + psi1) /
      (1.0 + negative_sigma_inverse * negative_sigma_inverse * psi2);
    const double next = sigma - step;
    if (std::isnan(next)) {
      return std::nullopt;
    }
    if (std::abs((sigma - next) / sigma) < tolerance) {
      sigma = next;
      break;
    }
    sigma = next;
  }

  return sigma;
}

}  // namespace

// Batched leave-one-out MLEs for designated starts. The algorithms are adapted
// from univariateML/cppunivariate and embedded here so kdensity does not depend
// on an unpublished C++ header ABI.
//
// [[Rcpp::export]]
NumericMatrix designated_loo_parameters_cpp(const NumericVector& x,
                                            const int start_id) {
  const R_xlen_t n = x.size();
  if (start_id < 1 || start_id > 6) {
    Rcpp::stop("Unknown compiled start id.");
  }

  const int parameter_count = start_parameter_count(start_id);
  NumericMatrix out(n, parameter_count);

  // Leave-one-out MLEs need at least two retained observations. Report the
  // degenerate case with NaN rather than throwing, so the R caller falls back
  // to the univariateML path instead of aborting kdensity().
  if (n < 3) {
    std::fill(out.begin(), out.end(), NA_REAL);
    return out;
  }
  const double n_loo = static_cast<double>(n - 1);
  const double tolerance =
    std::pow(std::numeric_limits<double>::epsilon(), 0.25);
  constexpr int kIterationLimit = 100;

  if (start_id <= 4 || start_id == 6) {
    NumericVector transformed = Rcpp::no_init(n);
    NumericVector transformed2 = Rcpp::no_init(n);
    double sum = 0.0;

    for (R_xlen_t i = 0; i < n; ++i) {
      const double value = start_id == 3
        ? std::log(x[i])
        : x[i];
      transformed[i] = value;
      sum += value;
      if (start_id == 4) {
        transformed2[i] = std::log(x[i]);
      } else if (start_id == 6) {
        transformed[i] = std::log(x[i]);
        transformed2[i] = std::log1p(-x[i]);
      }
    }

    double centered_sum_squares = 0.0;
    const double full_mean = start_id <= 3
      ? sum / static_cast<double>(n)
      : 0.0;
    if (start_id == 1 || start_id == 3) {
      for (R_xlen_t i = 0; i < n; ++i) {
        const double difference = transformed[i] - full_mean;
        centered_sum_squares += difference * difference;
      }
    }

    double sum_transformed = 0.0;
    double sum_transformed2 = 0.0;
    if (start_id == 4 || start_id == 6) {
      for (R_xlen_t i = 0; i < n; ++i) {
        sum_transformed += transformed[i];
        sum_transformed2 += transformed2[i];
      }
    }

    for (R_xlen_t i = 0; i < n; ++i) {
      if ((i & 255) == 0) {
        Rcpp::checkUserInterrupt();
      }

      if (start_id <= 3) {
        const double value = transformed[i];
        double mean_loo = (sum - value) / n_loo;

        if (start_id == 2) {
          // Subtracting the omitted observation from a much larger total can
          // erase the retained observations completely. Recompute this rare
          // pathological case after scaling by the largest retained value.
          if (!std::isfinite(mean_loo) || mean_loo <= 0.0) {
            double maximum_loo = 0.0;
            for (R_xlen_t j = 0; j < n; ++j) {
              if (j != i) {
                maximum_loo = std::max(maximum_loo, x[j]);
              }
            }
            if (!std::isfinite(maximum_loo) || maximum_loo <= 0.0) {
              out(i, 0) = NA_REAL;
              continue;
            }

            double scaled_sum = 0.0;
            for (R_xlen_t j = 0; j < n; ++j) {
              if (j != i) {
                scaled_sum += x[j] / maximum_loo;
              }
            }
            mean_loo = maximum_loo * (scaled_sum / n_loo);
          }

          const double rate = 1.0 / mean_loo;
          out(i, 0) = std::isfinite(rate) && rate > 0.0
            ? rate
            : NA_REAL;
          continue;
        }

        double variance_loo = (
          centered_sum_squares -
          (value - full_mean) * (value - mean_loo)
        ) / n_loo;
        if (variance_loo < 0.0 &&
            variance_loo > -16.0 * std::numeric_limits<double>::epsilon() *
              std::abs(centered_sum_squares / n_loo)) {
          variance_loo = 0.0;
        }
        out(i, 0) = mean_loo;
        out(i, 1) = std::sqrt(variance_loo);
        continue;
      }

      if (start_id == 4) {
        const double mean_loo = (sum - x[i]) / n_loo;
        const double mean_log_loo =
          (sum_transformed2 - transformed2[i]) / n_loo;
        const double s = std::log(mean_loo) - mean_log_loo;
        const std::optional<double> shape =
          fit_gamma_shape_loo(s, tolerance, kIterationLimit);
        if (!shape) {
          out(i, 0) = NA_REAL;
          out(i, 1) = NA_REAL;
          continue;
        }

        out(i, 0) = *shape;
        out(i, 1) = *shape / mean_loo;
        continue;
      }

      const double s = (sum_transformed - transformed[i]) / n_loo;
      const double r = (sum_transformed2 - transformed2[i]) / n_loo;
      const std::optional<BetaShapes> beta_fit =
        fit_beta_shapes_loo(s, r, tolerance, kIterationLimit);
      if (!beta_fit) {
        out(i, 0) = NA_REAL;
        out(i, 1) = NA_REAL;
        continue;
      }

      out(i, 0) = beta_fit->alpha;
      out(i, 1) = beta_fit->exp_b;
    }

    return out;
  }

  NumericVector negative_log_x = Rcpp::no_init(n);
  double sum_negative_log = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    negative_log_x[i] = -std::log(x[i]);
    sum_negative_log += negative_log_x[i];
  }
  const double full_mean_z =
    sum_negative_log / static_cast<double>(n);
  double full_centered_sum_squares = 0.0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double difference = negative_log_x[i] - full_mean_z;
    full_centered_sum_squares += difference * difference;
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    if ((i & 63) == 0) {
      Rcpp::checkUserInterrupt();
    }

    const double mean_z =
      (sum_negative_log - negative_log_x[i]) / n_loo;
    const double centered_sum_squares =
      full_centered_sum_squares -
      (negative_log_x[i] - full_mean_z) *
        (negative_log_x[i] - mean_z);
    const double sigma_init = std::sqrt(centered_sum_squares / n_loo) *
      2.45 / kPi * (1.0 - 1.0 / n_loo);

    const std::optional<double> sigma_opt = fit_weibull_sigma_loo(
      negative_log_x, mean_z, i, n, n_loo, sigma_init, tolerance,
      kIterationLimit
    );
    if (!sigma_opt) {
      out(i, 0) = NA_REAL;
      out(i, 1) = NA_REAL;
      continue;
    }
    const double sigma = *sigma_opt;

    double psi0 = 0.0;
    for (R_xlen_t j = 0; j < n; ++j) {
      if (j != i) {
        psi0 += std::exp(
          -(negative_log_x[j] - mean_z) / sigma
        );
      }
    }
    const double mu = -sigma * std::log(psi0 / n_loo);
    const double s = -mean_z;
    out(i, 0) = 1.0 / sigma;
    out(i, 1) = std::exp(-mu + s);
  }

  return out;
}

// The second term of unbiased cross-validation for a compiled kernel and a
// designated parametric start. Each row of parameters is fitted without x[i].
// Log-sum-exp accumulation avoids overflow in start-density ratios.
//
// [[Rcpp::export]]
double ucv_loo_term_cpp(const NumericVector& x, const double h,
                        const int kernel_id, const int start_id,
                        const NumericMatrix& parameters) {
  const R_xlen_t n = x.size();
  if (n < 2) {
    Rcpp::stop("At least two observations are needed for cross-validation.");
  }
  if (!(h > 0.0)) {
    Rcpp::stop("'h' must be positive.");
  }
  if (kernel_id < 1 || kernel_id > 15) {
    Rcpp::stop("Unknown compiled kernel id.");
  }
  const int parameter_count = start_parameter_count(start_id);
  if (parameter_count == 0) {
    Rcpp::stop("Unknown compiled start id.");
  }
  if (parameters.nrow() != n || parameters.ncol() != parameter_count) {
    Rcpp::stop("The leave-one-out parameter matrix has the wrong dimensions.");
  }

  NumericVector transformed1(n);
  NumericVector log_data(n);
  NumericVector log1m_data(n);
  if (kernel_id == 11) {
    for (R_xlen_t i = 0; i < n; ++i) {
      transformed1[i] = R::qnorm5(x[i], 0.0, 1.0, true, false);
    }
  }
  if (start_id >= 3 || kernel_id >= 12) {
    for (R_xlen_t i = 0; i < n; ++i) {
      log_data[i] = std::log(x[i]);
    }
  }
  if (start_id == 6 || kernel_id >= 14) {
    for (R_xlen_t i = 0; i < n; ++i) {
      log1m_data[i] = std::log1p(-x[i]);
    }
  }

  double mean_term = 0.0;
  const double log_average_scale =
    -std::log(h) - std::log(static_cast<double>(n - 1));

  for (R_xlen_t i = 0; i < n; ++i) {
    if ((i & 63) == 0) {
      Rcpp::checkUserInterrupt();
    }

    const double start_constant =
      start_log_constant(start_id, parameters, i);
    const double log_start_i = start_log_density(
      x[i], log_data[i], log1m_data[i], start_id, parameters, i,
      start_constant
    );
    const PairKernelConstants constants = pair_kernel_constants(x[i], h, kernel_id);
    double maximum = R_NegInf;
    double scaled_sum = 0.0;

    for (R_xlen_t j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }

      // Chen's asymmetric kernels take their shape from the FIRST argument
      // (the evaluation point) and are evaluated at the SECOND. The leave-one-
      // out density f_hat_{-i}(x_i) is a sum over j != i of kernel(x_i, x_j),
      // shape from the eval point x[i], density at the data point x[j] -- the
      // same convention as term1 and the estimator. (Symmetric and gcopula
      // kernels are symmetric in the two arguments, so the order is moot.)
      const double qy = kernel_id == 11 ? transformed1[i] : 0.0;
      const double qx = kernel_id == 11 ? transformed1[j] : 0.0;
      const double log_x = kernel_id >= 12 ? log_data[j] : 0.0;
      const double log1m_x = kernel_id >= 14 ? log1m_data[j] : 0.0;
      const double log_kernel = pair_log_kernel(
        x[i], x[j], h, kernel_id, qy, qx, log_x, log1m_x, constants
      );
      const double log_start_j = start_log_density(
        x[j], log_data[j], log1m_data[j], start_id, parameters, i,
        start_constant
      );
      const double log_term = log_kernel + log_start_i - log_start_j;

      if (std::isnan(log_term)) {
        return NA_REAL;
      }
      if (log_term == R_PosInf) {
        return R_PosInf;
      }
      if (log_term == R_NegInf) {
        continue;
      }

      if (log_term > maximum) {
        scaled_sum = maximum == R_NegInf
          ? 1.0
          : scaled_sum * std::exp(maximum - log_term) + 1.0;
        maximum = log_term;
      } else {
        scaled_sum += std::exp(log_term - maximum);
      }
    }

    if (maximum != R_NegInf) {
      mean_term += std::exp(maximum + log_average_scale) * scaled_sum;
    }
  }

  return 2.0 * mean_term / static_cast<double>(n);
}
