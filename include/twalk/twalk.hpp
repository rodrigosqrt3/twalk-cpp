// twalk.hpp -- header-only C++ implementation of the t-walk
// (Christen & Fox, 2010, "A general purpose sampling algorithm for
// continuous distributions (the t-walk)", Bayesian Analysis 5(2), 263-282).
//
// This is the kernel of the R package Rtwalk (src/twalk_engine.cpp) without
// any dependency on R. The random number generator is a template parameter:
// any type with member functions
//
//     double unif();   // uniform on (0, 1)
//     double norm();   // standard normal
//
// can be used. twalk/gsl_rng.hpp provides one based on the GNU Scientific
// Library. With an adapter to R's generator, the chains are identical, draw by
// draw, to those of Rtwalk (see tests/).
//
// The log posterior is any callable `double f(const std::vector<double>&)`.
// It returns -infinity outside the support; NaN or +infinity are errors.

#ifndef TWALK_TWALK_HPP
#define TWALK_TWALK_HPP

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace twalk {

enum Kernel { TRAVERSE = 0, WALK = 1, BLOW = 2, HOP = 3 };

// Constants of the algorithm, with the defaults of Christen & Fox (2010,
// Section 2.5).
struct Settings {
  double n1 = 4.0;   // expected number of coordinates moved per iteration
  double aw = 1.5;   // walk move parameter
  double at = 6.0;   // traverse move parameter
  double kernel_probs[4] = {0.4918, 0.4918, 0.0082, 0.0082};
};

// Current state of the chain: the two points and their log posteriors.
struct State {
  std::vector<double> x, xp;
  double lp_x, lp_xp;
};

// What happened in one transition.
struct StepInfo {
  Kernel kernel;
  bool moved_xp;   // whether xp (rather than x) was the point proposed
  bool accepted;
  bool changed;    // accepted and the point actually moved
  bool evaluated;  // whether the log posterior was evaluated
};

// Output of run(). Matrices are stored by rows (iteration-major).
struct Output {
  std::size_t n = 0;
  std::vector<double> samples, companion;
  std::vector<double> log_post, companion_log_post;
  long proposed[4] = {0, 0, 0, 0};
  long accepted[4] = {0, 0, 0, 0};
  long moved[4] = {0, 0, 0, 0};
  long n_evals = 0;

  std::size_t n_stored() const { return n == 0 ? 0 : samples.size() / n; }
  double acceptance_rate() const {
    long p = 0, a = 0;
    for (int k = 0; k < 4; ++k) { p += proposed[k]; a += accepted[k]; }
    return p == 0 ? 0.0 : static_cast<double>(a) / p;
  }
};

namespace detail {

template <class Rng>
double simulate_beta(Rng& rng, double at) {
  if (rng.unif() < (at - 1.0) / (2.0 * at)) {
    return std::exp(1.0 / (at + 1.0) * std::log(rng.unif()));
  }
  return std::exp(1.0 / (1.0 - at) * std::log(rng.unif()));
}

// -log g(h | a, b) for the Blow (centre b) and Hop (centre a) proposals.
// Only called with n_phi > 0.
inline double neg_log_density(const std::vector<double>& h,
                              const std::vector<double>& a,
                              const std::vector<double>& b,
                              const std::vector<double>& phi,
                              int n_phi, bool hop) {
  double sigma = 0.0;
  for (std::size_t j = 0; j < h.size(); ++j) {
    double s = phi[j] * std::fabs(b[j] - a[j]);
    if (s > sigma) sigma = s;
  }
  if (hop) sigma = sigma / 3.0;
  if (sigma < DBL_EPSILON) sigma = DBL_EPSILON;

  const std::vector<double>& centre = hop ? a : b;
  long double ss = 0.0;
  for (std::size_t j = 0; j < h.size(); ++j) {
    if (phi[j] != 0.0) {
      double d = h[j] - centre[j];
      ss += d * d;
    }
  }

  const double log_2pi = 1.8378770664093454836;
  return (n_phi / 2.0) * log_2pi + n_phi * std::log(sigma) +
    0.5 * static_cast<double>(ss) / (sigma * sigma);
}

template <class LogPost>
double checked_log_post(LogPost& log_post, const std::vector<double>& h) {
  double lp = log_post(h);
  if (std::isnan(lp) || lp == std::numeric_limits<double>::infinity()) {
    throw std::domain_error(
      "log posterior must return a finite value or -inf outside the support");
  }
  return lp;
}

}  // namespace detail

// One t-walk transition. Random numbers are drawn in the same order as in
// Rtwalk, so the same generator gives the same chain.
template <class Rng, class LogPost>
StepInfo step(State& s, LogPost& log_post, Rng& rng,
              const Settings& set = Settings()) {
  const int n = static_cast<int>(s.x.size());
  if (n == 0 || s.xp.size() != s.x.size()) {
    throw std::invalid_argument("x and xp must have the same, positive, length");
  }
  const double p_phi = std::min(static_cast<double>(n), set.n1) / n;

  // Cumulative kernel probabilities, accumulated as R's cumsum() does.
  double cum[3];
  long double acc = 0.0;
  for (int k = 0; k < 3; ++k) {
    acc += set.kernel_probs[k];
    cum[k] = static_cast<double>(acc);
  }

  StepInfo info;
  double u_kernel = rng.unif();
  int kernel = HOP;
  for (int k = 0; k < 3; ++k) {
    if (u_kernel < cum[k]) {
      kernel = k;
      break;
    }
  }
  info.kernel = static_cast<Kernel>(kernel);
  info.moved_xp = rng.unif() < 0.5;

  // `a` is the point being moved and `b` the other point.
  std::vector<double>& a = info.moved_xp ? s.xp : s.x;
  const std::vector<double>& b = info.moved_xp ? s.x : s.xp;
  double& lp_a = info.moved_xp ? s.lp_xp : s.lp_x;

  std::vector<double> h(n), phi(n);
  int n_phi = 0;
  double beta = 0.0;

  if (kernel == TRAVERSE) {
    beta = detail::simulate_beta(rng, set.at);
    for (int j = 0; j < n; ++j) {
      phi[j] = (rng.unif() < p_phi) ? 1.0 : 0.0;
      n_phi += (phi[j] != 0.0);
      h[j] = (phi[j] != 0.0) ? b[j] + beta * (b[j] - a[j]) : a[j];
    }
  } else if (kernel == WALK) {
    std::vector<double> u(n);
    for (int j = 0; j < n; ++j) u[j] = rng.unif();
    for (int j = 0; j < n; ++j) {
      phi[j] = (rng.unif() < p_phi) ? 1.0 : 0.0;
      n_phi += (phi[j] != 0.0);
    }
    const double aw = set.aw;
    for (int j = 0; j < n; ++j) {
      double z = (aw / (1.0 + aw)) * (aw * (u[j] * u[j]) + 2.0 * u[j] - 1.0);
      z = z * phi[j];
      h[j] = a[j] + (a[j] - b[j]) * z;
    }
  } else {
    for (int j = 0; j < n; ++j) {
      phi[j] = (rng.unif() < p_phi) ? 1.0 : 0.0;
      n_phi += (phi[j] != 0.0);
    }
    double sigma = 0.0;
    for (int j = 0; j < n; ++j) {
      double sj = phi[j] * std::fabs(b[j] - a[j]);
      if (sj > sigma) sigma = sj;
    }
    if (kernel == HOP) sigma = sigma / 3.0;
    if (sigma < DBL_EPSILON) sigma = DBL_EPSILON;

    std::vector<double> noise(n);
    for (int j = 0; j < n; ++j) noise[j] = rng.norm() * phi[j];

    for (int j = 0; j < n; ++j) {
      if (kernel == BLOW) {
        h[j] = b[j] * phi[j] + sigma * noise[j] + a[j] * (1.0 - phi[j]);
      } else {
        h[j] = (phi[j] != 0.0) ? a[j] + sigma * noise[j] : a[j];
      }
    }
  }

  double alpha = 0.0;
  double lp_h = lp_a;
  info.evaluated = false;

  if (n_phi == 0) {
    // Identity proposal: accepted with probability one, no evaluation needed.
    alpha = 1.0;
  } else {
    bool degenerate = false;
    if (kernel != TRAVERSE) {
      degenerate = std::equal(h.begin(), h.end(), b.begin());
    }

    if (!degenerate) {
      lp_h = detail::checked_log_post(log_post, h);
      info.evaluated = true;

      if (std::isfinite(lp_h)) {
        double log_alpha = lp_h - lp_a;
        if (kernel == TRAVERSE) {
          log_alpha = log_alpha + (n_phi - 2) * std::log(beta);
        } else if (kernel == BLOW || kernel == HOP) {
          bool hop = (kernel == HOP);
          double w1 = detail::neg_log_density(h, a, b, phi, n_phi, hop);
          double w2 = detail::neg_log_density(a, h, b, phi, n_phi, hop);
          log_alpha = log_alpha + (w1 - w2);
        }
        alpha = (log_alpha >= 0.0) ? 1.0 : std::exp(log_alpha);
      }
    }
  }

  info.accepted = rng.unif() < alpha;
  info.changed = false;
  if (info.accepted) {
    info.changed = !std::equal(h.begin(), h.end(), a.begin());
    a = h;
    lp_a = lp_h;
  }
  return info;
}

// Starting state, evaluating the log posterior at both points. Both must be
// in the support and must differ in every coordinate.
template <class LogPost>
State initial_state(LogPost& log_post, const std::vector<double>& x0,
                    const std::vector<double>& xp0) {
  if (x0.empty() || x0.size() != xp0.size()) {
    throw std::invalid_argument("x0 and xp0 must have the same, positive, length");
  }
  for (std::size_t j = 0; j < x0.size(); ++j) {
    if (x0[j] == xp0[j]) {
      throw std::invalid_argument("x0 and xp0 must differ in every coordinate");
    }
  }
  State s{x0, xp0, detail::checked_log_post(log_post, x0),
          detail::checked_log_post(log_post, xp0)};
  if (!std::isfinite(s.lp_x) || !std::isfinite(s.lp_xp)) {
    throw std::invalid_argument("x0 and xp0 must be in the support");
  }
  return s;
}

// Run n_iter transitions from `s` (which is updated, so a later call
// continues the chain), storing every `thin`-th state.
template <class Rng, class LogPost>
Output run(State& s, LogPost& log_post, Rng& rng, long n_iter,
           const Settings& set = Settings(), long thin = 1,
           bool store_companion = true) {
  if (n_iter < 0 || thin < 1) {
    throw std::invalid_argument("n_iter must be >= 0 and thin >= 1");
  }
  Output out;
  out.n = s.x.size();
  const std::size_t n_store = static_cast<std::size_t>(n_iter / thin);
  out.samples.reserve(n_store * out.n);
  out.log_post.reserve(n_store);
  if (store_companion) {
    out.companion.reserve(n_store * out.n);
    out.companion_log_post.reserve(n_store);
  }

  for (long i = 1; i <= n_iter; ++i) {
    StepInfo info = step(s, log_post, rng, set);
    out.proposed[info.kernel] += 1;
    if (info.accepted) out.accepted[info.kernel] += 1;
    if (info.changed) out.moved[info.kernel] += 1;
    if (info.evaluated) out.n_evals += 1;

    if (i % thin == 0) {
      out.samples.insert(out.samples.end(), s.x.begin(), s.x.end());
      out.log_post.push_back(s.lp_x);
      if (store_companion) {
        out.companion.insert(out.companion.end(), s.xp.begin(), s.xp.end());
        out.companion_log_post.push_back(s.lp_xp);
      }
    }
  }
  return out;
}

// Write the primary trajectory as CSV: columns x1..xn and log_post.
inline void write_csv(std::ostream& os, const Output& out) {
  os.precision(17);
  for (std::size_t j = 0; j < out.n; ++j) os << 'x' << (j + 1) << ',';
  os << "log_post\n";
  for (std::size_t i = 0; i < out.n_stored(); ++i) {
    for (std::size_t j = 0; j < out.n; ++j) os << out.samples[i * out.n + j] << ',';
    os << out.log_post[i] << '\n';
  }
}

}  // namespace twalk

#endif  // TWALK_TWALK_HPP
