// Runs twalk.hpp with R's random number generator, so that its chains can be
// compared draw by draw with those of the R package Rtwalk.

// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <twalk/twalk.hpp>

struct RRng {
  double unif() { return unif_rand(); }
  double norm() { return norm_rand(); }
};

// [[Rcpp::export]]
Rcpp::List twalk_hpp_run(Rcpp::Function log_posterior,
                         std::vector<double> x0, std::vector<double> xp0,
                         long n_iter, long thin = 1) {
  auto log_post = [&](const std::vector<double>& x) {
    return Rcpp::as<double>(log_posterior(Rcpp::wrap(x)));
  };
  twalk::State s = twalk::initial_state(log_post, x0, xp0);

  Rcpp::RNGScope scope;
  RRng rng;
  twalk::Output out = twalk::run(s, log_post, rng, n_iter, twalk::Settings(), thin);

  const int m = out.n_stored(), n = out.n;
  Rcpp::NumericMatrix samples(m, n), companion(m, n);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) {
      samples(i, j) = out.samples[i * n + j];
      companion(i, j) = out.companion[i * n + j];
    }
  Rcpp::IntegerVector proposed(4), accepted(4);
  for (int k = 0; k < 4; ++k) {
    proposed[k] = out.proposed[k];
    accepted[k] = out.accepted[k];
  }
  return Rcpp::List::create(
    Rcpp::_["samples"] = samples,
    Rcpp::_["companion_samples"] = companion,
    Rcpp::_["log_posterior"] = out.log_post,
    Rcpp::_["companion_log_posterior"] = out.companion_log_post,
    Rcpp::_["proposed"] = proposed,
    Rcpp::_["accepted"] = accepted,
    Rcpp::_["n_evals"] = out.n_evals);
}
