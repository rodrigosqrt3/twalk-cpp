// Example: sample a d-dimensional standard normal with the t-walk and write
// the chain to a CSV file.
//
//   g++ -O2 -std=c++17 -Iinclude examples/normal.cpp -o normal -lgsl -lgslcblas
//   ./normal 10 100000 chain.csv

#include <twalk/gsl_rng.hpp>
#include <twalk/twalk.hpp>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  const std::size_t d = argc > 1 ? std::atoi(argv[1]) : 10;
  const long n_iter = argc > 2 ? std::atol(argv[2]) : 100000;
  const char* path = argc > 3 ? argv[3] : "chain.csv";

  // The objective: log density of N(0, I), up to a constant.
  auto log_post = [](const std::vector<double>& x) {
    double ss = 0.0;
    for (double xi : x) ss += xi * xi;
    return -0.5 * ss;
  };

  twalk::GslRng rng(2026);
  std::vector<double> x0(d), xp0(d);
  for (std::size_t j = 0; j < d; ++j) {
    x0[j] = rng.norm();
    xp0[j] = rng.norm();
  }

  twalk::State state = twalk::initial_state(log_post, x0, xp0);
  twalk::Output out = twalk::run(state, log_post, rng, n_iter);

  std::vector<double> mean(d, 0.0), var(d, 0.0);
  const std::size_t m = out.n_stored();
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < d; ++j) mean[j] += out.samples[i * d + j] / m;
  for (std::size_t i = 0; i < m; ++i)
    for (std::size_t j = 0; j < d; ++j) {
      double z = out.samples[i * d + j] - mean[j];
      var[j] += z * z / (m - 1);
    }

  std::cout << "dimension " << d << ", " << n_iter << " iterations\n"
            << "acceptance rate: " << out.acceptance_rate() << '\n'
            << "x1: mean " << mean[0] << ", variance " << var[0]
            << " (target 0 and 1)\n";

  std::ofstream file(path);
  twalk::write_csv(file, out);
  std::cout << "chain written to " << path << '\n';
  return 0;
}
