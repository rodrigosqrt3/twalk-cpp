# twalk-cpp

A header-only C++ implementation of the t-walk (Christen & Fox, 2010), with
no dependency on R. The objective function is written in C++ and the program
is compiled directly with g++.

This is the same kernel as the R package
[Rtwalk](https://cran.r-project.org/package=Rtwalk) (version 2.2.0). When it
is driven by R's random number generator, it produces exactly the same chains
as Rtwalk, draw by draw (see [Validation](#validation)). For standalone use it
draws random numbers from the GNU Scientific Library (GSL).

## Files

| File | Contents |
|---|---|
| `include/twalk/twalk.hpp` | the t-walk; no dependencies beyond the C++ standard library |
| `include/twalk/gsl_rng.hpp` | a random number generator based on GSL |
| `examples/normal.cpp` | a complete example |
| `tests/` | comparison with Rtwalk |

## Quick start

```cpp
#include <twalk/gsl_rng.hpp>
#include <twalk/twalk.hpp>
#include <fstream>

int main() {
  // log of the objective, up to a constant; -infinity outside the support
  auto log_post = [](const std::vector<double>& x) {
    double ss = 0.0;
    for (double xi : x) ss += xi * xi;
    return -0.5 * ss;
  };

  twalk::GslRng rng(2026);                      // seed
  std::vector<double> x0(5, -1.0), xp0(5, 1.0);  // two starting points
  twalk::State s = twalk::initial_state(log_post, x0, xp0);
  twalk::Output out = twalk::run(s, log_post, rng, 100000);

  std::ofstream f("chain.csv");
  twalk::write_csv(f, out);                     // columns x1..xn, log_post
}
```

Compile with

```sh
g++ -O2 -std=c++17 -Iinclude my_model.cpp -o my_model -lgsl -lgslcblas
```

or run `make` to build `examples/normal.cpp`.

## Interface

- `twalk::initial_state(log_post, x0, xp0)` evaluates the objective at the
  two starting points. They must be in the support and must differ in every
  coordinate.
- `twalk::run(state, log_post, rng, n_iter, settings, thin, store_companion)`
  runs `n_iter` iterations and stores every `thin`-th state. `state` is
  updated in place, so a second call continues the chain.
- `twalk::step(state, log_post, rng, settings)` performs a single transition
  and reports the kernel used and whether it was accepted. It can be combined
  with other updates, for example Gibbs steps, as in Section 5 of the paper.
- `twalk::Settings` holds the constants of the algorithm (`n1`, `aw`, `at`
  and the probabilities of the four moves). The defaults are those of the
  paper (Section 2.5).
- `Output` holds the primary trajectory `samples` (stored by rows), the
  companion trajectory, both log objectives, and the number of proposals,
  acceptances and moves for each kernel.

The random number generator is a template parameter: any type with
`double unif()` (uniform on (0, 1)) and `double norm()` (standard normal)
works. `twalk::GslRng` uses the Mersenne Twister by default. Any other GSL
generator can be passed to its constructor, for example
`twalk::GslRng rng(1, gsl_rng_ranlxd2)`.

## Validation

`tests/compare_with_rtwalk.R` compiles `twalk.hpp` with an adapter to R's
generator and compares it with Rtwalk on five targets: normal in 1 and 10
dimensions, Rosenbrock, a target with a restricted support, and a Student t
in 50 dimensions. The samples, the companion samples and both log objectives
are identical to those of Rtwalk's C++ and R engines.

```sh
make check   # or: Rscript tests/compare_with_rtwalk.R
```

## Requirements

A C++17 compiler and GSL (`libgsl-dev` on Debian/Ubuntu, `gsl` on Homebrew;
on Windows, Rtools 4.3 includes it). On macOS, point the compiler to
Homebrew with `make CPPFLAGS=-I$(brew --prefix)/include
LDFLAGS=-L$(brew --prefix)/lib`. The validation also requires R with the
packages Rtwalk and Rcpp.

## Reference

Christen, J. A. and Fox, C. (2010). A general purpose sampling algorithm for
continuous distributions (the t-walk). *Bayesian Analysis*, 5(2), 263–282.
