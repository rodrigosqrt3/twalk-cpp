// gsl_rng.hpp -- random number generator for twalk.hpp based on the GNU
// Scientific Library. Link with -lgsl -lgslcblas.

#ifndef TWALK_GSL_RNG_HPP
#define TWALK_GSL_RNG_HPP

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <stdexcept>

namespace twalk {

class GslRng {
 public:
  // The default generator is the Mersenne Twister (gsl_rng_mt19937).
  explicit GslRng(unsigned long seed,
                  const gsl_rng_type* type = gsl_rng_mt19937)
      : r_(gsl_rng_alloc(type)) {
    if (r_ == nullptr) throw std::bad_alloc();
    gsl_rng_set(r_, seed);
  }
  ~GslRng() { gsl_rng_free(r_); }
  GslRng(const GslRng&) = delete;
  GslRng& operator=(const GslRng&) = delete;

  // Uniform on the open interval (0, 1): the traverse move takes log(u).
  double unif() { return gsl_rng_uniform_pos(r_); }
  double norm() { return gsl_ran_gaussian_ziggurat(r_, 1.0); }

  gsl_rng* get() { return r_; }

 private:
  gsl_rng* r_;
};

}  // namespace twalk

#endif  // TWALK_GSL_RNG_HPP
