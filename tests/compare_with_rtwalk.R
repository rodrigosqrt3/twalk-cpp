# Checks that twalk.hpp, driven by R's random number generator, produces
# exactly the same chains as the R package Rtwalk (C++ and R engines).
#
# Run from the twalk-cpp folder:  Rscript tests/compare_with_rtwalk.R

library(Rtwalk)
Sys.setenv(PKG_CPPFLAGS = paste0("-I", normalizePath("include", "/")))
Rcpp::sourceCpp("tests/rtwalk_adapter.cpp")

targets <- list(
  normal_d1 = list(f = function(x) -sum(x^2) / 2, d = 1),
  normal_d10 = list(f = function(x) -sum(x^2) / 2, d = 10),
  rosenbrock_d2 = list(
    f = function(x) -(100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2) / 20, d = 2),
  # support restricted to the positive orthant: tests rejections at -Inf
  gamma_d5 = list(
    f = function(x) if (any(x <= 0)) -Inf else sum(2 * log(x) - x), d = 5),
  t3_d50 = list(f = function(x) -(3 + 50) / 2 * log1p(sum(x^2) / 3), d = 50)
)

n_iter <- 20000
failures <- 0
for (nm in names(targets)) {
  tg <- targets[[nm]]
  set.seed(1)
  x0 <- runif(tg$d, 0.5, 2)
  xp0 <- runif(tg$d, 0.5, 2)

  for (engine in c("cpp", "R")) {
    set.seed(42)
    ref <- twalk(tg$f, n_iter, x0 = x0, xp0 = xp0, show_progress = FALSE,
                 control = twalk_control(engine = engine))
    set.seed(42)
    hpp <- twalk_hpp_run(tg$f, x0, xp0, n_iter)

    same <- identical(unname(ref$samples), hpp$samples) &&
      identical(unname(ref$companion_samples), hpp$companion_samples) &&
      identical(as.numeric(ref$log_posterior), hpp$log_posterior) &&
      identical(as.numeric(ref$companion_log_posterior),
                hpp$companion_log_posterior)
    cat(sprintf("%-14s engine %-3s: %s (acceptance %.3f)\n", nm, engine,
                if (same) "identical" else "DIFFERENT",
                sum(hpp$accepted) / n_iter))
    if (!same) failures <- failures + 1
  }
}

if (failures > 0) stop(failures, " comparison(s) failed")
cat("\nAll chains identical to Rtwalk.\n")
