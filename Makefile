# Build the examples:  make
# Check against Rtwalk: make check   (needs R with Rtwalk and Rcpp)

CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra
LDLIBS    = -lgsl -lgslcblas -lm

all: build/normal

build/normal: examples/normal.cpp include/twalk/twalk.hpp include/twalk/gsl_rng.hpp
	mkdir -p build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -Iinclude $< -o $@ $(LDFLAGS) $(LDLIBS)

check:
	Rscript tests/compare_with_rtwalk.R

clean:
	rm -rf build

.PHONY: all check clean
