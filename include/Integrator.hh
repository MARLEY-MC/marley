#pragma once
#include <functional>
#include <vector>

#include "marley_utils.hh"

namespace marley {

  class Integrator {
    public:
      inline Integrator(size_t n = N_DEFAULT) {
        // TODO: add error check for when n == 0
        // TODO: add error check for when n is ridiculously large
        N = n;
        // Precompute the 2N - 1 weights for speed
        size_t two_N = 2*N;
        double npi_over_two = -marley_utils::half_pi;
        for (size_t k = 0; k < two_N; ++k) {
          npi_over_two += marley_utils::half_pi;
          weights.push_back(std::cos(npi_over_two / N));
        }
      }
      // Numerically integrate a given function f (that takes a
      // double argument to integrate over and returns a double)
      // over the interval [a,b] using Clenshaw-Curtis quadrature
      // at 2N sampling points.
      // (see http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature)
      double num_integrate(const std::function<double(double)> &f,
        double a, double b) const;
  
    private:
      // Use 2*N sampling points to perform numerical integration
      size_t N;
      // Vector of precomputed weights to use to speed up integration
      std::vector<double> weights;
      // If the user doesn't specify a number of sampling points to use for an
      // integration, then use this value
      static constexpr size_t N_DEFAULT = 20;
  };

}
