#pragma once
#include <functional>
#include <vector>

class TMarleyIntegrator {
  public:
    // Numerically integrate a given function f (that takes a
    // double argument to integrate over and returns a double)
    // over the interval [a,b] using Clenshaw-Curtis quadrature
    // at 2N sampling points.
    // (see http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature)
    double num_integrate(const std::function<double(double)> &f,
      double a, double b, int N = N_DEFAULT);

  private:
    // Vector of precomputed weights to use to speed up integration
    // TODO: implement this
    //const std::vector<double> weights;
    // If the user doesn't specify a number of sampling points to use for an
    // integration, then use this value
    static constexpr int N_DEFAULT = 100;
};
