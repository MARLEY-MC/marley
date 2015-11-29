#include <cmath>

#include "marley_utils.hh"
#include "TMarleyIntegrator.hh"

// Numerically integrate a given function f (that takes a double argument to
// integrate over and returns a double) over the interval [a,b] using
// Clenshaw-Curtis quadrature at 2N sampling points.  (see
// http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature)
double TMarleyIntegrator::num_integrate(const std::function<double(double)> &f,
  double a, double b) const
{
  int two_N = 2*N;
  double A = (b - a)/2;
  double B = (b + a)/2;
  double C = (f(b) + f(a))/(two_N);
  double D = f(B)/N;

  double integral = 0.;
  int sign = 1;
  for(size_t k = 0; k <= N; sign *= -1, ++k) {

    double term = C + sign*D;
    for (size_t n = 1; n < N; ++n) {
      double epoint = A * weights[n];

      double nk_weight;
      size_t nk_mod = n*k % two_N;
      if (nk_mod >= N) nk_weight = -weights[2*(nk_mod % N)];
      else nk_weight = weights[2*nk_mod];
      
      term += (f(B + epoint) + f(B - epoint)) * nk_weight/N;
    }

    if (k != 0 && k != N) term *= 2;

    term /= 1 - 4*std::pow(k, 2);

    integral += term;
  }

  return A*integral;
}
