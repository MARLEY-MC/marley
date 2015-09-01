#include <cmath>

#include "marley_utils.hh"
#include "TMarleyIntegrator.hh"

// Numerically integrate a given function f (that takes a double argument to
// integrate over and returns a double) over the interval [a,b] using
// Clenshaw-Curtis quadrature at 2N sampling points.  (see
// http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature)
double TMarleyIntegrator::num_integrate(const std::function<double(double)> &f,
  double a, double b, int N)
{
  double twoN = 2*N;
  double A = (b - a)/2;
  double B = (b + a)/2;
  double C = (f(b) + f(a))/twoN;
  double D = f(B)/N;

  std::vector<double> eval_points;
  double npi = 0;
  for (int n = 1; n < N; ++n) {
    npi += marley_utils::pi;
    eval_points.push_back(A*std::cos(npi/twoN));
  }

  double integral = 0;

  for(int k = 0; k <= N; ++k) {

    npi = 0;
    int sign = 1;
    if (k % 2 == 1) sign = -1;
    double term = C + sign*D;
    for (int n = 1; n < N; ++n) {
      npi += marley_utils::pi;
      double epoint = eval_points[n-1];

      term += (f(B + epoint) + f(B - epoint))
        * std::cos(npi*k/N)/N;
    }

    if (k != 0 && k != N) term *= 2;

    term /= 1 - 4*std::pow(k, 2);

    integral += term;
  }

  return A*integral;
}
