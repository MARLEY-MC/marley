#include <cmath>

#include "marley_utils.hh"
#include "Integrator.hh"

marley::Integrator::Integrator(size_t n) : N_(n) {
  /// @todo add error check for when n == 0
  /// @todo add error check for when n is ridiculously large

  // Precompute the 2N - 1 weights for speed
  size_t two_N = 2*N_;
  double npi_over_two = -marley_utils::half_pi;
  for (size_t k = 0; k < two_N; ++k) {
    npi_over_two += marley_utils::half_pi;
    weights_.push_back(std::cos(npi_over_two / N_));
  }
}

double marley::Integrator::num_integrate(const std::function<double(double)>& f,
  double a, double b) const
{
  size_t two_N = 2*N_;
  double A = (b - a)/2;
  double B = (b + a)/2;
  double C = (f(b) + f(a))/(two_N);
  double D = f(B)/N_;

  double integral = 0.;
  int sign = 1;
  for(size_t k = 0; k <= N_; sign *= -1, ++k) {

    double term = C + sign*D;
    for (size_t n = 1; n < N_; ++n) {
      double epoint = A * weights_[n];

      double nk_weight;
      size_t nk_mod = n*k % two_N;
      if (nk_mod >= N_) nk_weight = -weights_[2*(nk_mod % N_)];
      else nk_weight = weights_[2*nk_mod];

      term += (f(B + epoint) + f(B - epoint)) * nk_weight/N_;
    }

    if (k != 0 && k != N_) term *= 2;

    term /= 1 - 4*std::pow(k, 2);

    integral += term;
  }

  return A*integral;
}
