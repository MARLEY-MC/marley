#include <cmath>
#include <iomanip>
#include <iostream>

#include "TMarleyIntegrator.hh"

double f(double x) {
  return std::exp(-5.1*x) + std::sin(0.1*x) + 3.5*x - 2.;
}

int main() {
  std::cout << std::setprecision(16) << std::scientific;
  TMarleyIntegrator integrator;
  for (size_t j = 1; j < 1e5; ++j)
    std::cout << "j = " << j << ", answer = "
      << integrator.num_integrate(&f, -5., 5.1) << std::endl;
  return 0;
}
