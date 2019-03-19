#pragma once
#include <complex>

std::complex<double> coulomb_H_plus(int l, double eta, double rho);

void marley_gsl_error_handler(const char* reason, const char* file, int line,
  int gsl_errno);
