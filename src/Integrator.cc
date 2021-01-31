/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#include <cmath>

#include "marley/marley_utils.hh"
#include "marley/Integrator.hh"

marley::Integrator::Integrator(size_t num) : N_(num), weights_(num + 1, 0.),
  offsets_(num - 1)
{
  /// @todo add error check for when num == 0 or too small
  /// @todo add error check for when num is ridiculously large

  // Precompute the N_ - 1 offsets for speed
  double arg = 0.;
  for (size_t n = 0; n < N_ - 1; ++n) {
    arg += marley_utils::half_pi / N_;
    offsets_[n] = std::cos(arg);
  }

  // Also precompute the N_ + 1 weights
  for (size_t n = 0; n <= N_; ++n) {
    for (size_t k = 0; k <= N_; ++k) {
      double weight_piece = std::cos(n * k * marley_utils::pi / N_)
        / (1. - std::pow(2*k, 2));
      if (k != 0 && k != N_) weight_piece *= 2.;
      weights_[n] += weight_piece;
    }
    weights_[n] /= N_;
  }
}

double marley::Integrator::num_integrate(const std::function<double(double)>& f,
  double a, double b) const
{
  double A = (b - a) / 2.;
  double B = (b + a) / 2.;
  double C = (f(a) + f(b)) / 2.;

  double integral = weights_[0] * C; // n = 0 term
  integral += weights_[N_] * f(B); // n = N_ term

  // n = 1 to n = N_ - 1 terms
  for (size_t n = 1; n < N_; ++n) {
    double epoint = A * offsets_[n - 1];
    integral += weights_[n] * (f(B + epoint) + f(B - epoint));
  }

  return A * integral;
}
