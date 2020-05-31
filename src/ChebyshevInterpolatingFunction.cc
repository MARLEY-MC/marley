/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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

// Standard library includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

// FFTPACK4 includes
#include "fftpack4/fftpack4.h"
#include "fftpack4/fftpack4_precision.h"

// MARLEY includes
#include "marley/ChebyshevInterpolatingFunction.hh"

namespace {
  constexpr double MY_EPSILON = std::numeric_limits<double>::epsilon();
}

marley::ChebyshevInterpolatingFunction::ChebyshevInterpolatingFunction(
  const std::function<double(double)>& func, double x_min, double x_max,
  size_t N)
  : x_min_( x_min ), x_max_( x_max )
{
  bool ok;

  if ( N != 0 ) {
    ok = true;
    N_ = N;
  }
  else {
    ok = false;
    N_ = 1;
  }

  // Adaptively find a good grid size
  // @todo add more explanation
  do {

    if ( !ok ) N_ *= 2;

    //std::cout << "N = " << N_ << '\n';

    // @todo optimize this
    Xs_.clear();
    Fs_.clear();

    for ( size_t j = 0; j <= N_; ++j ) {
      double x = chebyshev_point( j );
      double f = func(x);
      Xs_.push_back( x );
      Fs_.push_back( f );
    }

    chebyshev_coeffs_ = Fs_;
    wsave_ = std::vector<double>(3*Fs_.size() + 15, 0.);
    ifac_ = std::vector<int>(Fs_.size() / 2, 0);

    int my_size = Fs_.size();
    costi( &my_size, wsave_.data(), ifac_.data() );
    cost( &my_size, chebyshev_coeffs_.data(), wsave_.data(), ifac_.data() );

    double biggest_coeff = *std::max_element(chebyshev_coeffs_.cbegin(),
      chebyshev_coeffs_.cend(), [](double left, double right) -> double
      { return std::abs(left) < std::abs(right); });

    double upper_limit = 2. * std::abs( biggest_coeff ) * MY_EPSILON;

    double last_coeff_mag = std::abs( chebyshev_coeffs_.back() );
    double next_to_last_coeff_mag = std::abs(
      chebyshev_coeffs_.at( chebyshev_coeffs_.size() - 2 ));
    if ( last_coeff_mag < upper_limit && next_to_last_coeff_mag < upper_limit )
    {
      ok = true;
    }

    if ( N_ >= N_MAX_ ) {
      /// @todo PRINT WARNING MESSAGE
      ok = true;
    }

  } while ( !ok );

  // Normalize the Chebyshev coefficients by dividing by N
  for (double& d : chebyshev_coeffs_) d /= N_;

  compute_integral();
}

void marley::ChebyshevInterpolatingFunction::compute_integral() {
  // Compute the integral over [x_min_, x_max_] via Clenshaw-Curtis
  // quadrature
  integral_ = 0.;
  for (size_t k = 0; k <= N_; k += 2) {
    // Only even Chebyshev polynomials contribute to this integral
    double term = chebyshev_coeffs_.at( k ) * 2.0 / ( 1. - std::pow(k, 2) );
    if ( k == 0 || k == N_ ) term /= 2.;
    integral_ += term;
  }

  // Extra factor to account for integration over [x_min_, x_max_] instead
  // of the standard interval [-1, 1] for the Chebyshev polynomials
  integral_ *= (x_max_ - x_min_) / 2.;
}


marley::ChebyshevInterpolatingFunction
  marley::ChebyshevInterpolatingFunction::cdf() const
{
  marley::ChebyshevInterpolatingFunction result;
  result.N_ = this->N_ + 1;
  result.x_min_ = this->x_min_;
  result.x_max_ = this->x_max_;

  double beta_0 = this->chebyshev_coeffs_.at(0);
  beta_0 += -0.5 * this->chebyshev_coeffs_.at(1);
  for ( size_t j = 2; j <= this->N_; ++j ) {
    beta_0 += 2.*this->chebyshev_coeffs_.at(j)
      * std::pow(-1.0, j + 1) / (j*j - 1.0);
  }

  result.chebyshev_coeffs_.push_back( beta_0 );
  for ( size_t k = 1; k < this->N_; ++k ) {
    double ak_minus_one = this->chebyshev_coeffs_.at( k - 1 );
    double ak_plus_one = this->chebyshev_coeffs_.at( k + 1 );
    result.chebyshev_coeffs_.push_back(
      (ak_minus_one - ak_plus_one) / (2 * k) );
  }
  result.chebyshev_coeffs_.push_back( this->chebyshev_coeffs_.at( N_ - 1 )
    / (2 * N_) );
  result.chebyshev_coeffs_.push_back( this->chebyshev_coeffs_.at( N_ )
    / (2 * (N_ + 1)) );

  for (double& d : result.chebyshev_coeffs_) d *= (x_max_ - x_min_) / 2.;

  result.compute_integral();

  for ( size_t j = 0; j <= result.N_; ++j ) {
    double x = result.chebyshev_point( j );
    result.Xs_.push_back( x );
  }

  result.wsave_ = std::vector<double>(
    3*result.chebyshev_coeffs_.size() + 15, 0.);
  result.ifac_ = std::vector<int>(result.chebyshev_coeffs_.size() / 2, 0);

  int my_size = result.chebyshev_coeffs_.size();
  result.Fs_ = result.chebyshev_coeffs_;
  costi( &my_size, result.wsave_.data(), result.ifac_.data() );
  cost( &my_size, result.Fs_.data(), result.wsave_.data(), result.ifac_.data() );
  for (double& d : result.Fs_) d *= 0.5;

  return result;
}
