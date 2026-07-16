/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
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

// MARLEY includes
#include "marley/LinearInterpolatingFunction.hh"

marley::LinearInterpolatingFunction::LinearInterpolatingFunction(
  const std::function<double(double)>& func, double x_min, double x_max,
  size_t N )
  : marley::InterpolatingFunction( x_min, x_max )
{
  // Force at least two points to be used
  // TODO: add error handling here
  if ( N < 2u ) N = 2u;

  std::vector< double > x_vec;
  std::vector< double > y_vec;

  double dx = ( x_max - x_min ) / ( N - 1 );
  for ( size_t s = 0u; s < N; ++s ) {
    double x = x_min + s*dx;
    double y = func( x );

    x_vec.push_back( x );
    y_vec.push_back( y );
  }

  interp_grid_ = std::make_shared<
    marley::InterpolationGrid< double, double > >( x_vec, y_vec );
}
