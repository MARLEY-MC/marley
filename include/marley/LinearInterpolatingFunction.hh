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

#pragma once
#include <functional>
#include <vector>

#include "marley/InterpolationGrid.hh"
#include "marley/InterpolatingFunction.hh"

namespace marley {

  /// @brief Default number of grid points to use
  constexpr size_t DEFAULT_N_LINEAR = 200u;

  /// @brief Approximates a 1D function using linear interpolation on a grid
  /// of regularly-spaced points
  class LinearInterpolatingFunction : public InterpolatingFunction {

    public:

      LinearInterpolatingFunction( const std::function< double(double) >& func,
        double x_min, double x_max, size_t N = DEFAULT_N_LINEAR );

      inline double evaluate( double x ) const override {
        return interp_grid_->interpolate( x );
      }

      inline int N() const { return interp_grid_->size(); }

    protected:

      std::shared_ptr< InterpolationGrid<double, double> > interp_grid_;
  };

}
