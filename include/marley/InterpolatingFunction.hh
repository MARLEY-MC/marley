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

#include "marley/marley_utils.hh"

namespace marley {

  /// @brief Abstract base class for an approximate representation of a 1D
  /// continuous function
  class InterpolatingFunction {

    public:

      InterpolatingFunction( double x_min, double x_max ) : x_min_( x_min ),
        x_max_( x_max ) {}

      virtual ~InterpolatingFunction() = default;

      /// @brief Returns an approximate value of the represented function
      virtual double evaluate( double x ) const = 0;

    protected:

      /// @brief Allows construction by derived classes without immediately
      /// defining the bounds of the interpolation range
      InterpolatingFunction() {}

      /// @brief Lower bound of the interpolation range
      double x_min_;

      /// @brief Upper bound of the interpolation range
      double x_max_;
  };

}
