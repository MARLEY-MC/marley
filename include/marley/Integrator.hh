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

#pragma once
#include <functional>
#include <vector>

#include "marley/marley_utils.hh"

namespace marley {

  /// @brief Numerical integrator that uses <a
  /// href="http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature">Clenshaw-Curtis
  /// quadrature</a>
  class Integrator {
    public:

      /// @brief Create a Clenshaw-Curtis quadrature integrator that uses
      /// 2*num sampling points
      /// @param num half the number of sampling points to use
      Integrator(size_t num = N_DEFAULT_);

      /// @brief Numerically integrate an arbitrary 1D function
      /// @details Numerically integrate a std::function<double(double)> over the
      /// interval [a,b] using <a
      /// href="http://en.wikipedia.org/wiki/Clenshaw-Curtis_quadrature">
      /// Clenshaw-Curtis quadrature</a> at 2N_ sampling points.
      double num_integrate(const std::function<double(double)>& f,
        double a, double b) const;

    private:

      /// @brief use 2N_ sampling points to perform numerical integration
      size_t N_;

      /// @brief precomputed weights to use to speed up integration
      std::vector<double> weights_;

      /// @brief storage for function evaluation point offsets
      std::vector<double> offsets_;

      static constexpr size_t N_DEFAULT_ = 20; ///< default value of N_
  };

}
