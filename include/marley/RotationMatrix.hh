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

// Standard library includes
#include <array>

// MARLEY includes
#include "marley/Particle.hh"

namespace marley {

  /// @brief Simple rotation matrix implementation used to reorient
  /// Particle objects based on the incident neutrino direction
  class RotationMatrix {

    using ThreeVector = std::array<double, 3>;
    using ThreeThreeMatrix = std::array<ThreeVector, 3>;

    public:

      /// @brief Creates a 3&times;3 identity matrix
      RotationMatrix();

      /// @brief Create a 3&times;3 rotation matrix that rotates the 3-vector
      /// from_vec into the 3-vector to_vec
      RotationMatrix(const ThreeVector& from_vec, const ThreeVector& to_vec);

      /// @brief Create a rotated copy of the 3-vector v
      ThreeVector rotate_copy(const ThreeVector& v);

      // Returns a copy of the 3-vector v normalized to have unit magnitude
      static ThreeVector normalize(const ThreeVector& v);

      /// @brief Rotate a 3-vector v in place
      void rotate_inplace(ThreeVector& v);

      /// @brief Rotate the 3-momentum of a marley::Particle in place
      void rotate_particle_inplace(marley::Particle& p);

    protected:

      /// @brief 3&times;3 rotation matrix
      ThreeThreeMatrix matrix_;
  };

}
