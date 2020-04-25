/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/COPYING or
// visit http://opensource.org/licenses/GPL-3.0

#pragma once

// MARLEY includes
#include "marley/EventProcessor.hh"
#include "marley/RotationMatrix.hh"

namespace marley {

  /// @brief If needed, rotates the coordinate system of an Event so that the
  /// projectile 3-momentum lies along a desired direction
  class ProjectileDirectionRotator : public EventProcessor {

    using ThreeVector = std::array<double, 3>;

    public:

      /// @param dir A 3-vector pointing in the desired direction
      /// of the projectile in the rotated coordinate system
      ProjectileDirectionRotator( const ThreeVector& dir = {0., 0., 1.} );

      inline virtual ~ProjectileDirectionRotator() = default;

      /// @brief Rotates all 3-momenta in the input event so that
      /// the projectile 3-momentum lies along dir_vec_ in the
      /// new coordinate system
      virtual void process_event( marley::Event& ev,
        marley::Generator& gen ) override;

      inline const ThreeVector& projectile_direction() const
        { return dir_vec_; }

      inline void set_projectile_direction( const ThreeVector& dir ) {
        dir_vec_ = marley::RotationMatrix::normalize( dir );
        last_pdir_ = dir_vec_;
      }

      ThreeVector sample_isotropic_direction( marley::Generator& gen ) const;

      void set_randomize_directions( bool do_sampling );

    protected:

      /// @brief 3-vector that points in the desired direction of the
      /// projectile
      ThreeVector dir_vec_ = { 0., 0., 1. };

      /// @brief Stores the direction 3-vector for the last projectile
      /// that triggered a recalculation of the rotation matrix.
      /// @details Using this information avoids unnecessary recalculations
      ThreeVector last_pdir_ = { 0., 0., 1. };

      /// @brief RotationMatrix used to rotate the coordinate system
      /// of the input Event
      marley::RotationMatrix rot_matrix_;

      /// @brief Helper function that does the coordinate system rotation
      /// @param[in,out] ev Event whose Particle 3-vectors will be rotated
      void rotate_event( marley::Event& ev );

      /// @brief Flag that indicates whether the (rotated) projectile
      /// direction should be sampled isotropically for each event
      /// (true) or kept fixed across all events (false)
      bool randomize_projectile_direction_ = false;
  };

}
