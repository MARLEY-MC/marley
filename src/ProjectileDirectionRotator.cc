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

// Standard library includes

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Event.hh"
#include "marley/ProjectileDirectionRotator.hh"

marley::ProjectileDirectionRotator::ProjectileDirectionRotator(
  const std::array<double, 3>& dir ) : marley::EventProcessor(),
  dir_vec_( dir )
{
  constexpr ThreeVector null_three_vector = { 0., 0., 0. };
  if ( dir_vec_ == null_three_vector ) throw marley::Error( "Null 3-vector"
    " passed to constructor of marley::ProjectileDirectionRotator" );

  dir_vec_ = marley::RotationMatrix::normalize( dir_vec_ );
}

void marley::ProjectileDirectionRotator::process_event(marley::Event& ev,
  marley::Generator& /*gen*/)
{
  // First check that the projectile 3-momentum is not a null vector.
  // If it is, don't bother to rotate coordinates. Also don't bother if
  // (somehow) the magnitude of the projectile momentum is negative.
  const auto& projectile = ev.projectile();
  double pmom = projectile.momentum_magnitude();
  if ( pmom <= 0. ) return;

  ThreeVector pdir = { projectile.px() / pmom, projectile.py() / pmom,
    projectile.pz() / pmom };

  // If the projectile direction from the event matches the desired direction,
  // then we can skip the coordinate rotation.
  if ( pdir == dir_vec_ ) return;

  if ( pdir != last_pdir_ ) {
    last_pdir_ = pdir;
    rot_matrix_ = marley::RotationMatrix( pdir, dir_vec_ );
  }

  this->rotate_event( ev );
}

void marley::ProjectileDirectionRotator::rotate_event( marley::Event& ev ) {

  // Rotate the initial particles
  for ( auto* p : ev.get_initial_particles() ) {
    rot_matrix_.rotate_particle_inplace(*p);
  }

  // Rotate the final particles
  for ( auto* p : ev.get_final_particles() ) {
    rot_matrix_.rotate_particle_inplace(*p);
  }

}
