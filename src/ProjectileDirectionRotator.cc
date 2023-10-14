/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
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

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Event.hh"
#include "marley/Generator.hh"
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
  marley::Generator& gen)
{
  // First check that the projectile 3-momentum is not a null vector.
  // If it is, don't bother to rotate coordinates. Also don't bother if
  // (somehow) the magnitude of the projectile momentum is negative.
  const auto& projectile = ev.projectile();
  double pmom = projectile.momentum_magnitude();
  if ( pmom <= 0. ) return;

  ThreeVector pdir = { projectile.px() / pmom, projectile.py() / pmom,
    projectile.pz() / pmom };

  // If random projectile directions have been requested, then sample
  // a new one isotropically for this event
  if ( randomize_projectile_direction_ ) {
    ThreeVector random_dir = sample_isotropic_direction( gen );
    this->set_projectile_direction( random_dir );
  }

  // If the (unrotated) projectile direction from the event exactly matches the
  // desired direction, then we can skip the coordinate rotation.
  if ( pdir == dir_vec_ ) return;

  // If the (unrotated) projectile direction from the event differs from
  // the last one that was processed, then we need to recompute the
  // rotation matrix. Do so and save the unrotated projectile direction
  // to repeat this check next time.
  if ( pdir != last_pdir_ ) {
    last_pdir_ = pdir;
    rot_matrix_ = marley::RotationMatrix( pdir, dir_vec_ );
  }

  // Do the actual rotation of the Particle 3-momenta in the event
  this->rotate_event( ev );
}

void marley::ProjectileDirectionRotator::rotate_event( marley::Event& ev ) {

  // Rotate the initial particles
  for ( auto* p : ev.get_initial_particles() ) {
    rot_matrix_.rotate_particle_inplace( *p );
  }

  // Rotate the final particles
  for ( auto* p : ev.get_final_particles() ) {
    rot_matrix_.rotate_particle_inplace( *p );
  }

}

marley::ProjectileDirectionRotator::ThreeVector
  marley::ProjectileDirectionRotator::sample_isotropic_direction(
  marley::Generator& gen ) const
{
  // Sample a polar cosine on the interval [-1, 1]
  double cos_theta = gen.uniform_random_double( -1., 1., true );

  // Sample an azimuthal angle on the interval [0, 2*pi)
  double phi = gen.uniform_random_double( 0., marley_utils::two_pi, false );

  // Compute direction unit vector components
  double sin_theta = marley_utils::real_sqrt( 1. - std::pow(cos_theta, 2) );
  double ux = sin_theta * std::cos( phi );
  double uy = sin_theta * std::sin( phi );
  double uz = cos_theta;

  ThreeVector direction = { ux, uy, uz };
  return direction;
}

void marley::ProjectileDirectionRotator::set_randomize_directions(
  bool do_sampling )
{
  randomize_projectile_direction_ = do_sampling;
}
