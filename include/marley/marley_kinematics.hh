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

#pragma once
#include <cmath>

#include "marley/Error.hh"

namespace marley {
  class Particle;
}

namespace marley_kinematics {

  // Lorentz boost a particle, replacing its energy and momentum with the
  // boosted versions
  void lorentz_boost(double beta_x, double beta_y, double beta_z,
    marley::Particle& particle_to_boost);

  /// @brief Handles kinematic calculations needed to decay an initial
  /// particle into two final particles
  /// @details This function load two product Particle objects with the
  /// appropriate lab-frame energies and momenta (assuming they are produced on
  /// the mass shell) based on a decay of an initial particle (whose energy and
  /// momentum are assumed to be given in the lab frame). The direction cosine
  /// for the polar angle of the first particle and the azimuthal angle of the
  /// first particle (both as measured in the rest frame of the initial
  /// particle) must also be given as input. The two product Particle
  /// objects are expected to have their mass_ member variables set before
  /// this function is called.
  /// @param[in] initial_particle The mother particle that will decay
  /// into the two daughter particles
  /// @param[inout] first_product The first of the two daughter particles
  /// @param[inout] second_product The second of the two daughter particles
  /// @param[in] cos_theta_first The cosine of the polar emission angle for
  /// the first daughter particle (as measured in the rest frame of the
  /// mother particle)
  /// @param[in phi_first The azimuthal emission angle for the first daughter
  /// particle (as measured in the rest frame of the mother particle)
  void two_body_decay(const marley::Particle& initial_particle,
    marley::Particle& first_product, marley::Particle& second_product,
    double cos_theta_first, double phi_first);

  // Rotates a particle's 3-momentum so that it points in the (x, y, z)
  // direction
  void rotate_momentum_vector(double x, double y, double z,
    marley::Particle& particle_to_rotate);

  // Gets the square of the total energy in the center of momentum frame for
  // two particles
  double get_mandelstam_s(const marley::Particle& p1,
    const marley::Particle& p2);

  // Boost two particles into their mutual center-of-momentum frame
  void boost_to_cm_frame(marley::Particle& p1, marley::Particle& p2);
}
