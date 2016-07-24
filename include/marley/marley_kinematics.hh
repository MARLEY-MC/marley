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

  // Load two product particles with the appropriate lab-frame energies and
  // momenta based on a decay of an initial particle (whose energy and momentum
  // are assumed to be given in the lab frame). The direction cosine for the
  // polar angle of the first particle and the azimuthal angle of the first
  // particle (both as measured in the rest frame of the initial particle) must
  // also be given as input.
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
}
