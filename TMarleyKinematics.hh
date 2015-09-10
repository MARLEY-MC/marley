#pragma once

class TMarleyKinematics {
  public:
    // Lorentz boost a particle, replacing its energy and momentum
    // with the boosted versions
    static void lorentz_boost(double beta_x, double beta_y, double beta_z,
      TMarleyParticle& particle_to_boost);

    // Load two product particles with the appropriate lab-frame energies and
    // momenta based on a decay of an initial particle (whose energy and
    // momentum are assumed to be given in the lab frame). The direction
    // cosine for the polar angle of the first particle and the azimuthal
    // angle of the first particle (both as measured in the rest frame of the
    // initial particle) must also be given as input.
    static void two_body_decay(const TMarleyParticle& initial_particle,
      TMarleyParticle& first_product, TMarleyParticle& second_product,
      double cos_theta_first, double phi_first);

    // Rotates a particle's 3-momentum so that it points in the (x, y, z) direction
    static void rotate_momentum_vector(double x, double y, double z,
      TMarleyParticle& particle_to_rotate);
};
