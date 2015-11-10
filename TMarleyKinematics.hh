#pragma once
#include <cmath>
#include <stdexcept>

class TMarleyParticle;

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

    // Gets the square of the total energy in the center of momentum frame for
    // two particles
    static double get_mandelstam_s(const TMarleyParticle& p1,
      const TMarleyParticle& p2);

  private:
    inline static double get_beta2(double beta_x, double beta_y,
      double beta_z)
    {
      double beta2 = std::pow(beta_x, 2) + std::pow(beta_y, 2)
        + std::pow(beta_z, 2);
      if (beta2 == 1) throw std::runtime_error(std::string("Cannot perform")
        + " Lorentz boost because \u03B2^2 = 1 and therefore the Lorentz factor"
        + "\u03B3 is infinite.");
      else if (beta2 > 1) throw std::runtime_error(std::string("Cannot perform")
        + " Lorentz boost because \u03B2^2 = " + std::to_string(beta2) + " > 1,"
        + " which is unphysical.");
      return beta2;
    }
};
