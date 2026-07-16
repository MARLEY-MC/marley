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

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenParticle.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/marley_kinematics.hh"

namespace {

  double get_beta2(double beta_x, double beta_y, double beta_z)
  {
    double beta2 = std::pow(beta_x, 2) + std::pow(beta_y, 2)
      + std::pow(beta_z, 2);
    if (beta2 == 1) throw marley::Error(std::string("Cannot perform")
      + " Lorentz boost because \u03B2^2 = 1 and therefore the Lorentz factor"
      + "\u03B3 is infinite.");
    else if (beta2 > 1) throw marley::Error(std::string("Cannot perform")
      + " Lorentz boost because \u03B2^2 = " + std::to_string(beta2) + " > 1,"
      + " which is unphysical.");
    return beta2;
  }

}

// Rotates a particle's 3-momentum so that it points in the (x, y, z) direction
void marley_kinematics::rotate_momentum_vector( double x, double y, double z,
  HepMC3::GenParticle& particle_to_rotate )
{
  // Get the magnitude of the vector pointing in the desired direction
  double r = std::sqrt( std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2) );

  // Get magnitude of the particle's 3-momentum
  HepMC3::FourVector mom4 = particle_to_rotate.momentum();
  double rp = mom4.p3mod();

  // Rotate the particle's momentum vector into the desired direction
  double ratio = rp / r;

  double new_px = x * ratio;
  double new_py = y * ratio;
  double new_pz = z * ratio;

  mom4.set_px( new_px );
  mom4.set_py( new_py );
  mom4.set_pz( new_pz );

  particle_to_rotate.set_momentum( mom4 );
}

// Lorentz boost a particle, replacing its energy and momentum with the boosted
// versions
void marley_kinematics::lorentz_boost( double beta_x, double beta_y,
  double beta_z, HepMC3::GenParticle& particle_to_boost )
{
  double beta2 = get_beta2(beta_x, beta_y, beta_z);

  // If beta is zero in all directions, then we don't need to do the boost at
  // all
  if ( beta2 == 0. ) return;

  // Calculate the Lorentz factor based on the boost velocity
  double gamma = 1. / std::sqrt( 1. - beta2 );

  HepMC3::FourVector mom4 = particle_to_boost.momentum();

  double E = mom4.e();
  double px = mom4.px();
  double py = mom4.py();
  double pz = mom4.pz();
  double m = particle_to_boost.generated_mass();

  // Compute the boosted energy and 3-momentum for the particle. The
  // expressions we use here are based on
  // https://en.wikipedia.org/wiki/Lorentz_transformation#Boost_in_any_direction
  double beta_dot_p = beta_x * px + beta_y * py + beta_z * pz;
  double factor = ( gamma - 1. ) * beta_dot_p / beta2;

  double new_E = gamma * ( E - beta_dot_p );

  // The new energy could conceivably dip slightly below the mass
  // of the particle due to roundoff errors. If this is the case,
  // set it to the particle mass.
  if (new_E < m) new_E = m;

  double new_px = ( -gamma * E + factor ) * beta_x + px;
  double new_py = ( -gamma * E + factor ) * beta_y + py;
  double new_pz = ( -gamma * E + factor ) * beta_z + pz;

  mom4.set_e( new_E );
  mom4.set_px( new_px );
  mom4.set_py( new_py );
  mom4.set_pz( new_pz );

  particle_to_boost.set_momentum( mom4 );
}

// Right now, this function assumes that the coordinate axes in the lab and
// initial particle rest frames are coincident (i.e., the coordinate axes
// aren't rotated between the two frames at all). Since all interactions in the
// current version of the code are isotropic (except for the initial 2-body
// neutrino scattering reaction, which occurs for an incident neutrino
// traveling in the positive z direction in the lab frame [target nucleus's
// rest frame]), this is all we need for now.
// TODO: Expand this to allow the two frames to be rotated with respect to each
// other.
void marley_kinematics::two_body_decay(
  const std::shared_ptr< HepMC3::GenParticle >& initial_particle,
  std::shared_ptr< HepMC3::GenParticle >& first_product,
  std::shared_ptr< HepMC3::GenParticle >& second_product,
  double cos_theta_first, double phi_first)
{
  // Get the masses of all three particles
  double M = initial_particle->generated_mass();
  double mfirst = first_product->generated_mass();
  double msecond = second_product->generated_mass();

  // Check to make sure the decay is kinematically allowed
  if ( M < mfirst + msecond ) throw marley::Error( "A two-body decay was"
    " requested that is not kinematically allowed." );

  double M2 = std::pow( M, 2 );
  double mfirst2 = std::pow( mfirst, 2 );
  double msecond2 = std::pow( msecond, 2 );

  // Compute the energies for the decay products in the rest
  // frame of the initial particle
  double Efirst = ( M2 - msecond2 + mfirst2 ) / ( 2 * M );
  // M - Efirst == (M2 - mfirst2 + msecond2) / (2 * M)
  double Esecond = M - Efirst;

  // Avoid roundoff issues by not allowing the energies to dip below
  // the particle masses
  if ( Efirst < mfirst ) Efirst = mfirst;
  if ( Esecond < msecond ) Esecond = msecond;

  // Compute the 3-momenta for the decay products, still in the initial
  // particle's rest frame
  double pfirst = marley_utils::real_sqrt( std::pow(Efirst, 2) - mfirst2 );
  double sin_theta_first = marley_utils::real_sqrt( 1.
    - std::pow(cos_theta_first, 2) );
  double p1x = pfirst * sin_theta_first * std::cos( phi_first );
  double p1y = pfirst * sin_theta_first * std::sin( phi_first );
  double p1z = pfirst * cos_theta_first;

  // Conservation of 3-momenta tells us that, in the rest frame of the
  // initial particle, p1 = -p2. We can use this as a shortcut.
  double p2x = -p1x;
  double p2y = -p1y;
  double p2z = -p1z;

  // Now that we have this information, load it into the product
  // particles.
  HepMC3::FourVector first_mom4( p1x, p1y, p1z, Efirst );
  first_product->set_momentum( first_mom4 );

  HepMC3::FourVector second_mom4( p2x, p2y, p2z, Esecond );
  second_product->set_momentum( second_mom4 );

  // Compute the parameters needed to boost these particles
  // from the initial particle's rest frame into the lab frame
  const HepMC3::FourVector& initial_mom4 = initial_particle->momentum();
  double E_i = initial_mom4.e();
  double px_i = initial_mom4.px();
  double py_i = initial_mom4.py();
  double pz_i = initial_mom4.pz();

  // Boost in the opposite direction (this gives us the minus signs below) from
  // the 3-momentum of rest_frame_particle. This takes the null 3-vector to the
  // initial particle's 3-momentum.
  double beta_x = -px_i / E_i;
  double beta_y = -py_i / E_i;
  double beta_z = -pz_i / E_i;

  // Boost both products to the lab frame by replacing their
  // energies and momenta with the boosted versions
  lorentz_boost( beta_x, beta_y, beta_z, *first_product );
  lorentz_boost( beta_x, beta_y, beta_z, *second_product );
}

// Get the square of the total energy of two particles in their center of
// momentum frame
double marley_kinematics::get_mandelstam_s(
  const HepMC3::GenParticle& p1, const HepMC3::GenParticle& p2)
{
  const HepMC3::FourVector p1_mom4 = p1.momentum();
  const HepMC3::FourVector p2_mom4 = p2.momentum();

  double E1 = p1_mom4.e();
  double m1 = p1.generated_mass();

  double E2 = p2_mom4.e();
  double m2 = p2.generated_mass();

  // If one of the particles is at rest, use a shortcut. Otherwise,
  // Lorentz transform to the center of momentum frame to determine
  // the total cm frame energy
  if (E1 == m1) return std::pow(m1, 2) + std::pow(m2, 2) + 2 * m1 * E2;
  else if (E2 == m2) return std::pow(m1, 2) + std::pow(m2, 2) + 2 * m2 * E1;
  else {
    // Get total energy, momentum, and mass values for the two particles
    double E_tot = E1 + E2;
    double px_tot = p1_mom4.px() + p2_mom4.px();
    double py_tot = p1_mom4.py() + p2_mom4.py();
    double pz_tot = p1_mom4.pz() + p2_mom4.pz();
    double m_tot = m1 + m2;

    // Get boost parameters for a Lorentz transform to the center of momentum
    // frame
    double beta_x = px_tot / E_tot;
    double beta_y = py_tot / E_tot;
    double beta_z = pz_tot / E_tot;

    // Calculate the Lorentz factor based on the boost velocity
    double beta2 = get_beta2( beta_x, beta_y, beta_z );
    double gamma = 1.0 / std::sqrt( 1.0 - beta2 );

    // Compute the total boosted energy in the cm frame
    double beta_dot_p_tot = beta_x * px_tot + beta_y * py_tot + beta_z * pz_tot;
    double E_tot_cm = gamma * ( E_tot - beta_dot_p_tot );

    // The new energy could conceivably dip slightly below the total rest mass
    // of the particles due to roundoff errors. If this is the case,
    // set it to the total rest mass.
    if ( E_tot_cm < m_tot ) E_tot_cm = m_tot;

    return std::pow( E_tot_cm, 2 );
  }

}


// Boost two particles into their mutual CM frame
void marley_kinematics::boost_to_cm_frame( HepMC3::GenParticle& p1,
  HepMC3::GenParticle& p2)
{
  // Get parameters for a Lorentz boost to the center of momentum frame
  const HepMC3::FourVector& p1_mom4 = p1.momentum();
  const HepMC3::FourVector& p2_mom4 = p2.momentum();

  double E_tot = p1_mom4.e() + p2_mom4.e();
  double px_tot = p1_mom4.px() + p2_mom4.px();
  double py_tot = p1_mom4.py() + p2_mom4.py();
  double pz_tot = p1_mom4.pz() + p2_mom4.pz();

  double beta_x = px_tot / E_tot;
  double beta_y = py_tot / E_tot;
  double beta_z = pz_tot / E_tot;

  lorentz_boost( beta_x, beta_y, beta_z, p1 );
  lorentz_boost( beta_x, beta_y, beta_z, p2 );
}
