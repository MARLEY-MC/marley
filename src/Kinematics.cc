#include <stdexcept>

#include "marley_utils.hh"
#include "Particle.hh"
#include "Kinematics.hh"

// Rotates a particle's 3-momentum so that it points in the (x, y, z) direction
void marley::Kinematics::rotate_momentum_vector(double x, double y, double z,
  marley::Particle& particle_to_rotate)
{
  // Get the magnitude of the vector pointing in the desired direction
  double r = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));

  // Get magnitude of the particle's momentum
  double rp = particle_to_rotate.get_momentum_magnitude();

  // Rotate the particle's momentum vector into the desired direction
  double ratio = rp / r;

  double new_px = x * ratio;
  double new_py = y * ratio;
  double new_pz = z * ratio;

  particle_to_rotate.set_px(new_px);
  particle_to_rotate.set_py(new_py);
  particle_to_rotate.set_pz(new_pz);
}

// Lorentz boost a particle, replacing its energy and momentum with the boosted
// versions
void marley::Kinematics::lorentz_boost(double beta_x, double beta_y,
  double beta_z, marley::Particle& particle_to_boost)
{
  double beta2 = get_beta2(beta_x, beta_y, beta_z);

  // If beta is zero in all directions, then we don't need to do the boost at all
  if (beta2 == 0) return;

  // Calculate the Lorentz factor based on the boost velocity
  double gamma = 1 / std::sqrt(1 - beta2);

  double E = particle_to_boost.get_total_energy();
  double px = particle_to_boost.get_px();
  double py = particle_to_boost.get_py();
  double pz = particle_to_boost.get_pz();
  double m = particle_to_boost.get_mass();

  // Compute the boosted energy and 3-momentum for the particle (the expressions
  // we use here are based on
  // https://en.wikipedia.org/wiki/Lorentz_transformation#Boost_in_any_direction)
  double beta_dot_p = beta_x * px + beta_y * py + beta_z * pz;
  double factor = (gamma - 1) * beta_dot_p / beta2;

  double new_E = gamma * (E - beta_dot_p);

  // The new energy could conceivably dip slightly below the mass
  // of the particle due to roundoff errors. If this is the case,
  // set it to the particle mass.
  if (new_E < m) new_E = m;

  double new_px = (-gamma * E + factor) * beta_x + px;
  double new_py = (-gamma * E + factor) * beta_y + py;
  double new_pz = (-gamma * E + factor) * beta_z + pz;

  particle_to_boost.set_total_energy(new_E);
  particle_to_boost.set_px(new_px);
  particle_to_boost.set_py(new_py);
  particle_to_boost.set_pz(new_pz);
}

// Right now, this function assumes that the coordinate axes in the lab
// and initial particle rest frames are coincident (i.e., the coordinate axes
// aren't rotated between the two frames at all). Since all interactions in
// the current version of the code are isotropic (except for the initial
// 2-body neutrino scattering reaction, which occurs for an incident neutrino
// traveling in the positive z direction in the lab frame [target nucleus's
// rest frame]), this is all we need for now.
// TODO: Expand this to allow the two frames to be rotated with respect to each other.
void marley::Kinematics::two_body_decay(const marley::Particle& initial_particle,
  marley::Particle& first_product, marley::Particle& second_product,
  double cos_theta_first, double phi_first)
{
  // Get the masses of all three particles
  double M = initial_particle.get_mass();
  double mfirst = first_product.get_mass();
  double msecond = second_product.get_mass();

  // Check to make sure the decay is kinematically allowed
  if (M < mfirst + msecond) throw std::runtime_error(std::string("A two-body")
    + " decay was requested that is not kinematically allowed.");

  double M2 = std::pow(M, 2);
  double mfirst2 = std::pow(mfirst, 2);
  double msecond2 = std::pow(msecond, 2);

  // Compute the energies for the decay products in the rest
  // frame of the initial particle
  double Efirst = (M2 - msecond2 + mfirst2) / (2 * M);
  double Esecond = M - Efirst; // M - Efirst == (M2 - mfirst2 + msecond2) / (2 * M)

  // Avoid roundoff issues by not allowing the energies to dip below
  // the particle masses
  if (Efirst < mfirst) Efirst = mfirst;
  if (Esecond < msecond) Esecond = msecond;

  // Compute the 3-momenta for the decay products, still in the initial
  // particle's rest frame
  double pfirst = marley_utils::real_sqrt(std::pow(Efirst, 2) - mfirst2);
  double sin_theta_first = marley_utils::real_sqrt(1
    - std::pow(cos_theta_first, 2));
  double p1x = pfirst * sin_theta_first * std::cos(phi_first);
  double p1y = pfirst * sin_theta_first * std::sin(phi_first);
  double p1z = pfirst * cos_theta_first;

  // Conservation of 3-momenta tells us that, in the rest frame of the
  // initial particle, p1 = -p2. We can use this as a shortcut.
  double p2x = -p1x;
  double p2y = -p1y;
  double p2z = -p1z;

  // Now that we have this information, load it into the product
  // particles.
  first_product.set_total_energy(Efirst);
  first_product.set_px(p1x);
  first_product.set_py(p1y);
  first_product.set_pz(p1z);

  second_product.set_total_energy(Esecond);
  second_product.set_px(p2x);
  second_product.set_py(p2y);
  second_product.set_pz(p2z);

  // Compute the parameters needed to boost these particles
  // from the initial particle's rest frame into the lab frame
  double E_i = initial_particle.get_total_energy();
  double px_i = initial_particle.get_px();
  double py_i = initial_particle.get_py();
  double pz_i = initial_particle.get_pz();

  // Boost in the opposite direction (this gives us the minus signs below) from
  // the 3-momentum of rest_frame_particle. This takes the null 3-vector to the
  // initial particle's 3-momentum.
  double beta_x = -px_i / E_i;
  double beta_y = -py_i / E_i;
  double beta_z = -pz_i / E_i;

  // Boost both products to the lab frame by replacing their
  // energies and momenta with the boosted versions
  lorentz_boost(beta_x, beta_y, beta_z, first_product);
  lorentz_boost(beta_x, beta_y, beta_z, second_product);
}

// Get the square of the total energy of two particles in their center of
// momentum frame
double marley::Kinematics::get_mandelstam_s(const marley::Particle& p1,
  const marley::Particle& p2)
{
  double E1 = p1.get_total_energy();
  double m1 = p1.get_mass();
  double E2 = p2.get_total_energy();
  double m2 = p2.get_mass();

  // If one of the particles is at rest, use a shortcut. Otherwise,
  // Lorentz transform to the center of momentum frame to determine
  // the total cm frame energy
  if (E1 == m1) return std::pow(m1, 2) + std::pow(m2, 2) + 2 * m1 * E2;
  else if (E2 == m2) return std::pow(m1, 2) + std::pow(m2, 2) + 2 * m2 * E1;
  else {
    // Get total energy, momentum, and mass values for the two particles
    double E_tot = E1 + E2;
    double px_tot = p1.get_px() + p2.get_px();
    double py_tot = p1.get_py() + p2.get_py();
    double pz_tot = p1.get_pz() + p2.get_pz();
    double m_tot = m1 + m2;

    // Get boost parameters for a Lorentz transform to the center of momentum
    // frame
    double beta_x = px_tot / E_tot;
    double beta_y = py_tot / E_tot;
    double beta_z = pz_tot / E_tot;

    // Calculate the Lorentz factor based on the boost velocity
    double beta2 = get_beta2(beta_x, beta_y, beta_z);
    double gamma = 1 / std::sqrt(1 - beta2);

    // Compute the total boosted energy in the cm frame
    double beta_dot_p_tot = beta_x * px_tot + beta_y * py_tot + beta_z * pz_tot;
    double E_tot_cm = gamma * (E_tot - beta_dot_p_tot);

    // The new energy could conceivably dip slightly below the total rest mass
    // of the particles due to roundoff errors. If this is the case,
    // set it to the total rest mass.
    if (E_tot_cm < m_tot) E_tot_cm = m_tot;

    return std::pow(E_tot_cm, 2);
  }
}

//// Compute final particle energies and momenta for the two-two scattering
//// reaction b(a,c)d and store the results in particles c and d. Note that
//// particles a and b are passed by value so that we can Lorentz boost
//// them to the center-of-momentum frame without altering the originals.
//void marley::Kinematics::two_two_scatter(marley::Particle a, marley::Particle b,
//  marley::Particle &c, marley::Particle &d, double cos_theta_c, double phi_c)
//{
//  // Get total energy and momentum values for the two initial particles
//  double E_tot = a.get_total_energy() + b.get_total_energy();
//  double px_tot = a.get_px() + b.get_px();
//  double py_tot = a.get_py() + b.get_py();
//  double pz_tot = a.get_pz() + b.get_pz();
//
//  // Boost both initial particles into the center of momentum frame
//  double beta_x = px_tot / E_tot;
//  double beta_y = py_tot / E_tot;
//  double beta_z = pz_tot / E_tot;
//  lorentz_boost(beta_x, beta_y, beta_z, a);
//  lorentz_boost(beta_x, beta_y, beta_z, b);
//
//  // Get the total center of mass energy
//  double E_tot_cm = a.get_total_energy() + b.get_total_energy();
//}
