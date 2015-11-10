#include "TMarleyElectronReaction.hh"
#include "TMarleyGenerator.hh"

TMarleyElectronReaction::TMarleyElectronReaction(size_t Z) {
  Zatom = Z;

  // The target is an electron. Since this is an elastic scattering
  // reaction, the residue is also an electron.
  pid_b = marley_utils::ELECTRON;
  mb = TMarleyMassTable::get_particle_mass(pid_b);
  mb2 = std::pow(mb, 2);

  pid_d = pid_b;
  md = mb;
  md2 = mb2;

  // The particle ID for the projectile (and ejectile) will be updated
  // as this reaction is used, but use an electron neutrino as the
  // default.
  pid_a = marley_utils::ELECTRON_NEUTRINO;
  pid_c = pid_a;

  // All of the recognized projectiles for this reaction (neutrinos)
  // are approximately massless, so set the values of ma and mc accordingly
  ma = 0.;
  ma2 = 0.;
  mc = ma;
  mc2 = ma2;
}

bool TMarleyElectronReaction::determine_coupling_constants(int particle_id_a) {
  if (particle_id_a == marley_utils::ELECTRON_NEUTRINO) {
    g1 = 0.5 + marley_utils::sin2thetaw;
    g2 = marley_utils::sin2thetaw;
  }
  else if (particle_id_a == marley_utils::ELECTRON_ANTINEUTRINO) {
    g1 = marley_utils::sin2thetaw;
    g2 = 0.5 + marley_utils::sin2thetaw;
  }
  else if (particle_id_a == marley_utils::MUON_NEUTRINO ||
    particle_id_a == marley_utils::TAU_NEUTRINO)
  {
    g1 = -0.5 + marley_utils::sin2thetaw;
    g2 = marley_utils::sin2thetaw;
  }
  else if (particle_id_a == marley_utils::MUON_ANTINEUTRINO ||
    particle_id_a == marley_utils::TAU_ANTINEUTRINO)
  {
    g1 = marley_utils::sin2thetaw;
    g2 = -0.5 + marley_utils::sin2thetaw;
  }
  // The projectile particle ID was unrecognized, so don't bother
  // to set the coupling constants, and return false.
  else return false;

  // Update this helper variable for g2
  g2_squared_over_three = std::pow(g2, 2) / 3.0;

  // Update the projectile and ejectile particle IDs
  pid_a = particle_id_a;
  pid_c = pid_a;

  // The projectile particle ID was recognized, so return true
  // now that the coupling constants have been updated.
  return true;
}

double TMarleyElectronReaction::total_xs(int particle_id_a, double Ea) {

  // Determine the coupling constants g1 and g2 based on the projectile
  bool projectile_recognized = determine_coupling_constants(particle_id_a);
  // If the projectile's particle ID was unrecognized, set the cross section to
  // zero.
  // TODO: consider throwing an exception instead.
  if (!projectile_recognized) return 0.;

  // Mandelstam s (square of the total center of mass energy)
  double s = ma2 + mb2 + 2 * mb * Ea;

  // CM frame projectile energy
  double Ea_cm = (s - mb2) / (2 * std::sqrt(s));

  // Helper variable
  double me2_over_s = mb2 / s;
  // Total cross section in natural units (MeV^(-2))
  double xs = (4 / marley_utils::pi) * std::pow(marley_utils::GF * Ea_cm, 2)
    * (std::pow(g1, 2) + (g2_squared_over_three - g1*g2)*me2_over_s
    + g2_squared_over_three*(1 + std::pow(me2_over_s, 2)));
  // Multiply the single-electron cross section by the number of electrons
  // present in this atom
  xs *= Zatom;
  return xs;
}

// Differential cross section dsigma/dcos(theta_c_CM) per electron in units
// convenient for sampling a CM frame ejectile scattering cosine. To convert
// into a physical CM frame differential cross section (in units of MeV^(-2))
// multiply the output of this function by (2 / pi) * (GF * Ea_CM)^2, where
// Ea_CM is the projectile's total energy in the CM frame.
double TMarleyElectronReaction::diff_xs_cm_for_sampling(double s,
  double cos_theta_c_cm)
{
  return std::pow(g1, 2) + g1*g2*mb2*(cos_theta_c_cm - 1)/s
    + std::pow(g2 * (1 + (0.5 - mb2/(2*s)) * (cos_theta_c_cm - 1)), 2);
}

// Sample an ejectile scattering cosine in the CM frame. Make sure to update
// the coupling constants g1 and g2 before calling this function. The first
// argument of this function is Mandelstam s.
double TMarleyElectronReaction::sample_cos_theta_c_cm(double s,
  TMarleyGenerator& gen)
{
  // Generate a forwarding call wrapper for the differential cross section that
  // will be used during rejection sampling Note that our rejection sampling
  // technique does not require that we normalize the differential cross
  // section before sampling from it.
  std::function<double(double)> diff_xs = std::bind(
    &TMarleyElectronReaction::diff_xs_cm_for_sampling, this, s,
    std::placeholders::_1);

  // Sample a CM frame scattering cosine using the differential cross section
  return gen.rejection_sample(diff_xs, -1, 1);
}


// Creates an event object by sampling the appropriate quantities and
// performing kinematic calculations
TMarleyEvent TMarleyElectronReaction::create_event(int particle_id_a, double Ea,
  TMarleyGenerator& gen)
{
  // Determine the coupling constants g1 and g2 based on the projectile
  bool projectile_recognized = determine_coupling_constants(particle_id_a);
  // If the projectile's particle ID was unrecognized, complain and refuse to
  // create an event.
  if (!projectile_recognized) throw std::runtime_error(std::string("Could")
    + "not create this event. The projectile particle ID "
    + std::to_string(particle_id_a) + " was not recognized as a valid"
    + " projectile for the current reaction");

  double s, Ec_cm, pc_cm, Ed_cm;
  two_two_scatter(Ea, s, Ec_cm, pc_cm, Ed_cm);

  // Sample a CM frame scattering cosine for the ejectile.
  double cos_theta_c_cm = sample_cos_theta_c_cm(s, gen);

  // Sample a CM frame azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the differential cross section is independent of
  // the azimuthal angle.
  double phi_c_cm = gen.uniform_random_double(0, 2*marley_utils::pi, false);

  // Create and return the completed event object
  return make_event_object(Ea, pc_cm, cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm);
}
