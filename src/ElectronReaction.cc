#include "marley/marley_utils.hh"
#include "marley/ElectronReaction.hh"
#include "marley/Generator.hh"

namespace {

  // Allowed range of cos(x)
  constexpr double COS_MIN = -1.;
  constexpr double COS_MAX =  1.;

}

marley::ElectronReaction::ElectronReaction(int pdg_a, unsigned int target_Z)
{
  process_type_ = ProcessType::NuElectronElastic;

  pdg_a_ = pdg_a;
  pdg_b_ = marley_utils::ELECTRON;
  pdg_c_ = get_ejectile_pdg(pdg_a_, process_type_);
  pdg_d_ = marley_utils::ELECTRON;

  // Set the description string based on the particle PDG codes
  description_ = marley_utils::get_particle_symbol( pdg_a_ )
    + " + " + marley_utils::get_particle_symbol( pdg_b_ );
  description_ += " --> " + marley_utils::get_particle_symbol( pdg_c_ );
  description_ += " + " + marley_utils::get_particle_symbol( pdg_d_ );

  // How many electrons are in the target atom?
  Z_atom_ = target_Z;

  const auto& mt = marley::MassTable::Instance();
  ma_ = mt.get_particle_mass( pdg_a_ );
  mb_ = mt.get_particle_mass( pdg_b_ );
  mc_ = mt.get_particle_mass( pdg_c_ );
  md_ = mt.get_particle_mass( pdg_d_ );

  this->set_coupling_constants();

  // Compute the kinetic energy of the projectile needed
  // for this reaction to proceed at threshold
  KEa_threshold_ = ( std::pow(mc_ + md_, 2)
    - std::pow(ma_ + mb_, 2) ) / ( 2.*mb_ );
}

void marley::ElectronReaction::set_coupling_constants()
{
  if (pdg_a_ == marley_utils::ELECTRON_NEUTRINO) {
    g1_ = marley_utils::ONE_HALF + marley_utils::sin2thetaw;
    g2_ = marley_utils::sin2thetaw;
  }
  else if (pdg_a_ == marley_utils::ELECTRON_ANTINEUTRINO) {
    g1_ = marley_utils::sin2thetaw;
    g2_ = marley_utils::ONE_HALF + marley_utils::sin2thetaw;
  }
  else if (pdg_a_ == marley_utils::MUON_NEUTRINO ||
    pdg_a_ == marley_utils::TAU_NEUTRINO)
  {
    g1_ = -marley_utils::ONE_HALF + marley_utils::sin2thetaw;
    g2_ = marley_utils::sin2thetaw;
  }
  else if (pdg_a_ == marley_utils::MUON_ANTINEUTRINO ||
    pdg_a_ == marley_utils::TAU_ANTINEUTRINO)
  {
    g1_ = marley_utils::sin2thetaw;
    g2_ = -marley_utils::ONE_HALF + marley_utils::sin2thetaw;
  }
  else throw marley::Error("Unrecognized projectile PDG code "
    + std::to_string(pdg_a_) + " encountered in marley::Electron"
    "Reaction::set_coupling_constants()");
}

double marley::ElectronReaction::total_xs(int pdg_a, double KEa) {

  // If the cross section was requested for a different projectile,
  // then just return zero.
  if ( pdg_a != pdg_a_ ) return 0.;

  // If we're below threshold, then just return zero
  if ( KEa < KEa_threshold_ ) return 0.;

  // Mandelstam s (square of the total center of momentum frame energy)
  double s = std::pow(ma_ + mb_, 2) + 2*mb_*KEa;

  // CM frame ejectile total energy
  double Ec_cm = (s + mc_*mc_ - md_*md_) / ( 2. * std::sqrt(s) );

  // Helper variables
  double me2_over_s = md_*md_ / s;
  double g2_squared_over_three = g2_*g2_ / 3.;

  // Total cross section in natural units (MeV^(-2))
  double xs = (4. / marley_utils::pi) * std::pow(marley_utils::GF * Ec_cm, 2)
    * (std::pow(g1_, 2) + (g2_squared_over_three - g1_*g2_)*me2_over_s
    + g2_squared_over_three*(1. + std::pow(me2_over_s, 2)));

  // Multiply the single-electron cross section by the number of electrons
  // present in this atom
  /// @todo Take into account effects of electron binding energy
  xs *= Z_atom_;
  return xs;
}

double marley::ElectronReaction::diff_xs(int pdg_a, double KEa,
  double cos_theta_c_cm)
{
  // If we're asked for the wrong projectile, then just return zero
  if ( pdg_a != pdg_a_ ) return 0.;

  // If we're below threshold, then just return zero
  if ( KEa < KEa_threshold_ ) return 0.;

  // Mandelstam s (square of the total center of momentum frame energy)
  double s = std::pow(ma_ + mb_, 2) + 2.*mb_*KEa;

  // CM frame ejectile total energy
  double Ec_cm = (s + mc_*mc_ - md_*md_) / ( 2. * std::sqrt(s) );

  // Helper variables
  double me2_over_s = md_*md_ / s;
  double overall_factor = (2. / marley_utils::pi)
    * std::pow(marley_utils::GF * Ec_cm, 2);
  double terms = std::pow(g1_, 2) + g1_*g2_*me2_over_s*(cos_theta_c_cm - 1.)
    + std::pow(g2_ * (1. + marley_utils::ONE_HALF*(1. - me2_over_s)
    * (cos_theta_c_cm - 1.)), 2);

  // Compute and return the cross section
  double diff_xsec = overall_factor * terms;
  return diff_xsec;
}

// Creates an event object by sampling the appropriate quantities and
// performing kinematic calculations
marley::Event marley::ElectronReaction::create_event(int pdg_a, double KEa,
  marley::Generator& gen)
{
  // If the projectile's PDG code doesn't match that stored in this object,
  // complain and refuse to create an event.
  if ( pdg_a != pdg_a_ ) throw marley::Error("Could"
    " not create this event. The requested projectile PDG code "
    + std::to_string(pdg_a) + " does not match the value "
    + std::to_string(pdg_a_) + " stored in the Reaction object.");

  // Also complain if we're below threshold
  if ( KEa < KEa_threshold_ ) throw marley::Error("Could not create"
    " this event. The kinetic energy (" + std::to_string(KEa) + " MeV)"
    " of the projectile is below the reaction threshold of "
    + std::to_string(KEa_threshold_) + " MeV.");

  double s, Ec_cm, pc_cm, Ed_cm;
  two_two_scatter(KEa, s, Ec_cm, pc_cm, Ed_cm);

  // Compute the maximum differential cross section to use for rejection
  // sampling. To do this, we analytically solve for the value of
  // cos_theta_c_cm (labeled cth below) for which the derivative of the
  // differential cross section vanishes. This is an extremum of the function
  // and might correspond to the maximum. If cth is within the allowed angular
  // range, then it is considered alongside the two endpoints (COS_MIN and
  // COS_MAX), and the largest of the three values (or just the endpoint values
  // if cth is outside the allowed range) is chosen as the maximum.
  double me2_over_s = md_*md_ / s;
  double B = marley_utils::ONE_HALF * std::pow(g2_*(1. - me2_over_s), 2);
  double A = g1_*g2_*me2_over_s + g2_*g2_*(1. - me2_over_s) - B;
  double cth = -A / B;

  // Set the differential cross section to a huge negative value
  // at cth. This value will be compared to those at the angular
  // endpoints if cth does not lie in the allowed range. Otherwise,
  // the correct value will replace this one.
  double dxs_at_cth = std::numeric_limits<double>::lowest();

  double max = 0.;
  if ( cth >= COS_MIN && cth <= COS_MAX ) {
    dxs_at_cth = this->diff_xs(pdg_a, KEa, cth);
  }

  double dxs_at_min = this->diff_xs(pdg_a, KEa, COS_MIN);
  double dxs_at_max = this->diff_xs(pdg_a, KEa, COS_MAX);

  // Find the maximum value of the differential cross section
  max = std::max( { dxs_at_min, dxs_at_max, dxs_at_cth } );

  // Sample a CM frame scattering cosine for the ejectile.
  double cos_theta_c_cm = gen.rejection_sample(
    [this, pdg_a, KEa](double ctheta) -> double
    { return this->diff_xs(pdg_a, KEa, ctheta); }, COS_MIN, COS_MAX, max);

  // Sample a CM frame azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the differential cross section is independent of
  // the azimuthal angle.
  double phi_c_cm = gen.uniform_random_double(0., 2.*marley_utils::pi, false);

  // Create and return the completed event object
  return make_event_object(KEa, pc_cm, cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm);
}
