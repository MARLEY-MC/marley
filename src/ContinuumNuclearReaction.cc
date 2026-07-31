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
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/MassTable.hh"
#include "marley/ContinuumNuclearReaction.hh"
#include "marley/marley_utils.hh"
#include "marley/hepmc3_utils.hh"

using SubContinuumMode = marley::ContinuumNuclearReaction::SubContinuumMode;

// Initialize the static class member defining how to deal with sub-continuum
// cross-section strength. By default, we pick the "accumulate" option.
SubContinuumMode marley::ContinuumNuclearReaction::sc_mode_
  = SubContinuumMode::ACCUMULATE;

// Initialize the static class member specifying a mapping between strings
// and SubContinuumMode enum values
const std::map< SubContinuumMode, std::string >
  marley::ContinuumNuclearReaction::sc_mode_string_map_ = {
    { SubContinuumMode::IGNORE, "ignore" },
    { SubContinuumMode::MIRROR, "mirror" },
    { SubContinuumMode::ACCUMULATE, "accumulate" }
};

marley::ContinuumNuclearReaction::ContinuumNuclearReaction(
  Reaction::ProcessType pt, int pdg_a, int pdg_b, int pdg_c, int pdg_d,
  int q_d, const std::shared_ptr<TabulatedXSec>& txsec,
  const std::string& source_file )
  : marley::NuclearReaction( pt, pdg_a, pdg_b, pdg_c, pdg_d, q_d,
      source_file ),
  xsec_( txsec )
{
}

double marley::ContinuumNuclearReaction::total_xs( int pdg_a,
  double KEa ) const
{
  if ( pdg_a != pdg_a_ ) return 0.;
  return xsec_->integral( pdg_a, KEa );
}

std::shared_ptr< HepMC3::GenEvent > marley::ContinuumNuclearReaction
  ::create_event( int pdg_a, double KEa, marley::Generator& gen ) const
{
  // TODO: reduce code duplication here with DiscreteNuclearReaction using the
  // common base class NuclearReaction

  // Check that the projectile supplied to this event is correct. If not, alert
  // the user that this event does not use the requested projectile.
  if ( pdg_a != pdg_a_ ) throw marley::Error( "Could not create this event."
    " The requested projectile particle ID, " + std::to_string( pdg_a )
    + ", does not match the projectile particle ID, "
    + std::to_string( pdg_a_ ) + ", in the reaction dataset." );

  // Sample a final residue energy level. First, check to make sure the given
  // projectile energy is above threshold for this reaction.
  if ( KEa < KEa_threshold_ ) throw std::range_error( "Could"
    " not create this event. Projectile kinetic energy "
    + std::to_string( KEa ) + " MeV is below the threshold value "
    + std::to_string( KEa_threshold_ ) + " MeV." );

  // Select a specific multipole to use for the current event using the
  // individual total cross sections
  std::vector< double > multipole_weights;
  std::vector< marley::TabulatedXSec::MultipoleLabel > multipoles;
  std::vector< double > diff_max_values;
  const auto& table_map = xsec_->get_table_map();
  double sum_of_xsecs = 0.;
  for ( const auto& pair : table_map ) {
    const auto& ml = pair.first;
    double diff_max;
    double total_xsec = xsec_->integral( pdg_a_, KEa, ml, diff_max );

    sum_of_xsecs += total_xsec;

    multipole_weights.push_back( total_xsec );
    multipoles.push_back( ml );
    diff_max_values.push_back( diff_max );
  }

  // If there are no multipole weights, we can't go on. Complain if this
  // is the case.
  if ( multipole_weights.empty() ) {
    throw marley::Error( "Could not create this event. The TabulatedXSec"
      " object associated with this reaction does not own any nuclear response"
      " tables." );
  }

  // Complain if the total cross section (the sum of all partial cross
  // sections) is zero or negative (the latter is just to cover all
  // possibilities).
  if ( sum_of_xsecs <= 0. ) {
    throw marley::Error( "Could not create this event. All multipole total"
      " cross sections are nonpositive." );
  }

  // Create a discrete distribution based on the weights. This will be
  // used to choose a single multipole for the current event.
  std::discrete_distribution< size_t > multipole_dist(
    multipole_weights.begin(), multipole_weights.end() );

  // Sample a matrix_element using our discrete distribution and the
  // current set of weights
  size_t multipole_index = gen.sample_from_distribution( multipole_dist );

  // Label of the multipole chosen for this event
  const auto& sampled_ml = multipoles.at( multipole_index );
  // Maximum value of the differential cross section for this multipole.
  // This will be used for rejection sampling of inclusive kinematics below.
  double diff_max = diff_max_values.at( multipole_index );

  // ResponseTable object to use for computing the differential cross section
  // during kinematic sampling below
  const auto& rt = xsec_->get_table( sampled_ml );

  // Get the values of the energy transfer that correspond to the edges
  // of the table of nuclear responses. Note that the table is actually given
  // in terms of the effective energy transfer, which differs by delta_ias
  // from the actual energy transfer. We therefore apply a shift here
  // to correct for this.
  double table_wmin = rt.w_min() - xsec_->delta_ias();
  double table_wmax = rt.w_max() - xsec_->delta_ias();

  // Choose a reasonable sampling interval for the energy transfer
  double Ea = KEa + ma_; // Projectile total energy

  // Set the lower bound for the energy transfer to be either zero or
  // the lowest tabulated value (whichever is larger)
  double wmin = std::max( 0., table_wmin );

  // Set the upper bound for the energy transfer to be either the projectile
  // energy minus the final lepton mass or the highest tabulated value
  // (whichever is lower)
  double wmax = std::min( Ea - mc_, table_wmax );

  // Sample values for the energy transfer and scattering cosine using the
  // differential cross section for the chosen multipole.
  // Use a simple rejection sampling technique.
  double w, ctl, diff, y;
  int sampling_attempts = 0;
  bool recomputed_diff_max = false;
  do {
    // Occasionally the value of diff_max retrieved above can be extremely
    // overestimated when using an optimized version of the total cross
    // section calculation (which relies on interpolation). This typically
    // happens very close to threshold and leads to sampling getting stuck
    // due to a very low acceptance efficiency.
    //
    // To guard against this situation, when the number of sampling attempts
    // exceeds a large value, the estimate of the maximum differential cross
    // section diff_max is recalculated at exactly the input projectile kinetic
    // energy rather than relying on the precomputed value. The updated
    // estimate is then used in a new set of sampling attempts.
    if ( sampling_attempts > marley_utils::LARGE_NUMBER_OF_ITERATIONS ) {

      if ( recomputed_diff_max ) {
        // If we make it here, then we are still hitting a huge number of
        // iterations in this sampling loop despite recalculating the maximum
        // differential cross section. This suggests that we are stuck in an
        // infinite loop, so abort with an exception indicating the problem.
        throw marley::Error( "Reached maximum iteration count for rejection"
          " sampling in marley::ContinuumNuclearReaction::create_event()" );
      }

      // The value of diff_max is updated by this call to
      // marley::TabulatedXSec::compute_integral()
      xsec_->compute_integral( pdg_a_, KEa, sampled_ml, diff_max );
      sampling_attempts = 0;
      recomputed_diff_max = true;
    }

    w = gen.uniform_random_double( wmin, wmax, true );
    ctl = gen.uniform_random_double( -1., 1., true );
    diff = xsec_->diff_xsec( pdg_a_, KEa, w, ctl, sampled_ml );
    y = gen.uniform_random_double( 0., diff_max, true );
    ++sampling_attempts;

    // If reassignment is enabled and the excitation energy corresponding to the
    // sampled energy transfer w falls below the continuum, reassign the value
    // of the energy transfer to lie within the continuum. The scattering cosine
    // ctl and projectile kinetic energy are needed for the calculation but are
    // never altered. If reassignment is attempted but fails due to kinematic
    // limits, then force another iteration of this rejection sampling loop even
    // if the event would have been accepted without reassignment.
    //
    // NOTE: The reassignment operation is included in the condition of the
    // do-while loop for efficiency. There is no need to perform the
    // reassignment for w values that would be rejected anyway. Use of the
    // logical OR operation (||) will only evaluate the first condition if it
    // is false, thus skipping the reassignment when it is obviously
    // unnecessary.
  } while ( y > diff || !this->reassign_sub_continuum(w, ctl, KEa) );

  // Sample a lab-frame azimuthal scattering angle uniformly
  double phi_c = gen.uniform_random_double( 0., marley_utils::two_pi, false );

  // Load the initial residue twoJ and parity values into twoJ and P. These
  // variables are included in the event record and used by NucleusDecayer to
  // start the Hauser-Feshbach decay cascade.
  // NOTE: right now, these are taken directly from the multipole involved in
  // the current event. This is only valid for scattering on a 0+ target
  // nucleus
  // TODO: revisit this assumption and do something better
  int twoJ = 2 * sampled_ml.J_;
  marley::Parity P = sampled_ml.Pi_; // defaults to positive parity

  // Sine of the ejectile scattering angle
  double stl = marley_utils::real_sqrt( 1. - std::pow(ctl, 2) );

  // Calculate the full kinematics of the primary interaction based on the
  // lepton scattering cosine (ctl) and energy transfer (w) sampled above.

  // Determine the components of the ejectile's lab-frame 4-momentum
  double Ec = Ea - w;
  double pc = marley_utils::real_sqrt( Ec*Ec - mc_*mc_ );
  double pc_x = stl * std::cos( phi_c ) * pc;
  double pc_y = stl * std::sin( phi_c ) * pc;
  double pc_z = ctl * pc;

  // Determine the magnitude of the lab-frame 3-momentum of the projectile
  double pa = marley_utils::real_sqrt( Ea*Ea - ma_*ma_ );

  // Construct the lab-frame 4-momenta of the projectile, target, and ejectile
  HepMC3::FourVector pro_mom4( 0., 0., pa, Ea );
  HepMC3::FourVector tar_mom4( 0., 0., 0., mb_ );
  HepMC3::FourVector eje_mom4( pc_x, pc_y, pc_z, Ec );

  // Get the 4-momentum of the residue in the lab frame using conservation
  double Ed = pro_mom4.e() + tar_mom4.e() - eje_mom4.e();
  double pd_x = pro_mom4.px() + tar_mom4.px() - eje_mom4.px();
  double pd_y = pro_mom4.py() + tar_mom4.py() - eje_mom4.py();
  double pd_z = pro_mom4.pz() + tar_mom4.pz() - eje_mom4.pz();

  // Determine the residue mass from its 4-momentum
  md_ = marley_utils::real_sqrt( Ed*Ed - pd_x*pd_x - pd_y*pd_y - pd_z*pd_z );

  // The excitation energy is the mass difference between this mass and
  // the residue's ground-state mass
  double Ex = md_ - md_gs_;

  // Create particle objects representing the ejectile and residue
  auto ejectile = marley_hepmc3::make_particle( eje_mom4, pdg_c_,
    marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, mc_ );

  auto residue = marley_hepmc3::make_particle( pdg_d_, pd_x, pd_y, pd_z, Ed,
    marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, md_ );

  // Make the event object (this also sets the charge and nuclear level
  // attributes)
  auto event = this->make_nuclear_event_object( KEa, ejectile, residue, Ex,
    twoJ, P );

  return event;
}

// Adds an indication to the description that the daughter nucleus will always
// be left in an excited state.
// TODO: revisit this as needed. I assume here that the
// ContinuumNuclearReaction class will always be used for calculations in the
// unbound continuum of nuclear levels
void marley::ContinuumNuclearReaction::set_description() {
  marley::NuclearReaction::set_description();
  description_ += '*';
}

// Implements a requirement that the energy transfer sampled in create_event()
// lies within the excitation energy continuum. According to the user
// configuration, this requirement can be turned on or off and defined in
// different ways.
bool marley::ContinuumNuclearReaction::reassign_sub_continuum( double& w,
  const double ctl, const double KEa ) const
{
  // If reassignment of the sub-continuum strength is disabled, then this
  // function returns immediately without doing anything
  if ( sc_mode_ == SubContinuumMode::IGNORE ) return true;

  // Lab-frame total energy and 3-momentum of the projectile
  double Ea = KEa + ma_;
  double pa = marley_utils::real_sqrt( Ea*Ea - ma_*ma_ );

  // Lab-frame total energy and 3-momentum of the ejectile
  double Ec = Ea - w;
  double pc = marley_utils::real_sqrt( Ec*Ec - mc_*mc_ );

  // Squared magnitude of the 3-momentum transfer
  double kappa2 = pa*pa + pc*pc - 2.*pa*pc*ctl;

  // Total energy and excitation energy of the residue
  double Ed = mb_ + w;
  double Ex = marley_utils::real_sqrt( Ed*Ed - kappa2 ) - md_gs_;

  // Get the "unbound threshold" used to determine the start of the continuum
  const auto& mt = marley::MassTable::Instance();
  double unbound_threshold = mt.unbound_threshold( pdg_d_ );

  // If the excitation energy corresponding to the sampled energy transfer
  // is already within the continuum, no special action is needed. Just return
  // without making any modifications.
  if ( Ex >= unbound_threshold ) return true;

  // If we've made it here, then we need to reassign the excitation energy
  // and compute a new value of the energy transfer. First choose the new
  // excitation energy based on the recipe selected by the user configuration.
  if ( sc_mode_ == SubContinuumMode::ACCUMULATE ) {
    // For the "accumulate" option, just update the excitation energy to be
    // exactly at the unbound threshold
    Ex = unbound_threshold;

    MARLEY_LOG( DEBUG, "physics.reaction" ) << "Excitation energy " << Ex
      << " MeV is below the unbound threshold " << unbound_threshold
      << " MeV. Sampling exactly at the unbound threshold.";
  }
  else if ( sc_mode_ == SubContinuumMode::MIRROR ) {
    // For the "mirror" option, "reflect" the original excitation energy to the
    // upper side of the unbound threshold so it is the same distance above
    // as it was originally below.
    Ex = 2.*unbound_threshold - Ex;

    MARLEY_LOG( DEBUG, "physics.reaction" ) << "Excitation energy " << Ex
      << " MeV is below the unbound threshold " << unbound_threshold
      << " MeV. Mirroring the energy transfer around the unbound"
      << " threshold.";
  }
  else {
    throw marley::Error( "Unrecognized sub-continuum mode encountered"
      " in marley::ContinuumNuclearReaction::reassign_sub_continuum()" );
    return false;
  }

  // Solve for the new outgoing lepton total energy given the updated
  // excitation energy value. Update the value of Ec with the solution.
  Ec = this->get_Ec_from_Ex( Ex, ctl, KEa );

  // Now update the energy transfer accordingly
  w = Ea - Ec;

  // If the total energy of the outgoing lepton is now below its rest mass,
  // then the reassignment procedure failed due to the kinematic threshold.
  // Indicate this failure in the return value.
  if ( Ec < mc_ ) return false;

  // Otherwise, everything worked out, so indicate success.
  return true;
}

double marley::ContinuumNuclearReaction::get_Ec_from_Ex( const double Ex,
  const double cos_theta, const double KEa, double* jacobian ) const
{
  // Mass of the final-state ion (including excitation energy)
  double md = md_gs_ + Ex;

  // Total energy of the projectile
  double Ea = KEa + ma_;

  // Total energy of the two-body system in the lab frame
  double Etot = Ea + mb_;

  // Construct helper variables
  double pa = marley_utils::real_sqrt( Ea*Ea - ma_*ma_ );
  double help = md*md - mc_*mc_ + pa*pa - Etot*Etot;
  double other_help = 4.*pa*pa*cos_theta*cos_theta;

  // Quadratic coefficients (a*Ec^2 + b*Ec + c == 0)
  double a = 4.*Etot*Etot - other_help;
  double b = 4.*Etot*help;
  double c = help*help + other_help*mc_*mc_;

  // Get both solutions to the quadratic equation
  double sol_plus, sol_minus;
  marley_utils::solve_quadratic_equation( a, b, c,
    sol_plus, sol_minus );

  // Now for a trick: due to the way we derived the results above,
  // the two solutions correspond to positive (sol_plus) and
  // negative (sol_minus) values of cos_theta, with the two
  // solutions exactly equal when cos_theta == 0. Choose the
  // appropriate one to return here.
  double Ec = sol_plus;
  if ( cos_theta < 0. ) Ec = sol_minus;

  if ( jacobian ) {
    double pc = marley_utils::real_sqrt( Ec*Ec - mc_*mc_ );
    *jacobian = md / ( Ea + mb_ - pa*Ec*cos_theta/pc );
  }

  return Ec;
}

// Convert a string to a SubContinuumMode value
SubContinuumMode marley::ContinuumNuclearReaction
  ::sub_continuum_mode_from_string( const std::string& str )
{
  for ( const auto& pair : sc_mode_string_map_ ) {
    if ( str == pair.second ) return pair.first;
  }
  throw marley::Error( "The string \"" + str + "\" was not recognized"
    " as a valid sub-continuum mode setting" );
}

// Convert a SubContinuumMode value to a string
std::string marley::ContinuumNuclearReaction::string_from_sub_continuum_mode(
  SubContinuumMode mode )
{
  auto it = sc_mode_string_map_.find( mode );
  if ( it != sc_mode_string_map_.end() ) return it->second;
  else throw marley::Error( "Unrecognized sub-continuum mode value encountered"
    " in marley::ContinuumNuclearReaction::string_from_sub_continuum_mode()" );
}
