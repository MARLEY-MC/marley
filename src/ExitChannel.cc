/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/COPYING or
// visit http://opensource.org/licenses/GPL-3.0

#include "marley/marley_utils.hh"
#include "marley/ExitChannel.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/Logger.hh"
#include "marley/OpticalModel.hh"
#include "marley/marley_kinematics.hh"

namespace {

  // Helper function that issues a warning when the final nuclear excitation
  // energy exceeds the energetically accessible maximum value
  void issue_Exf_warning( double Exf, double Exf_max ) {
    MARLEY_LOG_WARNING() << "Final nuclear excitation energy Exf = "
      << Exf << " MeV exceeds the maximum accessible value of " << Exf_max
      << " MeV. The decay width for this exit channel will be set to zero.";
  }

  // Helper function that issues a warning when the final nuclear excitation
  // energy lies outside the accessible continuum
  void issue_Exf_continuum_warning( double Exf, double E_c_min,
    double E_c_max )
  {
    MARLEY_LOG_WARNING() << "Final nuclear excitation energy Exf = "
      << Exf << " MeV lies outside the accessible continuum [ " << E_c_min
      << " MeV, " << E_c_max << " MeV ]. The differential decay width"
      << " will be set to zero.";
  }

  // Helper function that throws an error in response to invalid excitation
  // energy continuum bounds
  void throw_continuum_bounds_error( double E_c_min, double E_c_max ) {
    throw marley::Error( "Invalid continuum excitation energy bounds "
      " E_c_min = " + std::to_string(E_c_min) + " MeV and E_c_max = "
      + std::to_string(E_c_max) + " MeV encountered" );
  }

}

using TrType = marley::GammaStrengthFunctionModel::TransitionType;

int marley::FragmentExitChannel::final_nucleus_pdg() const {
  int Zi = marley_utils::get_particle_Z( pdgi_ );
  int Ai = marley_utils::get_particle_A( pdgi_ );

  int Zfrag = marley_utils::get_particle_Z( fragment_pdg_ );
  int Afrag = marley_utils::get_particle_A( fragment_pdg_ );

  int Zf = Zi - Zfrag;
  int Af = Ai - Afrag;

  int pdgf = marley_utils::get_nucleus_pid( Zf, Af );
  return pdgf;
}

double marley::FragmentExitChannel::max_Exf() const {
  const auto& mt = marley::MassTable::Instance();
  double Sa = mt.get_fragment_separation_energy( pdgi_, fragment_pdg_ );
  double Exf_max = Exi_ - Sa;
  return Exf_max;
}

void marley::FragmentDiscreteExitChannel::compute_total_width() {

  int remnant_pdg = this->final_nucleus_pdg();
  marley::OpticalModel& om = sdb_.get_optical_model( remnant_pdg );

  // Get information about the emitted fragment
  const marley::Fragment& f = *sdb_.get_fragment( fragment_pdg_ );

  int two_s = f.get_two_s(); // two times the fragment spin
  marley::Parity Pa = f.get_parity(); // intrinsic parity

  // Maximum possible excitation energy in the daughter nucleus after the
  // fragment is emitted
  double Exf_max = this->max_Exf();

  // Excitation energy of the final nuclear level
  double Exf = final_level_.energy();

  // Initialize the total decay width to zero, just in case
  width_ = 0.;

  if ( Exf >= Exf_max ) {
    issue_Exf_warning( Exf, Exf_max );
    return;
  }

  // Total kinetic energy in the CM frame immediately after the binary decay
  double total_KE_CM_frame = Exf_max - Exf;

  int twoJf = final_level_.twoJ();
  marley::Parity Pf = final_level_.parity();

  for (int two_j = std::abs(twoJi_ - twoJf); two_j <= twoJi_ + twoJf;
    two_j += 2)
  {
    int j_plus_s = (two_j + two_s) / 2;
    // TODO: consider adding a check that l does not exceed its maximum
    // value l_max that is also used as a cutoff for transmission
    // coefficient calculations in the continuum.
    // For each new iteration, increment l and flip the overall final
    // state parity
    int l = std::abs(two_j - two_s) / 2;
    bool l_is_odd = l % 2;
    marley::Parity P_final_state = Pf * Pa * marley::Parity( !l_is_odd );
    for ( ; l <= j_plus_s; ++l, !P_final_state )
    {
      // The current term in the sum only contributes to the total decay
      // width if parity is conserved
      if ( Pi_ == P_final_state ) {

        double Tlj = om.transmission_coefficient( total_KE_CM_frame,
          fragment_pdg_, two_j, l, two_s );

        double partial_width = one_over_two_pi_rho_i_ * Tlj;

        width_ += partial_width;

        // TODO: cache term indexed by l, two_j
      }
    }
  }
}

double marley::GammaExitChannel::gamma_energy( double Exf ) const {
  // Approximate the gamma energy by the excitation energy difference between
  // the initial and final nuclear states
  // TODO: consider adding a nuclear recoil correction here
  double E_gamma = Exi_ - Exf;
  return E_gamma;
}


TrType marley::GammaExitChannel::get_transition_type( int mpol,
  marley::Parity Pf ) const
{
  // Electric transitions satisfy the parity relation Pi = (-1)^{\ell} * Pf,
  // where Pi (Pf) are the initial (final) nuclear parities and \ell is the
  // multipolarity of the transition. Magnetic transitions satisfy
  // Pi = (-1)^{\ell + 1} * Pf. We choose the appropriate transition type
  // here based on these rules.
  bool mpol_is_odd = mpol % 2;
  marley::Parity P_final_state = Pf * marley::Parity( !mpol_is_odd );
  TrType type = ( Pi_ == P_final_state ) ? TrType::electric : TrType::magnetic;
  return type;
}

void marley::GammaDiscreteExitChannel::compute_total_width() {

  // Retrieve the gamma strength function model used to compute transmission
  // coefficients
  marley::GammaStrengthFunctionModel& gsfm
    = sdb_.get_gamma_strength_function_model( pdgi_ );

  // Get properties of the final nuclear level
  double Exf = final_level_.energy();
  int twoJf = final_level_.twoJ();
  marley::Parity Pf = final_level_.parity();

  // Initialize the total width to zero, just in case.
  width_ = 0.;

  if ( Exf >= Exi_ ) {
    issue_Exf_warning( Exf, Exi_ );
    return;
  }

  // Compute the energy of the emitted gamma-ray
  double E_gamma = this->gamma_energy( Exf );

  // Initialize the multipolarity for the start of the loop.
  // There is no monopole radiation, so the minimum allowed value is one.
  int start_mpol = std::max( 1, std::abs(twoJi_ - twoJf) / 2 );
  int end_mpol = ( twoJi_ + twoJf ) / 2;

  // Sum contributions to the decay width from different multipolarities
  for ( int mpol = start_mpol; mpol <= end_mpol; ++mpol ) {

    // Use the multipolarity and final-state nuclear parity to determine
    // whether the current partial width represents an electric or magnetic
    // transition
    TrType type = this->get_transition_type( mpol, Pf );

    double TXl = gsfm.transmission_coefficient( type, mpol, E_gamma );

    double term = one_over_two_pi_rho_i_ * TXl;

    // TODO: add caching of terms by Xl pair

    width_ += term;
  }
}

double marley::FragmentContinuumExitChannel::differential_width( double Exf,
  bool store_jpi_widths ) const
{
  if ( store_jpi_widths ) jpi_widths_table_.clear();

  int remnant_pdg = this->final_nucleus_pdg();
  marley::OpticalModel& om = sdb_.get_optical_model( remnant_pdg );
  marley::LevelDensityModel& ldm = sdb_.get_level_density_model( remnant_pdg );

  // Get the maximum accessible final excitation energy
  double Exf_max = this->max_Exf();

  // Initialize the return value to zero
  double diff_width = 0.;

  if ( Exf_max < E_c_min_ ) throw_continuum_bounds_error( E_c_min_, Exf_max );

  // Check that Exf lies within the continuum
  if ( Exf < E_c_min_ || Exf > Exf_max ) {
    // If it doesn't, complain and return zero
    issue_Exf_continuum_warning( Exf, E_c_min_, Exf_max );
    return 0.;
  }

  double total_KE_CM_frame = Exf_max - Exf;

  // Get information about the emitted fragment
  const marley::Fragment& f = *sdb_.get_fragment( fragment_pdg_ );

  int two_s = f.get_two_s(); // two times the fragment spin
  marley::Parity Pa = f.get_parity(); // intrinsic parity

  // Final nuclear parity
  marley::Parity Pf;
  // The orbital parity starts as (-1)^0 = 1. Rather than applying parity
  // conservation each time, just find the final state parity Pf for l = 0.
  // Then we can safely flip Pf without further thought for each new l value in
  // the loop.
  if (Pi_ == Pa) Pf = 1;
  else Pf = -1;
  // For each new iteration, increment l and flip the final-state parity
  for (int l = 0; l <= marley::HauserFeshbachDecay::l_max_; ++l, !Pf) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      for (int twoJf = std::abs(twoJi_ - two_j);
        twoJf <= twoJi_ + two_j; twoJf += 2)
      {
        double Tlj = om.transmission_coefficient( total_KE_CM_frame,
          fragment_pdg_, two_j, l, two_s );

        double rho_f = ldm.level_density( Exf, twoJf, Pf );

        double term = one_over_two_pi_rho_i_ * Tlj * rho_f;

        diff_width += term;

        if ( store_jpi_widths ) {
          jpi_widths_table_.emplace_back( twoJf, Pf, term );
          // TODO: include (l, two_j) in cached term
        }
      }
    }
  }
  return diff_width;
}

void marley::ContinuumExitChannel::compute_total_width() {

  // Initialize the total width to zero
  width_ = 0.;

  double Ec_max = this->E_c_max();

  if ( Ec_max < E_c_min_ ) {
    throw_continuum_bounds_error( E_c_min_, Ec_max );
    return;
  }

  // Create a function object to use for integration of the differential decay
  // width
  std::function<double(double)> dw = [this](double Exf) -> double {
    return this->differential_width( Exf );
  };

  // Numerically integrate over the bounds of the continuum using the
  // function object prepared above
  width_ = marley_utils::num_integrate( dw, E_c_min_, Ec_max );

  // TODO: consider switching to doing the integration with a
  // ChebyshevInterpolatingFunction object. This avoids needing to create one
  // later (and may be comparable in terms of computational cost)
}

double marley::GammaContinuumExitChannel::differential_width( double Exf,
  bool store_jpi_widths ) const
{
  if ( store_jpi_widths ) jpi_widths_table_.clear();

  auto& ldm = sdb_.get_level_density_model( pdgi_ );
  auto& gsfm = sdb_.get_gamma_strength_function_model( pdgi_ );

  // Initialize the return value to zero
  double diff_width = 0.;

  // Check that the continuum bounds make sense
  if ( Exi_ < E_c_min_ ) throw_continuum_bounds_error( E_c_min_, Exi_ );

  // Check that Exf lies within the continuum
  if ( Exf < E_c_min_ || Exf > Exi_ ) {
    // If it doesn't, complain and return zero
    issue_Exf_continuum_warning( Exf, E_c_min_, Exi_ );
    return 0.;
  }

  // Compute the energy of the emitted gamma-ray
  double E_gamma = this->gamma_energy( Exf );

  // Array containing both possible parity values. It is used in
  // the loop below.
  constexpr std::array<marley::Parity, 2>
    parities = { marley::Parity(true), marley::Parity(false) };

  // Sum over multipolarities. There is no monopole radiation, so
  // the sum begins at mpol = 1.
  for ( int mpol = 1; mpol <= marley::HauserFeshbachDecay::l_max_; ++mpol ) {

    int two_mpol = 2 * mpol;

    for ( int twoJf = std::abs(twoJi_ - two_mpol); twoJf <= twoJi_ + two_mpol;
      twoJf += 2 )
    {
      for ( const auto& Pf : parities ) {

        // Use the multipolarity and final-state nuclear parity to determine
        // whether the current partial differential width represents an
        // electric or magnetic transition
        TrType type = this->get_transition_type( mpol, Pf );

        double Txl = gsfm.transmission_coefficient( type, mpol, E_gamma );
        double rho_f = ldm.level_density(Exf, twoJf, Pf);

        double term = one_over_two_pi_rho_i_ * Txl * rho_f;

        if ( store_jpi_widths ) {
          jpi_widths_table_.emplace_back( twoJf, Pf, term );
          // TODO: include Xl in cached values
        }

        diff_width += term;
      }
    }
  }
  return diff_width;
}

void marley::DiscreteExitChannel::do_decay(double& Exf, int& two_Jf,
  marley::Parity& Pf, const marley::Particle& compound_nucleus,
  marley::Particle& emitted_particle, marley::Particle& residual_nucleus,
  marley::Generator& gen) const
{
  Exf = final_level_.energy();
  two_Jf = final_level_.twoJ();
  Pf = final_level_.parity();
  this->prepare_products( compound_nucleus, emitted_particle,
    residual_nucleus, Exf, gen );
}

void marley::ExitChannel::prepare_products(
  const marley::Particle& compound_nucleus, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, double Exf, marley::Generator& gen) const
{
  const auto& mt = marley::MassTable::Instance();
  int ep_pdg = this->emitted_particle_pdg();
  double ep_mass = mt.get_particle_mass( ep_pdg );

  emitted_particle = marley::Particle( ep_pdg, ep_mass );

  // Proton number of the emitted particle
  int ep_Z = marley_utils::get_particle_Z( ep_pdg );

  // Final ion charge after particle emission
  int qf = qi_ - ep_Z;

  double me = mt.get_particle_mass( marley_utils::ELECTRON );
  int remnant_pdg = this->final_nucleus_pdg();

  // Approximate the ground state mass of the ion formed when the fragment is
  // emitted by subtracting qf electron masses from the atomic mass for the
  // final nuclide.
  double Mfgs_ion = mt.get_atomic_mass( remnant_pdg ) - qf*me;

  residual_nucleus = marley::Particle( remnant_pdg, Mfgs_ion + Exf, qf );

  // Now that the PDG codes, masses, and net charges of the binary decay
  // products have been set, choose a direction for the emitted particle.
  // TODO: Consider changing this to a more realistic model instead of
  // isotropic emissions.
  double cos_theta_emitted_particle
    = gen.uniform_random_double( -1., 1., true );
  double phi_emitted_particle
    = gen.uniform_random_double( 0., marley_utils::two_pi, false );

  // Handle the kinematics calculations for this decay. Load the
  // final-state particle objects with their full 4-momenta.
  marley_kinematics::two_body_decay( compound_nucleus, emitted_particle,
    residual_nucleus, cos_theta_emitted_particle, phi_emitted_particle );
}

double marley::ContinuumExitChannel::sample_Exf(marley::Generator& gen) const
{
  // The maximum accessible excitation energy for this exit channel. It
  // will be used when creating the ChebyshevInterpolatingFunction below
  double Emax = this->E_c_max();

  // If we haven't built a CDF for sampling the final nuclear excitation
  // energy yet, then build it before continuing
  if ( !Exf_cdf_ ) {
    // Build a polynomial approximant (at Chebyshev points) to the PDF for the
    // final nuclear excitation energy
    marley::ChebyshevInterpolatingFunction pdf_cheb( [this](double Exf)
      -> double { return this->differential_width(Exf); }, E_c_min_, Emax,
      marley::DEFAULT_N_CHEBYSHEV );

    // Store the cumulative density function for possible re-use
    Exf_cdf_ = std::make_unique<marley::ChebyshevInterpolatingFunction>(
      pdf_cheb.cdf() );
  }

  // Sample a final nuclear excitation energy using the Chebyshev polynomial
  // approximant to the CDF
  double Exf = gen.inverse_transform_sample( *Exf_cdf_, E_c_min_, Emax );
  return Exf;
}

void marley::ContinuumExitChannel::do_decay(double& Exf, int& two_Jf,
  marley::Parity& Pf, const marley::Particle& compound_nucleus,
  marley::Particle& emitted_particle, marley::Particle& residual_nucleus,
  marley::Generator& gen) const
{
  Exf = this->sample_Exf( gen );

  // Sample a final nuclear spin-parity, unless the user has
  // explicitly turned this off (presumably in unit tests
  // where we only care about the final excitation energy).
  if ( !skip_jpi_sampling_ ) sample_spin_parity(Exf, two_Jf, Pf, gen);

  // TODO: sample a sub-Jpi variable set (e.g., l + j)

  this->prepare_products( compound_nucleus, emitted_particle,
    residual_nucleus, Exf, gen );
}

void marley::ContinuumExitChannel::sample_spin_parity(double Exf, int& twoJ,
  marley::Parity& Pi, marley::Generator& gen) const
{
  // Clear any previous table entries of spin-parities and decay widths
  jpi_widths_table_.clear();

  // Load table of partial differential widths via a call to
  // differential_width()
  double diff_width = this->differential_width( Exf, true );

  // Throw an error if all decays are impossible
  if ( diff_width <= 0. ) throw marley::Error( "Cannot "
    "continue Hauser-Feshbach decay. All partial differential decay widths"
    "are zero." );

  // Sample a final spin and parity
  const auto begin = marley::IteratorToMember<
    std::vector<SpinParityWidth>::const_iterator,
    const double>( jpi_widths_table_.cbegin(),
    &SpinParityWidth::diff_width );

  const auto end = marley::IteratorToMember<
    std::vector<SpinParityWidth>::const_iterator,
    const double>( jpi_widths_table_.cend(),
    &SpinParityWidth::diff_width );

  std::discrete_distribution<size_t> jpi_dist( begin, end );
  size_t jpi_index = gen.sample_from_distribution( jpi_dist );

  // Store the results
  const SpinParityWidth& Jpi = jpi_widths_table_.at( jpi_index );
  twoJ = Jpi.twoJf;
  Pi = Jpi.Pf;
}
