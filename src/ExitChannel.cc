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
#include "marley/OpticalModel.hh"

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
  const marley::Fragment& f = *marley::StructureDatabase
    ::get_fragment( fragment_pdg_ );

  int two_s = f.get_two_s(); // two times the fragment spin
  marley::Parity Pa = f.get_parity(); // intrinsic parity

  // Maximum possible excitation energy in the daughter nucleus after the
  // fragment is emitted
  double Exf_max = this->max_Exf();

  // Excitation energy of the final nuclear level
  double Exf = final_level_.energy();

  if ( Exf >= Exf_max ) {
    // TODO: do warning
    width_ = 0.;
    return;
  }

  // Initialize the total decay width to zero, just in case
  width_ = 0.;

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
        double partial_width = one_over_two_pi_rho_i_
          * om.transmission_coefficient( total_KE_CM_frame, fragment_pdg_,
          two_j, l, two_s );

        width_ += partial_width;

        // TODO: cache term indexed by l, two_j
      }
    }
  }
}

void marley::GammaDiscreteExitChannel::compute_total_width() {

  marley::GammaStrengthFunctionModel& gsfm
    = sdb_.get_gamma_strength_function_model( pdgi_ );

  double Exf = final_level_.energy();
  if ( Exf >= Exi_ ) {
    // TODO: do warning
    width_ = 0.;
    return;
  }

  // TODO: refactor this to do the multipolarity sum right
  // TODO: add caching of terms by Xl pair
  width_ = one_over_two_pi_rho_i_ * gsfm.transmission_coefficient( Exi_,
    twoJi_, Pi_, final_level_ );
}

double marley::FragmentContinuumExitChannel::differential_width( double Exf,
  bool store_jpi_widths ) const
{
  if ( store_jpi_widths ) jpi_widths_table_.clear();

  int remnant_pdg = this->final_nucleus_pdg();
  marley::OpticalModel& om = sdb_.get_optical_model( remnant_pdg );
  marley::LevelDensityModel& ldm = sdb_.get_level_density_model( remnant_pdg );

  double Exf_max = this->max_Exf();

  // Check that Exf lies within the continuum
  if ( Exf < E_c_min_ || Exf > Exf_max ) {
    // If it doesn't just return zero
    // TODO: do warning
    // TODO: check that Exf_max > E_c_min
    return 0.;
  }

  double total_KE_CM_frame = Exf_max - Exf;

  // Initialize the return value to zero
  double diff_width = 0.;

  // Get information about the emitted fragment
  const marley::Fragment& f = *marley::StructureDatabase
    ::get_fragment( fragment_pdg_ );

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
        // TODO: since only the spin distribution changes in the twoJf loop,
        // you can optimize this by precomputing most of the level density
        // and the multiplying here by the appropriate spin distribution
        double term = one_over_two_pi_rho_i_ * om.transmission_coefficient(
          total_KE_CM_frame, fragment_pdg_, two_j, l, two_s )
          * ldm.level_density( Exf, twoJf, Pf );

        diff_width += term;

        if ( store_jpi_widths ) {
          jpi_widths_table_.emplace_back( twoJf, Pf, term );
          // TODO: optionally cache term by (l, two_j, twoJf)
        }
      }
    }
  }
  return diff_width;
}

void marley::ContinuumExitChannel::compute_total_width() {

  double Ec_max = this->E_c_max();
  if ( Ec_max <= E_c_min_ ) {
    // TODO: do warning
    width_ = 0.;
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
  using TrType = marley::GammaStrengthFunctionModel::TransitionType;

  if ( store_jpi_widths ) jpi_widths_table_.clear();

  marley::LevelDensityModel& ldm = sdb_.get_level_density_model( pdgi_ );
  auto& gsfm = sdb_.get_gamma_strength_function_model( pdgi_ );

  // Initialize the return value to zero
  double diff_width = 0.;

  // Check that Exf lies within the continuum. If it
  // doesn't, then just return zero.
  if ( Exf < E_c_min_ || Exf > Exi_ ) {
    // TODO: do warning
    return diff_width;
  }

  // Approximate the gamma energy by the energy difference between the two
  // levels
  // TODO: consider adding a nuclear recoil correction here
  double E_gamma = Exi_ - Exf;
  bool initial_spin_is_zero = twoJi_ == 0;

  // TODO: revisit and correct the sums here
  for (int l = 0; l <= marley::HauserFeshbachDecay::l_max_; ++l) {
    int two_l = 2*l;
    int twoJf = twoJi_ + two_l;
    // Determine the multipolarity being considered in this trip through the
    // loop based on the value of the orbital angular momentum l.
    int mpol = ( (l == 0) ? 1 : l );
    // Consider both transition types (equivalently, both final parities). Check
    // for Ji = 0 to Jf = 0 transitions, which are not allowed by angular
    // momentum conservation.
    if ( !initial_spin_is_zero || (twoJf > 0) ) {

      double tcE = gsfm.transmission_coefficient(TrType::electric, mpol,
        E_gamma);
      double tcM = gsfm.transmission_coefficient(TrType::magnetic, mpol,
        E_gamma);

      if ( !initial_spin_is_zero ) {

        diff_width += store_gamma_jpi_width(Exf, twoJf, Pi_, tcE, tcM,
          mpol, ldm, store_jpi_widths);

	// Consider the other possible final spin value for this multipolarity
	// if it is allowed (Jf = Ji - l is positive)
        twoJf = twoJi_ - two_l;
        if ( twoJf >= 0 ) diff_width += store_gamma_jpi_width(Exf, twoJf,
          Pi_, tcE, tcM, mpol, ldm, store_jpi_widths);
      }
      // twoJf > 0 but twoJi == 0
      else diff_width += store_gamma_jpi_width(Exf, twoJf, Pi_, tcE,
        tcM, mpol, ldm, store_jpi_widths);
    }
  }

  return diff_width;
}

void marley::DiscreteExitChannel::do_decay(double& Exf, int& two_Jf,
  marley::Parity& Pf, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
  const
{
  Exf = final_level_.energy();
  two_Jf = final_level_.twoJ();
  Pf = final_level_.parity();
  this->prepare_products( emitted_particle, residual_nucleus, Exf );
}

void marley::ExitChannel::prepare_products(
  marley::Particle& emitted_particle, marley::Particle& residual_nucleus,
  double Exf) const
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
  marley::Parity& Pf, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& gen)
  const
{
  Exf = this->sample_Exf( gen );

  // Sample a final nuclear spin-parity, unless the user has
  // explicitly turned this off (presumably in unit tests
  // where we only care about the final excitation energy).
  if ( !skip_jpi_sampling_ ) sample_spin_parity(Exf, two_Jf, Pf, gen);

  // TODO: sample a sub-jpi variable set (e.g., l + j)

  this->prepare_products( emitted_particle, residual_nucleus, Exf );
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

// Helper function used when sampling continuum spin-parities for gamma-ray
// transitions
double marley::GammaContinuumExitChannel::store_gamma_jpi_width(double Exf,
  int twoJf, marley::Parity Pi, double tcE, double tcM, int mpol,
  marley::LevelDensityModel& ldm, bool store_jpi_widths) const
{
  marley::Parity Pf;
  double combined_width = 0.;

  // Electric transitions represent a parity flip for odd multipolarities
  if (mpol % 2) Pf = -Pi;
  else Pf = Pi;
  double rho = ldm.level_density(Exf, twoJf, Pf);

  // Compute and store information for the electric transition
  double width = one_over_two_pi_rho_i_ * tcE * rho;
  combined_width += width;

  if ( store_jpi_widths ) jpi_widths_table_.emplace_back(twoJf, Pf, width);

  // Magnetic transitions have opposite final parity from electric
  // transitions of the same multipolarity
  !Pf;
  // Compute and store information for the magnetic transition
  rho = ldm.level_density(Exf, twoJf, Pf);
  width = one_over_two_pi_rho_i_ * tcM * rho;
  combined_width += width;

  if ( store_jpi_widths ) jpi_widths_table_.emplace_back(twoJf, Pf, width);

  return combined_width;
}
