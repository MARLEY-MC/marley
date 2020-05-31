/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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

#include <memory>

#include "marley/ExitChannel.hh"
#include "marley/Generator.hh"
#include "marley/MassTable.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/marley_utils.hh"

marley::HauserFeshbachDecay::HauserFeshbachDecay(const marley::Particle&
  compound_nucleus, double Exi, int twoJi, marley::Parity Pi,
  marley::StructureDatabase& sdb) : compound_nucleus_(compound_nucleus),
  Exi_(Exi), twoJi_(twoJi), Pi_(Pi)
{
  build_exit_channels( sdb );
}

void marley::HauserFeshbachDecay::build_exit_channels(
  marley::StructureDatabase& sdb)
{
  // Remove any pre-existing ExitChannel objects, just in case
  exit_channels_.clear();

  int pdgi = compound_nucleus_.pdg_code();
  int Zi = marley_utils::get_particle_Z( pdgi );
  int Ai = marley_utils::get_particle_A( pdgi );
  int qi = compound_nucleus_.charge(); // Get net charge of initial ion

  // Get the initial nuclear level density (MeV^{-1}) in the vicinity of the
  // initial nuclear level. This will be used to apply an overall normalization
  // factor when computing exit channel decay widths. This isn't strictly
  // needed for MC sampling, but it's helpful to work with physically
  // meaningful units when possible.
  marley::LevelDensityModel& ldm = sdb.get_level_density_model( Zi, Ai );
  double rho_i = ldm.level_density( Exi_, twoJi_, Pi_ );

  total_width_ = 0.; // total compound nucleus decay width

  for ( const auto& pair : sdb.fragments() ) {

    const marley::Fragment& f = pair.second;

    // Get information about the current fragment
    int two_s = f.get_two_s(); // spin
    marley::Parity Pa = f.get_parity();
    int fragment_pid = f.get_pid();
    int Za = f.get_Z(); // atomic number
    double Ma = f.get_mass(); // mass

    // Get information about the final-state nucleus
    int Zf = Zi - f.get_Z(); // atomic number
    int Af = Ai - f.get_A(); // mass number
    int pdg_final = marley_utils::get_nucleus_pid(Zf, Af);

    // Approximate the ground state mass of the ion formed when the fragment f
    // is emitted by adding (Za - qi) electron masses to the atomic mass for
    // the final nucleus.
    const auto& mt = marley::MassTable::Instance();
    double Sa = mt.get_fragment_separation_energy( Zi, Ai, fragment_pid );

    // Get discrete level data (if any) and models for the final nucleus
    marley::DecayScheme* ds = sdb.get_decay_scheme(pdg_final);

    // Determine the maximum excitation energy available after fragment
    // emission in the final nucleus. This is simply the difference
    // between the initial excitation energy and the fragment separation
    // energy.
    double Exf_max = Exi_ - Sa;

    // Check if emission of this fragment is energetically allowed. If we're
    // exactly at threshold, still refuse to emit the fragment to avoid
    // numerical problems.
    if ( Exf_max <= 0. ) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0.;

    // If discrete level data are available for the final nuclide, get decay
    // widths for each accessible level
    if (ds) {

      // Get a vector of pointers to levels in the decay scheme. The levels are
      // sorted in order of increasing excitation energy.
      const auto& levels = ds->get_levels();

      // Use the maximum discrete level energy from the decay scheme object as
      // the lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (levels.size() > 0) E_c_min = levels.back()->energy();

      // Loop over the final discrete nuclear levels in order of increasing
      // energy until the new level energy exceeds the maximum value. For each
      // energetically allowed level, if a transition to it for a given
      // fragment orbital angular momentum l and total angular momentum j
      // conserves parity, then compute an optical model transmission
      // coefficient and add it to the total.
      for (const auto& level : levels) {
        double Exf = level->energy();
        if (Exf < Exf_max)  {

          // Store information for this decay channel
	  auto ec = std::make_unique<marley::FragmentDiscreteExitChannel>(
            pdgi, qi, Exi_, twoJi_, Pi_, rho_i, sdb, *level, f );

          total_width_ += ec->width();

          exit_channels_.push_back( std::move(ec) );
        }
        else break;
      }
    }

    // If transitions to the energy continuum are possible, include the
    // continuum in the decay channels
    if ( Exf_max > E_c_min ) {

      // Create an ExitChannel object to handle decays to the continuum
      auto ec = std::make_unique<marley::FragmentContinuumExitChannel>(
        pdgi, qi, Exi_, twoJi_, Pi_, rho_i, sdb, E_c_min, f );

      total_width_ += ec->width();

      exit_channels_.push_back( std::move(ec) );
    }
  }

  marley::DecayScheme* ds = sdb.get_decay_scheme( pdgi );

  // For gamma-ray emission, let the continuum go down to Ex = 0 MeV unless
  // there is a decay scheme object available (we'll check this in a second).
  double E_c_min = 0.;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if ( ds ) {

    // Loop over the final discrete nuclear levels in order of increasing
    // energy until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission
    // coefficient for it
    const auto& levels = ds->get_levels();

    // Use the maximum discrete level energy from the decay scheme object as
    // the lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if ( levels.size() > 0 ) E_c_min = levels.back()->energy();

    for (const auto& level_f : levels) {
      double Exf = level_f->energy();
      if (Exf < Exi_) {
        auto ec = std::make_unique<marley::GammaDiscreteExitChannel>( pdgi, qi,
          Exi_, twoJi_, Pi_, rho_i, sdb, *level_f );

        total_width_ += ec->width();

        exit_channels_.push_back( std::move(ec) );
      }
      else break;
    }
  }

  // If gamma transitions to the energy continuum are possible, include them
  // in the possible decay channels
  if ( Exi_ > E_c_min ) {

    // Create an exit channel object to handle gamma-ray emission into the
    // continuum
    auto ec = std::make_unique<marley::GammaContinuumExitChannel>( pdgi, qi,
      Exi_, twoJi_, Pi_, rho_i, sdb, E_c_min );

    total_width_ += ec->width();

    exit_channels_.push_back( std::move(ec) );
  }
}

bool marley::HauserFeshbachDecay::do_decay(double& Exf, int& twoJf,
  marley::Parity& Pf, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& gen)
{
  const auto& ec = this->sample_exit_channel( gen );

  ec->do_decay( Exf, twoJf, Pf, compound_nucleus_, emitted_particle,
    residual_nucleus, gen );

  bool discrete_level = !ec->is_continuum();
  return !discrete_level;
}

void marley::HauserFeshbachDecay::print(std::ostream& out) const {

  // Needed to print results in conventional units
  constexpr double hbar = 6.58211951e-22; // MeV * s

  out << "Compound nucleus " << compound_nucleus_.pdg_code()
    << " with Ex = " << Exi_ << ", spin = " << twoJi_ / 2;
  if (twoJi_ % 2) out << ".5";
  out << ", and parity = " << Pi_ << '\n';
  out << "Total width = " << total_width_ << " MeV\n";
  out << "Mean lifetime = " << hbar / total_width_ << " s\n";
  for (const auto& ec : exit_channels_) {
    double width = ec->width();
    bool continuum = ec->is_continuum();
    bool frag = ec->emits_fragment();
    int pdg = ec->emitted_particle_pdg();
    std::string symbol = marley_utils::particle_symbols.at( pdg );
    out << "  ";
    if ( frag ) out << symbol;
    else out << "gamma-ray";
    if ( continuum ) out << " emission to the continuum width = ";
    else {
      auto* dec = dynamic_cast< marley::DiscreteExitChannel* >( ec.get() );
      if ( !dec ) throw marley::Error( "Dynamic cast failed in marley::"
        "HauserFeshbachDecay::print()" );
      out << " emission to level at " << dec->get_final_level().energy()
        << " MeV width = ";
    }
    out << width << " MeV\n";
  }
}

const std::unique_ptr<marley::ExitChannel>&
  marley::HauserFeshbachDecay::sample_exit_channel(
  marley::Generator& gen) const
{
  // Throw an error if all decays are impossible
  if ( total_width_ <= 0. ) throw marley::Error("Cannot sample an exit channel"
    " for a Hauser-Feshbach decay. All partial decay widths are zero.");

  // Sample an exit channel using a std::discrete distribution and a table of
  // partial decay widths
  const auto widths_begin
    = marley::ExitChannel::make_width_iterator( exit_channels_.cbegin() );
  const auto widths_end
    = marley::ExitChannel::make_width_iterator( exit_channels_.cend() );

  std::discrete_distribution<size_t> exit_channel_dist(widths_begin,
    widths_end);
  size_t exit_channel_index = gen.sample_from_distribution( exit_channel_dist );

  const auto& ec = exit_channels_.at( exit_channel_index );
  return ec;
}
