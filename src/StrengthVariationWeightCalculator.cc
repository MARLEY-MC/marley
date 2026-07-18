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
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// MARLEY includes
#include "marley/DiscreteNuclearReaction.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/StrengthVariationWeightCalculator.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"

marley::StrengthVariationWeightCalculator
  ::StrengthVariationWeightCalculator( const marley::JSON& config,
    long instance_index, std::shared_ptr< std::mt19937_64 > rng,
    const std::string& resolved_reaction_file )
  : WeightCalculator( config.at( "name" ).to_string() + '_'
    + std::to_string( instance_index ) ),
  rng_( std::move( rng ) ),
  resolved_reaction_file_( resolved_reaction_file )
{}

void marley::StrengthVariationWeightCalculator::ensure_initialized(
  marley::Generator& gen ) const
{
  if ( initialized_ ) return;

  // Scan the Generator's reactions to find the one matching our
  // resolved reaction file path and a discrete nuclear process type
  for ( const auto& reaction_ptr : gen.get_reactions() ) {

    if ( reaction_ptr->source_file() != resolved_reaction_file_ ) continue;

    // Only nuclear reactions populating discrete nuclear levels may be
    // handled by this weight calculator
    auto pt = reaction_ptr->process_type();
    if ( pt != Reaction::ProcessType::NeutrinoCC_Discrete
      && pt != Reaction::ProcessType::AntiNeutrinoCC_Discrete
      && pt != Reaction::ProcessType::NC_Discrete )
    {
      continue;
    }

    auto* dnr = dynamic_cast< const DiscreteNuclearReaction* >(
      reaction_ptr.get() );
    if ( !dnr ) throw marley::Error( "Reaction data file \""
      + resolved_reaction_file_ + "\" was loaded as a discrete nuclear"
      " reaction type but the corresponding Reaction object is not a"
      " DiscreteNuclearReaction. This should not happen and likely"
      " indicates a bug in MARLEY." );

    dnr_ = dnr;
    process_type_ = pt;
    target_pdg_ = dnr->pdg_b();

    // Pre-generate one dimidiated Gaussian throw per matrix element
    // using the shared RNG
    const auto& matrix_els = dnr->matrix_elements();
    varied_.reserve( matrix_els.size() );
    for ( const auto& me : matrix_els ) {
      double nom = me.strength();
      double err_low = me.strength_err_low();
      double err_high = me.strength_err_high();

      double varied;
      if ( err_low == 0. && err_high == 0. ) {
        varied = nom;
      }
      else {
        double u = normal_dist_( *rng_ );
        double sigma = ( u >= 0. ) ? err_high : err_low;
        varied = nom + u * sigma;
        if ( varied < 0. ) varied = 0.;
      }
      varied_.push_back( varied );
    }

    initialized_ = true;
    return;
  }

  throw marley::Error( "Could not find a discrete nuclear reaction"
    " with the source file \"" + resolved_reaction_file_
    + "\" in the Generator." );
}

double marley::StrengthVariationWeightCalculator::weight(
  HepMC3::GenEvent& event, marley::Generator& gen ) const
{
  ensure_initialized( gen );

  // Check that the event's process type matches this calculator's
  auto sp_attr = event.attribute< HepMC3::IntAttribute >(
    "signal_process_id" );
  if ( !sp_attr ) return 1.;
  auto event_pt = marley_hepmc3::from_nuhepmc_proc_id(
    sp_attr->value() );
  if ( event_pt != process_type_ ) return 1.;

  // Verify that the event's target nucleus PDG code matches
  auto target = marley_hepmc3::get_target( event );
  if ( !target || target->pdg_id() != target_pdg_ ) return 1.;

  // Read the matrix element index from the event
  auto mi_attr = event.attribute< HepMC3::IntAttribute >(
    "me_index" );
  if ( !mi_attr ) return 1.;
  size_t mi = static_cast< size_t >( mi_attr->value() );
  if ( mi >= varied_.size() ) return 1.;

  // Compute the weight as the ratio of the varied strength to the
  // nominal strength
  double nom = dnr_->matrix_elements().at( mi ).strength();
  if ( nom == 0. ) return 1.;
  return varied_.at( mi ) / nom;
}
