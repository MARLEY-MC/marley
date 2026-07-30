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

// Standard library includes
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_set>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// MARLEY includes
#include "marley/DiscreteNuclearReaction.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/JSON.hh"
#include "marley/Generator.hh"
#include "marley/StrengthVariationWeightCalculator.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"

namespace {

  /// @brief Generates the shortest string label for each double that
  /// roundtrips via stod back to the original value. Since numerically
  /// distinct inputs produce distinct roundtripping labels, uniqueness
  /// is guaranteed without extra collision resolution.
  std::vector< std::string > shortest_roundtrip_labels(
    const std::vector< double >& values )
  {
    constexpr int MAX_PREC = std::numeric_limits< double >::max_digits10;
    std::vector< std::string > labels;
    labels.reserve( values.size() );

    for ( double v : values ) {
      std::string best;
      for ( int prec = 0; prec <= MAX_PREC; ++prec ) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision( prec ) << v;
        std::string s = ss.str();

        auto dot = s.find( '.' );
        if ( dot != std::string::npos ) {
          auto last = s.find_last_not_of( '0' );
          if ( last == dot ) s.erase( dot );
          else s.erase( last + 1 );
        }

        if ( std::stod( s ) == v ) { best = s; break; }
      }
      if ( best.empty() ) best = std::to_string( v );
      labels.push_back( std::move( best ) );
    }

    std::unordered_set< std::string > seen;
    for ( const auto& l : labels ) {
      if ( !seen.insert( l ).second ) {
        throw marley::Error( "Unexpected label collision in"
          " shortest_roundtrip_labels" );
      }
    }

    return labels;
  }

} // anonymous namespace

// -------------------------------------------------------------------
// Multisim constructor
// -------------------------------------------------------------------
marley::StrengthVariationWeightCalculator
  ::StrengthVariationWeightCalculator( const std::string& name,
    std::shared_ptr< std::mt19937_64 > rng,
    const std::string& resolved_reaction_file )
  : WeightCalculator( name ), rng_( std::move( rng ) ),
  resolved_reaction_file_( resolved_reaction_file ),
  mode_( VariationMode::multisim )
{}

// -------------------------------------------------------------------
// Sigma_shift constructor
// -------------------------------------------------------------------
marley::StrengthVariationWeightCalculator
  ::StrengthVariationWeightCalculator( const std::string& name,
    double sigma_factor, const std::string& resolved_reaction_file )
  : WeightCalculator( name ), rng_( nullptr ),
  resolved_reaction_file_( resolved_reaction_file ),
  mode_( VariationMode::sigma_shift ),
  sigma_factor_( sigma_factor )
{}

// -------------------------------------------------------------------
// Unisim constructor
// -------------------------------------------------------------------
marley::StrengthVariationWeightCalculator
  ::StrengthVariationWeightCalculator( const std::string& name,
    double sigma_factor, size_t matrix_element_index,
    const std::string& resolved_reaction_file )
  : WeightCalculator( name ), rng_( nullptr ),
  resolved_reaction_file_( resolved_reaction_file ),
  mode_( VariationMode::unisim ),
  sigma_factor_( sigma_factor ),
  me_idx_( matrix_element_index )
{}

// -------------------------------------------------------------------
// Static factory: create_instances
// -------------------------------------------------------------------
std::vector< std::shared_ptr<
  marley::StrengthVariationWeightCalculator > >
  marley::StrengthVariationWeightCalculator::create_instances(
    const marley::JSON& config, marley::Generator& gen )
{
  // Resolve the reaction file
  if ( !config.has_key( "reaction_file" ) ) {
    throw marley::Error( "Missing \"reaction_file\" key in a"
      " strength_variation weight calculator JSON configuration" );
  }
  std::string reaction_file = config.at( "reaction_file" ).to_string();
  std::string resolved_reaction_file = marley::FileManager::Instance()
    .find_file( reaction_file );
  if ( resolved_reaction_file.empty() ) {
    throw marley::Error( "Could not find reaction data file \""
      + reaction_file + "\" requested by a strength_variation"
      " weight calculator" );
  }

  // Read the optional variation mode (default "multisim")
  std::string mode_str = "multisim";
  if ( config.has_key( "mode" ) ) {
    mode_str = config.at( "mode" ).to_string();
  }

  // Extract the base name
  if ( !config.has_key( "name" ) ) {
    throw marley::Error( "Missing \"name\" key in a"
      " strength_variation weight calculator JSON configuration" );
  }
  std::string base_name = config.at( "name" ).to_string();

  std::vector< std::shared_ptr<
    StrengthVariationWeightCalculator > > instances;

  if ( mode_str == "multisim" ) {

    // Reject sigma_factor key in multisim mode
    if ( config.has_key( "sigma_factor" ) ) {
      throw marley::Error( "The \"sigma_factor\" key is not allowed"
        " in multisim mode for strength_variation weight calculators" );
    }

    // Read and validate the number of variations
    if ( !config.has_key( "num_variations" ) ) {
      throw marley::Error( "Missing \"num_variations\" key in a"
        " strength_variation weight calculator JSON configuration" );
    }
    const auto& nv = config.at( "num_variations" );
    if ( !nv.is_integer() ) {
      throw marley::Error( "The \"num_variations\" value must be a"
        " positive integer" );
    }
    long num_instances = nv.to_long();
    if ( num_instances <= 0 ) {
      throw marley::Error( "The \"num_variations\" value must be a"
        " positive integer" );
    }

    // Read the optional seed (default 0)
    long seed = 0;
    if ( config.has_key( "seed" ) ) {
      seed = config.at( "seed" ).to_long();
    }

    // Create one shared RNG for all instances
    auto rng = std::make_shared< std::mt19937_64 >(
      static_cast< std::mt19937_64::result_type >( seed ) );

    // Create num_instances weight calculators, each sharing the RNG
    for ( long idx = 0; idx < num_instances; ++idx ) {
      instances.push_back( std::shared_ptr<
        StrengthVariationWeightCalculator >(
          new StrengthVariationWeightCalculator(
            base_name + '_' + std::to_string( idx ),
            rng, resolved_reaction_file ) ) );
    }
  }
  else if ( mode_str == "sigma_shift" ) {

    // Reject keys inappropriate for sigma_shift
    if ( config.has_key( "num_variations" ) ) {
      throw marley::Error( "The \"num_variations\" key is not allowed"
        " in sigma_shift mode for strength_variation weight"
        " calculators" );
    }
    if ( config.has_key( "seed" ) ) {
      throw marley::Error( "The \"seed\" key is not allowed"
        " in sigma_shift mode for strength_variation weight"
        " calculators" );
    }

    // Read and validate sigma_factor
    if ( !config.has_key( "sigma_factor" ) ) {
      throw marley::Error( "Missing \"sigma_factor\" key in a"
        " strength_variation weight calculator JSON configuration"
        " with sigma_shift mode" );
    }
    const auto& sf = config.at( "sigma_factor" );

    std::vector< double > factors;
    if ( sf.is_array() ) {
      for ( const auto& elem : sf.array_range() ) {
        factors.push_back( elem.to_double_or_throw() );
      }
      if ( factors.empty() ) {
        throw marley::Error( "The \"sigma_factor\" array must have"
          " at least one element" );
      }
    }
    else {
      factors.push_back( sf.to_double_or_throw() );
    }

    // Reject duplicate sigma_factor values (exact equality)
    {
      std::set< double > seen;
      for ( double f : factors ) {
        if ( !seen.insert( f ).second ) {
          throw marley::Error( "Duplicate sigma_factor value \""
            + std::to_string( f ) + "\" in strength_variation"
            " weight calculator configuration" );
        }
      }
    }

    // Generate shortest roundtrip labels for all factors
    auto labels = shortest_roundtrip_labels( factors );

    // Create two instances per factor (+k and -k)
    for ( size_t i = 0; i < factors.size(); ++i ) {
      double k = factors[ i ];

      // +k "up" instance
      instances.push_back( std::shared_ptr<
        StrengthVariationWeightCalculator >(
          new StrengthVariationWeightCalculator(
            base_name + "-up@" + labels[ i ],
            +k, resolved_reaction_file ) ) );

      // -k "down" instance
      instances.push_back( std::shared_ptr<
        StrengthVariationWeightCalculator >(
          new StrengthVariationWeightCalculator(
            base_name + "-down@" + labels[ i ],
            -k, resolved_reaction_file ) ) );
    }
  }
  else if ( mode_str == "unisim" ) {

    // Reject keys inappropriate for unisim
    if ( config.has_key( "num_variations" ) ) {
      throw marley::Error( "The \"num_variations\" key is not allowed"
        " in unisim mode for strength_variation weight calculators" );
    }
    if ( config.has_key( "seed" ) ) {
      throw marley::Error( "The \"seed\" key is not allowed"
        " in unisim mode for strength_variation weight calculators" );
    }

    // Read and validate sigma_factor
    if ( !config.has_key( "sigma_factor" ) ) {
      throw marley::Error( "Missing \"sigma_factor\" key in a"
        " strength_variation weight calculator JSON configuration"
        " with unisim mode" );
    }
    const auto& sf = config.at( "sigma_factor" );

    std::vector< double > factors;
    if ( sf.is_array() ) {
      for ( const auto& elem : sf.array_range() ) {
        factors.push_back( elem.to_double_or_throw() );
      }
      if ( factors.empty() ) {
        throw marley::Error( "The \"sigma_factor\" array must have"
          " at least one element" );
      }
    }
    else {
      factors.push_back( sf.to_double_or_throw() );
    }

    // Reject duplicate sigma_factor values (exact equality)
    {
      std::set< double > seen;
      for ( double f : factors ) {
        if ( !seen.insert( f ).second ) {
          throw marley::Error( "Duplicate sigma_factor value \""
            + std::to_string( f ) + "\" in strength_variation"
            " weight calculator configuration" );
        }
      }
    }

    // Find the matching DiscreteNuclearReaction from the Generator
    const DiscreteNuclearReaction* dnr = nullptr;
    for ( const auto& rptr : gen.get_reactions() ) {
      if ( rptr->source_file() != resolved_reaction_file ) continue;
      auto pt = rptr->process_type();
      if ( pt != Reaction::ProcessType::NeutrinoCC_Discrete
        && pt != Reaction::ProcessType::AntiNeutrinoCC_Discrete
        && pt != Reaction::ProcessType::NC_Discrete )
      {
        continue;
      }
      dnr = dynamic_cast< const DiscreteNuclearReaction* >(
        rptr.get() );
      if ( !dnr ) throw marley::Error( "Reaction data file \""
        + resolved_reaction_file + "\" was loaded as a discrete nuclear"
        " reaction type but the corresponding Reaction object is not a"
        " DiscreteNuclearReaction. This should not happen and likely"
        " indicates a bug in MARLEY." );
      break;
    }
    if ( !dnr ) throw marley::Error( "Could not find a discrete nuclear"
      " reaction with the source file \"" + resolved_reaction_file
      + "\" in the Generator." );

    size_t M = dnr->matrix_elements().size();

    // Generate shortest roundtrip labels for all factors
    auto labels = shortest_roundtrip_labels( factors );

    // Create two instances per factor per matrix element
    for ( size_t i = 0; i < factors.size(); ++i ) {
      double k = factors[ i ];
      for ( size_t m = 0; m < M; ++m ) {
        // +k "up" instance for this matrix element
        instances.push_back( std::shared_ptr<
          StrengthVariationWeightCalculator >(
            new StrengthVariationWeightCalculator(
              base_name + "-me" + std::to_string( m )
                + "_up@" + labels[ i ],
              +k, m, resolved_reaction_file ) ) );

        // -k "down" instance for this matrix element
        instances.push_back( std::shared_ptr<
          StrengthVariationWeightCalculator >(
            new StrengthVariationWeightCalculator(
              base_name + "-me" + std::to_string( m )
                + "_down@" + labels[ i ],
              -k, m, resolved_reaction_file ) ) );
      }
    }
  }
  else {
    throw marley::Error( "Unrecognized variation mode \""
      + mode_str + "\" for strength_variation weight calculator."
      " Allowed values are \"multisim\", \"sigma_shift\","
      " and \"unisim\"" );
  }

  return instances;
}

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

    // Pre-generate varied strengths depending on the variation mode
    const auto& matrix_els = dnr->matrix_elements();
    varied_.reserve( matrix_els.size() );

    if ( mode_ == VariationMode::multisim ) {
      std::normal_distribution< double > normal_dist;
      for ( const auto& me : matrix_els ) {
        double nom = me.strength();
        double err_low = me.strength_err_low();
        double err_high = me.strength_err_high();

        double varied;
        if ( err_low == 0. && err_high == 0. ) {
          varied = nom;
        }
        else {
          double u = normal_dist( *rng_ );
          double sigma = ( u >= 0. ) ? err_high : err_low;
          varied = nom + u * sigma;
          if ( varied < 0. ) varied = 0.;
        }
        varied_.push_back( varied );
      }
    }
    else if ( mode_ == VariationMode::sigma_shift ) {
      for ( const auto& me : matrix_els ) {
        double nom = me.strength();
        double err_low = me.strength_err_low();
        double err_high = me.strength_err_high();

        double varied;
        if ( err_low == 0. && err_high == 0. ) {
          varied = nom;
        }
        else {
          double sigma = ( sigma_factor_ >= 0. ) ? err_high : err_low;
          varied = nom + sigma_factor_ * sigma;
          if ( varied < 0. ) varied = 0.;
        }
        varied_.push_back( varied );
      }
    }
    else { // unisim
      for ( const auto& me : matrix_els ) {
        varied_.push_back( me.strength() );
      }
      const auto& target_me = matrix_els[ me_idx_ ];
      double nom = target_me.strength();
      double err_low = target_me.strength_err_low();
      double err_high = target_me.strength_err_high();
      if ( err_low != 0. || err_high != 0. ) {
        double sigma = ( sigma_factor_ >= 0. ) ? err_high : err_low;
        double varied = nom + sigma_factor_ * sigma;
        if ( varied < 0. ) varied = 0.;
        varied_[ me_idx_ ] = varied;
      }
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
