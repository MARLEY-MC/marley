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

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/OMPWeightCalculator.hh"
#include "marley/WeightCalculator.hh"
#include "marley/Weighter.hh"

// The fixed name of the central-value weight is defined by the NuHepMC
// standard
const std::string marley::Weighter::CV_WEIGHT_NAME = "CV";

marley::Weighter::Weighter( const marley::JSON& config ) {

  // Always create a trivial weight calculator to handle the central-value
  // weights, and add it to the vector of owned calculators
  auto cv_wgt = make_cv_weight_calc();
  calc_vec_.push_back( cv_wgt );

  // Now parse the JSON configuration to instantiate any other requested
  // weight calculators
  if ( !config.is_array() ) {
    throw marley::Error( "Non-array JSON configuration passed to constructor"
      " of marley::Weighter" );
  }

  for ( const auto& obj : config.array_range() ) {

    if ( !obj.is_object() ) {
      throw marley::Error( "Each weight calculator configuration must be"
        " specified as a JSON object" );
    }

    if ( !obj.has_key("type") ) {
      throw marley::Error( "Missing \"type\" key in a weight calculator"
        " JSON configuration" );
    }

    // Based on the type key included with each weight calculator
    // configuration, build the appropriate kind of WeightCalculator object
    std::shared_ptr< WeightCalculator > wc;
    auto type = obj.at( "type" ).to_string();
    if ( type == "trivial" ) {
      auto twc = std::make_shared< marley::TrivialWeightCalculator >( obj );
      wc = std::static_pointer_cast< marley::WeightCalculator >( twc );
    }
    else if ( type == "optical_model" ) {
      auto omp_wc = std::make_shared< marley::OMPWeightCalculator >( obj );
      wc = std::static_pointer_cast< marley::WeightCalculator >( omp_wc );
    }
    // TODO: add more options here
    else {
      throw marley::Error( "Unrecognized weight calculator type"
        " specification \"" + type + '\"' );
    }

    // Double-check that we don't have a requested weight calculator with
    // the same name as the required central-value one
    if ( wc->name() == CV_WEIGHT_NAME ) {
      throw marley::Error( "The weight calculator name \"" + CV_WEIGHT_NAME
        + "\" is reserved for internal MARLEY use" );
    }

    // Add the completed WeightCalculator object to the owned vector
    calc_vec_.push_back( wc );

  } // loop over weight calculator configuration JSON objects
}

void marley::Weighter::process_event( HepMC3::GenEvent& event,
  marley::Generator& gen )
{
  // Get non-const access to the vector of weights in the input event
  auto& weights_vec = event.weights();
  size_t num_weights = weights_vec.size();

  // Double-check that we have exactly the right number of weight
  // calculators configured
  if ( num_weights != calc_vec_.size() ) {
    throw marley::Error( "The number of configured weights is not the same"
      " as the number of configured weight calculators" );
  }

  // Evaluate all configured weights and store the results in the event
  auto temp_weights = compute_weights( event, gen );

  for ( size_t w = 0u; w < num_weights; ++w ) {
    double wgt = temp_weights.at( w );

    // Preserve any pre-existing event weights by multiplying the
    // current values in the event by the calculated result
    double& evw = weights_vec.at( w );
    evw *= wgt;
  }

}

std::vector< double > marley::Weighter::compute_weights(
  HepMC3::GenEvent& event, marley::Generator& gen ) const
{
  // Evaluate all configured weights and store the results in the output vector
  std::vector< double > weights;
  for ( const auto& weight_calc : calc_vec_ ) {
    double wgt = weight_calc->weight( event, gen );
    weights.push_back( wgt );
  }
  return weights;
}

std::vector< std::string > marley::Weighter::get_weight_names() const {
  // Look up the name of each owned weight calculator and store the results in
  // the output vector
  std::vector< std::string > names;
  for ( const auto& wgt_calc : calc_vec_ ) {
    names.push_back( wgt_calc->name() );
  }
  return names;
}

void marley::Weighter::set_use_cv_weight( bool use_it ) {
  // Determine whether or not the CV weight is currently enabled by
  // checking the first element of the vector of weight calculators
  bool currently_using = false;
  if ( !calc_vec_.empty() ) {
    auto first_calc_name = calc_vec_.front()->name();
    currently_using = ( first_calc_name == marley::Weighter::CV_WEIGHT_NAME );
  }

  // If the requested setting matches the current state, return without
  // doing anything
  if ( currently_using == use_it ) return;

  // If we need to add it, then put a TrivialWeightCalculator with the
  // correct name at the start of the vector
  else if ( use_it ) {
    auto cv_wgt_calc = make_cv_weight_calc();
    calc_vec_.insert( calc_vec_.begin(), cv_wgt_calc );
    return;
  }

  // If we need to remove it, then do so
  else {
    calc_vec_.erase( calc_vec_.begin() );
  }
}

// Helper function to construct a trivial weight calculator to manage
// the central-value weight assignment
std::shared_ptr< marley::WeightCalculator > marley::Weighter
  ::make_cv_weight_calc()
{
  marley::JSON temp_js;
  temp_js[ "name" ] = marley::Weighter::CV_WEIGHT_NAME;
  auto cv_wgt = std::make_shared< marley
    ::TrivialWeightCalculator >( temp_js );
  return std::static_pointer_cast< marley::WeightCalculator >( cv_wgt );
}
