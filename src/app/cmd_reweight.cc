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
#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "cmd_helpers.hh"
#include "marley/CommandHandler.hh"
#include "marley/EventFileReader.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"
#include "marley/OutputFile.hh"
#include "marley/WeightCalculator.hh"
#include "marley/Weighter.hh"
#include "marley/marley_utils.hh"

bool marley::CommandHandler::cmd_reweight( std::deque< std::string >& args ) {

  // If we have fewer than two arguments, decide whether the
  // user intended to request help with this command
  if ( args.size() < 2u ) {
    std::string first_arg;
    if ( !args.empty() ) first_arg = args.front();

    // Print the help message either way
    args.clear();
    args.push_front( "reweight" );
    marley::CommandHandler::cmd_help( args );

    // Return a boolean status based on whether the help message was
    // explicitly requested (normal behavior) or not (an error condition)
    if ( first_arg == "-h" || first_arg == "--help" ) return true;
    return false;
  }

  // Extract the configuration file name and collect all input files
  std::string config_file_name( args.front() );
  args.pop_front();
  std::vector< std::string > input_files( args.begin(), args.end() );

  // Load the reweight configuration
  marley::JSON rw_config = marley::JSON::load_file( config_file_name );
  if ( !rw_config.has_key("weights") ) throw marley::Error( "Missing"
    " \"weights\" key in marley reweight configuration file \""
    + config_file_name + "\"" );

  const auto& json_weights = rw_config.at( "weights" );

  // Read the first event from the first input file to extract the run
  // information needed to reconstruct the Generator and check weight names.
  marley::EventFileReader first_reader( input_files[0] );
  HepMC3::GenEvent first_ev;
  if ( !(first_reader >> first_ev) ) {
    throw marley::Error( "Failed to read the first event from input file \""
      + input_files[0] + "\". The file may be empty or corrupt." );
  }

  auto first_run_info = first_ev.run_info();
  const std::vector< std::string > wgt_names = first_run_info->weight_names();

  // Reconstruct the original Generator from the saved configuration
  auto prior_config_str = first_run_info->attribute< HepMC3::StringAttribute >(
    "MARLEY.JSONconfig" );

  if ( !prior_config_str ) {
    throw marley::Error( "Failed to retrieve previous generator"
      " configuration from the input file \"" + input_files[0] + "\"" );
  }

  auto prior_json_config = marley::JSON::load( prior_config_str->value() );
  marley::JSONConfig jc( prior_json_config );
  auto gen = std::make_unique< marley::Generator >( jc.create_generator() );

  // Create the Weighter and check for name conflicts
  marley::Weighter weighter( json_weights, *gen );
  weighter.set_use_cv_weight( false );

  auto& calc_vec = weighter.get_weight_calculators();

  for ( const auto& wc : calc_vec ) {
    if ( std::find( wgt_names.cbegin(), wgt_names.cend(), wc->name() )
      != wgt_names.cend() )
    {
      throw marley::Error( "Weight name \"" + wc->name()
        + "\" from the reweight configuration file \"" + config_file_name
        + "\" conflicts with an existing weight in the input file" );
    }
  }

  // Prepend TrivialWeightCalculators for the existing weight names so that
  // the Weighter preserves them in the output. Iterate in reverse order
  // and insert at the beginning to maintain the original ordering.
  for ( auto riter = wgt_names.crbegin();
    riter != wgt_names.crend(); ++riter )
  {
    const auto& w_name = *riter;
    marley::JSON temp_json;
    temp_json[ "name" ] = w_name;
    auto w_calc = std::make_shared< marley
      ::TrivialWeightCalculator >( temp_json );
    calc_vec.insert( calc_vec.begin(), w_calc );
  }

  auto full_name_vec = weighter.get_weight_names();

  // Read output settings from the optional "reweight" section (if present)
  marley::JSON rw_section;
  bool has_rw_section = false;
  if ( rw_config.has_key("reweight") ) {
    const marley::JSON& rw_section_ref = rw_config.at( "reweight" );
    if ( !rw_section_ref.is_object() ) throw marley::Error(
      "The \"reweight\" section in the marley reweight configuration"
      " file \"" + config_file_name + "\" must be a JSON object" );
    rw_section = rw_section_ref;
    has_rw_section = true;
  }

  std::vector< std::shared_ptr<marley::OutputFile> > output_files;

  if ( has_rw_section && rw_section.has_key("output") ) {
    marley::JSON output_set = rw_section.at( "output" );
    if ( !output_set.is_array() ) throw marley::Error( "The"
      " \"output\" key in the reweighting configuration must have a value"
      " that is a JSON array." );
    else for ( const auto& el : output_set.array_range() ) {
      if ( el.has_key("mode") ) {
        std::string mode_str = el.at( "mode" ).to_string();
        if ( mode_str != "overwrite" ) throw marley::Error( "Only the"
          " \"overwrite\" output file mode is allowed for a reweighting"
          " job." );
      }
      output_files.push_back( marley::OutputFile::make_OutputFile(el) );
    }
  }
  else {
    std::string out_config_str = "{ format: \"ascii\","
      " file: \"reweighted_events.hepmc3\", mode: \"overwrite\" }";
    auto out_config = marley::JSON::load( out_config_str );

    output_files.push_back( marley::OutputFile::make_OutputFile(out_config) );
  }

  // Build the reweighted GenRunInfo from a copy of the first file's run info
  bool multi_file = ( input_files.size() > 1 );
  auto reweighted_run_info = std::make_shared< HepMC3::GenRunInfo >(
    *first_run_info );
  reweighted_run_info->set_weight_names( full_name_vec );

  // For multi-file reweight, strip the RNG seed to prevent unsafe resume.
  // Single-file reweight preserves the seed so that resume remains possible
  // (the accumulated Weighter will be reconstructed from the saved reweight
  // provenance attributes when needed).
  if ( multi_file ) reweighted_run_info->remove_attribute(
    "MARLEY.RNGseed" );

  // Save the reweight configuration as run info provenance attributes
  {
    int rw_index = 0;
    auto count_attr = first_run_info->attribute< HepMC3::IntAttribute >(
      "MARLEY.ReweightConfig.count" );
    if ( count_attr ) rw_index = count_attr->value();

    // Build a combined provenance object with the weights array and
    // (if present) the reweight section containing output settings
    marley::JSON prov_obj = marley::JSON::object();
    prov_obj["weights"] = json_weights;
    if ( has_rw_section ) {
      prov_obj["reweight"] = rw_section;
    }

    reweighted_run_info->add_attribute(
      "MARLEY.ReweightConfig." + std::to_string( rw_index ),
      std::make_shared< HepMC3::StringAttribute >(
        prov_obj.dump_string() ) );

    reweighted_run_info->add_attribute(
      "MARLEY.ReweightConfig.count",
      std::make_shared< HepMC3::IntAttribute >( rw_index + 1 ) );
  }

  // Process all events across all input files
  int event_count = 0;
  for_each_event( input_files,
    [ & ]( HepMC3::GenEvent& ev, bool /*first_event*/,
      double /*flux_avg_xsec*/, const auto& /*first_info*/ )
    {
      std::cout << "Event " << event_count << '\n';

      // Save the original weight values before set_run_info resizes
      // the event's weight vector to match the reweighted run info.
      auto orig_weights = ev.weights();

      // Apply the reweighted GenRunInfo. HepMC3's set_run_info resizes
      // m_weights to match the weight_names count, filling with 1.0.
      ev.set_run_info( reweighted_run_info );

      // Restore the original weight values into the first slots
      for ( size_t i = 0; i < orig_weights.size(); ++i )
        ev.weights()[ i ] = orig_weights[ i ];

      // Compute the new weight values
      weighter.process_event( ev, *gen );

      // Write the event to all output files
      for ( const auto& file : output_files )
        file->write_event( &ev );

      ++event_count;
    } );

  return true;
}
