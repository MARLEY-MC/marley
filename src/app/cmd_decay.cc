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

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>
#include <vector>

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "marley/CommandHandler.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/OutputFile.hh"
#include "marley/Parity.hh"
#include "marley/Reaction.hh"
#include "marley/hepmc3_utils.hh"

using ProcType = marley::Reaction::ProcessType;

constexpr int TARGET_NET_CHARGE = 0;

constexpr std::array< ProcType, 6 > nuclear_proc_types
  = { ProcType::NeutrinoCC_Discrete, ProcType::AntiNeutrinoCC_Discrete,
      ProcType::NC_Discrete, ProcType::NeutrinoCC_Continuum,
      ProcType::AntiNeutrinoCC_Continuum, ProcType::NC_Continuum };

namespace {

  std::shared_ptr< HepMC3::GenEvent > make_event_object(
    int pdg_a, int pdg_b, int pdg_c, ProcType process_type, double Ex,
    int twoJ, const marley::Parity& P,
    std::shared_ptr< HepMC3::GenParticle >& residue, int residue_net_charge )
  {
    auto event = std::make_shared< HepMC3::GenEvent >( HepMC3::Units::MEV,
      HepMC3::Units::CM );

    int signal_process_id = marley_hepmc3::get_nuhepmc_proc_id( process_type );
    event->add_attribute( "signal_process_id",
      std::make_shared< HepMC3::IntAttribute >( signal_process_id )
    );

    auto prim_vtx = std::make_shared< HepMC3::GenVertex >();
    prim_vtx->set_status( marley_hepmc3::NUHEPMC_PRIMARY_VERTEX );

    event->add_vertex( prim_vtx );

    auto projectile = marley_hepmc3::make_particle( pdg_a, 0., 0., 0., 0.,
      marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, 0. );

    auto target = marley_hepmc3::make_particle( pdg_b,
      marley_hepmc3::NUHEPMC_TARGET_STATUS, 0. );

    auto ejectile = marley_hepmc3::make_particle( pdg_c, 0., 0.,
      0., 0., marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, 0. );

    prim_vtx->add_particle_in( projectile );
    prim_vtx->add_particle_in( target );

    prim_vtx->add_particle_out( ejectile );
    prim_vtx->add_particle_out( residue );

    residue->add_attribute( "Ex",
      std::make_shared< HepMC3::DoubleAttribute >(Ex) );
    residue->add_attribute( "twoJ",
      std::make_shared< HepMC3::IntAttribute >(twoJ) );
    residue->add_attribute( "parity",
      std::make_shared< HepMC3::IntAttribute >(static_cast<int>( P )) );

    marley_hepmc3::set_particle_charge( *target, TARGET_NET_CHARGE );
    marley_hepmc3::set_particle_charge( *residue, residue_net_charge );

    return event;
  }

}

bool marley::CommandHandler::cmd_decay( std::deque< std::string >& args ) {

  std::string config_file_name;
  if ( !args.empty() ) config_file_name = args.front();

  if ( config_file_name.empty() || config_file_name == "-h"
    || config_file_name == "--help" )
  {
    args.clear();
    args.push_front( "decay" );
    return marley::CommandHandler::cmd_help( args );
  }

  marley::JSONConfig jc( config_file_name );

  marley::Generator gen = jc.create_generator();
  gen.set_up_run_info();

  const marley::JSON& json = jc.get_json();
  bool ok;

  const std::string decay_config_label( "decay" );

  marley::JSON decay_config;
  ok = get_from_json< marley::JSON >( decay_config_label, json, decay_config );
  if ( !ok ) throw marley::Error( "Missing key '" + decay_config_label
    + "' in job configuration file" );

  long num_events = assign_from_json< long >( "events", decay_config, ok, 1000 );

  int projectile_pdg = assign_from_json< int >( "projectile", decay_config,
    ok, marley_utils::ELECTRON_NEUTRINO );

  int Zi = assign_from_json< int >( "target_Z", decay_config, ok );
  int Ai = assign_from_json< int >( "target_A", decay_config, ok );

  auto proc_type = static_cast< ProcType >(
    assign_from_json< int >( "proc_type", decay_config, ok,
      static_cast<int>(ProcType::NC_Continuum) )
  );

  auto iter = std::find( nuclear_proc_types.cbegin(),
    nuclear_proc_types.cend(), proc_type );

  bool process_is_a_nuclear_reaction = ( iter != nuclear_proc_types.cend() );

  if ( !process_is_a_nuclear_reaction ) {
    std::cerr << "The 'decay' command handles nuclear reactions only\n";
    return false;
  }

  double Ex = 0.;
  double Ex_min = 0.;
  double Ex_max = 0.;
  bool sample_Ex = decay_config.has_key( "Ex_max" );
  if ( !sample_Ex ) {
    Ex = assign_from_json< double >( "Ex", decay_config, ok, -1.0 );
  }
  else {
    Ex_min = assign_from_json< double >( "Ex_min", decay_config, ok, -1.0 );
    Ex_max = assign_from_json< double >( "Ex_max", decay_config, ok, -1.0 );
  }

  if ( Ex < 0. || Ex_min < 0. || Ex_max < 0. ) {
    throw marley::Error( "Negative excitation energy encountered" );
  }
  if ( Ex_max < Ex_min ) {
    throw marley::Error( "Upper Ex bound is less than lower Ex bound" );
  }

  std::vector< int > twoJ_vec;
  if ( !decay_config.has_key("twoJ") ) {
    throw marley::Error( "Missing \"twoJ\" key specifying the nuclear spin" );
  }
  else {
    const marley::JSON& twoJ_obj = decay_config.at( "twoJ" );
    if ( !twoJ_obj.is_array() ) {
     throw marley::Error( "The \"twoJ\" key must have a value that"
       " is a JSON array." );
    }
    else {
      convert_json< std::vector<int> >( twoJ_obj, twoJ_vec );
    }
  }

  for ( const auto& twoJ : twoJ_vec ) {
    if ( twoJ < 0 ) throw marley::Error( "Negative nuclear spin"
      " value encountered" );
  }

  std::vector< double > twoJ_sampling_weights( twoJ_vec.size(), 1. );
  const auto twoJ_begin = twoJ_sampling_weights.cbegin();
  const auto twoJ_end = twoJ_sampling_weights.cend();
  std::discrete_distribution< size_t > twoJ_dist( twoJ_begin, twoJ_end );

  auto parity_str = assign_from_json< std::string >( "parity", decay_config,
    ok, "+" );
  if ( parity_str != "+" && parity_str != "-" && parity_str != "random" ) {
    throw marley::Error( "Invalid parity setting \"" + parity_str + '\"' );
  }

  const std::vector< std::string > parity_strings = { "+", "-" };
  const std::vector< double > parity_sampling_weights = { 1., 1. };
  const auto par_begin = parity_sampling_weights.cbegin();
  const auto par_end = parity_sampling_weights.cend();
  std::discrete_distribution< size_t > par_dist( par_begin, par_end );

  int ejectile_pdg = marley::Reaction
    ::get_ejectile_pdg( projectile_pdg, proc_type );

  int projectile_charge = marley_utils::get_particle_charge( projectile_pdg );
  int ejectile_charge = marley_utils::get_particle_charge( ejectile_pdg );

  int Af = Ai;
  int Delta_Z = projectile_charge - ejectile_charge;
  int Zf = Zi + Delta_Z;

  int residue_net_charge = Delta_Z + TARGET_NET_CHARGE;

  std::vector< std::shared_ptr<marley::OutputFile> > output_files;

  if ( decay_config.has_key("output") ) {
    marley::JSON output_set = decay_config.at( "output" );
    if ( !output_set.is_array() ) throw marley::Error( "The"
      " \"output\" key must have a value that is a JSON array." );
    else for ( const auto& el : output_set.array_range() ) {
      output_files.push_back( marley::OutputFile::make_OutputFile(el) );
    }
  }
  else {
    std::string out_config_str = "{ format: \"ascii\","
      " file: \"decay_events.hepmc3\", mode: \"overwrite\" }";
    auto out_config = marley::JSON::load( out_config_str );

    output_files.push_back( marley::OutputFile::make_OutputFile(out_config) );
  }

  for ( long evnum = 0; evnum < num_events; ++evnum ) {

    const auto& mt = marley::MassTable::Instance();

    if ( sample_Ex ) {
      Ex = gen.uniform_random_double( Ex_min, Ex_max, true );
    }

    size_t twoJ_index = gen.sample_from_distribution( twoJ_dist );
    int twoJ = twoJ_vec.at( twoJ_index );

    std::string par_str = parity_str;
    if ( par_str == "random" ) {
      size_t par_index = gen.sample_from_distribution( par_dist );
      par_str = parity_strings.at( par_index );
    }
    marley::Parity parity;
    std::istringstream temp_iss( par_str );
    temp_iss >> parity;

    if ( Ex <= mt.unbound_threshold(Zf, Af) ) {
      auto* ds = gen.get_structure_db().get_decay_scheme( Zf, Af );

      if ( ds ) {
        auto* lev = ds->get_pointer_to_closest_level( Ex );
        Ex = lev->energy();
      }
    }

    int target_pdg = marley_utils::get_nucleus_pid( Zi, Ai );
    int residue_pdg = marley_utils::get_nucleus_pid( Zf, Af );

    double residue_gs_mass = mt.get_atomic_mass( residue_pdg )
      - ( residue_net_charge * mt.get_particle_mass(marley_utils::ELECTRON) );

    double m_residue = residue_gs_mass + Ex;

    auto residue = marley_hepmc3::make_particle( residue_pdg, 0., 0., 0.,
      m_residue, marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, m_residue );

    auto event = make_event_object( projectile_pdg, target_pdg, ejectile_pdg,
      proc_type, Ex, twoJ, parity, residue, residue_net_charge );

    marley::NucleusDecayer nd;
    nd.process_event( *event, gen );

    gen.finish_event_metadata( *event );

    for ( const auto& file : output_files ) {
      file->write_event( event.get() );
    }

    std::cout << "Event " << evnum << "\n";
  }

  return true;
}
