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

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "marley/CommandHandler.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/Logger.hh"
#include "marley/MassTable.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/OutputFile.hh"
#include "marley/Parity.hh"
#include "marley/Reaction.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"

using ProcType = marley::Reaction::ProcessType;

// ---------------------------------------------------------------------------
// Helper: derive Z and A from a nuclear PDG code.
// Delegates to the canonical marley_utils helpers, which already handle the
// neutron (2112), proton (2212), and general nuclear ion cases.
// Returns false if the PDG code is not recognized as a nucleus or nucleon.
// ---------------------------------------------------------------------------
static bool pdg_to_ZA( int pdg, int& Z, int& A ) {
  if ( pdg == marley_utils::NEUTRON || pdg == marley_utils::PROTON
    || marley_utils::is_ion(pdg) )
  {
    Z = marley_utils::get_particle_Z( pdg );
    A = marley_utils::get_particle_A( pdg );
    return true;
  }
  return false;
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

  // -------------------------------------------------------------------------
  // Load the JSON config from file, then patch in null entries for "reactions"
  // and "source" if those keys are absent.  This allows create_generator() to
  // take its early-return path and skip cross-section normalisation, which is
  // irrelevant for standalone de-excitation events.
  // -------------------------------------------------------------------------
  marley::JSON json = marley::JSON::load_file( config_file_name );

  if ( !json.has_key("reactions") ) json["reactions"] = marley::JSON( nullptr );
  if ( !json.has_key("source")    ) json["source"]    = marley::JSON( nullptr );

  marley::JSONConfig jc( json );
  marley::Generator gen = jc.create_generator();
  gen.set_up_run_info();

  bool ok;

  // -------------------------------------------------------------------------
  // Parse the required top-level "decay" block
  // -------------------------------------------------------------------------
  const std::string decay_config_label( "decay" );
  marley::JSON decay_config;
  ok = get_from_json< marley::JSON >( decay_config_label, json, decay_config );
  if ( !ok ) throw marley::Error( "Missing key '" + decay_config_label
    + "' in job configuration file" );

  long num_events = assign_from_json< long >( "events", decay_config, ok, 1000 );

  // -------------------------------------------------------------------------
  // Parse the required "nucleus" sub-object
  // -------------------------------------------------------------------------
  if ( !decay_config.has_key("nucleus") ) {
    throw marley::Error( "Missing required \"nucleus\" key in the"
      " \"decay\" configuration block" );
  }
  const marley::JSON& nuc_config = decay_config.at( "nucleus" );

  // Read the nucleus identity.  The user may supply either (or both) of:
  //   Option A: Z and A integer keys
  //   Option B: pdg integer key (nuclear PDG code)
  // If both are supplied they must be consistent.

  bool has_pdg = nuc_config.has_key("pdg");
  bool has_Z   = nuc_config.has_key("Z");
  bool has_A   = nuc_config.has_key("A");

  if ( !has_pdg && !( has_Z && has_A ) ) {
    throw marley::Error( "The \"nucleus\" block must specify the nuclide"
      " using either the \"pdg\" key or both the \"Z\" and \"A\" keys" );
  }

  int nucleus_pdg = 0;
  int Z = 0, A = 0;

  if ( has_pdg ) {
    nucleus_pdg = assign_from_json< int >( "pdg", nuc_config, ok );
    int Z_from_pdg = 0, A_from_pdg = 0;
    if ( !pdg_to_ZA( nucleus_pdg, Z_from_pdg, A_from_pdg ) ) {
      throw marley::Error( "The value " + std::to_string(nucleus_pdg)
        + " given for \"nucleus.pdg\" is not a recognized nuclear"
        " or nucleon PDG code" );
    }
    Z = Z_from_pdg;
    A = A_from_pdg;
  }

  if ( has_Z && has_A ) {
    int Z_cfg = assign_from_json< int >( "Z", nuc_config, ok );
    int A_cfg = assign_from_json< int >( "A", nuc_config, ok );
    if ( has_pdg ) {
      // Both representations present — verify consistency
      if ( Z_cfg != Z || A_cfg != A ) {
        throw marley::Error( "Inconsistent nucleus specification: pdg="
          + std::to_string(nucleus_pdg) + " implies Z=" + std::to_string(Z)
          + ", A=" + std::to_string(A) + ", but Z=" + std::to_string(Z_cfg)
          + ", A=" + std::to_string(A_cfg) + " were also given" );
      }
    }
    else {
      Z = Z_cfg;
      A = A_cfg;
      nucleus_pdg = marley_utils::get_nucleus_pid( Z, A );
    }
  }

  if ( Z < 0 ) throw marley::Error( "Negative Z encountered" );
  if ( A < 1 ) throw marley::Error( "A < 1 encountered" );

  // Net ionic charge of the nucleus (protons minus electrons).
  // Default 0 = neutral atom.
  int net_charge = assign_from_json< int >( "net_charge", nuc_config, ok, 0 );

  // -------------------------------------------------------------------------
  // Parse excitation energy: fixed value or uniform sampling interval
  // -------------------------------------------------------------------------
  double Ex = 0.;
  double Ex_min = 0.;
  double Ex_max = 0.;
  bool sample_Ex = nuc_config.has_key( "Ex_max" );
  if ( !sample_Ex ) {
    Ex = assign_from_json< double >( "Ex", nuc_config, ok, -1.0 );
    if ( !ok ) throw marley::Error( "Missing \"Ex\" key in \"nucleus\" block"
      " (or use \"Ex_min\" + \"Ex_max\" for uniform sampling)" );
    if ( Ex < 0. ) throw marley::Error( "Negative excitation energy"
      " encountered (Ex = " + std::to_string(Ex) + " MeV)" );
  }
  else {
    Ex_min = assign_from_json< double >( "Ex_min", nuc_config, ok, -1.0 );
    Ex_max = assign_from_json< double >( "Ex_max", nuc_config, ok, -1.0 );
    if ( Ex_min < 0. ) throw marley::Error( "Negative lower excitation"
      " energy bound (Ex_min = " + std::to_string(Ex_min) + " MeV)" );
    if ( Ex_max < 0. ) throw marley::Error( "Negative upper excitation"
      " energy bound (Ex_max = " + std::to_string(Ex_max) + " MeV)" );
    if ( Ex_max < Ex_min ) throw marley::Error( "Upper Ex bound"
      " (" + std::to_string(Ex_max) + " MeV) is less than lower Ex bound"
      " (" + std::to_string(Ex_min) + " MeV)" );
  }

  // -------------------------------------------------------------------------
  // Parse nuclear spin: required array of 2J values
  // -------------------------------------------------------------------------
  std::vector< int > twoJ_vec;
  if ( !nuc_config.has_key("twoJ") ) {
    throw marley::Error( "Missing \"twoJ\" key in \"nucleus\" block" );
  }
  else {
    const marley::JSON& twoJ_obj = nuc_config.at( "twoJ" );
    if ( !twoJ_obj.is_array() ) {
      throw marley::Error( "The \"twoJ\" key in \"nucleus\" must have"
        " a value that is a JSON array" );
    }
    convert_json< std::vector<int> >( twoJ_obj, twoJ_vec );
  }
  if ( twoJ_vec.empty() ) {
    throw marley::Error( "The \"twoJ\" array in \"nucleus\" must not"
      " be empty" );
  }
  for ( const auto& tJ : twoJ_vec ) {
    if ( tJ < 0 ) throw marley::Error( "Negative 2J value encountered"
      " in \"twoJ\" array" );
  }

  std::vector< double > twoJ_weights( twoJ_vec.size(), 1. );
  std::discrete_distribution< size_t > twoJ_dist(
    twoJ_weights.cbegin(), twoJ_weights.cend() );

  // -------------------------------------------------------------------------
  // Parse parity: optional, default "+"
  // -------------------------------------------------------------------------
  auto parity_str = assign_from_json< std::string >( "parity", nuc_config,
    ok, "+" );
  if ( parity_str != "+" && parity_str != "-" && parity_str != "random" ) {
    throw marley::Error( "Invalid \"parity\" setting \""
      + parity_str + "\" in \"nucleus\" block."
      " Expected \"+\", \"-\", or \"random\"" );
  }

  const std::vector< std::string > parity_strings = { "+", "-" };
  const std::vector< double > parity_weights = { 1., 1. };
  std::discrete_distribution< size_t > parity_dist(
    parity_weights.cbegin(), parity_weights.cend() );

  // -------------------------------------------------------------------------
  // Parse output configuration
  // -------------------------------------------------------------------------
  std::vector< std::shared_ptr<marley::OutputFile> > output_files;

  if ( decay_config.has_key("output") ) {
    marley::JSON output_set = decay_config.at( "output" );
    if ( !output_set.is_array() ) throw marley::Error( "The"
      " \"output\" key must have a value that is a JSON array." );
    for ( const auto& el : output_set.array_range() ) {
      output_files.push_back( marley::OutputFile::make_OutputFile(el) );
    }
  }
  else {
    std::string out_cfg_str = "{ format: \"ascii\","
      " file: \"decay_events.hepmc3\", mode: \"overwrite\" }";
    auto out_cfg = marley::JSON::load( out_cfg_str );
    output_files.push_back( marley::OutputFile::make_OutputFile(out_cfg) );
  }

  for ( const auto& file : output_files ) {
    if ( file->mode_is_resume() ) {
      throw marley::Error( "The \"resume\" output mode is not supported"
        " by the \"marley decay\" command" );
    }
  }

  // -------------------------------------------------------------------------
  // Pre-compute fixed quantities used each event
  // -------------------------------------------------------------------------
  const auto& mt = marley::MassTable::Instance();

  double gs_mass = mt.get_atomic_mass( nucleus_pdg )
    - net_charge * mt.get_particle_mass( marley_utils::ELECTRON );

  double unbound_threshold = mt.unbound_threshold( nucleus_pdg );

  int signal_proc_id = marley_hepmc3::get_nuhepmc_proc_id(
    ProcType::StandaloneDecay );

  // -------------------------------------------------------------------------
  // Event loop
  // -------------------------------------------------------------------------
  // One-time warning flags for the level-snap behaviour
  bool warned_snap        = false; // suppress after first snap warning
  bool warned_snap_sample = false; // suppress after first "will continue" note

  for ( long evnum = 0; evnum < num_events; ++evnum ) {

    // --- Sample Ex, twoJ, and parity for this event -----------------------
    if ( sample_Ex ) {
      Ex = gen.uniform_random_double( Ex_min, Ex_max, true );
    }

    size_t twoJ_index = gen.sample_from_distribution( twoJ_dist );
    int twoJ = twoJ_vec.at( twoJ_index );

    std::string par_str = parity_str;
    if ( par_str == "random" ) {
      size_t par_index = gen.sample_from_distribution( parity_dist );
      par_str = parity_strings.at( par_index );
    }
    marley::Parity parity;
    std::istringstream temp_iss( par_str );
    temp_iss >> parity;

    // --- Snap to nearest discrete level if below the unbound threshold ----
    if ( Ex <= unbound_threshold ) {
      auto* ds = gen.get_structure_db().get_decay_scheme( nucleus_pdg );
      if ( ds ) {
        auto* lev = ds->get_pointer_to_closest_level( Ex );

        // Check what changes on the snap
        double Ex_lev   = lev->energy();
        int    twoJ_lev = lev->twoJ();
        marley::Parity P_lev = lev->parity();

        bool Ex_changed   = ( std::abs(Ex_lev - Ex) > 1e-5 );
        bool twoJ_changed = ( twoJ_lev != twoJ );
        bool P_changed    = ( static_cast<int>(P_lev)
                              != static_cast<int>(parity) );

        if ( ( Ex_changed || twoJ_changed || P_changed ) && !warned_snap ) {
          std::ostringstream warn_msg;
          warn_msg << "User-specified nuclear state (Ex=" << Ex
            << " MeV, 2J=" << twoJ
            << ", P=" << parity
            << ") was snapped to nearest discrete level (Ex="
            << Ex_lev << " MeV, 2J=" << twoJ_lev
            << ", P=" << P_lev << ").";
          if ( !Ex_changed   ) warn_msg << " Ex unchanged.";
          if ( !twoJ_changed ) warn_msg << " 2J unchanged.";
          if ( !P_changed    ) warn_msg << " Parity unchanged.";
          warn_msg << " This warning will not be repeated.";
          MARLEY_LOG( WARN, "cmd.decay" ) << warn_msg.str();
          warned_snap = true;

          if ( sample_Ex && !warned_snap_sample ) {
            MARLEY_LOG( WARN, "cmd.decay" ) << "Discrete level matching"
              " will continue for subsequent events as Ex values below the"
              " unbound threshold (" << unbound_threshold << " MeV) are"
              " sampled. This message will not be repeated.";
            warned_snap_sample = true;
          }
        }

        // Apply the snap — use the level's quantum numbers
        Ex    = Ex_lev;
        twoJ  = twoJ_lev;
        parity = P_lev;
      }
    }

    // --- Build the HepMC3 event object ------------------------------------
    //
    // Primary vertex layout:
    //
    //   IN:  dummy projectile (PDG 0, status 4, zero 4-momentum)
    //        target nucleus   (status 20, at rest, gs mass)
    //
    //   OUT: dummy ejectile   (PDG 0, status 1, zero 4-momentum)
    //        residue nucleus  (status 27, mass = gs + Ex, with Ex/twoJ/parity
    //                          attributes for NucleusDecayer)

    auto event = std::make_shared< HepMC3::GenEvent >(
      HepMC3::Units::MEV, HepMC3::Units::CM );

    event->add_attribute( "signal_process_id",
      std::make_shared< HepMC3::IntAttribute >( signal_proc_id ) );

    auto prim_vtx = std::make_shared< HepMC3::GenVertex >();
    prim_vtx->set_status( marley_hepmc3::NUHEPMC_PRIMARY_VERTEX );
    event->add_vertex( prim_vtx );

    // Dummy projectile — PDG 0, zero 4-momentum
    auto projectile = marley_hepmc3::make_particle( 0, 0., 0., 0., 0.,
      marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, 0. );

    // Target nucleus — at rest, ground-state mass
    auto target = marley_hepmc3::make_particle( nucleus_pdg,
      marley_hepmc3::NUHEPMC_TARGET_STATUS, gs_mass );

    // Dummy ejectile — clone of projectile (PDG 0, zero 4-momentum), final-state
    auto ejectile = marley_hepmc3::make_particle( 0, 0., 0., 0., 0.,
      marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, 0. );

    // Residue nucleus — excited state
    double m_residue = gs_mass + Ex;
    auto residue = marley_hepmc3::make_particle( nucleus_pdg, 0., 0., 0.,
      m_residue, marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, m_residue );

    prim_vtx->add_particle_in( projectile );
    prim_vtx->add_particle_in( target );
    prim_vtx->add_particle_out( ejectile );
    prim_vtx->add_particle_out( residue );

    // Set attributes for the target + residue now that these particles
    // are attached to the event
    marley_hepmc3::set_particle_charge( *target, net_charge );
    marley_hepmc3::set_particle_charge( *residue, net_charge );

    residue->add_attribute( "Ex",
      std::make_shared< HepMC3::DoubleAttribute >( Ex ) );
    residue->add_attribute( "twoJ",
      std::make_shared< HepMC3::IntAttribute >( twoJ ) );
    residue->add_attribute( "parity",
      std::make_shared< HepMC3::IntAttribute >( static_cast<int>(parity) ) );

    // --- Run de-excitation cascade ----------------------------------------
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
