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
#include <iostream>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/PrintStreams.h"

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/EventFileReader.hh"
#include "marley/hepmc3_utils.hh"

namespace {

  enum class PrintFormat {
    Pretty,
    HepMC3,
    Legacy
  };

  void print_particle_info( HepMC3::GenParticle& p ) {
    const auto& p4 = p.momentum();
    std::cout << "  particle with PDG code = " << p.pid()
      << " has total energy " << p4.e() << " MeV,"
      << '\n' << "    3-momentum = (" << p4.px() << " MeV, " << p4.py()
      << " MeV, " << p4.pz() << " MeV)," << '\n'
      << "    mass = " << p.generated_mass() << " MeV, and charge = "
      << marley_hepmc3::get_particle_charge( p )
      << " times the proton charge." << '\n';
  }

  void print_event_info( HepMC3::GenEvent& e, const size_t num ) {

    const auto initials_proj = marley_hepmc3::get_particles_with_status(
      marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, e );
    const auto initials_targ = marley_hepmc3::get_particles_with_status(
      marley_hepmc3::NUHEPMC_TARGET_STATUS, e );
    const auto finals = marley_hepmc3::get_particles_with_status(
      marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, e );
    auto residue = marley_hepmc3::get_residue( e );

    size_t num_initial = initials_proj.size() + initials_targ.size();
    size_t num_final = finals.size();

    std::cout << "\n*** Event " << num << " has "
      << num_initial << " initial particles and "
      << num_final << " final particles. ***" << '\n';

    int twoJ = 0;
    double Ex = 0.;
    int parity = 0;
    if ( residue ) {
      auto Ex_attr = residue->attribute< HepMC3::DoubleAttribute >( "Ex" );
      if ( Ex_attr ) Ex = Ex_attr->value();
      auto twoJ_attr = residue->attribute< HepMC3::IntAttribute >( "twoJ" );
      if ( twoJ_attr ) twoJ = twoJ_attr->value();
      auto parity_attr = residue->attribute< HepMC3::IntAttribute >( "parity" );
      if ( parity_attr ) parity = parity_attr->value();
    }
    bool twoJ_is_odd = ( twoJ % 2 == 1 );

    std::cout << "The residual nucleus initially had excitation energy "
      << Ex << " MeV and spin-parity ";
    if ( twoJ_is_odd ) std::cout << twoJ << "/2";
    else std::cout << twoJ / 2;
    std::cout << ( parity >= 0 ? '+' : '-' ) << '\n';

    std::cout << "Initial particles" << '\n';
    for ( const auto& particle_i : initials_proj ) {
      print_particle_info( *particle_i );
    }
    for ( const auto& particle_i : initials_targ ) {
      print_particle_info( *particle_i );
    }

    std::cout << "Final particles" << '\n';
    for ( const auto& particle_f : finals ) {
      print_particle_info( *particle_f );
    }
  }

}

bool marley::CommandHandler::cmd_print( std::deque< std::string >& args ) {

  std::string first_arg;
  if ( !args.empty() ) first_arg = args.front();

  if ( first_arg.empty() || first_arg == "-h" || first_arg == "--help" ) {
    args.clear();
    args.push_front( "print" );
    return marley::CommandHandler::cmd_help( args );
  }

  PrintFormat format = PrintFormat::Pretty;
  if ( first_arg == "pretty" ) {
    // Pretty is the default format (set above), so just
    // drop the format specifier from the input arguments
    args.pop_front();
  }
  else if ( first_arg == "hepmc3" ) {
    format = PrintFormat::HepMC3;
    args.pop_front();
  }
  else if ( first_arg == "legacy" ) {
    format = PrintFormat::Legacy;
    args.pop_front();
  }

  // If there are no input files listed, print the help message
  // and signal that an error condition was encountered
  if ( args.empty() ) {
    args.push_front( "print" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  // Now do the actual printing by iterating over the events in each
  // input file specified in the remaining command-line arguments
  for ( const auto& file_name : args ) {
    marley::EventFileReader reader( file_name );
    HepMC3::GenEvent ev;
    int event_number = 0;
    while ( reader >> ev ) {
      if ( format == PrintFormat::Pretty ) {
        marley_hepmc3::print_event( ev );
      }
      else if ( format == PrintFormat::HepMC3 ) {
        std::cout << ev;
      }
      else {
        // Legacy format
        print_event_info( ev, event_number );
        ++event_number;
      }
    }
  }
  return true;
}
