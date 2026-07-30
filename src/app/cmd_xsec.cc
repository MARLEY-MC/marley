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
#include <fstream>
#include <iostream>
#include <string>

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/Logger.hh"
#include "marley/marley_utils.hh"

namespace {

  constexpr double DEFAULT_KE_MIN = 0.;
  constexpr double DEFAULT_KE_MAX = 100.;
  constexpr int DEFAULT_NUM_STEPS = 10000;
  constexpr int DEFAULT_PDG = marley_utils::ELECTRON_NEUTRINO;

}

bool marley::CommandHandler::cmd_xsec( std::deque< std::string >& args ) {

  std::string output_path;
  std::string config_file_path;
  bool force = false;

  while ( !args.empty() ) {
    std::string arg = args.front();
    args.pop_front();

    if ( arg == "-o" || arg == "--output" ) {
      if ( args.empty() ) {
        std::cerr << "marley xsec: missing argument after '" << arg << "'\n";
        return false;
      }
      output_path = args.front();
      args.pop_front();
    }
    else if ( arg == "-f" || arg == "--force" ) {
      force = true;
    }
    else if ( arg == "-h" || arg == "--help" ) {
      args.clear();
      args.push_front( "xsec" );
      return marley::CommandHandler::cmd_help( args );
    }
    else if ( arg.front() == '-' ) {
      std::cerr << "marley xsec: unrecognized option '" << arg << "'\n";
      return false;
    }
    else if ( config_file_path.empty() ) {
      config_file_path = arg;
    }
    else {
      std::cerr << "marley xsec: unexpected extra argument '"
        << arg << "'\n";
      return false;
    }
  }

  if ( output_path.empty() ) {
    std::cerr << "marley xsec: missing required output file\n";
    args.push_front( "xsec" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  if ( config_file_path.empty() ) {
    std::cerr << "marley xsec: missing required configuration file\n";
    args.push_front( "xsec" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  if ( !force ) {
    std::ifstream temp_stream( output_path );
    if ( temp_stream ) {
      bool overwrite = marley_utils::prompt_yes_no(
        "Really overwrite " + output_path + '?');
      if ( !overwrite ) {
        std::cout << "Total cross section dump aborted.\n";
        return true;
      }
    }
  }

  std::ofstream out_file( output_path );

  marley::JSONConfig config( config_file_path );
  marley::Generator gen = config.create_generator();

  double KEmin = DEFAULT_KE_MIN;
  double KEmax = DEFAULT_KE_MAX;
  int num_steps = DEFAULT_NUM_STEPS;
  int projectile_pdg = DEFAULT_PDG;

  const marley::JSON& json = config.get_json();

  if ( json.has_key("xsec") ) {
    const marley::JSON& xsec_settings = json.at( "xsec" );

    if ( xsec_settings.has_key("KEmin") ) {
      bool ok = false;
      KEmin = xsec_settings.at("KEmin").to_double( ok );
      if ( !ok ) throw marley::Error("Unrecognized KEmin value "
        + xsec_settings.at("KEmin").to_string() + " encountered in the"
        " \"xsec\" section of the job configuration file.");
    }

    if ( xsec_settings.has_key("KEmax") ) {
      bool ok = false;
      KEmax = xsec_settings.at("KEmax").to_double( ok );
      if ( !ok ) throw marley::Error("Unrecognized KEmax value "
        + xsec_settings.at("KEmax").to_string() + " encountered in the"
        " \"xsec\" section of the job configuration file.");
    }

    if ( xsec_settings.has_key("steps") ) {
      bool ok = false;
      num_steps = xsec_settings.at("steps").to_long( ok );
      if ( !ok ) throw marley::Error("Unrecognized steps value "
        + xsec_settings.at("steps").to_string() + " encountered in the"
        " \"xsec\" section of the job configuration file.");
    }

    if ( xsec_settings.has_key("pdg") ) {
      bool ok = false;
      projectile_pdg = xsec_settings.at("pdg").to_long( ok );
      if ( !ok ) throw marley::Error("Unrecognized pdg value "
        + xsec_settings.at("pdg").to_string() + " encountered in the"
        " \"xsec\" section of the job configuration file.");
    }
  }

  double KE = KEmin;
  int steps = num_steps;

  if ( steps <= 1 ) {
    double xsec = gen.total_xs( projectile_pdg, KE );
    xsec *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2 * 1e2;

    out_file << KE << ' ' << xsec << '\n';

    MARLEY_LOG( INFO, "app" ) << "KE = " << KE
      << " MeV, abundance-weighted total xsec = "
      << xsec << " × 10^{-42} cm^2 / atom";
  }
  else {
    double delta = ( KEmax - KEmin ) / ( steps - 1 );
    for ( int s = 0; s < steps; ++s ) {
      KE = KEmin + s * delta;
      double xsec = gen.total_xs( projectile_pdg, KE );
      xsec *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2 * 1e2;

      out_file << KE << ' ' << xsec << '\n';

      MARLEY_LOG( INFO, "app" ) << "KE = " << KE
        << " MeV, abundance-weighted total xsec = "
        << xsec << " × 10^{-42} cm^2 / atom";
    }
  }

  return true;
}
