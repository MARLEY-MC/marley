/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
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
#include "marley/DecayScheme.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/NeutrinoSource.hh"
#include "marley/NuclearReaction.hh"
#include "marley/StructureDatabase.hh"
#include "marley/marley_utils.hh"

#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#else
  #include "marley/JSONConfig.hh"
#endif

namespace {
  // Default settings for the projectile kinetic energy range and number of
  // steps for the total cross section dump
  constexpr double DEFAULT_KE_MIN = 0.;
  constexpr double DEFAULT_KE_MAX = 100.;
  constexpr int DEFAULT_NUM_STEPS = 10000;

  // Default projectile PDG code
  constexpr int DEFAULT_PDG = marley_utils::ELECTRON_NEUTRINO;

  // Helper functions for loading custom dump parameters from the job
  // configuration file
  // TODO: reduce code duplication between these two functions
  void get_double_dump_param( const marley::JSON& json,
    const std::string& param_key, double& value )
  {
    // If the key doesn't appear in the job configuration file, just return
    // without doing anything (all cross section dump parameters are optional)
    if ( !json.has_key(param_key) ) return;

    const auto& temp_js = json.at( param_key );
    bool ok = false;
    value = temp_js.to_double( ok );
    if ( !ok ) throw marley::Error("Unrecognized " + param_key
      + " value " + temp_js.to_string() + " encountered in the"
      " job configuration file.");
  }

  void get_int_dump_param( const marley::JSON& json,
    const std::string& param_key, int& value )
  {
    // If the key doesn't appear in the job configuration file, just return
    // without doing anything (all cross section dump parameters are optional)
    if ( !json.has_key(param_key) ) return;

    const auto& temp_js = json.at( param_key );
    bool ok = false;
    value = temp_js.to_long( ok );
    if ( !ok ) throw marley::Error("Unrecognized " + param_key
      + " value " + temp_js.to_string() + " encountered in the"
      " job configuration file.");
  }

}

int main(int argc, char* argv[]) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if (argc <= 2) {
    std::cout << "Usage: " << argv[0] << " OUTPUT_FILE CONFIG_FILE\n";
    return 1;
  }

  // Get the output and config file names from the command line
  std::string output_file_name( argv[1] );
  std::string config_file_name( argv[2] );

  // Check whether the output file exists and warn the user before
  // overwriting it if it does
  std::ifstream temp_stream( output_file_name );
  if ( temp_stream ) {
    bool overwrite = marley_utils::prompt_yes_no(
      "Really overwrite " + output_file_name + '?');
    if ( !overwrite ) {
      std::cout << "Total cross section dump aborted.\n";
      return 0;
    }
  }

  // Open the output file for writing
  std::ofstream out_file( output_file_name );

  // Configure a new Generator object
  #ifdef USE_ROOT
    marley::RootJSONConfig config( config_file_name );
  #else
    marley::JSONConfig config( config_file_name );
  #endif
  marley::Generator gen = config.create_generator();

  // Initialize the projectile kinetic energy range and step size to use for
  // the total cross section dump
  double KEmin = DEFAULT_KE_MIN;
  double KEmax = DEFAULT_KE_MAX;
  int num_steps = DEFAULT_NUM_STEPS;
  int projectile_pdg = DEFAULT_PDG;

  // If the user has specified a non-default energy range or step size, then
  // retrieve that information from the config file.
  const marley::JSON& json = config.get_json();

  get_double_dump_param( json, "xsec_dump_KEmin", KEmin );
  get_double_dump_param( json, "xsec_dump_KEmax", KEmax );

  get_int_dump_param( json, "xsec_dump_steps", num_steps );
  get_int_dump_param( json, "xsec_dump_pdg", projectile_pdg );

  // Compute the (linear) step size based on the number of steps and the
  // desired projectile energy range
  double delta_KE_step = ( KEmax - KEmin ) / num_steps;

  // Initial projectile energy to use in the loop
  double KE = KEmin;

  // Ready for the dump now
  // Energies are dumped in units of MeV
  // Cross section values are dumped in units of 10^(-42) cm^2 / atom
  for ( int s = 0; s < num_steps; ++s ) {

    // Update the projectile energy for this step
    KE += delta_KE_step;

    // Compute the total cross section at this projectile
    // energy (summed over all active reactions and weighted
    // by nuclide abundance in the target material)
    double xsec = gen.total_xs( projectile_pdg, KE );

    // Convert the total cross section from MARLEY natural units (MeV^{-2}) to
    // conventional units (10^{-42} cm^2)
    xsec *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2 * 1e2;

    // Write the current kinetic energy and total cross section to the output
    // file
    out_file << KE << ' ' << xsec << '\n';

    // Also write these quantities to the Logger
    MARLEY_LOG_INFO() << "KE = " << KE << " MeV, abundance-weighted total xsec = "
      << xsec << " \u00D7 10^{-42} cm^2 / atom";
  }

  return 0;
}
