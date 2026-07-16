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
#include "cmd_helpers.hh"
#include "marley/Error.hh"
#include "marley/EventFileReader.hh"
#include "marley/hepmc3_utils.hh"

void for_each_event(
  const std::vector< std::string >& input_files,
  std::function< void( HepMC3::GenEvent&, bool,
    double, const std::shared_ptr< HepMC3::GenRunInfo >& ) > callback )
{
  std::shared_ptr< HepMC3::GenRunInfo > first_run_info;
  double flux_avg_xsec = 0.;
  bool first_event = true;
  bool first_file = true;

  for ( const auto& input_file : input_files ) {
    marley::EventFileReader reader( input_file );
    HepMC3::GenEvent ev;
    while ( reader >> ev ) {
      if ( !first_run_info ) {
        first_run_info = ev.run_info();
        flux_avg_xsec = reader.flux_averaged_xsec( true );
      }
      else if ( !first_file ) {
        std::string issue = marley_hepmc3::check_run_info_compatibility(
          *first_run_info, *ev.run_info() );
        if ( !issue.empty() ) {
          throw marley::Error( "File '" + input_file
            + "' has incompatible run information: " + issue );
        }
      }

      callback( ev, first_event, flux_avg_xsec, first_run_info );
      first_event = false;
    }
    first_file = false;
  }
}
