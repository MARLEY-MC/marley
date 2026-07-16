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
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "cmd_helpers.hh"
#include "marley/CommandHandler.hh"
#include "marley/Error.hh"
#include "marley/EventFileReader.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/JSON.hh"
#include "marley/OutputFile.hh"
#include "marley/marley_utils.hh"

namespace {

  std::shared_ptr<HepMC3::GenRunInfo> make_cleaned_run_info(
    std::shared_ptr<HepMC3::GenRunInfo> original)
  {
    auto cleaned = std::make_shared<HepMC3::GenRunInfo>(*original);
    cleaned->remove_attribute("MARLEY.RNGseed");
    cleaned->remove_attribute("MARLEY.JSONconfig");
    return cleaned;
  }

  std::shared_ptr<HepMC3::GenParticle> find_legacy_residue(
    std::shared_ptr<HepMC3::GenParticle> residue )
  {
    auto current = residue;
    while ( auto vtx = current->end_vertex() ) {
      std::shared_ptr<HepMC3::GenParticle> daughter;
      for ( const auto& out : vtx->particles_out() ) {
        if ( !marley_utils::is_ion( out->pid() ) ) continue;
        if ( out->status()
          != marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS )
        {
          daughter = out;
          break;
        }
        if ( !daughter ) daughter = out;
      }
      if ( !daughter ) break;
      current = daughter;
    }
    return current;
  }

  struct LegacyEventView {
    std::shared_ptr<HepMC3::GenParticle> projectile;
    std::shared_ptr<HepMC3::GenParticle> target;
    std::shared_ptr<HepMC3::GenParticle> ejectile;
    std::shared_ptr<HepMC3::GenParticle> residue;
    std::shared_ptr<HepMC3::GenParticle> legacy_residue;
    double Ex = 0.;
    int twoJ = 0;
    int parity = 1;
    std::vector< std::shared_ptr< HepMC3::GenParticle > > final_state_particles;
  };

  LegacyEventView extract_legacy_event_view( HepMC3::GenEvent& ev ) {
    LegacyEventView v;
    v.projectile = marley_hepmc3::get_projectile( ev );
    v.target = marley_hepmc3::get_target( ev );
    v.ejectile = marley_hepmc3::get_ejectile( ev );
    v.residue = marley_hepmc3::get_residue( ev );

    if ( v.residue ) {
      auto Ex_a
        = v.residue->attribute< HepMC3::DoubleAttribute >( "Ex" );
      if ( Ex_a ) v.Ex = Ex_a->value();
      auto twoJ_a
        = v.residue->attribute< HepMC3::IntAttribute >( "twoJ" );
      if ( twoJ_a ) v.twoJ = twoJ_a->value();
      auto par_a
        = v.residue->attribute< HepMC3::IntAttribute >( "parity" );
      if ( par_a ) v.parity = par_a->value();
    }

    if ( v.residue ) v.legacy_residue = find_legacy_residue( v.residue );

    v.final_state_particles
      = marley_hepmc3::get_particles_with_status(
        marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, ev );

    return v;
  }

  void write_old_ascii_event( std::ostream& os,
    const LegacyEventView& v )
  {
    int num_final = 2;
    for ( const auto& p : v.final_state_particles ) {
      if ( v.ejectile && p->id() == v.ejectile->id() ) continue;
      if ( v.legacy_residue && p->id() == v.legacy_residue->id() ) continue;
      ++num_final;
    }

    std::ostringstream tmp;
    tmp << std::scientific
      << std::setprecision( std::numeric_limits< double >::max_digits10 );

    tmp << "2 " << num_final << ' '
      << v.Ex << ' ' << v.twoJ << ' ' << ( v.parity >= 0 ? '+' : '-' )
      << '\n';

    auto write_part = [ &tmp ]( const HepMC3::GenParticle& p ) {
      const auto& m = p.momentum();
      int q = marley_hepmc3::get_particle_charge(
        const_cast< HepMC3::GenParticle& >( p ) );
      tmp << p.pid() << ' ' << m.e() << ' ' << m.px() << ' '
        << m.py() << ' ' << m.pz() << ' '
        << p.generated_mass() << ' ' << q << '\n';
    };

    if ( v.projectile ) write_part( *v.projectile );
    if ( v.target ) write_part( *v.target );
    if ( v.ejectile ) write_part( *v.ejectile );
    if ( v.legacy_residue ) write_part( *v.legacy_residue );
    for ( const auto& p : v.final_state_particles ) {
      if ( v.ejectile && p->id() == v.ejectile->id() ) continue;
      if ( v.legacy_residue && p->id() == v.legacy_residue->id() ) continue;
      write_part( *p );
    }

    os << tmp.str();
  }

  void write_hepevt_event( std::ostream& os,
    const LegacyEventView& v, unsigned long event_num,
    double flux_avg_xsec_natural )
  {
    int nhep = 5;
    for ( const auto& p : v.final_state_particles ) {
      if ( v.ejectile && p->id() == v.ejectile->id() ) continue;
      if ( v.legacy_residue && p->id() == v.legacy_residue->id() ) continue;
      ++nhep;
    }

    constexpr double MEV2GEV = 0.001;

    std::ostringstream tmp;
    tmp << std::scientific
      << std::setprecision( std::numeric_limits< double >::max_digits10 );

    tmp << event_num << ' ' << nhep << '\n';

    auto dump_line = [ &tmp ]( const HepMC3::GenParticle& p,
      int status, int jmo1 = 0, int jmo2 = 0 )
    {
      const auto& m = p.momentum();
      tmp << status << ' ' << p.pid() << ' ' << jmo1 << ' ' << jmo2
        << " 0 0 "
        << m.px() * MEV2GEV << ' ' << m.py() * MEV2GEV << ' '
        << m.pz() * MEV2GEV << ' ' << m.e() * MEV2GEV << ' '
        << p.generated_mass() * MEV2GEV
        << " 0. 0. 0. 0." << '\n';
    };

    if ( v.projectile ) dump_line( *v.projectile, 3 );
    if ( v.target ) dump_line( *v.target, 3 );

    tmp << "11 0 " << v.twoJ << ' ' << v.parity << " 0 0 "
      << "0. 0. 0. " << v.Ex << ' ' << flux_avg_xsec_natural
      << " 0. 0. 0. 0." << '\n';

    if ( v.ejectile ) dump_line( *v.ejectile, 1 );
    if ( v.legacy_residue ) dump_line( *v.legacy_residue, 1 );
    for ( const auto& p : v.final_state_particles ) {
      if ( v.ejectile && p->id() == v.ejectile->id() ) continue;
      if ( v.legacy_residue && p->id() == v.legacy_residue->id() ) continue;
      dump_line( *p, 1 );
    }

    os << tmp.str();
  }

  void convert_to_legacy(
    const std::vector< std::string >& input_files,
    const std::string& output_path )
  {
    std::ofstream out( output_path );
    if ( !out ) throw marley::Error( "Could not open output file \""
      + output_path + "\" for writing" );

    out << std::scientific
      << std::setprecision( std::numeric_limits< double >::max_digits10 );

    bool header_written = false;
    for_each_event( input_files,
      [ & ]( HepMC3::GenEvent& ev, bool, double xsec, const auto& ) {
        if ( !header_written ) {
          out << xsec << '\n';
          header_written = true;
        }
        write_old_ascii_event( out,
          extract_legacy_event_view( ev ) );
      } );
  }

  void convert_to_hepevt(
    const std::vector< std::string >& input_files,
    const std::string& output_path )
  {
    std::ofstream out( output_path );
    if ( !out ) throw marley::Error( "Could not open output file \""
      + output_path + "\" for writing" );

    unsigned long ev_num = 0;
    for_each_event( input_files,
      [ & ]( HepMC3::GenEvent& ev, bool, double xsec, const auto& ) {
        write_hepevt_event( out,
          extract_legacy_event_view( ev ), ev_num, xsec );
        ++ev_num;
      } );
  }

}

bool marley::CommandHandler::cmd_convert( std::deque< std::string >& args ) {

  std::string output_path;
  std::string output_format;
  bool force = false;
  std::vector< std::string > input_files;

  while ( !args.empty() ) {
    std::string arg = args.front();
    args.pop_front();

    if ( arg == "-o" || arg == "--output" ) {
      if ( args.empty() ) {
        std::cerr << "marley convert: missing argument after '" << arg
          << "'\n";
        return false;
      }
      output_path = args.front();
      args.pop_front();
    }
    else if ( arg == "--output-format" ) {
      if ( args.empty() ) {
        std::cerr << "marley convert: missing argument after '"
          << "--output-format'\n";
        return false;
      }
      output_format = args.front();
      args.pop_front();
    }
    else if ( arg == "--force" || arg == "-f" ) {
      force = true;
    }
    else if ( arg == "--help" || arg == "-h" ) {
      args.clear();
      args.push_front( "convert" );
      return marley::CommandHandler::cmd_help( args );
    }
    else if ( arg.front() == '-' ) {
      std::cerr << "marley convert: unrecognized option '" << arg << "'\n";
      return false;
    }
    else {
      input_files.push_back( arg );
    }
  }

  if ( output_path.empty() ) {
    std::cerr << "marley convert: missing required option -o OUTPUT_FILE\n";
    args.push_front( "convert" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  if ( input_files.empty() ) {
    std::cerr << "marley convert: no input files specified\n";
    args.push_front( "convert" );
    marley::CommandHandler::cmd_help( args );
    return false;
  }

  if ( output_format.empty() ) {
    if ( output_path.size() >= 5
      && output_path.substr( output_path.size() - 5 ) == ".root" )
    {
      output_format = "root";
    }
    else {
      output_format = "ascii";
    }
  }

  if ( output_format != "ascii" && output_format != "root"
    && output_format != "legacy" && output_format != "hepevt" )
  {
    std::cerr << "marley convert: invalid output format '"
      << output_format << "'. Supported formats:"
      " ascii, root, legacy, hepevt\n";
    return false;
  }

#ifndef USE_ROOT
  if ( output_format == "root" ) {
    std::cerr << "marley convert: ROOT output format requires a"
      " ROOT-enabled build of MARLEY.\n";
    return false;
  }
#endif

  if ( !force ) {
    std::ifstream test( output_path );
    if ( test ) {
      bool overwrite = marley_utils::prompt_yes_no(
        "Really overwrite " + output_path + "?" );
      if ( !overwrite ) {
        std::cout << "Action aborted.\n";
        return true;
      }
    }
  }

  if ( output_format == "legacy" ) {
    convert_to_legacy( input_files, output_path );
    return true;
  }

  if ( output_format == "hepevt" ) {
    convert_to_hepevt( input_files, output_path );
    return true;
  }

  std::string out_config_str = "{ format: \"" + output_format
    + "\", file: \"" + output_path
    + "\", mode: \"overwrite\", force: true }";
  auto out_config = marley::JSON::load( out_config_str );
  auto output_file = marley::OutputFile::make_OutputFile( out_config );

  bool multi_file = (input_files.size() > 1);
  std::shared_ptr< HepMC3::GenRunInfo > cleaned_run_info;
  for_each_event( input_files,
    [ & ]( HepMC3::GenEvent& ev, bool first_event, double,
      const auto& first_info )
    {
      if ( first_event && multi_file ) {
        cleaned_run_info = make_cleaned_run_info( first_info );
      }
      if ( cleaned_run_info ) {
        ev.set_run_info( cleaned_run_info );
      }
      output_file->write_event( &ev );
    } );

  return true;
}
