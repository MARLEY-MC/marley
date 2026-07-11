// Standard library includes
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/Error.hh"
#include "marley/EventFileReader.hh"
#include "marley/JSON.hh"
#include "marley/OutputFile.hh"
#include "marley/marley_utils.hh"

namespace {

  std::string check_compatibility(
    const HepMC3::GenRunInfo& ref,
    const HepMC3::GenRunInfo& candidate)
  {
    if (ref.weight_names() != candidate.weight_names()) {
      return "weight names differ between files";
    }

    auto ref_attr = ref.attribute<HepMC3::StringAttribute>(
      "MARLEY.JSONconfig");
    auto cand_attr = candidate.attribute<HepMC3::StringAttribute>(
      "MARLEY.JSONconfig");

    if (ref_attr && cand_attr) {
      try {
        marley::JSON ref_json = marley::JSON::load(ref_attr->value());
        marley::JSON cand_json = marley::JSON::load(cand_attr->value());

        auto strip_run_keys = [](marley::JSON& j) -> marley::JSON {
          marley::JSON result = marley::JSON::object();
          for (const auto& [key, value] : j.object_range()) {
            if (key != "seed" && key != "executable_settings") {
              result[key] = value;
            }
          }
          return result;
        };

        marley::JSON ref_stripped = strip_run_keys(ref_json);
        marley::JSON cand_stripped = strip_run_keys(cand_json);

        if (ref_stripped.dump_string() != cand_stripped.dump_string()) {
          return "MARLEY JSON configuration differs between files";
        }
      } catch (const std::exception& e) {
        return "failed to parse MARLEY JSON configuration: "
          + std::string(e.what());
      }
    } else if (static_cast<bool>(ref_attr)
      != static_cast<bool>(cand_attr))
    {
      return "one file has MARLEY JSON configuration and the other does not";
    }

    return {};
  }

  std::shared_ptr<HepMC3::GenRunInfo> make_cleaned_run_info(
    std::shared_ptr<HepMC3::GenRunInfo> original)
  {
    auto cleaned = std::make_shared<HepMC3::GenRunInfo>(*original);
    cleaned->remove_attribute("MARLEY.RNGseed");
    cleaned->remove_attribute("MARLEY.JSONconfig");
    return cleaned;
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

  if ( output_format != "ascii" && output_format != "root" ) {
    std::cerr << "marley convert: invalid output format '"
      << output_format << "'. Supported formats: ascii, root\n";
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

  std::string out_config_str = "{ format: \"" + output_format
    + "\", file: \"" + output_path
    + "\", mode: \"overwrite\", force: true }";
  auto out_config = marley::JSON::load( out_config_str );
  auto output_file = marley::OutputFile::make_OutputFile( out_config );

  bool multi_file = (input_files.size() > 1);
  std::shared_ptr<HepMC3::GenRunInfo> first_raw_run_info;
  std::shared_ptr<HepMC3::GenRunInfo> cleaned_run_info;
  bool first_file = true;

  for ( const auto& input_file : input_files ) {
    marley::EventFileReader reader( input_file );
    HepMC3::GenEvent ev;
    while ( reader >> ev ) {
      if ( !first_raw_run_info ) {
        first_raw_run_info = ev.run_info();
        if ( multi_file ) {
          cleaned_run_info = make_cleaned_run_info( first_raw_run_info );
        }
      } else if ( !first_file ) {
        std::string issue = check_compatibility(
          *first_raw_run_info, *ev.run_info() );
        if ( !issue.empty() ) {
          throw marley::Error( "File '" + input_file
            + "' has incompatible run information: " + issue );
        }
      }

      if ( cleaned_run_info ) {
        ev.set_run_info( cleaned_run_info );
      }
      output_file->write_event( &ev );
    }
    first_file = false;
  }

  return true;
}
