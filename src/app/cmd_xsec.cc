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

  void get_double_dump_param( const marley::JSON& json,
    const std::string& param_key, double& value )
  {
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
    if ( !json.has_key(param_key) ) return;

    const auto& temp_js = json.at( param_key );
    bool ok = false;
    value = temp_js.to_long( ok );
    if ( !ok ) throw marley::Error("Unrecognized " + param_key
      + " value " + temp_js.to_string() + " encountered in the"
      " job configuration file.");
  }

}

bool marley::CommandHandler::cmd_xsec( std::deque< std::string >& args ) {

  // If we have an unexpected number of arguments, decide whether the
  // user intended to request help with this command
  if ( args.size() != 2u ) {
    std::string first_arg;
    if ( !args.empty() ) first_arg = args.front();

    // Print the help message either way
    args.clear();
    args.push_front( "xsec" );
    marley::CommandHandler::cmd_help( args );

    // Return a boolean status based on whether the help message was
    // explicitly requested (normal behavior) or not (an error condition)
    if ( first_arg == "-h" || first_arg == "--help" ) return true;
    return false;
  }

  // We know that args has exactly two elements if we make it here
  std::string output_file_name( args.front() );
  std::string config_file_name( args.back() );

  std::ifstream temp_stream( output_file_name );
  if ( temp_stream ) {
    bool overwrite = marley_utils::prompt_yes_no(
      "Really overwrite " + output_file_name + '?');
    if ( !overwrite ) {
      std::cout << "Total cross section dump aborted.\n";
      return true;
    }
  }

  std::ofstream out_file( output_file_name );

  marley::JSONConfig config( config_file_name );
  marley::Generator gen = config.create_generator();

  double KEmin = DEFAULT_KE_MIN;
  double KEmax = DEFAULT_KE_MAX;
  int num_steps = DEFAULT_NUM_STEPS;
  int projectile_pdg = DEFAULT_PDG;

  const marley::JSON& json = config.get_json();

  get_double_dump_param( json, "xsec_dump_KEmin", KEmin );
  get_double_dump_param( json, "xsec_dump_KEmax", KEmax );

  get_int_dump_param( json, "xsec_dump_steps", num_steps );
  get_int_dump_param( json, "xsec_dump_pdg", projectile_pdg );

  double delta_KE_step = ( KEmax - KEmin ) / num_steps;
  double KE = KEmin;

  for ( int s = 0; s < num_steps; ++s ) {

    KE += delta_KE_step;
    double xsec = gen.total_xs( projectile_pdg, KE );
    xsec *= marley_utils::hbar_c2 * marley_utils::fm2_to_minus40_cm2 * 1e2;

    out_file << KE << ' ' << xsec << '\n';

    MARLEY_LOG( INFO, "app" ) << "KE = " << KE
      << " MeV, abundance-weighted total xsec = "
      << xsec << " × 10^{-42} cm^2 / atom";

  }

  return true;
}
