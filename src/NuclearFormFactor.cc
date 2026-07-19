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

// MARLEY includes
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/NuclearFormFactor.hh"

std::shared_ptr< marley::NuclearFormFactor > marley::NuclearFormFactor::create(
  int Z, int A, const marley::JSON& ff_config )
{
  static bool need_to_log_nuc_ff_info = true;

  // If the user has requested use of the allowed approximation, then configure
  // the trivial nuclear form factor model and return without logging this
  // choice (it will be handled elsewhere)
  if ( marley::JSONConfig::check_for_allowed_approximation(ff_config) ) {
    auto tnff = std::make_shared< marley::TrivialNuclearFormFactor >( Z, A );
    need_to_log_nuc_ff_info = false;
    return tnff;
  }

  if ( !ff_config.has_key("nuclear_model") ) {
    throw marley::Error( "Missing nuclear form factor model configuration" );
  }
  const auto& nucl_json = ff_config.at( "nuclear_model" );
  if ( !nucl_json.is_string() ) {
    throw marley::Error( "Invalid nuclear form factor model "
      + nucl_json.dump_string() );
  }
  std::string nucl_ff_model = nucl_json.to_string();

  std::shared_ptr< marley::NuclearFormFactor > nuclear_ff;

  // Choose the model to use for the nuclear form factor
  if ( nucl_ff_model == "trivial" ) {
    nuclear_ff = std::make_shared< marley::TrivialNuclearFormFactor >( Z, A );
    if ( need_to_log_nuc_ff_info ) {
      MARLEY_LOG( INFO, "physics.formfactor" ) << "Using trivial nuclear form factor";
    }
  }
  else if ( nucl_ff_model == "helm" ) {
    nuclear_ff = std::make_shared< marley::HelmNuclearFormFactor >( Z, A );
    if ( need_to_log_nuc_ff_info ) {
      MARLEY_LOG( INFO, "physics.formfactor" ) << "Using Helm nuclear form factor";
    }
  }
  else if ( nucl_ff_model == "klein" ) {

    // Retrieve the optional "nucl_options" JSON object if it is present
    bool ok; // dummy flag used when parsing the JSON configuration
    marley::JSON nucl_opt = assign_from_json< marley::JSON >(
      "nucl_options", ff_config, ok, marley::JSON::object() );

    // Set the "adapted" option from the parameters if it is present
    bool adapted = assign_from_json< bool >( "adapted", nucl_opt, ok, true );

    nuclear_ff = std::make_shared< marley
      ::KleinNystrandNuclearFormFactor >( Z, A, adapted );

    if ( need_to_log_nuc_ff_info ) {
      std::string kn_name;
      if ( adapted ) kn_name += "adapted ";
      kn_name += "Klein-Nystrand";
      MARLEY_LOG( INFO, "physics.formfactor" ) << "Using " << kn_name << " nuclear form factor";
    }
  }
  else throw marley::Error( "Unrecognized nuclear form factor model name \""
    + nucl_ff_model + "\" in marley::NuclearFormFactor::create()" );

  need_to_log_nuc_ff_info = false;
  return nuclear_ff;
}

double marley::HelmNuclearFormFactor::F( double kappa ) const {
  if ( kappa == 0. ) return 1.;
  double kappa_in_inverse_fm = kappa / marley_utils::hbar_c;
  double x = kappa_in_inverse_fm * R_;
  double j1 = ( std::sin( x ) / x - std::cos( x ) ) / x;
  double F = 3. * j1 / x * std::exp( -kappa_in_inverse_fm
    * kappa_in_inverse_fm * s_ * s_ / 2. );
  return F;
}

double marley::KleinNystrandNuclearFormFactor::F( double kappa ) const {
  if ( kappa == 0. ) return 1.;
  double kappa_in_inverse_fm = kappa / marley_utils::hbar_c;
  double x = kappa_in_inverse_fm * R_;
  double j1 = ( std::sin( x ) / x - std::cos( x ) ) / x;
  double F = 3. * j1 / x * ( 1. / ( 1.
    + kappa_in_inverse_fm*kappa_in_inverse_fm*a_*a_ ) );
  return F;
}

void marley::KleinNystrandNuclearFormFactor::initialize_r0_table() {

  // Instantiate the file manager and use it to find the data file containing
  // the rms charge radii for many nuclei
  const auto& fm = marley::FileManager::Instance();
  std::string full_r0_file_name = fm.find_file( r0_data_file_name_ );

  if ( full_r0_file_name.empty() ) {
    throw marley::Error( "Could not find the MARLEY nuclear rms charge radii"
      " data file " + r0_data_file_name_ + ". Please ensure that"
      " the folder containing it is on the MARLEY search path."
      " If needed, the folder can be appended to the MARLEY_SEARCH_PATH"
      " environment variable." );
  }

  MARLEY_LOG( INFO, "init.structure" ) << "Loading ground-state nuclear rms charge radii from "
    << full_r0_file_name;

  marley::JSON r0_json_obj = marley::JSON::load_file( full_r0_file_name );
  if ( !r0_json_obj.has_key("nuclear_charge_radii") ) {
    throw marley::Error( "Missing \"nuclear_charge_radii\" key in "
      + full_r0_file_name );
  }
  const marley::JSON& r0_json_array = r0_json_obj.at( "nuclear_charge_radii" );
  if ( !r0_json_array.is_array() ) {
    throw marley::Error( "Invalid \"nuclear_charge_radii\" array in "
      + full_r0_file_name );
  }

  r0_table_ = std::make_unique< std::map< int, double > >();

  bool ok; // Helper flag to use when reading the JSON array elements
  for ( const auto& r0_js : r0_json_array.array_range() ) {
    int Z = assign_from_json< int >( "Z", r0_js, ok );
    int A = assign_from_json< int >( "A", r0_js, ok );
    double R = assign_from_json< double >( "R", r0_js, ok );

    // TODO: load and use the uncertainty on the rms charge radius
    //double R_unc = assign_from_json< double >( "R_unc", r0_js, ok );

    // Add the completed entry to the table of nuclear radii
    int pdg = marley_utils::get_nucleus_pid( Z, A );
    r0_table_->operator[]( pdg ) = R;

    MARLEY_LOG( TRACE, "init.structure.masstable" ) << "Nucleus with PDG code "
      << pdg << " has rms charge radius " << R << " fm";
  }

}
