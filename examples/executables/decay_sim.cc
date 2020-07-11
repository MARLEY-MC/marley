/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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
#include <array>
#include <iostream>

#ifdef USE_ROOT
  // ROOT includes
  #include "TFile.h"
  #include "TParameter.h"
  #include "TTree.h"
#endif

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/Parity.hh"
#include "marley/Reaction.hh"

#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#else
  #include "marley/JSONConfig.hh"
#endif

using ProcType = marley::Reaction::ProcessType;

// Assume that the target is an atom and thus has zero net charge
constexpr int TARGET_NET_CHARGE = 0;

// TODO: Refactor the marley::JSONConfig class to use something like
// these templates

// Extract a parameter of an arbitrary type from a JSON object with
// error handling
template <typename T> T get_param_from_object(const marley::JSON& obj,
  const std::string& obj_key, const std::string& param_key)
{
  if ( !obj.has_key(param_key) ) {
    std::string msg( "Missing \"" );
    msg += param_key + "\" key in the";
    if ( obj_key.empty() ) msg += " JSON configuration file";
    else msg += " \"" + obj_key + "\" object";

    throw marley::Error( msg );
  }

  bool ok = false;
  const marley::JSON& param = obj.at( param_key );

  // Helper variables
  T result;
  int proc_int;
  std::string parity_str;

  // TODO: "if constexpr" is a C++17 feature. Keep this in mind if you
  // want to use these tricks in the main MARLEY code base.
  // Based on https://tinyurl.com/type-traits
  if constexpr ( std::is_floating_point<T>::value ) {
    result = param.to_double( ok );
  }
  else if constexpr ( std::is_integral<T>::value ) {
    result = param.to_long( ok );
  }
  else if constexpr ( std::is_same<T, ProcType>::value ) {
    proc_int = param.to_long( ok );
  }
  else if constexpr ( std::is_same<T, std::string>::value )
  {
    result = param.to_string( ok );
  }
  else if constexpr ( std::is_same<T, marley::Parity>::value )
  {
    parity_str = param.to_string( ok );
  }
  else if constexpr ( std::is_same<T, bool>::value ) {
    result = param.to_bool( ok );
  }
  else if constexpr ( std::is_same<T, marley::JSON>::value ) {
    return param;
  }

  // No else statement with a static_assert needed here. See discussion
  // at https://stackoverflow.com/q/38304847

  if ( !ok ) {
    std::string combined_label = obj_key + '/' + param_key;
    marley::JSONConfig::handle_json_error( combined_label, param );
  }

  if constexpr ( std::is_same<T, ProcType>::value ) {
    result = static_cast<ProcType>( proc_int );
  }
  else if constexpr ( std::is_same<T, marley::Parity>::value ) {
    std::stringstream temp_ss;
    temp_ss >> result;
  }

  return result;
}

// Similar to get_param_from_object_def, but return a default value if
// the requested parameter key is not present in the JSON object
template <typename T> T get_param_from_object_def(const marley::JSON& obj,
  const std::string& obj_key, const std::string& param_key,
  const T& default_value)
{
  if ( !obj.has_key(param_key) ) return default_value;
  return get_param_from_object<T>( obj, obj_key, param_key );
}

// List of nuclear processes to consider. Other processes will not be
// considered by this program.
constexpr std::array< ProcType, 3 > nuclear_proc_types = { ProcType::NeutrinoCC,
  ProcType::AntiNeutrinoCC, ProcType::NC };


int main( int argc, char* argv[] ) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if ( argc < 2 ) {
    std::cout << "Usage: " << argv[0] << " CONFIG_FILE\n";
    return 1;
  }

  // Get the configuration file name from the command line
  std::string config_file_name( argv[1] );

  // Set up the generator using the job configuration file
  #ifdef USE_ROOT
    marley::RootJSONConfig jc( config_file_name );
  #else
    marley::JSONConfig jc( config_file_name );
  #endif

  marley::Generator gen = jc.create_generator();

  // Parse the extra "decay sim" configuration parameters
  // TODO: This is hacky. You can make this better by refactoring some
  // other parts of the code that read similar configurations (e.g., the
  // output specification in the "executable_settings" object used by
  // the main MARLEY executable). Lots of repetition here too.
  const marley::JSON& json = jc.get_json();

  // Check that there is an object with the right label in the configuration
  // file
  const std::string decay_config_label( "decays" );

  const auto& decays = get_param_from_object<marley::JSON>(
    json, "", decay_config_label );

  // Get the total number of events to be simulated
  long num_events = get_param_from_object<long>( decays, decay_config_label,
    "events" );

  // Get the projectile PDG code. If it is missing, assume it's an electron
  // neutrino.
  int projectile_pdg = get_param_from_object_def<int>( decays, decay_config_label,
    "projectile", marley_utils::ELECTRON_NEUTRINO );

  // Get the target proton and nucleon numbers
  int Zi = get_param_from_object<int>( decays, decay_config_label,
    "target_Z" );

  int Ai = get_param_from_object<int>( decays, decay_config_label,
    "target_A" );

  // Get the process type for the primary 2 --> 2 reaction to assume.
  // If it is absent, assume a NC process.
  auto proc_type = get_param_from_object<ProcType>( decays, decay_config_label,
    "proc_type" );

  auto iter = std::find( nuclear_proc_types.cbegin(),
    nuclear_proc_types.cend(), proc_type );

  bool process_is_a_nuclear_reaction = ( iter != nuclear_proc_types.cend() );

  if ( !process_is_a_nuclear_reaction ) {
    std::cerr << "This program handles nuclear reactions only\n";
    return 2;
  }

  // Get the initial excitation energy
  double Ex = get_param_from_object<double>( decays, decay_config_label, "Ex" );
  if ( Ex < 0. ) {
    throw marley::Error( "Negative excitation energy encountered" );
  }

  // Get two times the initial nuclear spin
  int twoJ = get_param_from_object<int>( decays, decay_config_label, "twoJ" );
  if ( twoJ < 0 ) {
    throw marley::Error( "Negative nuclear spin value encountered" );
  }

  // Get the initial nuclear parity
  auto parity = get_param_from_object<marley::Parity>( decays,
    decay_config_label, "parity" );

  // Get the name to use for the output file
  auto output_file_name = get_param_from_object<std::string>( decays,
    decay_config_label, "output_file_name" );

///////////////////////////////////

  // Determine the final-state particle PDG codes based on the process type
  int ejectile_pdg = marley::Reaction
    ::get_ejectile_pdg( projectile_pdg, proc_type );

  // Use conservation of electric charge to obtain the residue net charge
  int projectile_charge = marley_utils::get_particle_charge( projectile_pdg );
  int ejectile_charge = marley_utils::get_particle_charge( ejectile_pdg );

  // Determine the residue nucleon and proton numbers
  int Af = Ai;
  int Delta_Z = projectile_charge - ejectile_charge;
  int Zf = Zi + Delta_Z;

  int residue_net_charge = Delta_Z + TARGET_NET_CHARGE;

  #ifdef USE_ROOT
    // Create a ROOT tree to store the events
    TFile out_file( output_file_name.c_str(), "recreate" );

    TTree* out_tree = new TTree( "MARLEY_event_tree",
      "Compound nucleus decay events generated by MARLEY" );

    // We create a branch to store the events here, but set the branch
    // address to nullptr. This will be fixed later when write_event()
    // is called.
    out_tree->Branch( "event", "marley::Event", nullptr );
  #else
    std::ofstream out_ascii_file( output_file_name );
  #endif

  // Write a dummy flux-averaged cross section value to the output file. This
  // keeps utilities like marsum happy.
  #ifdef USE_ROOT
    TParameter<double> dummy_xsec( "MARLEY_flux_avg_xsec", 0. );
    out_file.cd();
    dummy_xsec.Write();
  #else
    double dummy_xsec = 0.;
    out_ascii_file << dummy_xsec << '\n';
  #endif

  for ( long evnum = 0; evnum < num_events; ++evnum ) {

    // Get access to the MassTable
    const auto& mt = marley::MassTable::Instance();

    // If the requested excitation energy corresponds to a bound nuclear
    // state, match it to a known discrete level instead of taking
    // the input at face value
    // TODO: Revisit this. Maybe do something more sophisticated.
    if ( Ex <= mt.unbound_threshold(Zf, Af) ) {
      // Try to load a discrete level DecayScheme from the generator's owned
      // StructureDatabase. If we can find one, use it to assign a refined
      // excitation energy value. If not, give up and move on.
      auto* ds = gen.get_structure_db().get_decay_scheme( Zf, Af );

      if ( ds ) {
        // Replace our current excitation energy value with the closest
        // known discrete level
        auto* lev = ds->get_pointer_to_closest_level( Ex );
        // TODO: perhaps check that the spin-parity of the level matches
        // the expected one, issue a warning if it doesn't?
        Ex = lev->energy();
      }
    }

    // Create the event object with a dummy projectile, target, and ejectile.
    int target_pdg = marley_utils::get_nucleus_pid( Zi, Ai );
    int residue_pdg = marley_utils::get_nucleus_pid( Zf, Af );

    marley::Particle projectile( projectile_pdg, 0., 0., 0., 0., 0 );
    marley::Particle target( target_pdg, 0., 0., 0., 0., 0., 0 );
    marley::Particle ejectile( ejectile_pdg, 0., 0., 0., 0., 0., 0 );

    // Approximate the ground-state mass of the (possibly ionized) residue by
    // subtracting the appropriate number of electron masses from its atomic
    // (i.e., neutral) ground state mass.
    double residue_gs_mass = mt.get_atomic_mass( residue_pdg )
      - ( residue_net_charge * mt.get_particle_mass(marley_utils::ELECTRON) );

    // Create a Particle object representing the residue. Since we're
    // simulating its decays without worrying about a primary 2 --> 2
    // interaction, let it be at rest in the lab frame.
    double m_residue = residue_gs_mass + Ex;
    marley::Particle residue( residue_pdg, m_residue, 0., 0., 0.,
      m_residue, residue_net_charge );

    marley::Event event( projectile, target, ejectile,
      residue, Ex, twoJ, parity );

    // Pass the event to the nuclear de-excitation simulation
    marley::NucleusDecayer nd;
    nd.process_event( event, gen );

    // We're done
    #ifdef USE_ROOT
      // Write the new event to the output ROOT TTree
      marley::Event* temp_event_ptr = &event;
      out_tree->SetBranchAddress( "event", &temp_event_ptr );
      out_tree->Fill();
    #else
      out_ascii_file << event << '\n';
    #endif

    std::cout << "Event " << evnum << ":\n" << event << '\n';
  }

  #ifdef USE_ROOT
    out_tree->Write();
  #endif

}
