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
#include <iostream>
#include <string>

// MARLEY includes
#include "marley/Event.hh"
#ifdef USE_ROOT
  #include "marley/RootEventFileReader.hh"
#else
  #include "marley/EventFileReader.hh"
#endif

// These functions largely duplicate code that also exists
// in marley::Event::print_human_readable(). They are included
// here as an example of code that accesses all data members
// of the Particle and Event objects.
void print_particle_info( const marley::Particle& p ) {
  std::cout << "  particle with PDG code = " << p.pdg_code()
    << " has total energy " << p.total_energy() << " MeV,"
    << '\n' << "    3-momentum = (" << p.px() << " MeV, " << p.py()
    << " MeV, " << p.pz() << " MeV)," << '\n'
    << "    mass = " << p.mass() << " MeV, and charge = "
    << p.charge() << " times the proton charge." << '\n';
}

void print_event_info(const marley::Event& e, const size_t num) {

  const std::vector< marley::Particle* > initials = e.get_initial_particles();
  const std::vector< marley::Particle* > finals = e.get_final_particles();

  size_t num_initial = initials.size();
  size_t num_final = finals.size();

  std::cout << "\n*** Event " << num << " has "
    << num_initial << " initial particles and "
    << num_final << " final particles. ***" << '\n';

  int twoJ = e.twoJ();
  bool twoJ_is_odd = ( twoJ % 2 == 1 );
  std::cout << "The residual nucleus initially had excitation energy "
    << e.Ex() << " MeV and spin-parity ";
  if ( twoJ_is_odd ) std::cout << twoJ << "/2";
  else std::cout << twoJ / 2;
  std::cout << e.parity() << '\n';

  std::cout << "Initial particles" << '\n';
  for ( const auto* particle_i : initials ) {
    print_particle_info( *particle_i );
  }

  std::cout << "Final particles" << '\n';
  for ( const auto* particle_f : finals ) {
    print_particle_info( *particle_f );
  }
}

// Prints events from an input file in a human-readable format. This function
// also updates an overall event number as it works through the file.
void print_all_events(const std::string& file_name, int& ev_num) {

  #ifdef USE_ROOT
    marley::RootEventFileReader reader( file_name );
  #else
    marley::EventFileReader reader( file_name );
  #endif

  marley::Event ev;

  while ( reader >> ev ) {
    print_event_info( ev, ev_num );
    ++ev_num;
  }
}

int main(int argc, char* argv[]) {

  // If the user has not supplied any command-line arguments, display the
  // standard help message and exit
  if (argc <= 1) {
    std::cout << "Usage: " << argv[0] << " INPUT_FILE...\n";
    return 0;
  }

  std::vector<std::string> file_names;
  for ( int s = 1; s < argc; ++s ) file_names.push_back( argv[s] );

  int event_number = 0;
  for ( const auto& file_name : file_names ) {
    print_all_events( file_name, event_number );
  }

  return 0;
}
