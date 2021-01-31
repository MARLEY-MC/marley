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
#include <memory>

// Geant4 includes
#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"

// marg4 includes
#include "DetectorConstruction.hh"
#include "MarleyPrimaryGeneratorAction.hh"
#include "EventAction.hh"

namespace {

  // Retrieves the desired number of events from the first command
  // line argument. Based on https://stackoverflow.com/a/2797823/4081973
  // Returns true if everything went well, or false if there was a problem
  bool get_num_events( const std::string& arg, int& x ) {
    try {
      size_t pos;
      x = std::stoi( arg, &pos );
      if ( pos < arg.size() ) {
        std::cerr << "Trailing characters after number of events: " << arg
          << '\n';
        return false;
      }
      if ( x < 0 ) {
        std::cerr << "Negative number of events: " << arg << '\n';
        return false;
      }
    }
    catch (std::invalid_argument const &ex) {
      std::cerr << "Invalid number: " << arg << '\n';
      return false;
    }
    catch (std::out_of_range const &ex) {
      std::cerr << "Number out of range: " << arg << '\n';
      return false;
    }

    return true;
  }

}


int main( int argc, char* argv[] ) {

  if ( argc <= 2 ) {
    std::cout << "Usage: marg4 NUM_EVENTS MARLEY_CONFIG_FILE\n";
    return 1;
  }

  int num_events = 0;
  bool num_ok = get_num_events( argv[1], num_events );
  if ( !num_ok ) return 2;

  // Retrieve the configuration file name from the command-line argument
  std::string config_file_name( argv[2] );

  // Initialize the Geant4 run manager
  std::unique_ptr<G4RunManager> rm( new G4RunManager );

  // ** Set mandatory initialization classes **
  // Define the geometry for the simulation
  DetectorConstruction* det = new DetectorConstruction();
  rm->SetUserInitialization( det );

  // Set up the built-in QGSP_BIC_HP physics list
  // Note: High-precision tracking of neutrons is set up by standard
  // physics lists with the "_HP" suffix. Transport of neutrons generated
  // by MARLEY will not be trustworthy unless the high-precision tracking
  // is enabled.
  G4PhysListFactory factory;
  G4VModularPhysicsList* refList
    = factory.GetReferencePhysList( "QGSP_BIC_HP" );
  rm->SetUserInitialization( refList );

  // ** Set user actions **

  // The primary generator action interfaces with MARLEY
  MarleyPrimaryGeneratorAction* mpga
    = new MarleyPrimaryGeneratorAction( config_file_name );
  rm->SetUserAction( mpga );

  // The event action prints the current event number at the beginning of
  // every hundredth event without doing anything else.
  EventAction* eva = new EventAction;
  rm->SetUserAction( eva );

  rm->Initialize();

  std::cout << "Simulating " << num_events << " events . . .\n";

  rm->BeamOn( num_events );

  return 0;
}
