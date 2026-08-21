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
#include <iostream>

// Geant4 includes
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
//#include "globals.hh"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/JSONConfig.hh"
#include "marley/hepmc3_utils.hh"

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// marg4 includes
#include "MarleyPrimaryGeneratorAction.hh"

MarleyPrimaryGeneratorAction::MarleyPrimaryGeneratorAction(
  const std::string& config_file_name) : G4VUserPrimaryGeneratorAction()
{
  // Create a new marley::Generator object using the settings from the
  // configuration file.
  marley::JSONConfig config( config_file_name );

  marley_generator_= config.create_generator();
}

void MarleyPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Create a new primary vertex at the spacetime origin.
  G4PrimaryVertex* vertex = new G4PrimaryVertex(0., 0., 0., 0.); // x,y,z,t0

  // Generate a new MARLEY event using the owned marley::Generator object
  auto ev = marley_generator_.create_event();

  // Account for possibly different systems of units for the particle
  // 4-momenta by using this conversion factor
  auto mom4_conv_factor = MeV;

  // Query the event to determine what it uses for energy/momentum units
  auto ev_energy_unit = ev->momentum_unit();
  if ( ev_energy_unit == HepMC3::Units::GEV ) {
    mom4_conv_factor = GeV;
  }
  else if ( ev_energy_unit != HepMC3::Units::MEV ) {
    throw marley::Error( "Unrecognized momentum unit encountered in"
      " MarleyPrimaryGeneratorAction::GeneratePrimaries()" );
  }

  // Collect the final-state particles from the HepMC3 event
  auto finals = marley_hepmc3::get_particles_with_status(
    marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, *ev);

  // Loop over each of the final particles in the MARLEY event
  for ( const auto& fp : finals ) {

    // Access the 4-momentum from the HepMC3 particle
    const auto& mom = fp->momentum();

    // Convert each one from a HepMC3::GenParticle into a G4PrimaryParticle.
    // Do this by first setting the PDG code and the 4-momentum components.
    G4PrimaryParticle* particle = new G4PrimaryParticle( fp->pid(),
      mom.px() * mom4_conv_factor,
      mom.py() * mom4_conv_factor,
      mom.pz() * mom4_conv_factor,
      mom.e() * mom4_conv_factor );

    // Also set the charge of the G4PrimaryParticle appropriately.
    // Use the MARLEY utility function that checks for a stored charge
    // attribute before falling back to a PDG-code-based lookup.
    particle->SetCharge( marley_hepmc3::get_particle_charge( *fp ) );

    // Add the fully-initialized G4PrimaryParticle to the primary vertex
    vertex->SetPrimary( particle );
  }

  // The primary vertex has been fully populated with all final-state particles
  // from the MARLEY event. Add it to the G4Event object so that Geant4 can
  // begin tracking the particles through the simulated geometry.
  anEvent->AddPrimaryVertex( vertex );
}
