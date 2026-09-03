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
#include <cmath>

// Geant4 includes
#include "G4Event.hh"
#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RandomDirection.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

// marg4 includes
#include "NeutronSource.hh"

NeutronSource::NeutronSource()
{
}

NeutronSource::~NeutronSource()
{
}

void NeutronSource::GeneratePrimaries(G4Event* anEvent)
{
  G4double kinetic_energy = G4UniformRand() * (energy_max_ - energy_min_) + energy_min_;
  G4ThreeVector direction = G4RandomDirection();

  G4ParticleDefinition* neutron = G4Neutron::Neutron();
  G4double total_mass = neutron->GetPDGMass();

  G4double total_energy = kinetic_energy + total_mass;
  G4double momentum = std::sqrt(total_energy * total_energy - total_mass * total_mass);

  G4double px = direction.x() * momentum;
  G4double py = direction.y() * momentum;
  G4double pz = direction.z() * momentum;

  G4PrimaryVertex* vertex = new G4PrimaryVertex(0., 0., 0., 0.);
  G4PrimaryParticle* particle = new G4PrimaryParticle(neutron, px, py, pz);
  particle->SetMass(total_mass);

  vertex->SetPrimary(particle);
  anEvent->AddPrimaryVertex(vertex);
}
