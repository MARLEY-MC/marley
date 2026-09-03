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

// marg4 includes
#include "PhysicsList.hh"

// Geant4 includes
#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4CascadeInterface.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4Neutron.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4ProcessManager.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList()
{
  SetDefaultCutValue(10 * keV);
  SetVerboseLevel(1);
}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConstructParticle()
{
  // Register all particle families that the Bertini cascade and
  // G4Evaporation de-excitation may produce as secondaries.
  G4BaryonConstructor{}.ConstructParticle();
  G4BosonConstructor{}.ConstructParticle();
  G4LeptonConstructor{}.ConstructParticle();
  G4MesonConstructor{}.ConstructParticle();
  G4IonConstructor{}.ConstructParticle();
  G4ShortLivedConstructor{}.ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  // Bertini cascade (G4CascadeInterface) with PreCompound+G4Evaporation
  // de-excitation, covering the full energy range of the neutron source.
  G4CascadeInterface* bertini = new G4CascadeInterface();
  bertini->usePreCompoundDeexcitation();
  bertini->SetMinEnergy(0.);
  bertini->SetMaxEnergy(10. * GeV);

  // Attach a cross-section dataset and the cascade model to the process.
  G4NeutronInelasticProcess* neutron_inelastic = new G4NeutronInelasticProcess();
  neutron_inelastic->AddDataSet(new G4NeutronInelasticXS());
  neutron_inelastic->RegisterMe(bertini);

  // Register the process with the neutron's process manager.
  G4ProcessManager* pman = G4Neutron::Neutron()->GetProcessManager();
  pman->AddDiscreteProcess(neutron_inelastic);
}
