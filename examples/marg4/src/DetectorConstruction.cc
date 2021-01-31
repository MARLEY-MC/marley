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

// Geant4 includes
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Orb.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"

// marg4 includes
#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction()
{
}

// Make a ball of liquid argon with a radius of 10 m
G4VPhysicalVolume* DetectorConstruction::Construct() {

  G4Orb* orb = new G4Orb("orb", 10.*meter);

  G4Material* liquid_argon = new G4Material("liquidArgon", 18, 39.95*g/mole,
    1.390*g/cm3, kStateLiquid, 87.*kelvin);

  G4LogicalVolume* logical_world = new G4LogicalVolume( orb,
    liquid_argon, liquid_argon->GetName() );

  world_ = new G4PVPlacement(0, // no rotation
    G4ThreeVector(), // at (0, 0, 0)
    logical_world, // its logical volume
    "WorldOrb", // its name
    nullptr, // its mother logical volume
    false, // no boolean operation
    0); // copy number

  return world_;
}

const G4VPhysicalVolume* DetectorConstruction::GetWorld() const {
  return world_;
}
