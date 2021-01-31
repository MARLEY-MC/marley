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

// Simple detector construction class
#pragma once

#include "G4VUserDetectorConstruction.hh"

class G4PhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();

    virtual G4VPhysicalVolume* Construct() override;

    const G4VPhysicalVolume* GetWorld() const;

  protected:

    G4VPhysicalVolume* world_ = nullptr;
};
