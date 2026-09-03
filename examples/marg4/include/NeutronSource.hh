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

#pragma once

// Geant4 includes
#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;

class NeutronSource : public G4VUserPrimaryGeneratorAction
{
  public:
    NeutronSource();
    ~NeutronSource();

    virtual void GeneratePrimaries(G4Event*) override;

  private:
    constexpr static double energy_min_ = 0.5 * CLHEP::MeV;
    constexpr static double energy_max_ = 500 * CLHEP::MeV;
};
