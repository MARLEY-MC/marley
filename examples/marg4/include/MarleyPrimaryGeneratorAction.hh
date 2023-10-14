/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
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

// Standard library includes
#include <string>

// Geant4 includes
#include "G4VUserPrimaryGeneratorAction.hh"

// MARLEY includes
#include "marley/Generator.hh"

class G4Event;

class MarleyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    MarleyPrimaryGeneratorAction(const std::string& config_file_name);

    virtual void GeneratePrimaries(G4Event*) override;

  protected:
    // MARLEY event generator object
    marley::Generator marley_generator_;
};
