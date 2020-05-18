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
