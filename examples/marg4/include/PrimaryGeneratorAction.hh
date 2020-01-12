#pragma once

// Geant4 includes
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

// MARLEY includes
#include "marley/Generator.hh"

class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(const std::string& config_file_name);

    virtual void GeneratePrimaries(G4Event*) override;

  protected:
    // MARLEY event generator object
    marley::Generator marley_generator_;
};
