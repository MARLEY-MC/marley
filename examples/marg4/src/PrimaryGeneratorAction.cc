// standard library includes
#include <iostream>

// Geant4 includes
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"

// MARLEY includes
#include "marley/Event.hh"
#include "marley/JSON.hh"
#include "marley/Logger.hh"
#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#else
  #include "marley/JSONConfig.hh"
#endif

#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(
  const std::string& config_file_name) : G4VUserPrimaryGeneratorAction()
{
  // Create a new generator object using the settings from the configuration
  // file
  #ifdef USE_ROOT
  marley::RootJSONConfig config( config_file_name );
  #else
  marley::JSONConfig config( config_file_name );
  #endif

  marley_generator_= config.create_generator();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4PrimaryVertex* vertex = new G4PrimaryVertex(0., 0., 0., 0.); // x,y,z,t0

  marley::Event ev = marley_generator_.create_event();

  //std::cout << ev << '\n';

  for ( const auto& fp : ev.get_final_particles() ) {

    G4PrimaryParticle* particle = new G4PrimaryParticle( fp->pdg_code(),
      fp->px(), fp->py(), fp->pz(), fp->total_energy() );

    particle->SetCharge( fp->charge() );

    vertex->SetPrimary(particle);
  }

  anEvent->AddPrimaryVertex(vertex);
}
