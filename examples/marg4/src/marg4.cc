#include <memory>

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace {
  // Retrieves the desired number of events from the first command
  // line argument. Based on https://stackoverflow.com/a/2797823/4081973
  // Returns true if everything went well, or false if there was a problem
  bool get_num_events( const std::string& arg, int& x ) {
    try {
      size_t pos;
      x = std::stoi( arg, &pos );
      if ( pos < arg.size() ) {
        std::cerr << "Trailing characters after number of events: " << arg
          << '\n';
        return false;
      }
      if ( x < 0 ) {
        std::cerr << "Negative number of events: " << arg << '\n';
        return false;
      }
    }
    catch (std::invalid_argument const &ex) {
      std::cerr << "Invalid number: " << arg << '\n';
      return false;
    }
    catch (std::out_of_range const &ex) {
      std::cerr << "Number out of range: " << arg << '\n';
      return false;
    }

    return true;
  }
}


int main(int argc, char* argv[]) {

  if ( argc <= 2 ) {
    std::cout << "Usage: example NUM_EVENTS MARLEY_CONFIG_FILE\n";
    return 1;
  }

  int num_events = 0;
  bool num_ok = get_num_events( argv[1], num_events );
  if ( !num_ok ) return 2;

  // Retrieve the configuration file name from the command-line argument
  std::string config_file_name( argv[2] );

  // Initialize the Geant4 run manager
  std::unique_ptr<G4RunManager> rm( new G4RunManager );

  // ** Set mandatory initialization classes **
  // Define the geometry for the simulation
  DetectorConstruction* det = new DetectorConstruction();
  rm->SetUserInitialization(det);

  // Set up the QGSP_BIC_HP physics list
  // Note: High-precision tracking of neutrons is set up by standard
  // physics lists with the "_HP" suffix. Transport of neutrons generated
  // by MARLEY will not be trustworthy unless the high-precision tracking
  // is enabled.
  G4PhysListFactory factory;
  G4VModularPhysicsList* refList = factory.GetReferencePhysList("QGSP_BIC_HP");
  rm->SetUserInitialization(refList);

  // ** Set user actions **

  // The primary generator action interfaces with MARLEY
  PrimaryGeneratorAction* pga = new PrimaryGeneratorAction(config_file_name);
  rm->SetUserAction(pga);

  TrackingAction* ta = new TrackingAction;
  rm->SetUserAction(ta);

  SteppingAction* sa = new SteppingAction();
  rm->SetUserAction(sa);

  EventAction* eva = new EventAction;
  rm->SetUserAction(eva);

  rm->Initialize();

  std::cout << "Simulating " << num_events << " events . . .\n";

  rm->BeamOn(num_events);

  return 0;
}
