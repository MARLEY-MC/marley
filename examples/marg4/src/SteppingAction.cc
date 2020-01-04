#include <iostream>
#include <string>

#include "SteppingAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"

SteppingAction::SteppingAction() : G4UserSteppingAction()
{
}

void SteppingAction::UserSteppingAction(const G4Step* /*aStep*/)
{
}
