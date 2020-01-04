#pragma once

#include "G4UserSteppingAction.hh"

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();

    virtual void UserSteppingAction(const G4Step*) override;
};
