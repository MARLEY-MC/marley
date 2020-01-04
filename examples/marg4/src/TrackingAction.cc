#include <iostream>

#include "TrackingAction.hh"

#include "G4Neutron.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

TrackingAction::TrackingAction() : G4UserTrackingAction()
{
}

void TrackingAction::PreUserTrackingAction(const G4Track* /*aTrack*/)
{
}

void TrackingAction::PostUserTrackingAction(const G4Track* /*aTrack*/)
{
}
