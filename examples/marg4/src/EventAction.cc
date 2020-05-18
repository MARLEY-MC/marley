#include <iostream>

#include "EventAction.hh"

EventAction::EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event* /*anEvent*/)
{
  // Print the event number at the beginning of every hundredth event
  static unsigned long long event_count = 0;
  if ( event_count % 100 == 0 ) std::cout << "Beginning event #" << event_count
    << '\n';
  ++event_count;
}

void EventAction::EndOfEventAction(const G4Event*)
{
}
