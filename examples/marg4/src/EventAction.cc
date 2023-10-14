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
