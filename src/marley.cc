/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
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

#include "marley/CommandHandler.hh"

int main( int argc, char* argv[] ) {

  // Initalize the CommandHandler which will manage calling
  // the functions requested on the command line
  marley::CommandHandler ch( argc, argv );

  bool ok = ch.execute();
  if ( !ok ) return 1;
  return 0;
}
