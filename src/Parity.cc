/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

// Standard library includes
#include <string>

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Parity.hh"

marley::Parity::Parity( int i ) {
  if (i == 1) is_positive_ = true;
  else if (i == -1) is_positive_ = false;
  else throw marley::Error( "Invalid parity " + std::to_string(i)
    + " passed to constructor of marley::Parity" );
}

/// @brief Assigns a parity value using an integer
/// @details If an int value other than +1 or -1 is used, a marley::Error
/// will be thrown.
marley::Parity& marley::Parity::operator=( const int& i ) {
  // Do the assignment
  if (i == 1) is_positive_ = true;
  else if (i == -1) is_positive_ = false;
  else throw marley::Error( "Invalid parity " + std::to_string(i)
    + " assigned to variable of type marley::Parity" );

  // Return the existing object
  return *this;
}

std::istream& operator>>(std::istream& in, marley::Parity& p) {
  char c;
  in >> c;

  if (c == '+') p = marley::Parity(true);
  else if (c == '-') p = marley::Parity(false);
  else throw marley::Error( std::string("Invalid parity \"") + c + "\" assigned"
    " via the >> operator to a variable of type marley::Parity" );

  return in;
}
