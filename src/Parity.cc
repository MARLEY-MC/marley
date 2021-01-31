/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
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

void marley::Parity::from_char( const char c ) {
  if ( c == '+' ) is_positive_ = true;
  else if ( c == '-' ) is_positive_ = false;
  else throw marley::Error( std::string("Invalid parity character \"") + c
    + "\" encountered in marley::Parity::from_char()" );
}

std::istream& operator>>(std::istream& in, marley::Parity& p) {
  // Read in the next non-whitespace character
  char c;
  in >> std::ws >> c;

  // If we suceeded in reading one, then convert it to a parity value
  if ( in ) p.from_char( c );

  // Either way, return the input stream
  return in;
}
