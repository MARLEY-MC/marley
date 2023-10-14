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



#include "marley/Error.hh"
#include "marley/TargetAtom.hh"
#include "marley/marley_utils.hh"

marley::TargetAtom::TargetAtom( int pdg ) : pdg_( pdg ) {
  this->check_pdg_validity();
}

marley::TargetAtom::TargetAtom( int Z, int A ) {
  pdg_ = marley_utils::get_nucleus_pid( Z, A );
  this->check_pdg_validity();
}

int marley::TargetAtom::Z() const {
  return marley_utils::get_particle_Z( pdg_ );
}

int marley::TargetAtom::A() const {
  return marley_utils::get_particle_A( pdg_ );
}

std::string marley::TargetAtom::to_string() const {
  std::string result = std::to_string( this->A() );
  result += marley_utils::element_symbols.at( this->Z() );
  return result;
}

void marley::TargetAtom::check_pdg_validity() const {
  int A = marley_utils::get_particle_A( pdg_ );
  if ( A <= 0 ) throw marley::Error("Invalid nuclear PDG code "
    + std::to_string(pdg_) + " encountered upon construction of"
    + " a marley::TargetAtom object.");
}
