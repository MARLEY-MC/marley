/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/COPYING or
// visit http://opensource.org/licenses/GPL-3.0

#include "marley/Error.hh"
#include "marley/Target.hh"
#include "marley/marley_utils.hh"

marley::Target::Target( int pdg ) {
  this->initialize_single_nuclide( pdg );
}

marley::Target::Target( int Z, int A ) {
  // Get the nuclear PDG code corresponding to these Z and A values
  int pdg = marley_utils::get_nucleus_pid( Z, A );

  // Delegate the rest of the work
  this->initialize_single_nuclide( pdg );
}

void marley::Target::initialize_single_nuclide( int pdg ) {
  // Create trivial vectors of a single nuclide
  // whose atom fraction is one
  std::vector<double> atom_fracs = { 1. };
  std::vector<marley::TargetAtom> atoms;
  atoms.emplace_back( pdg );

  // Delegate the rest to the initialization method for the general case
  this->initialize( atoms, atom_fracs );
}

marley::Target::Target( const std::vector<marley::TargetAtom>& nuclides,
  const std::vector<double>& atom_fracs )
{
  this->initialize( nuclides, atom_fracs );
}

void marley::Target::initialize( const std::vector<marley::TargetAtom>& nuclides,
  const std::vector<double>& atom_fracs )
{
  // Check that we have the same number of target atoms and atom fractions
  if ( nuclides.size() != atom_fracs.size() ) {
    throw marley::Error("Different numbers of target atoms and atom fractions"
      " encountered in the constructor of marley::Target");
  }

  // Check that we have at least one target nuclide. Note that the vector's
  // size is unsigned and therefore can't go negative.
  if ( nuclides.size() == 0u ) throw marley::Error( "No target atoms"
    " specified when constructing a marley::Target object" );

  // Make sure the atom fractions sum to unity by explicitly renormalizing them
  double sum = 0.;
  for ( size_t j = 0u; j < atom_fracs.size(); ++j ) {
    double af = atom_fracs.at( j );
    const auto& nuc = nuclides.at( j );
    if ( af < 0. ) throw marley::Error("Invalid atom fraction "
      + std::to_string(af) + " encountered for a "
      + nuc.to_string() + " target atom");
    else if ( af == 0. ) MARLEY_LOG_WARNING() << "Atom fraction of zero"
      << " encountered for a " << nuc << " target atom";
    sum += af;
  }

  // If the sum of the fractions is non-positive, we'll have trouble.
  if ( sum <= 0. ) throw marley::Error( "Sum of atom fractions = "
    + std::to_string(sum) + " in the constructor of marley::Target" );

  std::vector<double> renorm_fracs;
  for ( const auto& af : atom_fracs ) renorm_fracs.push_back( af / sum );

  // We're ready. Initialize the target atom map.
  for ( size_t j = 0u; j < nuclides.size(); ++j ) {
    const auto& nuc = nuclides.at( j );
    const auto& frac = renorm_fracs.at( j );
    // If we have duplicate nuclides, throw an error (likely to be a mistake)
    if ( atom_fractions_.count(nuc) ) throw marley::Error("Duplicate atom"
      " fractions specified for the target atom " + nuc.to_string() + " in"
      " the constructor of marley::Target");
    atom_fractions_[ nuc ] = frac;
  }
}

double marley::Target::atom_fraction( const marley::TargetAtom& atom ) const {
  // If we have an atom fraction stored in the map for this nuclide,
  // then just return it
  auto it = atom_fractions_.find( atom );
  if ( it != atom_fractions_.end() ) return it->second;
  // If we couldn't find it, then just return zero
  return 0.;
}

bool marley::Target::contains( const marley::TargetAtom& atom ) const {
  return ( this->atom_fraction(atom) != 0. );
}

void marley::Target::print( std::ostream& out ) const {
  size_t num_atom_types = atom_fractions_.size();
  size_t count = 1;
  for ( const auto& pair : atom_fractions_ ) {
    const auto& nuc = pair.first;
    const auto& frac = pair.second;
    out << nuc << " = " << frac;
    // Add a newline to all but the last entry
    if ( count < num_atom_types ) out << '\n';
    ++count;
  }
}
