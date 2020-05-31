/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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

#pragma once

// Standard library includes
#include <ostream>
#include <map>
#include <vector>

// MARLEY includes
#include "TargetAtom.hh"

namespace marley {

  /// @brief Description of a macroscopic target for scattering reactions
  /// @details The target is composed of one or more nuclides and stores
  /// their relative abundances as atom fractions
  class Target {

    public:

      /// @brief Create a Target composed of a single nuclide
      /// @param pdg Particle Data Group code for the atomic species
      explicit Target( int pdg );

      /// @brief Create a Target composed of a single nuclide
      /// @param Z The proton number of the target atom
      /// @param A The mass number of the target atom
      explicit Target( int Z, int A );

      /// @brief Create a composite Target
      /// @details Atom fractions that do not sum to unity will automatically
      /// be renormalized to do so
      /// @param nuclides TargetAtom objects specifying each atomic species
      /// @param atom_fracs Atom fraction for each TargetAtom
      /// @note The nuclides and atom_fracs vectors must have the same length.
      /// If this is not the case, a marley::Error will be thrown
      Target( const std::vector<TargetAtom>& nuclides,
        const std::vector<double>& atom_fracs );

      /// @brief Returns the atom fraction for the requested nuclide
      double atom_fraction( const marley::TargetAtom& atom ) const;

      /// @brief Return true if the target contains the requested atom, or
      /// false otherwise
      /// @details Even if the target atom appears in the map, this function
      /// will return false if its atom fraction is exactly zero
      bool contains( const marley::TargetAtom& atom ) const;

      /// @brief Returns true if the target consists of a single kind of
      /// target atom, or false otherwise
      inline bool has_single_nuclide() const
        { return (atom_fractions_.size() == 1u); }

      /// @brief Print a textual representation of the Target to a std::ostream
      void print( std::ostream& out ) const;

      /// @brief Get read-only access to the map of atom fractions
      inline const std::map< marley::TargetAtom, double >& atom_fraction_map()
        const { return atom_fractions_; }

    protected:

      /// @brief Helper function for the constructors that use a single
      /// nuclide
      void initialize_single_nuclide( int pdg );

      /// @brief General helper function for the constructors
      void initialize( const std::vector<TargetAtom>& nuclides,
        const std::vector<double>& atom_fracs );

      /// @brief Map storing the atom fraction for each nuclide in the target
      /// material
      /// @details Keys are TargetAtom objects representing each
      /// nuclide, values are atom fractions. The Target class ensures
      /// that the atom fractions always sum to unity.
      std::map< marley::TargetAtom, double > atom_fractions_;
  };

}

inline std::ostream& operator<<(std::ostream& out, const marley::Target& t) {
  t.print( out );
  return out;
}
