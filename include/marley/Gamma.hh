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
#include "marley/IteratorToMember.hh"

namespace marley {

  class Level;

  /// @brief A gamma-ray transition between two nuclear levels
  class Gamma {

    public:

      /// @param energy Energy of the emitted &gamma;-ray (MeV)
      /// @param rel_intensity Relative intensity of this transition
      /// @param start_lev Pointer to the Level that de-excites
      /// by emitting this &gamma;-ray
      Gamma(double energy = 0., double rel_intensity = 0.,
        marley::Level* start_lev = nullptr, marley::Level* end_lev = nullptr);

      /// @brief Get a pointer to the Level that emits this &gamma;-ray
      inline marley::Level* start_level() const;

      /// @brief Set the Level that emits this &gamma;-ray
      inline void set_start_level(marley::Level* start_lev);

      /// @brief Get a pointer to the Level that absorbs this &gamma;-ray
      inline marley::Level* end_level() const;

      /// @brief Set the Level that abosrbs this &gamma;-ray
      inline void set_end_level(marley::Level* end_lev);

      /// @brief Get the energy of the emitted &gamma;-ray (MeV)
      inline double energy() const;

      /// @brief Get the relative intensity for this transition
      inline double relative_intensity() const;

      /// @brief Convert an iterator that points to a Gamma object into an
      /// iterator that points to the Gamma's relative_intensity_ member
      /// variable.
      /// @details This function is used to load the std::discrete_distribution
      /// in the starting Level object with the intensities of the gammas that
      /// it owns without redundant storage.
      template<typename It> inline static marley::IteratorToMember<It, double>
        make_intensity_iterator(It it);

    protected:

      double energy_; ///< Energy of the emitted gamma (MeV)
      double relative_intensity_; ///< Relative intensity of the transition

      /// @brief Pointer to the Level that emits this &gamma;-ray
      marley::Level* start_level_;

      /// @brief Pointer to the Level that absorbs this &gamma;-ray
      marley::Level* end_level_;
  };

  // Inline function definitions
  inline marley::Level* Gamma::start_level() const { return start_level_; }
  inline void Gamma::set_start_level(marley::Level* start_lev)
    { start_level_ = start_lev; }

  inline marley::Level* Gamma::end_level() const { return end_level_; }
  inline void Gamma::set_end_level(marley::Level* end_lev)
    { end_level_ = end_lev; }

  inline double Gamma::energy() const { return energy_; }
  inline double Gamma::relative_intensity() const
    { return relative_intensity_; }

  template<typename It> inline marley::IteratorToMember<It, double>
    Gamma::make_intensity_iterator(It it)
  {
    return marley::IteratorToMember<It, double>(it,
      &marley::Gamma::relative_intensity_);
  }
}
