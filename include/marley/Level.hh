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

#pragma once
#include <random>

#include "marley/Gamma.hh"
#include "marley/IteratorToPointerMember.hh"
#include "marley/Parity.hh"

namespace marley {

  class Generator;

  /// @brief A discrete nuclear energy level
  class Level {
    public:

      /// @param E excitation energy of this level (MeV)
      /// @param twoJ two times the level spin
      /// @param pi level parity
      Level(double E, int twoJ, marley::Parity pi);

      /// @brief Retrieve a const reference to the vector of gamma rays owned
      /// by this level
      inline const std::vector<marley::Gamma>& gammas() const;

      /// Retrieve a reference to the vector of gamma rays owned by this level
      inline std::vector<marley::Gamma>& gammas();

      /// Get the excitation energy of this level (MeV)
      inline double energy() const;
      /// Set the excitation energy of this level (MeV)
      inline void set_energy(double E);

      /// Get two times the level spin
      inline int twoJ() const;
      /// Set two times the level spin
      inline void set_twoJ(int twoJ);

      /// Get the level parity
      inline marley::Parity parity() const;
      /// Set the level parity
      inline void set_parity(marley::Parity pi);

      /// @return true if the level owns at least one Gamma object
      inline bool has_gammas() const;

      /// @brief Add a new gamma-ray transition to this level
      /// @return a reference to the newly added gamma
      marley::Gamma& add_gamma(const marley::Gamma& gamma);

      /// @brief Add a new gamma-ray transition to this level
      /// @param energy Gamma-ray energy (MeV)
      /// @param branching_ratio Branching ratio for the new
      /// gamma-ray transition
      /// @param end_lev Pointer to the level that absorbs this gamma
      /// @return a reference to the newly added gamma
      marley::Gamma& add_gamma(double energy, double branching_ratio,
        marley::Level* end_lev = nullptr);

      /// Remove all gamma ray information from this level
      void clear_gammas();

      /// @brief Choose a gamma owned by this level randomly based on the
      /// relative intensities of all of the gammas.
      /// @return a pointer to the selected Gamma object, or nullptr
      /// if the level doesn't have any gammas
      const marley::Gamma* sample_gamma(marley::Generator& gen);

      /// Returns the level spin-parity as a string
      std::string spin_parity_string() const;

      /// @brief Convert an iterator that points to a marley::Level* (or a
      /// smart pointer to a marley::Level) into an iterator that points to the
      /// Level's energy_ member variable.
      /// @details This function is used by the DecayScheme class to keep the
      /// vector of std::unique_ptr<marley::Level> objects that it owns sorted
      /// in order of increasing excitation energy.
      template<typename It> inline static marley::IteratorToPointerMember<It,
        double> make_energy_iterator(It it);

    private:

      double energy_; ///< excitation energy (MeV)

      /// @brief two times the level spin
      /// @note the factor of two allows us to represent half-integer spins as
      /// integers
      int twoJ_;

      marley::Parity parity_; ///< level parity

      /// @brief gamma-ray transitions owned by this level
      std::vector<marley::Gamma> gammas_;

      /// @brief discrete distribution object used to sample gamma-ray
      /// de-excitations
      std::discrete_distribution<size_t> gamma_dist_;

      /// @brief helper function that updates gamma-ray distribution when Gamma
      /// objects are added or removed from the level
      void update_gamma_distribution();
  };

  // Inline function definitions
  inline double Level::energy() const { return energy_; }
  inline void Level::set_energy(double E) { energy_ = E; }

  inline int Level::twoJ() const { return twoJ_; }
  inline void Level::set_twoJ(int twoJ) { twoJ_ = twoJ; }

  inline marley::Parity Level::parity() const { return parity_; }
  inline void Level::set_parity(marley::Parity pi) { parity_ = pi; }

  inline const std::vector<marley::Gamma>& Level::gammas() const
    { return gammas_; }
  inline std::vector<marley::Gamma>& Level::gammas() { return gammas_; }

  inline bool Level::has_gammas() const { return !gammas_.empty(); }

  template<typename It> inline marley::IteratorToPointerMember<It,
    double> Level::make_energy_iterator(It it)
  {
    return marley::IteratorToPointerMember<It, double>(it,
      &marley::Level::energy_);
  }
}
