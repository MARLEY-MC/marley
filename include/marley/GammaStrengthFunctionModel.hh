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

namespace marley {

  class Level;
  class Parity;

  /// @brief Abstract base class for models of gamma-ray strength functions
  /// @details Classes derived from GammaStrengthFunctionModel may be used to
  /// simulate gamma-ray emission in situations where nuclear level data are
  /// unavailable. In particular, the HauserFeshbachDecay class uses instances
  /// of GammaStrengthFunctionModel to model the competition between gamma-ray
  /// emission and particle evaporation for highly-excited nuclear states.
  class GammaStrengthFunctionModel {

    public:

      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      GammaStrengthFunctionModel(int Z, int A);

      virtual ~GammaStrengthFunctionModel() = default;

      /// @brief Electromagnetic transitions in nuclei may be classified by
      /// their <a href="http://tinyurl.com/hscpepq">multipolarity</a>
      /// (electric vs. magnetic multipole radiation)
      enum class TransitionType { electric, magnetic, unphysical };

      /// @brief Determines whether a given electromagnetic transition between
      /// two nuclear states corresponds to electric or magnetic multipole
      /// radiation.
      /// @param twoJi Two times the initial nuclear spin
      /// @param Pi Initial parity
      /// @param twoJf Two times the final nuclear spin
      /// @param Pf Final parity
      /// @param[out] l Multipolarity of the transition (@f$\ell@f$ = 1 &hArr;
      /// dipole, @f$\ell@f$ = 2 &hArr; quadrupole, &hellip;)
      /// @note This function returns TransitionType::unphysical if the
      /// requested transition is impossible (e.g., nuclear spin changes by
      /// half, initial and final spins are both zero)
      static TransitionType determine_transition_type(int twoJi,
        marley::Parity Pi, int twoJf, marley::Parity Pf, int& l);

      /// @brief Determines whether a given electromagnetic transition between
      /// two nuclear states corresponds to electric or magnetic multipole
      /// radiation.
      /// @param twoJi Two times the initial nuclear spin
      /// @param Pi Initial parity
      /// @param[in] level_f Reference to the final nuclear Level
      /// @param[out] l Multipolarity of the transition (@f$\ell@f$ = 1 &hArr;
      /// dipole, @f$\ell@f$ = 2 &hArr; quadrupole, &hellip;)
      /// @note This function returns TransitionType::unphysical if the
      /// requested transition is impossible (e.g., nuclear spin changes by
      /// half, initial and final spins are both zero)
      static TransitionType determine_transition_type(int twoJi,
        marley::Parity Pi, marley::Level& level_f, int& l);

      /// @brief Determines whether a given electromagnetic transition between
      /// two nuclear states corresponds to electric or magnetic multipole
      /// radiation.
      /// @param[in] level_i Reference to the initial nuclear Level
      /// @param[in] level_f Reference to the final nuclear Level
      /// @param[out] l Multipolarity of the transition (@f$\ell@f$ = 1 &hArr;
      /// dipole, @f$\ell@f$ = 2 &hArr; quadrupole, &hellip;)
      /// @note This function returns TransitionType::unphysical if the
      /// requested transition is impossible (e.g., nuclear spin changes by
      /// half, initial and final spins are both zero)
      static TransitionType determine_transition_type(marley::Level& level_i,
        marley::Level& level_f, int& l);

      /// @brief Returns the gamma-ray strength function
      /// (MeV<sup> &ndash;2@f$\ell@f$&ndash;1</sup>) for the requested gamma
      /// energy and multipolarity
      /// @param type Electric or magnetic transition
      /// @param l Multipolarity of the transition
      /// @param e_gamma Gamma-ray energy (MeV)
      virtual double strength_function(TransitionType type, int l,
        double e_gamma) = 0;

      /// @brief Returns the gamma-ray transmission coefficient (dimensionless)
      /// for the requested gamma energy and multipolarity
      /// @details The gamma-ray transmission coefficient and strength function
      /// are related via @f$\text{T}_{\text{X}\ell}(\text{E}_\gamma)
      /// = 2\pi f_{\text{X}\ell}(\text{E}_\gamma)\text{E}_\gamma^{(2\ell
      /// + 1)},@f$ where X is the type of transition (electric or magnetic),
      /// @f$\ell@f$ is the multipolarity, @f$\text{T}_{\text{X}\ell}@f$ is
      /// the transmission coefficient, @f$f_{\text{X}\ell}@f$ is the strength
      /// function, and @f$\text{E}_\gamma@f$ is the gamma-ray energy.
      /// @param type Electric or magnetic transition
      /// @param l Multipolarity of the transition
      /// @param e_gamma Gamma-ray energy (MeV)
      virtual double transmission_coefficient(TransitionType type, int l,
        double e_gamma) = 0;

      /// @brief Returns the gamma-ray transmission coefficient (dimensionless)
      /// for a transition from a given initial state to a final discrete Level
      /// @param Exi Initial nuclear excitation energy
      /// @param twoJ Two times the initial nuclear spin
      /// @param Pi Initial nuclear parity
      /// @param[in] level_f Reference to the final nuclear Level
      double transmission_coefficient(double Exi, int twoJ,
        marley::Parity Pi, marley::Level& level_f);

    protected:

      /// @brief Check that l > 0 and throw a marley::Error if it is not.
      static void check_multipolarity(int l);

      int Z_; ///< Atomic number
      int A_; ///< Mass number
  };

}
