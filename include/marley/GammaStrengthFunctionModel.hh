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

    protected:

      /// @brief Check that l > 0 and throw a marley::Error if it is not.
      static void check_multipolarity(int l);

      int Z_; ///< Atomic number
      int A_; ///< Mass number
  };

}
