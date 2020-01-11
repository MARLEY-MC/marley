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

#pragma once
#include "marley/Parity.hh"

namespace marley {

  /// @brief Abstract base class for models of nuclear level densities
  class LevelDensityModel {

    public:

      virtual ~LevelDensityModel() = default;

      /// Nuclear level density @f$ \rho(E_x) @f$ including all spins and
      /// parities.
      /// @param Ex Excitation energy in MeV
      /// @return %Level density in MeV<sup> -1</sup>
      virtual double level_density(double Ex) = 0;

      /// %Level density @f$ \rho(E_x, J) @f$ for a specific nuclear spin.
      /// @param Ex Excitation energy in MeV
      /// @param two_J Two times the nuclear spin
      /// @return %Level density in MeV<sup> -1</sup>
      virtual double level_density(double Ex, int two_J) = 0;

      /// %Level density @f$ \rho(E_x, J, \Pi) @f$ for a specific nuclear spin
      /// and parity.
      /// @param Ex Excitation energy in MeV
      /// @param two_J Two times the nuclear spin
      /// @param Pi The nuclear parity
      /// @return %Level density in MeV<sup> -1</sup>
      virtual double level_density(double Ex, int two_J, marley::Parity Pi) = 0;
  };
}
