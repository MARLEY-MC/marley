#pragma once
#include "marley/Parity.hh"

namespace marley {

  /// @brief Abstract base class for models of nuclear level densities
  class LevelDensityModel {

    public:

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
