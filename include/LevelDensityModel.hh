#pragma once
#include "Parity.hh"

namespace marley {

  // Abstract base class for models of nuclear level densities
  class LevelDensityModel {

    public:

      // rho(Ex)
      virtual double level_density(double Ex) = 0;

      // rho(Ex, J)
      virtual double level_density(double Ex, int two_J) = 0;

      // rho(Ex, J, Pi)
      virtual double level_density(double Ex, int two_J, marley::Parity Pi) = 0;
  };
}
