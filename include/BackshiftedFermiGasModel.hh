#pragma once
#include "LevelDensityModel.hh"

namespace marley {

  class BackshiftedFermiGasModel : public LevelDensityModel {

    public:

      // constructor
      BackshiftedFermiGasModel(int Z, int A);
      // TODO: add other constructors to allow setting the parameters based on
      // local fits

      // rho(Ex)
      virtual double level_density(double Ex) override;

      // rho(Ex, J)
      virtual double level_density(double Ex, int two_J) override;

      // rho(Ex, J, Pi) with the assumption of parity equipartition
      virtual double level_density(double Ex, int two_J, marley::Parity Pi)
        override;

    protected:

      int Z_; // atomic number for this nuclide
      int A_; // mass number for this nuclide

      double sigma_; // spin cut-off parameter

      double a_tilde_; // asymptotic level density parameter
      double gamma_; // damping parameter
      double delta_W_; // shell correction energy
      double Delta_BFM_; // energy shift
      double sigma_d_global_; // global fit for discrete-region spin cut-off parameter
      double Sn_; // neutron separation energy
  };
}
