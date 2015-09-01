#pragma once
// TODO: change this from an all-static class so that one can create objects
// with local fit parameters for a specific nuclide
class TMarleyBackshiftedFermiGasModel {
  public:
    static double level_density(int Z, int A, double Ex); // rho(Ex)
    // rho(Ex, J, Pi) with the assumption of parity equipartition
    static double level_density(int Z, int A, double Ex, int two_J);
  private:
    // Version of the level density function used internally
    static double level_density(int Z, int A, double Ex, double& sigma);

    // Parameters from global level density fit performed in A. J. Koning,
    // et al., Nucl. Phys. A810 (2008) pp. 13-76.
    static constexpr double alpha = 0.0722396;
    static constexpr double beta = 0.195267;
    static constexpr double gamma_1 = 0.410289;
    static constexpr double delta_global = 0.173015;
};
