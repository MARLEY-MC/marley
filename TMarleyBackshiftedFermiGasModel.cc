#include <cmath>

#include "marley_utils.hh"
#include "TMarleyMassTable.hh"
#include "TMarleyBackshiftedFermiGasModel.hh"

double TMarleyBackshiftedFermiGasModel::level_density(int Z, int A, double Ex) {
  double dummy;
  return level_density(Z, A, Ex, dummy);
}

double TMarleyBackshiftedFermiGasModel::level_density(int Z, int A, double Ex, int two_J) {
  double sigma;
  double rho = level_density(Z, A, Ex, sigma);
  double two_sigma2 = 2 * std::pow(sigma, 2);
  return 0.5 * ((two_J + 1) / two_sigma2) * std::exp(-0.25 * std::pow(two_J + 1, 2)
    / two_sigma2) * rho;
}

// Returns the total level density rho(Ex) and store the spin cut-off parameter in sigma
double TMarleyBackshiftedFermiGasModel::level_density(int Z, int A, double Ex, double& sigma) {
  int N = A - Z;
  double A_to_the_one_third = std::pow(A, 1.0/3.0);

  // Asymptotic level density parameter
  double a_tilde = alpha*A + beta*std::pow(A_to_the_one_third, 2);
  // Damping parameter
  double gamma = gamma_1 / A_to_the_one_third;
  // Shell correction energy
  double delta_W = TMarleyMassTable::get_mass_excess(Z, A)
    - TMarleyMassTable::liquid_drop_model_mass_excess(Z, A);
  // Energy shift
  double Delta_BFM = delta_global;
  bool z_odd = Z % 2;
  bool n_odd = N % 2;
  // There will be no change to the energy shift if the nucleus is odd-even
  if (z_odd && n_odd) Delta_BFM += -12/std::sqrt(A);
  else if (!z_odd && !n_odd) Delta_BFM += 12/std::sqrt(A);

  // Effective excitation energy
  double U = Ex - Delta_BFM;

  // Level density parameter
  double a;
  if (U <= 0) { // Equivalently, Ex <= Delta_BFM
    // Use first-order Taylor expansion for small energies
    a = a_tilde * (1 + gamma * delta_W);
  }
  else {
    a = a_tilde * (1 + (delta_W / U) * (1 - std::exp(-gamma * U)));
  }

  // Spin cut-off parameter
  double sigma_d_global = 0.83*std::pow(A, 0.26);
  double Sn = TMarleyMassTable::get_fragment_separation_energy(Z, A,
    marley_utils::NEUTRON);
  double sigma_F2;
  // TODO: Replace Ed here with database of local fits taken from RIPL-3 or TALYS
  double Ed = 0;

  // To avoid numerical problems, we will always use the discrete spin cutoff parameter for U <= Ed.
  // The TALYS manual suggests using this for Ex <= Ed, but this isn't a huge change. The actual TALYS
  // code may make this same choice.
  // TODO: Look into this more.
  if (U <= Ed) sigma = sigma_d_global;
  else {
    sigma_F2 = 0.01389 * std::pow(A_to_the_one_third, 5)
      * std::sqrt(a * U) / a_tilde;
    if (Ex >= Sn) sigma = std::sqrt(sigma_F2);
    else {
      // Ed < Ex < Sn
      double sigma_d2 = std::pow(sigma_d_global, 2);
      sigma = std::sqrt(sigma_d2 + (Ex - Ed)
        * (sigma_F2 - sigma_d2) / (Sn - Ed));
    }
  }

  // For very small excitation energies, take the limit of the total
  // level density as U -> 0 to prevent numerical issues.
  if (U <= 0) {
    return std::exp(1) * a / (12 * sigma);
  }

  double aU = a * U;
  double sqrt_aU = std::sqrt(aU);
  return std::pow(12 * sigma * (std::sqrt(2 * sqrt_aU)*U*std::exp(-2 * sqrt_aU)
    + std::exp(-aU - 1)/a), -1);
}
