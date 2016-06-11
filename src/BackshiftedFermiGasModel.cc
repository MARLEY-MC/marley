#include <cmath>

#include "marley_utils.hh"
#include "MassTable.hh"
#include "BackshiftedFermiGasModel.hh"

marley::BackshiftedFermiGasModel::BackshiftedFermiGasModel(int Z, int A)
  : Z_(Z), A_(A)
{
  int N = A_ - Z_;
  double A_third = std::pow(A_, 1.0/3.0);

  // Parameters from global level density fit performed in A. J. Koning,
  // et al., Nucl. Phys. A810 (2008) pp. 13-76.
  static constexpr double alpha = 0.0722396;
  static constexpr double beta = 0.195267;
  static constexpr double gamma_1 = 0.410289;
  static constexpr double delta_global = 0.173015;

  // Asymptotic level density parameter
  a_tilde_ = alpha*A_ + beta*std::pow(A_third, 2);

  // Damping parameter
  gamma_ = gamma_1 / A_third;

  // Shell correction energy
  delta_W_ = marley::MassTable::get_mass_excess(Z_, A_)
    - marley::MassTable::liquid_drop_model_mass_excess(Z_, A_);

  // Energy shift
  Delta_BFM_ = delta_global;
  bool z_odd = Z % 2;
  bool n_odd = N % 2;
  // There will be no change to the energy shift if the nucleus is odd-even
  if (z_odd && n_odd) Delta_BFM_ += -12/std::sqrt(A_);
  else if (!z_odd && !n_odd) Delta_BFM_ += 12/std::sqrt(A_);

  // Spin cut-off parameter
  sigma_d_global_ = 0.83*std::pow(A_, 0.26);
  Sn_ = marley::MassTable::get_fragment_separation_energy(Z_, A_,
    marley_utils::NEUTRON);
}

// rho(Ex, J, Pi) assuming equipartition of parity (the parameter Pi is unused)
double marley::BackshiftedFermiGasModel::level_density(double Ex, int two_J,
  marley::Parity /*Pi*/)
{
  return 0.5 * level_density(Ex, two_J);
}

double marley::BackshiftedFermiGasModel::level_density(double Ex, int two_J) {
  double rho = level_density(Ex);
  // Spin-cutoff parameter sigma_ is updated by previous call to
  // this->level_density(Ex)
  double two_sigma2 = 2 * std::pow(sigma_, 2);
  return ((two_J + 1) / two_sigma2) * std::exp(-0.25 * std::pow(two_J + 1, 2)
    / two_sigma2) * rho;
}

// Returns the total level density rho(Ex) and store the spin cut-off parameter in sigma
double marley::BackshiftedFermiGasModel::level_density(double Ex) {

  // Effective excitation energy
  double U = Ex - Delta_BFM_;

  // Level density parameter
  double a;
  if (U <= 0) { // Equivalently, Ex <= Delta_BFM_
    // Use first-order Taylor expansion for small energies
    a = a_tilde_ * (1 + gamma_ * delta_W_);
  }
  else {
    a = a_tilde_ * (1 + (delta_W_ / U) * (1 - std::exp(-gamma_ * U)));
  }

  // TODO: Replace Ed here with database of local fits taken from RIPL-3 or TALYS
  const double Ed = 0.;

  // Spin cut-off parameter
  double sigma_F2 = 0.;

  // To avoid numerical problems, we will always use the discrete spin cutoff parameter for U <= Ed.
  // The TALYS manual suggests using this for Ex <= Ed, but this isn't a huge change. The actual TALYS
  // code may make this same choice.
  // TODO: Look into this more.
  if (U <= Ed) sigma_ = sigma_d_global_;
  else {
    sigma_F2 = 0.01389 * std::pow(A_, 5.0/3.0)
      * std::sqrt(a * U) / a_tilde_;
    if (Ex >= Sn_) sigma_ = std::sqrt(sigma_F2);
    else {
      // Ed < Ex < Sn_
      double sigma_d2 = std::pow(sigma_d_global_, 2);
      sigma_ = std::sqrt(sigma_d2 + (Ex - Ed)
        * (sigma_F2 - sigma_d2) / (Sn_ - Ed));
    }
  }

  // For very small excitation energies, take the limit of the total
  // level density as U -> 0 to prevent numerical issues.
  if (U <= 0) {
    constexpr double exp1 = std::exp(1);
    return exp1 * a / (12 * sigma_);
  }

  double aU = a * U;
  double sqrt_aU = std::sqrt(aU);
  return std::pow(12 * sigma_ * (std::sqrt(2 * sqrt_aU)*U*std::exp(-2 * sqrt_aU)
    + std::exp(-aU - 1)/a), -1);
}
