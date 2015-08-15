#include <cmath>
#include <complex>
#include <iostream>

// TODO: fix this ugly include pattern that the
// cwfcomp people forced you to use
#include "complex_functions.hh"
#include "cwfcomp.hh"
#include "cwfcomp.cc"

#include "/home/sjg/Desktop/marley/TMarleyMassTable.hh"
#include "/home/sjg/Desktop/marley/marley_utils.hh"

// Constants
const double alpha = 7.2973525698e-3; // Fine structure constant
const double hbarc = 197.3269631; // MeV * fm
const double hbarc2 = std::pow(hbarc, 2); // MeV^2 * fm^2
// Constant to use when converting from mb to MeV^(-2)
const double mb = 1/3.89379338e5; // MeV^(-2) mb^(-1)
const double e2 = hbarc * alpha; // Elementary charge (in MeV*fm)
// Mass of a charged pion
const double mpiplus = 139.57018; // MeV
// Squared pion Compton wavelength in fm
const double lambda_piplus2 = std::pow(hbarc / mpiplus, 2);

// Fragment particle IDs
const int NEUTRON = 2112;
const int PROTON = 2212;
const int DEUTERON = 1000010020;
const int TRITON = 1000010030;
const int HELION = 1000020030;
const int ALPHA = 1000020040;

class SphericalOpticalModel {
  public:
    SphericalOpticalModel(int z, int a);
    std::complex<double> optical_model_potential(double r, double E,
      int fragment_pid, int two_j, int l, int two_s);
    double transmission_coefficient(double E, int fragment_pid, int two_j,
      int l, int two_s, double h);

    inline double get_fragment_reduced_mass(int fragment_pid) {
      return reduced_masses.at(fragment_pid);
    }

    // Non-derivative radial Schrödinger equation terms to use for computing transmission
    // coefficients via the Numerov method
    inline std::complex<double> a(double r, double E, int fragment_pid,
      int two_j, int l, int two_s)
    {
      return (-l*(l+1) / std::pow(r, 2)) +
        2 * reduced_masses.at(fragment_pid) * (E - optical_model_potential(r,
        E, fragment_pid, two_j, l, two_s)) / hbarc2;
    }

  private:
    // Woods-Saxon shape
    inline double f(double r, double R, double a) {
      return std::pow(1 + std::exp((r - R) / a), -1);
    }

    // Partial derivative with respect to r of the Woods-Saxon shape
    inline double dfdr(double r, double R, double a) {
      // In the limit as r -> +-infinity, this goes to zero.
      // We pick an upper limit for the exponent to avoid evaluating
      // the function explicitly when r gets too large (otherwise, C++
      // returns NaN because the function becomes indeterminate in double
      // precision (infinity/infinity or 0/0)
      double exponent = (r - R) / a;
      if (std::abs(exponent) > 100.) return 0;
      double temp = std::exp(exponent);
      return -temp / (a * std::pow(1 + temp, 2));
    }

    // Coulomb potential for a point particle with charge z*e interacting with a
    // uniformly charged sphere with radius R and charge Z*e
    inline double Vc(double r, double R, int Z, int z) {
      if (Z == 0 || z == 0) return 0;
      else if (r < R) return Z * z * e2 * (3 - std::pow(r / R, 2)) / (2 * R);
      else return Z * z * e2 / r;
    }

    // Nuclear atomic and mass numbers
    int Z, A;
    // Neutron parameters
    double v1n, v2n, v3n, v4n, w1n, w2n, d1n, d2n, d3n, vso1n, vso2n;
    double wso1n, wso2n, Efn, Rvn, avn, Rdn, adn, Rso_n, aso_n;
    // Proton parameters
    double v1p, v2p, v3p, v4p, w1p, w2p, d1p, d2p, d3p, vso1p, vso2p;
    double wso1p, wso2p, Efp, Vcbar_p, Rvp, avp, Rdp, adp, Rso_p, aso_p;
    // Radius for nuclear Coulomb potential
    double Rc;
    // Reduced mass lookup table
    std::unordered_map<int, double> reduced_masses;
};

std::complex<double> SphericalOpticalModel::optical_model_potential(double r, double E,
  int fragment_pid, int two_j, int l, int two_s)
{

  // Fragment atomic, mass, and neutron numbers
  int z = TMarleyMassTable::get_particle_Z(fragment_pid);
  int a = TMarleyMassTable::get_particle_A(fragment_pid);
  int n = a - z;

  // Geometrical parameters
  double Rv = 0;
  double av = 0;
  double Rd = 0;
  double ad = 0;
  double Rso = 0;
  double aso = 0;

  // Terms in the spherical optical model potential
  double Vv = 0;
  double Wv = 0;
  double Wd = 0;
  double Vso = 0;
  double Wso = 0;

  // Energy to use when computing folded potentials
  double E_eff = E / a;

  if (n > 0) {
    double Ediff_n = E_eff - Efn;
    double Ediff_n2 = std::pow(Ediff_n, 2);
    double Ediff_n3 = std::pow(Ediff_n, 3);

    Vv += n * v1n * (1 - v2n*Ediff_n + v3n*Ediff_n2 - v4n*Ediff_n3);
    Wv += n * w1n * Ediff_n2 / (Ediff_n2 + std::pow(w2n, 2));
    Wd += n * d1n * Ediff_n2 * std::exp(-d2n * Ediff_n)
      / (Ediff_n2 + std::pow(d3n, 2));

    Rv += n * Rvn;
    av += n * avn;
    Rd += n * Rdn;
    ad += n * adn;
    Rso += n * Rso_n;
    aso += n * aso_n;

    if (two_s != 0) {
      Vso += vso1n * std::exp(-vso2n * Ediff_n);
      Wso += wso1n * Ediff_n2 / (Ediff_n2 + std::pow(wso2n, 2));
    }
  }

  if (z > 0) {
    double Ediff_p = E_eff - Efp;
    double Ediff_p2 = std::pow(Ediff_p, 2);
    double Ediff_p3 = std::pow(Ediff_p, 3);

    Vv += z * v1p * (1 - v2p*Ediff_p + v3p*Ediff_p2 - v4p*Ediff_p3
      + Vcbar_p*(v2p - 2*v3p*Ediff_p + 3*v4p*Ediff_p2));
    Wv += z * w1p * Ediff_p2 / (Ediff_p2 + std::pow(w2p, 2));
    Wd += z * d1p * Ediff_p2 * std::exp(-d2p * Ediff_p)
      / (Ediff_p2 + std::pow(d3p, 2));

    Rv += z * Rvp;
    av += z * avp;
    Rd += z * Rdp;
    ad += z * adp;
    Rso += z * Rso_p;
    aso += z * aso_p;

    if (two_s != 0) {
      Vso += vso1p * std::exp(-vso2p * Ediff_p);
      Wso += wso1p * Ediff_p2 / (Ediff_p2 + std::pow(wso2p, 2));
    }
  }

  if (a != 1) {
    Rv /= a;
    av /= a;
    Rd /= a;
    ad /= a;
    Rso /= a;
    aso /= a;
  }

  double f_v = f(r, Rv, av);
  double dfdr_d = dfdr(r, Rd, ad);

  Vv *= f_v;
  Wv *= f_v;
  Wd *= -4 * ad * dfdr_d;

  if (two_s != 0) {
    // Eigenvalue of the spin-orbit operator
    // l.sigma = (j*(j + 1)  - l*(l + 1) -  s*(s + 1))/2
    double spin_orbit_eigenvalue = 0.5*(0.25*(two_j*(two_j + 2)
      - two_s*(two_s + 2)) - l*(l + 1));
    double factor_so = lambda_piplus2 * dfdr(r, Rso, aso)
      * spin_orbit_eigenvalue / r;

    Vso *= factor_so;
    Wso *= factor_so;

    // Adjust spin-orbit potentials for multi-nucleon fragments
    if (a > 1) {
      Vso /= a * std::max(z, n);
    }
  }
 
  return std::complex<double>(-Vv + Vso + Vc(r, Rc, z, Z), -Wv - Wd + Wso);
}


SphericalOpticalModel::SphericalOpticalModel(int z, int a) {
  Z = z;
  A = a;

  int N = A - Z; // Neutron number

  double A_to_the_one_third = std::pow(A, 1.0/3.0);

  // Initialize the spherical optical model parameters (see TALYS 1.6 manual)

  // Neutrons
  v1n = 59.30 - 21.0*(N - Z)/A - 0.024*A; // MeV
  v2n = 0.007228 - 1.48e-6*A; // MeV^(-1)
  v3n = 1.994e-5 - 2.0e-8*A; // MeV^(-2)
  v4n = 7e-9; // MeV^(-3)
  w1n = 12.195 + 0.0167*A; // MeV
  w2n = 73.55 + 0.0795*A; // MeV
  d1n = 16.0 - 16.0*(N - Z)/A; // MeV
  d2n = 0.0180 + 0.003802/(1 + std::exp((A - 156.)/8.0)); // MeV^(-1)
  d3n = 11.5; // MeV
  vso1n = 5.922 + 0.0030*A; // MeV
  vso2n = 0.0040; // MeV^(-1)
  wso1n = -3.1; // MeV
  wso2n = 160.; // MeV
  Efn = -11.2814 + 0.02646*A; // MeV
  Rvn = 1.3039*A_to_the_one_third - 0.4054; // fm
  avn = 0.6778 - 1.487e-4*A; // fm
  Rdn = 1.3424*A_to_the_one_third - 0.01585*std::pow(A_to_the_one_third, 2); // fm
  adn = 0.5446 - 1.656e-4*A; // fm
  Rso_n = 1.1854*A_to_the_one_third - 0.647; // fm
  aso_n = 0.59; // fm

  // Protons
  v1p = 59.30 + 21.0*(N - Z)/A - 0.024*A; // MeV
  v2p = 0.007067 + 4.23e-6*A; // MeV^(-1)
  v3p = 1.729e-5 + 1.136e-8*A; // MeV^(-2)
  v4p = v4n; // MeV^(-3)
  w1p = 14.667 + 0.009629*A; // MeV
  w2p = w2n; // MeV
  d1p = 16.0 + 16.0*(N - Z)/A; // MeV
  d2p = d2n; // MeV^(-1)
  d3p = d3n; // MeV
  vso1p = vso1n; // MeV
  vso2p = vso2n; // MeV^(-1)
  wso1p = wso1n; // MeV
  wso2p = wso2n; // MeV
  Efp = -8.4075 + 0.01378*A; // MeV
  Rvp = Rvn;
  avp = avn;
  Rdp = Rdn;
  adp = 0.5187 + 5.205e-4*A; // fm
  Rso_p = Rso_n;
  aso_p = aso_n;
  Rc = 1.198*A_to_the_one_third + 0.697/A_to_the_one_third
    + 12.994*std::pow(A_to_the_one_third, -4); // fm
  Vcbar_p = 1.73 * Z / Rc; // MeV

  // Compute and store reduced masses for all of the fragments of interest
  double m_nucleus = TMarleyMassTable::get_atomic_mass(Z, A);
  for (auto pid : TMarleyMassTable::get_fragment_pids()) {
    double m_fragment = TMarleyMassTable::get_particle_mass(pid);
    reduced_masses[pid] = (m_fragment * m_nucleus) / (m_fragment + m_nucleus);
  }
}

double SphericalOpticalModel::transmission_coefficient(double E,
  int fragment_pid, int two_j, int l, int two_s, double h)
{
  double h2_over_twelve = std::pow(h, 2) / 12.0;

  // TODO: adjust method for determining matching radius
  double r_max_1 = 12; // Matching radius (fm)
  double r_max_2 = 1.2*r_max_1;

  std::complex<double> u1 = 0, u2 = 0;

  std::complex<double> a_n_minus_two;
  // a(r) really blows up at the origin for the optical model potential, but we're saved
  // by the boundary condition that u(0) = 0. We just need something finite here, but we might
  // as well make it zero.
  std::complex<double> a_n_minus_one = 0;
  std::complex<double> a_n = a(h, E, fragment_pid, two_j, l, two_s);

  std::complex<double> u_n_minus_two;
  // Boundary condition that the wavefunction vanishes at the origin (the optical model
  // potential blows up at r = 0)
  std::complex<double> u_n_minus_one = 0;
  // Asymptotic approximation for a regular potential (see J. Thijssen, Computational
  // Physics, p. 20 for details). We really just need something finite and nonzero here, since
  // our specific choice only determines the overall normalization, which isn't important for
  // determining the transmission coefficients.
  std::complex<double> u_n = std::pow(h, l + 1);

  double r;

  for (r = 2*h; r < r_max_1; r += h) {
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;
    a_n = a(r, E, fragment_pid, two_j, l, two_s);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*h2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + h2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + h2_over_twelve*a_n);
  }

  u1 = u_n;

  for (; r < r_max_2; r += h) {
    a_n_minus_two = a_n_minus_one;
    a_n_minus_one = a_n;
    a_n = a(r, E, fragment_pid, two_j, l, two_s);

    u_n_minus_two = u_n_minus_one;
    u_n_minus_one = u_n;

    u_n = ((2.0 - 10*h2_over_twelve*a_n_minus_one)*u_n_minus_one
      - (1.0 + h2_over_twelve*a_n_minus_two)*u_n_minus_two)
      / (1.0 + h2_over_twelve*a_n);
  }

  u2 = u_n;

  //std::cout << "u1(" << r_max_1 << " fm) = " << u1 << std::endl;
  //std::cout << "u2(" << r_max_2 << " fm) = " << u2 << std::endl;

  // Coulomb parameter
  double mu = get_fragment_reduced_mass(fragment_pid);
  int z = TMarleyMassTable::get_particle_Z(fragment_pid);
  double k = std::sqrt(2.0 * mu * E) / hbarc;
  double eta = mu * Z * z * e2 / (hbarc2 * k);

  //std::cout << "k = " << k << std::endl;
  //std::cout << "eta = " << eta << std::endl;

  // Compute the Coulomb wavefunctions at the matching radii
  Coulomb_wave_functions cwf(true, l, eta);
  std::complex<double> dummy, Hplus1, Hminus1, Hplus2, Hminus2;
  cwf.H_dH(1, k*r_max_1, Hplus1, dummy);
  cwf.H_dH(-1, k*r_max_1, Hminus1, dummy);
  cwf.H_dH(1, k*r_max_2, Hplus2, dummy);
  cwf.H_dH(-1, k*r_max_2, Hminus2, dummy);

  //std::cout << "Hplus1 = " << Hplus1 << std::endl;
  //std::cout << "Hminus1 = " << Hminus1 << std::endl;
  //std::cout << "Hplus2 = " << Hplus2 << std::endl;
  //std::cout << "Hminus2 = " << Hminus2 << std::endl;

  // Compute the transmission coefficient using the radial wavefunction
  // evaluated at the two matching radii
  std::complex<double> S = (u1*Hminus2 - u2*Hminus1) / (u1*Hplus2 - u2*Hplus1);
  return 1.0 - std::norm(S);
}

class BackshiftedFermiGasModel {
  public:
    static double liquid_drop_model_mass_excess(int Z, int A);
    static double level_density(int Z, int A, double Ex); // rho(Ex)
    // rho(Ex, J, Pi) with the assumption of parity equipartition
    static double level_density(int Z, int A, double Ex, int two_J);
  private:
    // Version of the level density function used internally
    static double level_density(int Z, int A, double Ex, double& sigma);

    // Parameters from global level density fit performed in A. J. Koning,
    // et al., Nucl. Phys. A810 (2008) pp. 13-76.
    static const double alpha;
    static const double beta;
    static const double gamma_1;
    static const double delta_global;

    // Liquid drop model parameters
    static const double Mn;
    static const double MH;
    static const double a1;
    static const double a2;
    static const double kappa;
    static const double c3;
    static const double c4;
};

// Parameters from global level density fit performed in A. J. Koning,
// et al., Nucl. Phys. A810 (2008) pp. 13-76.
const double BackshiftedFermiGasModel::alpha = 0.0722396;
const double BackshiftedFermiGasModel::beta = 0.195267;
const double BackshiftedFermiGasModel::gamma_1 = 0.410289;
const double BackshiftedFermiGasModel::delta_global = 0.173015;

// Liquid drop model parameters
const double BackshiftedFermiGasModel::Mn = 8.07144; // MeV
const double BackshiftedFermiGasModel::MH = 7.28899; // MeV
const double BackshiftedFermiGasModel::a1 = 15.677; // MeV
const double BackshiftedFermiGasModel::a2 = 18.56; // MeV
const double BackshiftedFermiGasModel::kappa = 1.79;
const double BackshiftedFermiGasModel::c3 = 0.717; // MeV
const double BackshiftedFermiGasModel::c4 = 1.21129; // MeV

double BackshiftedFermiGasModel::level_density(int Z, int A, double Ex) {
  double dummy;
  return level_density(Z, A, Ex, dummy);
}

double BackshiftedFermiGasModel::level_density(int Z, int A, double Ex, int two_J) {
  double sigma;
  double rho = level_density(Z, A, Ex, sigma);
  double two_sigma2 = 2 * std::pow(sigma, 2);
  return 0.5 * ((two_J + 1) / two_sigma2) * std::exp(-0.25 * std::pow(two_J + 1, 2)
    / two_sigma2) * rho;
}

// Returns the total level density rho(Ex) and store the spin cut-off parameter in sigma
double BackshiftedFermiGasModel::level_density(int Z, int A, double Ex, double& sigma) {
  int N = A - Z;
  double A_to_the_one_third = std::pow(A, 1.0/3.0);

  // Asymptotic level density parameter
  double a_tilde = alpha*A + beta*std::pow(A_to_the_one_third, 2);
  // Damping parameter
  double gamma = gamma_1 / A_to_the_one_third;
  // Shell correction energy
  double delta_W = TMarleyMassTable::get_mass_excess(Z, A)
    - liquid_drop_model_mass_excess(Z, A);
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
  if (Ex <= Delta_BFM) {
    // Use first-order Taylor expansion for small energies
    a = a_tilde * (1 + gamma * delta_W);
  }
  else {
    a = a_tilde * (1 + (delta_W / U) * (1 - std::exp(-gamma * U)));
  }

  // Spin cut-off parameter
  double sigma_d_global = 0.83*std::pow(A, 0.26);
  double Sn = TMarleyMassTable::get_neutron_separation_energy(Z, A);
  double sigma_F2;
  double Ed = 0;
  // Use Delta_BFM as Ed if the excitation energy is too small to avoid
  // numerical problems when computing sigma_F.
  // TODO: Replace Ed here with database of local fits taken from RIPL-3 or TALYS
  //if (U > 0) Ed = 0;
  //else Ed = Delta_BFM;

  if (Ex <= Ed) sigma = sigma_d_global;
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

  double aU = a * U;
  double sqrt_aU = std::sqrt(aU);
  return std::pow(12 * sigma * (std::sqrt(2 * sqrt_aU)*U*std::exp(-2 * sqrt_aU)
    + std::exp(-aU - 1)/a), -1);
}

double BackshiftedFermiGasModel::liquid_drop_model_mass_excess(int Z, int A) {
  int N = A - Z;

  double kappa_term = kappa * std::pow((N - Z) / static_cast<double>(A), 2);
  double c1 = a1 * (1 - kappa_term);
  double c2 = a2 * (1 - kappa_term);

  double Evol = -c1 * A; 
  double Esur = c2 * std::pow(A, 2.0/3.0);
  double Ecoul = (c3 / std::pow(A, 1.0/3.0) - c4 / A) * std::pow(Z, 2);

  double delta_LDM = 0;
  bool z_odd = Z % 2;
  bool n_odd = N % 2;
  // delta_LDM will be zero if the nucleus is odd-even
  if (z_odd && n_odd) delta_LDM = 11/std::sqrt(A);
  else if (!z_odd && !n_odd) delta_LDM = -11/std::sqrt(A);

  return Mn * N + MH * Z + Evol + Esur + Ecoul + delta_LDM;
}

enum class TransitionType { electric, magnetic };
double gamma_strength_function(int Z, int A, TransitionType type,
  int l, double e_gamma)
{
  // TODO: improve this error message
  if (l < 1) throw std::runtime_error(std::string("Invalid multipolarity")
    + std::to_string(l) + " given for gamma ray strength"
    + " function calculation");

  // The strength, energy, and width of the giant resonance for a transition of
  // type x (E or M) and multipolarity l
  double sigma_xl, e_xl, gamma_xl;

  if (type == TransitionType::electric) {
    if (l == 1) {
      e_xl = 31.2*std::pow(A, -1.0/3.0) + 20.6*std::pow(A, -1.0/6.0);
      gamma_xl = 0.026*std::pow(e_xl, 1.91);
      sigma_xl = 1.2*120*(A-Z)*Z/(A*marley_utils::pi*gamma_xl)*mb;
    }
    if (l > 1) {
      // Values for E2 transitions
      e_xl = 63*std::pow(A, -1.0/3.0);
      gamma_xl = 6.11 - 0.012*A;
      sigma_xl = 0.00014 * std::pow(Z, 2) * e_xl
        / (std::pow(A, 1.0/3.0) * gamma_xl) * mb;
      // If this is an E2 transition, we're done. Otherwise,
      // compute the giant resonance strength iteratively
      for (int i = 2; i < l; ++i) {
        sigma_xl *= 8e-4;
      }
    }
  }
  else if (type == TransitionType::magnetic) {
    // Values for M1 transitions
    // The commented-out version is for RIPL-1
    //double factor_m1 = 1.58e-9*std::pow(A, 0.47);
    // RIPL-2 factor
    const double e_gamma_ref = 7.0; // MeV
    double factor_m1 = gamma_strength_function(Z, A,
      TransitionType::electric, 1, e_gamma_ref)
      / (0.0588 * std::pow(A, 0.878));
    gamma_xl = 4.0;
    e_xl = 41*std::pow(A, -1.0/3.0);
    sigma_xl = (std::pow(std::pow(e_gamma_ref, 2) - std::pow(e_xl, 2), 2)
      + std::pow(e_gamma_ref, 2) * std::pow(gamma_xl, 2))
      * (3 * std::pow(marley_utils::pi, 2) * factor_m1)
      / (e_gamma_ref * std::pow(gamma_xl, 2));
    // If this is an M1 transition, we're done. Otherwise,
    // compute the giant resonance strength iteratively
    for (int i = 1; i < l; ++i) {
      sigma_xl *= 8e-4;
    }
  }
  // TODO: improve this error message
  else throw std::runtime_error(std::string("Invalid transition type")
    + " given for gamma ray strength function calculation");

  // Now that we have the appropriate giant resonance parameters,
  // calculate the strength function using the Brink-Axel expression.
  // Note that the strength function has units of MeV^(-3)
  double f_xl = (sigma_xl * std::pow(e_gamma, 3 - 2*l)
    * std::pow(gamma_xl, 2)) / ((2*l + 1) * std::pow(marley_utils::pi, 2)
    * (std::pow(std::pow(e_gamma, 2) - std::pow(e_xl, 2), 2)
    + std::pow(e_gamma, 2) * std::pow(gamma_xl, 2)));

  //std::cout << "E_xl = " << e_xl
  //  << " gamma_xl = " << gamma_xl
  //  << " sigma_xl = " << sigma_xl
  //  << std::endl;

  return f_xl;
}

double gamma_transmission_coefficient(int Z, int A, TransitionType type,
  int l, double e_gamma)
{
  return 2 * marley_utils::pi * gamma_strength_function(Z, A, type,
    l, e_gamma) * std::pow(e_gamma, 2*l + 1);
}

//int main() {
//  int Z = 19;
//  int A = 40;
//  std::string s_type;
//  int l;
//  TransitionType type;
//  double e_gamma;
//
//  std::cout << "Input type, l, and E_gamma" << std::endl;
//  while(true) {
//    std::cin >> s_type >> l >> e_gamma;
//    if (s_type == "e") type = TransitionType::electric;
//    else type = TransitionType::magnetic;
//    std::cout << gamma_strength_function(Z, A, type, l, e_gamma)
//      << std::endl;
//  }
//}

int main() {

  int Z = 19;
  int A = 40;

  SphericalOpticalModel som(Z, A);

  for (double E = 0.1; E < 50.; E += 0.1) {
    double T = som.transmission_coefficient(E, PROTON, 1, 0, 1, 0.1);
    double rho = BackshiftedFermiGasModel::level_density(Z, A, E);
    std::cout << "E = " << E << ", T = " << T << ", rho = "
      << rho << std::endl;
  }
  std::cout << TMarleyMassTable::get_mass_excess(Z, A)
    - BackshiftedFermiGasModel::liquid_drop_model_mass_excess(Z, A) << std::endl;

}
