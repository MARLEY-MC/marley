#include <cmath>
#include <complex>
#include <iostream>

// TODO: fix this ugly include pattern that the
// cwfcomp people forced you to use
#include "complex_functions.hh"
#include "cwfcomp.hh"
#include "cwfcomp.cc"

#include "/home/sjg/Desktop/marley/TMarleyMassTable.hh"

// Constants
const double alpha = 7.2973525698e-3; // Fine structure constant
const double hbarc = 197.3269631; // MeV * fm
const double hbarc2 = std::pow(hbarc, 2); // MeV^2 * fm^2
const double e2 = hbarc * alpha; // Elementary charge (in MeV*fm)
// Mass of a charged pion
const double mpiplus = 139.57018; // MeV
// Squared pion Compton wavelength in fm
const double lambda_piplus2 = std::pow(hbarc / mpiplus, 2);

// Particle IDs
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

// Using the Numerov method (see https://en.wikipedia.org/wiki/Numerov%27s_method for details),
// compute the value of the radial wavefunction u(r) that satisfies the Schrödinger equation
// (d^2/dr^2 + a(r))u(r) = 0, where a(r) = -l*(l+1)/r^2 + k^2 - 2*mu*V(r)/hbar^2,
// V(r) is the potential, and the wavenumber k = sqrt(2*mu*E)/hbar, where E is
// the particle's (nonrelativistic) energy and mu is its reduced mass. In this function,
// we use 2*l rather than l in the input so that half-integer orbital angular momenta
// may be represented using integer values.
//template<typename numType>
//numType numerov_radial_wavefunction(std::function<numType(double)> V,
//  double R, int l, double mu, double E, numType u0, numType u1, double h)
//{
//  double h2_over_twelve = std::pow(h, 2) / 12.0;
//
//  double centrifugal_term = -l * (l + 1) / std::pow(r, 2);
//  double factor = 2 * mu / hbarc2;
//  double k2 = factor * E;
//  double constant_terms = centrifugal_term + k2;
//
//  numType u_n = u1;
//  numType u_n_minus_one = u0;
//
//  numType a_n_minus_one, a_n, a_n_plus_one, temp;
//
//  double r;
//
//  for (r = 2*h,
//    a_n = constant_terms - factor * V(h),
//    a_n_plus_one = constant_terms - factor * V(2*h); r < R; r += h)
//  {
//    a_n_minus_one = a_n;
//    a_n = a_n_plus_one;
//    a_n_plus_one = constant_terms - factor * V(r + h);
//
//    temp = u_n;
//
//    u_n = ((2.0 - 10*h2_over_twelve*a_n)*u_n
//      - (1.0 + h2_over_twelve*a_n_minus_one)*u_n_minus_one)
//      / (1.0 + h2_over_twelve*a_n_plus_one);
//
//    u_n_minus_one = temp;
//  }
//
//  return u_n;
//}

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

int main() {

  double Z = 19;
  int A = 40;

  SphericalOpticalModel som(Z, A);

  for (double E = 0.1; E < 50.; E += 0.1) {
    double T = som.transmission_coefficient(E, PROTON, 1, 0, 1, 0.1);
    std::cout << "E = " << E << ", T = " << T << std::endl;
  }

}
