#include <cmath>
#include <complex>
#include <iostream>

// TODO: fix this ugly include pattern that the
// cwfcomp people forced you to use
#include "complex_functions.hh"
#include "cwfcomp.hh"
#include "cwfcomp.cc"

// Constants
const double alpha = 7.2973525698e-3; // Fine structure constant
const double hbarc = 197.3269631; // MeV * fm
const double hbarc2 = std::pow(hbarc, 2); // MeV^2 * fm^2
const double e2 = hbarc * alpha; // Elementary charge (in MeV*fm)
// Mass of a charged pion
const double mpiplus = 139.57018; // MeV
// Squared pion Compton wavelength in fm
const double lambda_piplus2 = std::pow(hbarc / mpiplus, 2);

const double micro_amu = 0.000931494061;
const double m_fragment = 1008664.91585 * micro_amu; //TMarleyMassTable::get_particle_mass(2112);
const double m_nucleus = 39963998.166 * micro_amu; //TMarleyMassTable::get_atomic_mass(19, 40);
const double mu  = 1 / (1/m_fragment + 1/m_nucleus); // Reduced mass

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

std::complex<double> optical_model_potential(double r, double E,
  int z, int Z, int A, double j, double l)
{
  //double SnA = TMarleyMassTable::get_neutron_separation_energy(Z, A);
  //double SnB = TMarleyMassTable::get_neutron_separation_energy(Z, A + 1);
  //double Ef = -0.5*(SnA + SnB); // Fermi energy
  const double ONE_THIRD = 1.0/3.0;

  double N_minus_Z_over_A = (A - 2.0*Z)/A;
  double v1 = 59.30 - 21.0*N_minus_Z_over_A - 0.024*A; // MeV
  double v2 = 0.007228 - 1.48e-6*A; // MeV^(-1)
  double v3 = 1.994e-5 - 2.0e-8*A; // MeV^(-2)
  double v4 = 7e-9; // MeV^(-3)
  double w1 = 12.195 + 0.0167*A; // MeV
  double w2 = 73.55 + 0.0795*A; // MeV
  double d1 = 16.0 - 16.0*N_minus_Z_over_A; // MeV
  double d2 = 0.0180 + 0.003802/(1 + std::exp((A - 156.0)/8.0)); // MeV^(-1)
  double d3 = 11.5; // MeV
  double vso1 = 5.922 + 0.0030*A; // MeV
  double vso2 = 0.0040; // MeV^(-1)
  double wso1 = -3.1; // MeV
  double wso2 = 160.0; // MeV
  double Ef = -11.2814 + 0.02646*A; // MeV
  
  double A_to_the_one_third = std::pow(A, ONE_THIRD);
  double Rv = 1.3039*A_to_the_one_third - 0.4054; // fm
  double av = 0.6778 - 1.487e-4*A; // fm
  double Rd = 1.3424*A_to_the_one_third - 0.01585; // fm
  double ad = 0.5446 - 1.656e-4*A; // fm // Different for n,p
  double Rso = 1.1854*A_to_the_one_third - 0.647; // fm
  double aso = 0.59; // fm
  double Rc = 1.198*A_to_the_one_third + 0.697/A_to_the_one_third
    + 12.994*std::pow(A_to_the_one_third, -4);

  double Ediff = E - Ef;
  double Ediff2 = std::pow(Ediff, 2);
  double Ediff3 = std::pow(Ediff, 3);

  double f_v = f(r, Rv, av);
  double dfdr_d = dfdr(r, Rd, ad);
  const double THREE_EIGHTHS = 3.0/8.0;
  // Eigenvalue of the spin-orbit operator
  // l.sigma = j*(j + 1)  - l*(l + 1) -  s*(s + 1)/2
  // where the particle spin s = 1/2 for nucleons
  double spin_orbit_eigenvalue = 0.5*(j*(j + 1) - l*(l + 1) - THREE_EIGHTHS);
  double factor_so = lambda_piplus2 * dfdr(r, Rso, aso)
    * spin_orbit_eigenvalue / r;

  double Vv = v1 * (1 - v2*Ediff + v3*Ediff2 - v4*Ediff3) * f_v;
  double Wv = w1 * Ediff2 / (Ediff2 + std::pow(w2, 2)) * f_v;
  double Wd = d1 * Ediff2 * std::exp(-d2 * Ediff) * dfdr_d
    / (Ediff2 + std::pow(d3, 2));
  double Vso = vso1 * std::exp(-vso2 * Ediff) * factor_so;
  double Wso = wso1 * Ediff2 * factor_so / (Ediff2 + std::pow(wso2, 2));

  //std::cout << "V = " << Vv << std::endl;
  //std::cout << "W = " << Wv << std::endl;
  //std::cout << "Wd = " << Wd << std::endl;
  //std::cout << "Vso = " << Vso << std::endl;
  //std::cout << "Wso = " << Wso << std::endl;

  return std::complex<double>(-Vv + Vso + Vc(r, Rc, z, Z), -Wv - Wd + Wso);
}


inline std::complex<double> a(double r, double E, int z, int Z,
  int A, double j, double l)
{
  return -l*(l + 1)/std::pow(r, 2) + 2*mu*E/hbarc2
    - 2*mu*optical_model_potential(r, E, z, Z, A, j, l)/hbarc2;
}

int main() {

  double Z = 19;
  double z = 0;
  int A = 40;
  double j = 0.5;
  double l = 0.0;

  double h = 0.1; // Step size (fm)
  double h2 = std::pow(h, 2);
  std::complex<double> u1 = 0, u2 = 0, u3 = 0;
  std::complex<double> u_n = h;
  std::complex<double> u_n_minus_one = 0;
  std::complex<double> temp = 0;
  const double FIVE_SIXTHS = 5.0/6.0;
  const double ONE_TWELFTH = 1.0/12.0;

  double r;
  double r_max_1 = 12; // Matching radius (fm)
  double r_max_2 = 1.5*r_max_1;
  double r_max_3 = 2*r_max_1;

  double E = 0.0013; // MeV
  std::cout << "E = " << E << std::endl;

  for (r = 2*h; r < r_max_1; r += h) {
    temp = u_n;
    u_n = ((2.0 - FIVE_SIXTHS*h2*a(r, E, z, Z, A, j, l))*u_n
      - (1.0 + ONE_TWELFTH*h2*a(r - h, E, z, Z, A, j, l))*u_n_minus_one)
      / (1.0 + ONE_TWELFTH*h2*a(r + h, E, z, Z, A, j, l));
    u_n_minus_one = temp;
  }

  u1 = u_n;

  for (; r < r_max_2; r += h) {
    temp = u_n;
    u_n = ((2.0 - FIVE_SIXTHS*h2*a(r, E, z, Z, A, j, l))*u_n
      - (1.0 + ONE_TWELFTH*h2*a(r - h, E, z, Z, A, j, l))*u_n_minus_one)
      / (1.0 + ONE_TWELFTH*h2*a(r + h, E, z, Z, A, j, l));
    u_n_minus_one = temp;
  }

  u2 = u_n;

  for (; r < r_max_3; r += h) {
    temp = u_n;
    u_n = ((2.0 - FIVE_SIXTHS*h2*a(r, E, z, Z, A, j, l))*u_n
      - (1.0 + ONE_TWELFTH*h2*a(r - h, E, z, Z, A, j, l))*u_n_minus_one)
      / (1.0 + ONE_TWELFTH*h2*a(r + h, E, z, Z, A, j, l));
    u_n_minus_one = temp;
  }

  u3 = u_n;

  std::cout << "u1(" << r_max_1 << " fm) = " << u1 << std::endl;
  std::cout << "u2(" << r_max_2 << " fm) = " << u2 << std::endl;
  std::cout << "u3(" << r_max_3 << " fm) = " << u3 << std::endl;

  // Coulomb parameter
  double k = std::sqrt(2.0 * mu * E) / hbarc;
  double eta = mu * Z * z * e2 / (hbarc2 * k);

  std::cout << "k = " << k << std::endl;
  std::cout << "eta = " << eta << std::endl;

  // Compute the Coulomb wavefunctions at the matching radii
  Coulomb_wave_functions cwf(false, 0, eta);
  std::complex<double> dummy, Hplus1, Hminus1, Hplus2, Hminus2, Hplus3, Hminus3;
  //std::complex<double> dummy, F1, G1, F2, G2, F3, G3;
  //cwf.F_dF(k*r_max_1, F1, dummy);
  //cwf.G_dG(k*r_max_1, G1, dummy);
  //cwf.F_dF(k*r_max_2, F2, dummy);
  //cwf.G_dG(k*r_max_2, G2, dummy);
  //cwf.F_dF(k*r_max_3, F3, dummy);
  //cwf.G_dG(k*r_max_3, G3, dummy);
  cwf.H_dH(1, k*r_max_1, Hplus1, dummy);
  cwf.H_dH(-1, k*r_max_1, Hminus1, dummy);
  cwf.H_dH(1, k*r_max_2, Hplus2, dummy);
  cwf.H_dH(-1, k*r_max_2, Hminus2, dummy);
  cwf.H_dH(1, k*r_max_3, Hplus3, dummy);
  cwf.H_dH(-1, k*r_max_3, Hminus3, dummy);

  std::cout << "Hplus1 = " << Hplus1 << std::endl;
  std::cout << "Hminus1 = " << Hminus1 << std::endl;
  std::cout << "Hplus2 = " << Hplus2 << std::endl;
  std::cout << "Hminus2 = " << Hminus2 << std::endl;
  std::cout << "Hplus3 = " << Hplus3 << std::endl;
  std::cout << "Hminus3 = " << Hminus3 << std::endl;

  // Compute the transmission coefficient a few different ways
  // to compare different matching radii
  std::complex<double> S1 = (u1*Hminus2 - u2*Hminus1) / (u1*Hplus2 - u2*Hplus1);
  std::complex<double> S2 = (u1*Hminus3 - u3*Hminus1) / (u1*Hplus3 - u3*Hplus1);
  std::complex<double> S3 = (u3*Hminus2 - u2*Hminus3) / (u3*Hplus2 - u2*Hplus3);
  double T1 = 1.0 - std::norm(S1);
  double T2 = 1.0 - std::norm(S2);
  double T3 = 1.0 - std::norm(S3);
  std::cout << "S1 = " << S1 << std::endl;
  std::cout << "S2 = " << S2 << std::endl;
  std::cout << "S3 = " << S3 << std::endl;
  std::cout << "T1 = " << T1 << std::endl;
  std::cout << "T2 = " << T2 << std::endl;
  std::cout << "T3 = " << T3 << std::endl;
}
