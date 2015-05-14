#include <cmath>
#include <complex>
#include <iostream>

#include "marley_utils.hh"
#include "TMarleyReaction.hh"

TMarleyReaction::TMarleyReaction() {
  // TODO: add stuff here
}


// TODO: Fix both Fermi functions to use natural units and energies in MeV (this keeps
// everything consistent)

// Fermi function used in calculating cross-sections
// The equation was taken from http://en.wikipedia.org/wiki/Beta_decay
// Input: atomic number Z and electron kinetic energy T (in keV), positron or electron?
double TMarleyReaction::fermi_function(double z, double t, bool electron){
  
  double charge = z;
  double kin_e = t;
  double atomic_mass = 40.;
  
  // Constants
  const double alpha = 7.297e-3; // Fine-structure constant
  const double hbar_c = 1.97e+3; // hbar*c in keV*fm
  const double e = std::exp(1); // Euler's number
  const double pi = 4*std::atan(1);
  const double r_0 = 1.2; // in fm
  const double m = 510.9989; // in keV

  double eta;
  double p = std::sqrt(kin_e*kin_e + 2*m*kin_e); // Electron momentum
  double s = std::sqrt(1 - alpha*alpha*charge*charge);
  double energy = std::sqrt(p*p + m*m);
  double r_n = r_0*std::pow(atomic_mass, 1./3);
  double rho = r_n/hbar_c;
  
  if(electron)
    eta = alpha*charge*energy/p;
  else
    eta = -alpha*charge*energy/p;

  // Complex variables for the gamma function
  std::complex<double> a(s, eta);
  std::complex<double> b(1+2*s, 0);

  return 2*(1 + s)* std::pow(2*p*rho, 2*s-2)*std::pow(e, pi*eta)*std::norm(marley_utils::gamma(a))
    / (std::abs(marley_utils::gamma(b)) * std::abs(marley_utils::gamma(b)));
}

double TMarleyReaction::fermi_approx(int z, double t, bool electron){

  // This is a test for comparison with the "exact" calculation, which involves
  // complex gamma functions. This holds true for Q << mc^2
  
  int charge = z;
  double kin_e = t;

  const double alpha = 7.297e-3; // Fine-structure constant
  const double e = std::exp(1); // Euler's number
  const double pi = 4*std::atan(1);
  const double m = 510.9989; // in keV

  double eta;
  double p = std::sqrt(kin_e*kin_e + 2*kin_e*m);
  double energy = std::sqrt(p*p + m*m);
  
  if(electron)
    eta = alpha*charge*energy/p;
  else
    eta = -alpha*charge*energy/p;

  return 2*pi*eta/(1-std::pow(e, -2*pi*eta));
}

double TMarleyReaction::ejectile_energy(double E_level, double Ea, double cos_theta_c) {

  // All energies are in MeV
  double md = md_gs + E_level;

  double Etot = Ea + mb;

  // Squared 3-momentum of the projectile
  double p2a = Ea*Ea - ma*ma;

  // Quadratic formula coefficients from 4-momentum
  // conservation solution in the lab frame
  double A = 4*Etot*Etot - 4*p2a*cos_theta_c*cos_theta_c;

  double B = 4*Etot*(md*md + p2a - mc*mc - Etot*Etot);

  double C = std::pow((Etot*Etot + mc*mc - md*md - p2a), 2)
    + 4*mc*mc*p2a*cos_theta_c*cos_theta_c;

  // Restructure the calculation to avoid some potentially bad cancellations
  // See http://www.petebecker.com/js/js200010.html for details
  double c = C/A;
  double b = B/(2*A);

  // There are two solutions of the quadratic equation 
  // that correspond to forward and backward scattering
  // solutions
  double Ecplus, Ecminus;   if (b > 0) {
    Ecminus = -b - std::sqrt(b*b - c);
    Ecplus = c/Ecminus;
  }
  else {
    Ecplus = -b + std::sqrt(b*b - c);
    Ecminus = c/Ecplus;
  }

  // Assume we are above the backward threshold, and match the
  // correct solution to the physical ejectile energy.
  // TODO: generalize this to include cases where we're below
  // the backward threshold
  if (cos_theta_c >= 0) {
    return Ecplus; // Forward scattering
  }
  else {
    return Ecminus; // Backward scattering
  }

}

// TODO: change this to return an event object
void TMarleyReaction::create_event(double Ea) {

  // All hard-coded energies are in MeV
  double Etot = Ea + mb;

  // TODO: change this to sample a final 40K
  // level energy
  double E_level = 2.289871;
  double md = md_gs + E_level;

  // TODO: change this to sample a scattering angle
  // for the ejectile
  double theta_c = std::acos(0);
  double cos_theta_c = std::cos(theta_c);

  // Use conservation of 4-momentum to compute the ejectile energy
  // based on the sampled scattering angle
  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);

  // Compute the energy and scattering angle of the residue
  double Ed = Etot - Ec;
  double theta_d = std::asin(std::sqrt(Ec*Ec - mc*mc)
    * std::sin(theta_c) / std::sqrt(Ed*Ed - md*md));
  double cos_theta_d = std::cos(theta_d);

  // Print results to std::cout
  std::cout.precision(15);
  std::cout << std::scientific;
  std::cout << "Ec = " << Ec << std::endl;
  std::cout << "Ed = " << Ed << std::endl;
  std::cout << "cos_theta_d = " << cos_theta_d << std::endl;
}

// Compute the differential reaction cross section dsigma/dcos_theta
// for a given final residue level, projectile energy, and ejectile
// scattering cosine
double TMarleyReaction::differential_xs(double E_level, double Ea, double cos_theta_c) {

  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);
  double pc = std::sqrt(Ec*Ec - mc*mc);

  // TODO: compute the matrix element using B(F) + B(GT)
  // values and/or some other method
  double matrix_element = 1;
  
  // TODO: consider removing the constants GF, Vud, etc. since
  // for sampling they will not be needed (what matters is the
  // relative size of the differential cross section).
  // This may reduce numerical error.
  
  // TODO: adjust this differential cross section expression
  // as needed
  
  // This expression for the differential cross section
  // dsigma/dcos_theta comes from Kuramoto, et al.,
  // Nucl. Phys. A512 (1990) 711-736. It is modified
  // using a nuclear recoil correction factor taken from
  // J. D. Walecka, "Semileptonic Weak Interactions in Nuclei,"
  // In: Muon Physics, Volume II: Weak Interactions.
  // Ed. by V. W. Hughes and C. S. Wu.
  return (1.0/(2*std::acos(-1)))*GF*GF*Vud*Vud
    *pc*Ec*fermi_function(Zf, Ec - mc, true)*matrix_element
    /(1.0 + (Ea/mb)*(1 - (Ec/pc)*cos_theta_c));
}

// TODO: Consider adding this function to marley_utils
// rather than the TMarleyReaction class

// Numerically integrate a given function f (that takes a
// double argument to integrate over and returns a double)
// over the interval [a,b] using the composite trapezoidal
// rule over n subintervals.
// (see http://en.wikipedia.org/wiki/Numerical_integration)
double TMarleyReaction::num_integrate(const std::function<double(double)> &f,
  double a, double b, int n)
{
  double integral = 0; 
  for(int k = 1; k < n-1; k++) {
    integral += ((b - a)/n)*f(a + k*(b - a)/n);
  }
  integral += ((b - a)/n)*(f(a)/2 + f(b)/2);
  return integral;
}

// Numerically integrate the differential cross section over
// cos(theta_c) = [-1,1] to get the total reaction cross section
// for a given final residue level and projectile energy
double TMarleyReaction::total_xs(double E_level, double Ea) {
  // TODO: consider switching to a different integration method,
  // using tabulated data, etc.

  // Numerically integrate the differential cross section over
  // the interval cos(theta_c) = [-1,1] using the composite
  // trapezoidal rule over n subintervals
  static const int n = 1e3;

  // Create a forwarding call wrapper for the  differential_xs
  // member function that takes a single argument.
  std::function<double(double)> dxs = std::bind(
    &TMarleyReaction::differential_xs, this, E_level, Ea, std::placeholders::_1);

  // Numerically integrate using the call wrapper, the integration bounds,
  // and the number of subintervals
  return num_integrate(dxs, -1.0, 1.0, n); 
}
