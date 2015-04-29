// Fermi function using in calculating cross-sections
// The equation was taken from http://en.wikipedia.org/wiki/Beta_decay

#include <cmath>
#include <complex>
#include <iostream>

#include "marley_utils.hh"

// Input: atomic number Z and electron kinetic energy T (in keV), positron or electron?
double fermi_function(double z, double t, bool electron){
  
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

double fermi_approx(int z, double t, bool electron){

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

int main(){

  std::cout << "\"Exact:\" " << fermi_function(18., 20, true) << std::endl;
  std::cout << "Approximate: " << fermi_approx(18., 20, true) << std::endl; 
}
