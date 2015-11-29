#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "marley_utils.hh"
#include "meta_numerics.hh"
#include "TMarleyIntegrator.hh"
#include "TMarleyMassTable.hh"

class WoodsSaxon {
  public:
    inline WoodsSaxon(int z, int a) {
      if (z < 1) throw std::runtime_error(std::string("Cannot construct")
        + " a WoodsSaxon object for Z < 1");
      if (a < 1) throw std::runtime_error(std::string("Cannot construct")
        + " a WoodsSaxon object for A < 1");
      Z = z;
      A = a;
      double A_to_the_minus_one_third = std::pow(A, -1.0/3.0);
      hbar_omega = 45*A_to_the_minus_one_third
        - 25*std::pow(A_to_the_minus_one_third, 2);
      mN = 0.5 * (TMarleyMassTable::get_particle_mass(
        marley_utils::PROTON) + TMarleyMassTable::get_particle_mass(
        marley_utils::NEUTRON));
      b = marley_utils::hbar_c / std::sqrt(mN * hbar_omega);
      b2 = std::pow(b, 2);
      int N = A - Z;
      double factor = (static_cast<double>(N) - Z) / A;
      V0n = 51 - 33*factor; // MeV
      V0p = 51 + 33*factor; // MeV
      R = r0 / A_to_the_minus_one_third; // Nuclear radius (fm)
    }

    // The radial wavefunction for a 3D harmonic oscillator with principal
    // quantum number n (n >= 0), orbital angular momentum l (l >= 0), and
    // oscillator length b.
    double radial_3d_ho_wavefunction(int n, int l, double r) const;

    // Woods-Saxon shape
    inline double f(double r) const {
      return std::pow(1 + std::exp((r - R) / a), -1);
    }

    // Partial derivative with respect to r of the Woods-Saxon shape
    inline double dfdr(double r) const {
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
    inline double Vc(double r) const {
      if (r < R) return Z * marley_utils::e2
        * (3 - std::pow(r / R, 2)) / (2 * R);
      else return Z * marley_utils::e2 / r;
    }

    inline double matrix_element(const TMarleyIntegrator& integrator,
      int two_j, int l, int pid, int nu_1, int nu_2) const
    {
      std::function<double(double)> f = std::bind(
        &WoodsSaxon::transformed_integrand, this, std::placeholders::_1,
        two_j, l, pid, nu_1, nu_2);
      return integrator.num_integrate(f, 0, 1);
    }

    inline double transformed_integrand(double t, int two_j, int l,
      int pid, int nu_1, int nu_2) const
    {
      if (t == 1 || t == 0) return 0.;
      else return matrix_element_integrand(t / (1 - t), two_j, l, pid, nu_1,
        nu_2) / std::pow(1 - t, 2);
    }

    inline double matrix_element_integrand(double r, int two_j, int l,
      int pid, int nu_1, int nu_2) const
    {
      double V0, vc;
      if (pid == marley_utils::NEUTRON) {
        V0 = V0n;
        vc = 0.;
      }
      else if (pid == marley_utils::PROTON) {
        V0 = V0p;
        vc = Vc(r);
      }
      else throw std::runtime_error(std::string("Woods-Saxon matrix elements")
        + " may only be computed for either protons or neutrons.");
      double v0ls = 0.44*V0;
      double kinetic_and_centrifugal_terms = (marley_utils::hbar_c2 / (2*mN))
        * ((4*nu_1 + 2*l + 3)/b2 - std::pow(r / b2, 2));
      double vws = -V0 * f(r);
      double vls = v0ls * (r02 / r) * dfdr(r);
      vls *= (two_j*(two_j + 2) - 4*l*(l + 1) - 3)/8.0;
      double g_nu_1_l = radial_3d_ho_wavefunction(nu_1, l, r);
      double g_nu_2_l = radial_3d_ho_wavefunction(nu_2, l, r);
      return std::pow(r, 2) * g_nu_1_l * g_nu_2_l
        * (kinetic_and_centrifugal_terms + vws + vc + vls);
    }

  private:
    int Z, A; // Atomic number and mass number
    double b; // Oscillator length (fm)
    double b2; // Oscillator length squared (fm^2)
    double hbar_omega; // Energy quantum (MeV)
    double mN; // Nucleon mass
    // Woods-Saxon parameters for the phenomenological nuclear potential
    constexpr static double r0 = 1.27; // fm
    constexpr static double r02 = std::pow(r0, 2); // fm^2
    constexpr static double a = 0.67; // Surface diffuseness (fm)
    double V0n, V0p; // Woods Saxon potential depths for neutrons and protons
    double R; // Nuclear radius
};

double WoodsSaxon::radial_3d_ho_wavefunction(int n, int l, double r) const {
  if (n < 0) throw std::runtime_error(std::string("Cannot compute")
    + " harmonic oscillator wavefunction for principal quantum number"
    + " n = " + std::to_string(n));
  if (l < 0) throw std::runtime_error(std::string("Cannot compute")
    + " harmonic oscillator wavefunction for orbital angular momentum"
    + " l = " + std::to_string(l));
  double num = meta_numerics::LogFactorial(n);
  double denom = meta_numerics::LogGamma(n + l + 1.5);
  double norm_factor = std::sqrt(2*std::exp(num - denom) / std::pow(b, 3));
  double r_over_b = r / b;
  double r_over_b_squared = std::pow(r_over_b, 2);
  return norm_factor * std::pow(r_over_b, l) * std::exp(-r_over_b_squared / 2)
    * meta_numerics::LaguerreL(n, l + 0.5, r_over_b_squared);
}

#include <iomanip>
int main() {
  std::cout << std::setprecision(16);// << std::scientific;
  TMarleyIntegrator integrator(100);
  WoodsSaxon ws(19, 40);
  std::cout << ws.matrix_element(integrator, 1, 0,
    marley_utils::NEUTRON, 0, 0);
  //std::function<double(double)> f = [&ws](double t) -> double {
  //  if (t == 1) return 0.;
  //  double r = t / (1 - t);
  //  double g1 = ws.radial_3d_ho_wavefunction(1, 0, r);
  //  double g2 = ws.radial_3d_ho_wavefunction(1, 0, r);
  //  return std::pow(r, 2) * g1 * g2 / std::pow(1 - t, 2);
  //};
  //std::cout << integrator.num_integrate(f, 0, 1) << std::endl;
  return 0;
}
