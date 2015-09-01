#pragma once
#include <cmath>
#include <complex>
#include <vector>

#include "TMarleyDecayScheme.hh"
#include "TMarleyMassTable.hh"
#include "marley_utils.hh"

class TMarleySphericalOpticalModel {
  public:
    TMarleySphericalOpticalModel(int z, int a);
    std::complex<double> optical_model_potential(double r, double E,
      int fragment_pid, int two_j, int l, int two_s) const;
    double transmission_coefficient(double E, int fragment_pid, int two_j,
      int l, int two_s, double h) const;

    inline int get_Z() const {
      return Z;
    }

    inline int get_A() const {
      return A;
    }

  private:
    // Woods-Saxon shape
    inline double f(double r, double R, double a) const {
      return std::pow(1 + std::exp((r - R) / a), -1);
    }

    inline double get_fragment_reduced_mass(int fragment_pid) const {
      return reduced_masses.at(fragment_pid);
    }

    // Partial derivative with respect to r of the Woods-Saxon shape
    inline double dfdr(double r, double R, double a) const {
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
    inline double Vc(double r, double R, int Z, int z) const {
      if (Z == 0 || z == 0) return 0;
      else if (r < R) return Z * z * marley_utils::e2
        * (3 - std::pow(r / R, 2)) / (2 * R);
      else return Z * z * marley_utils::e2 / r;
    }

    // Non-derivative radial Schrödinger equation terms to use for computing transmission
    // coefficients via the Numerov method
    inline std::complex<double> a(double r, double E, int fragment_pid,
      int two_j, int l, int two_s) const
    {
      return (-l*(l+1) / std::pow(r, 2)) +
        2 * reduced_masses.at(fragment_pid) * (E - optical_model_potential(r,
        E, fragment_pid, two_j, l, two_s)) / marley_utils::hbar_c2;
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

    // Mass of a charged pion
    static constexpr double mpiplus = 139.57018; // MeV
    // Squared pion Compton wavelength
    static constexpr double lambda_piplus2 = std::pow(marley_utils::hbar_c
      / mpiplus, 2); // fm
};
