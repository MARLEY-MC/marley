#pragma once
#include <cmath>
#include <complex>
#include <vector>

#include "DecayScheme.hh"
#include "MassTable.hh"
#include "marley_utils.hh"

namespace marley {

  class SphericalOpticalModel {
    public:
      SphericalOpticalModel(int z, int a);
      std::complex<double> optical_model_potential(double r, double E,
        int fragment_pid, int two_j, int l, int two_s);
  
      double transmission_coefficient(double E, int fragment_pid, int two_j,
        int l, int two_s, double h);
  
      double total_cross_section(double E, int fragment_pid, int two_s,
        size_t l_max, double h);
  
      inline int get_Z() const {
        return Z_;
      }
  
      inline int get_A() const {
        return A_;
      }
  
    private:
  
      // Helper function for computing optical model transmission coefficients
      // and cross sections
      std::complex<double> s_matrix_element(double E, int fragment_pid, int two_j,
        int l, int two_s, double h);
  
      // Helper functions for computing the optical model potential
      void calculate_om_parameters(double E, int fragment_pid, int two_j,
        int l, int two_s);
  
      inline std::complex<double> omp(double r) const {
        return omp_minus_Vc(r) + Vc(r, Rc, z, Z_);
      }
  
      // Computes the optical model potential minus the Coulomb potential at
      // radius r
      std::complex<double> omp_minus_Vc(double r) const;
  
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
  
      // Coulomb potential for a point particle with charge q*e interacting with a
      // uniformly charged sphere with radius R and charge Q*e
      inline double Vc(double r, double R, int Q, int q) const {
        if (Q == 0 || q == 0) return 0;
        else if (r < R) return Q * q * marley_utils::e2
          * (3 - std::pow(r / R, 2)) / (2 * R);
        else return Q * q * marley_utils::e2 / r;
      }
  
      // Non-derivative radial Schrödinger equation terms to use for computing
      // transmission coefficients via the Numerov method
      inline std::complex<double> a(double r, double E, int fragment_pid, int l) //const
      {
        return (-l*(l+1) / std::pow(r, 2)) +
          2 * reduced_masses.at(fragment_pid) * (E - omp(r))
          / marley_utils::hbar_c2;
      }
  
      // Version of Schrodinger equation terms with the optical model potential U
      // pre-computed
      inline std::complex<double> a(double r, double E, int fragment_pid, int l,
        std::complex<double> U) const
      {
        return (-l*(l+1) / std::pow(r, 2)) +
          2 * reduced_masses.at(fragment_pid) * (E - U) / marley_utils::hbar_c2;
      }
  
      // Nuclear atomic and mass numbers
      int Z_, A_;
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
  
      // Temporary storage for optical model calculations
      double Rv, av, Rd, ad, Rso, aso; // Geometrical parameters
      double Vv, Wv, Wd, Vso, Wso; // Energy-dependent terms in the potential
      double spin_orbit_eigenvalue; // Eigenvalue of the spin-orbit operator
      int z; // Fragment atomic number
  
      // Threshold for abs(U - Vc) used to find a suitable matching radius for
      // computing transmission coefficients.
      // WARNING: This value should be chosen carefully.  Since the Numerov
      // method used for computing the fragment wavefunctions is only accurate to
      // order h^4 (where h is the step size used), choosing this threshold to be
      // comparable to or smaller than h^4 may cause numerical problems.
      // A value of h^3 seems to work pretty well. You may be able to get away
      // with a higher threshold, but check to make sure the transmission coefficients
      // aren't significantly affected before adopting a higher value.
      static constexpr double MATCHING_RADIUS_THRESHOLD = 1e-3;
  };

}
