/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#pragma once
#include <cmath>
#include <complex>
#include <vector>

#include "marley/DecayScheme.hh"
#include "marley/MassTable.hh"
#include "marley/OpticalModel.hh"
#include "marley/marley_utils.hh"

namespace marley {

  /// @brief Nuclear optical model for fragment emission calculations
  /// @details This class implements the global optical model potential
  /// of A. J. Koning and J. P. Delaroche, <a
  /// href="http://dx.doi.org/10.1016/S0375-9474(02)01321-0">
  /// Nucl. Phys. A 713 (2003) 231-310</a>.
  /// <a href="https://en.wikipedia.org/wiki/Numerov%27s_method">Numerov's
  /// method</a> is used to integrate the Schr&ouml;dinger equation during
  /// transmission coefficient and cross section calculations.
  class KoningDelarocheOpticalModel : public OpticalModel {

    public:

      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      /// @param step_size Step size (fm) to use for numerical integration of
      /// the Schr&ouml;dinger equation
      KoningDelarocheOpticalModel(int Z, int A, double step_size
        = DEFAULT_NUMEROV_STEP_SIZE_);

      virtual std::complex<double> optical_model_potential(double r,
        double fragment_KE_lab, int fragment_pdg, int two_j, int l, int two_s,
        int target_charge = 0) override;

      virtual double transmission_coefficient(double total_KE_CM,
        int fragment_pdg, int two_j, int l, int two_s, int target_charge = 0)
        override;

      virtual double total_cross_section(double fragment_KE_lab,
        int fragment_pdg, int two_s, size_t l_max, int target_charge = 0)
        override;

    private:

      /// @brief Request-local state used by optical-model calculations
      /// @details This removes mutable calculation state from the model object.
      /// External numerical-library calls require their own thread-safety
      /// guarantees before calculations can safely run concurrently.
      struct WorkingState {
        double total_CM_frame_KE = 0.;
        double fragment_mass = 0.;
        double fragment_KE_lab = 0.;
        double CM_frame_momentum_squared = 0.;
        double target_mass = 0.;

        double Rv = 0., av = 0., Rd = 0., ad = 0.;
        double Rso = 0., aso = 0.;
        double Vv = 0., Wv = 0., Wd = 0., Vso = 0., Wso = 0.;
        double spin_orbit_eigenvalue = 0.;
        int z = 0;
      };

      // Helper function for computing optical model transmission coefficients
      // and cross sections
      std::complex<double> s_matrix_element(int fragment_pdg, int two_j,
        int l, int two_s, WorkingState state) const;

      // Helper functions for computing the optical model potential
      void calculate_om_parameters(int fragment_pdg, int two_j, int l,
        int two_s, WorkingState& state) const;

      // Compute the optical model potential at radius r
      std::complex<double> omp(double r, const WorkingState& state) const;

      // Computes the optical model potential minus the Coulomb potential at
      // radius r
      std::complex<double> omp_minus_Vc(double r,
        const WorkingState& state) const;

      // Woods-Saxon shape
      double f(double r, double R, double a) const;

      // Partial derivative with respect to r of the Woods-Saxon shape
      double dfdr(double r, double R, double a) const;

      // Coulomb potential for a point particle with charge q*e interacting
      // with a uniformly charged sphere with radius R and charge Q*e
      double Vc(double r, double R, int Q, int q) const;

      // Non-derivative radial Schrödinger equation terms to use for computing
      // transmission coefficients via the Numerov method
      std::complex<double> a(double r, int l,
        const WorkingState& state) const;

      // Version of Schrodinger equation terms with the optical model potential
      // U pre-computed
      std::complex<double> a(double r, int l, std::complex<double> U,
        const WorkingState& state) const;

      // Neutron parameters
      double v1n, v2n, v3n, v4n, w1n, w2n, d1n, d2n, d3n, vso1n, vso2n;
      double wso1n, wso2n, Efn, Rvn, avn, Rdn, adn, Rso_n, aso_n;
      // Proton parameters
      double v1p, v2p, v3p, v4p, w1p, w2p, d1p, d2p, d3p, vso1p, vso2p;
      double wso1p, wso2p, Efp, Vcbar_p, Rvp, avp, Rdp, adp, Rso_p, aso_p;
      // Radius for nuclear Coulomb potential
      double Rc;

      // Mass of a charged pion
      static constexpr double mpiplus = 139.57018; // MeV
      // Squared pion Compton wavelength
      static constexpr double lambda_piplus2 = (marley_utils::hbar_c
        / mpiplus) * (marley_utils::hbar_c / mpiplus); // fm

      // Threshold for abs(U - Vc) used to find a suitable matching radius for
      // computing transmission coefficients.
      static constexpr double MATCHING_RADIUS_THRESHOLD = 1e-3;

      /// @brief Step size (fm) for integration of the Schr&ouml;dinger
      /// equation using the Numerov method
      double step_size_ = DEFAULT_NUMEROV_STEP_SIZE_;

      /// @brief Default step size (fm) for computing transmission coefficients
      /// via the
      /// <a href="https://en.wikipedia.org/wiki/Numerov%27s_method">Numerov
      /// method</a>
      /// @todo Make this a user-controlled value specified in the
      /// configuration file
      static constexpr double DEFAULT_NUMEROV_STEP_SIZE_ = 0.1;

      // More helper functions
      WorkingState make_working_state(double KE_tot_CM, int fragment_pdg,
        double target_mass) const;
      void calculate_kinematic_variables(double KE_tot_CM, int fragment_pdg,
        WorkingState& state) const;
      double target_mass_for_charge(int target_charge) const;
  };

}
