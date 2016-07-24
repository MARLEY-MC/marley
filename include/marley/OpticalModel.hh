#pragma once
#include <complex>

namespace marley {

  /// @brief Abstract base class for nuclear optical model implementations
  class OpticalModel {

    public:

      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      OpticalModel(int Z, int A) : Z_(Z), A_(A) {}

      /// @brief Calculate the optical model potential (including the Coulomb
      /// potential)
      /// @param r Distance from nuclear center (fm)
      /// @param E Fragment kinetic energy (MeV)
      /// @param fragment_pdg PDG code for the nuclear fragment in the
      /// potential
      /// @param two_j Two times the total angular momentum of the fragment
      /// @param l Orbital angular momentum of the fragment
      /// @param two_s Two times the spin of the fragment
      /// @returns Complex-valued optical model potential (MeV)
      virtual std::complex<double> optical_model_potential(double r, double E,
        int fragment_pdg, int two_j, int l, int two_s) = 0;

      /// @brief Calculate the transmission coefficient for a nuclear fragment
      virtual double transmission_coefficient(double E, int fragment_pdg,
        int two_j, int l, int two_s) = 0;

      /// @brief Calculate the total cross section for a nuclear fragment
      virtual double total_cross_section(double E, int fragment_pdg, int two_s,
        size_t l_max) = 0;

      /// @brief Get the atomic number
      inline int Z() const;

      /// @brief Get the mass number
      inline int A() const;

    protected:

      // Nuclear atomic and mass numbers
      int Z_, A_;
  };

  // Inline function definitions
  inline int OpticalModel::Z() const { return Z_; }

  inline int OpticalModel::A() const { return A_; }
}
