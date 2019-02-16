#pragma once
#include <complex>

namespace marley {

  /// @brief Abstract base class for nuclear optical model implementations
  class OpticalModel {

    public:

      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      OpticalModel(int Z, int A) : Z_(Z), A_(A) {}

      virtual ~OpticalModel() = default;

      /// @brief Calculate the optical model potential (including the Coulomb
      /// potential)
      /// @param r Distance from nuclear center (fm)
      /// @param fragment_KE_lab Fragment kinetic energy (MeV) in the lab frame
      /// @param fragment_pdg PDG code for the nuclear fragment
      /// @param two_j Two times the total angular momentum of the fragment
      /// @param l Orbital angular momentum of the fragment
      /// @param two_s Two times the spin of the fragment
      /// @param target_charge Net charge of the target atom
      /// @returns Complex-valued optical model potential (MeV)
      virtual std::complex<double> optical_model_potential(double r,
        double fragment_KE_lab, int fragment_pdg, int two_j, int l, int two_s,
        int target_charge = 0) = 0;

      /// @brief Calculate the transmission coefficient for a nuclear fragment
      /// @param total_KE_CM Total CM frame kinetic energy (MeV)
      /// @param fragment_pdg PDG code of the fragment
      /// @param two_j Two times the total angular momentum of the fragment
      /// @param l Orbital angular momentum of the fragment
      /// @param two_s Two times the spin of the fragment
      /// @param target_charge Net charge of the target atom
      virtual double transmission_coefficient(double total_KE_CM, int fragment_pdg,
        int two_j, int l, int two_s, int target_charge = 0) = 0;

      /// @brief Compute the energy-averaged total cross section
      /// (MeV<sup> -2</sup>) for a nuclear fragment projectile
      /// @details The total cross section given here by the optical model may
      /// be directly compared to experiment. Expressions exist for the optical
      /// model elastic and reaction (absorption) cross sections, but these are
      /// hard to directly compare to the data because the absorption cross
      /// section includes the compound elastic channel.
      ///
      /// For more details, see appendix A (especially equation A.12) of S.
      /// Gardiner, "Nuclear Effects in Neutrino Detection," PhD thesis,
      /// University of California, Davis, 2018. Note that we assume in this
      /// function that the target nucleus has zero spin (spherical optical
      /// model).
      /// @param fragment_KE_lab Lab-frame kinetic energy of the incident projectile (MeV)
      /// @param fragment_pdg Projectile's PDG code
      /// @param two_s Two times the spin of the projectile
      /// @param l_max The maximum value of the orbital angular momentum quantum
      /// number @f$\ell@f$ to include in the sum over S-matrix elements
      /// @param target_charge Net charge of the target atom (used to adjust the
      /// atomic mass by the appropriate number of electron masses if the target
      /// is ionized)
      /// @returns Energy-averaged total scattering cross section (MeV<sup> -2</sup>)
      virtual double total_cross_section(double fragment_KE_lab,
        int fragment_pdg, int two_s, size_t l_max, int target_charge = 0) = 0;

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
