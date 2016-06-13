#pragma once
#include <cmath>
#include <unordered_map>

#include "Fragment.hh"
#include "Level.hh"
#include "LevelDensityModel.hh"
#include "MassTable.hh"
#include "Parity.hh"
#include "SphericalOpticalModel.hh"
#include "StructureDatabase.hh"
#include "marley_utils.hh"

namespace marley {

  class NuclearPhysics {
    public:
      enum class TransitionType { electric, magnetic };

      static TransitionType determine_gamma_transition_type(int twoJi,
        marley::Parity Pi, int twoJf, marley::Parity Pf, int& l);

      static TransitionType determine_gamma_transition_type(int twoJi,
        marley::Parity Pi, marley::Level* level_f, int& l);

      static TransitionType determine_gamma_transition_type(marley::Level* level_i,
        marley::Level* level_f, int& l);

      static double weisskopf_partial_decay_width(int A, TransitionType type,
        int l, double e_gamma);

      static double gamma_strength_function_coefficient(int Z, int A,
        TransitionType type, int l, double e_gamma);

      static inline double gamma_strength_function(int Z, int A,
        TransitionType type, int l, double e_gamma)
      {
        return std::pow(e_gamma, 3 - 2*l) * gamma_strength_function_coefficient(Z,
          A, type, l, e_gamma);
      }


      static inline double gamma_transmission_coefficient(int Z, int A,
        TransitionType type, int l, double e_gamma)
      {
        return 2 * marley_utils::pi * gamma_strength_function_coefficient(Z, A,
          type, l, e_gamma) * std::pow(e_gamma, 4); // Eg^4 = (Eg^[2l + 1] * Eg^[3 - 2l]
      }

      // Approximates the electric potential energy (in MeV) of two nuclei (with
      // atomic numbers Z, z and mass numbers A, a) that are just touching by
      // modeling them as uniformly charged spheres with radii given by
      // R = (1.2 fm) * A^(1/3).
      static inline double coulomb_barrier(int Z, int A, int z, int a) {
        return Z * z * marley_utils::e2 / (marley_utils::r0 * (std::pow(A, 1.0/3.0)
          + std::pow(a, 1.0/3.0)));
      }

      // Computes the kinetic energy (in MeV and in the initial nucleus's rest
      // frame) of the recoiling residual ion (with mass Mf = Mgsf + Exf +
      // z*melectron, where Mgsf is its ground-state mass and Exf is its
      // excitation energy) from a decay in which the original nucleus (with mass
      // Mi = Mgsi + Exi) emits a fragment with mass m and atomic number z.
      static inline double recoil_ke(double Mi, double Mf, double m)
      {
        return (std::pow(Mf - Mi, 2) - std::pow(m, 2)) / (2 * Mi);
      }

      // Computes the recoil kinetic energy (in MeV and in the initial nucleus's
      // rest frame) of a residual negative ion with atomic number Zf, mass
      // number Af, and excitation energy Exf that was created when an initial
      // atom (Zi, Ai, Exi) emitted a nuclear fragment (z = Zi - Zf, a = Ai - Af)
      // which is assumed to have been emitted in its ground state.
      static double recoil_ke(int Zi, int Ai, double Exi, int Zf, int Af,
        double Exf);

      // Computes the partial width (multiplied by 2*pi and the initial nuclear level
      // density) predicted by the Hauser-Feshbach statistical model for decay of a
      // given nuclear level via emission of the fragment f.
      static double hf_fragment_partial_width(int Zi, int Ai, double Ex,
        int twoJi, marley::Parity Pi, const marley::Fragment& f,
        marley::SphericalOpticalModel& om, const marley::DecayScheme& ds,
        marley::LevelDensityModel& ldm);

      static double hf_gamma_partial_width(double Ex, int twoJi,
        marley::Parity Pi, const marley::DecayScheme& ds,
        marley::LevelDensityModel& ldm);

      // Table of nuclear fragments that will be considered when computing
      // branching ratios for nuclear de-excitations.
      inline static const std::vector<marley::Fragment>& get_fragments() {
        return fragments;
      }

      static bool hauser_feshbach_decay(int Zi, int Ai,
        const marley::Particle& initial_particle, marley::Particle& first_product,
        marley::Particle& second_product, double& Ex, int& twoJ, marley::Parity& Pi,
        marley::StructureDatabase& db, marley::Generator& gen);

      // Based on some initial parameters (including the nuclear 2J and parity)
      // determine the final nuclear spin and parity after a Hauser-Feshbach
      // decay and store them in twoJ and Pi
      static void sample_gamma_spin_parity(int Z, int A, int& twoJ,
        marley::Parity& Pi, double Exi, double Exf, marley::Generator& gen);
      static void sample_fragment_spin_parity(int& twoJ, marley::Parity& Pi,
        const marley::Fragment& f, marley::SphericalOpticalModel& om,
        marley::Generator& gen, double Exf, double Ea);

      static double get_fragment_emission_threshold(const int Zi, const int Ai,
        const marley::Fragment& f);

    private:

      // Mass of a charged pion
      static constexpr double mpiplus = 139.57018; // MeV
      // Squared pion Compton wavelength
      static constexpr double lambda_piplus2 = std::pow(marley_utils::hbar_c
        / mpiplus, 2); // fm

      // Table of nuclear fragments that will be considered when computing
      // branching ratios for nuclear de-excitations.
      static const std::vector<marley::Fragment> fragments;

      // Maximum value of the orbital angular momentum to use in Hauser-Feshbach
      // calculations
      // TODO: make this a user-controlled value specified in the configuration file
      static constexpr int l_max = 2;

      // Default step size for computing optical model transmission coefficients via
      // the Numerov method
      // TODO: make this a user-controlled value specified in the configuration file
      static constexpr double DEFAULT_NUMEROV_STEP_SIZE = 0.1; // fm

      // Default number of subintervals to use when integrating fragment and gamma
      // decay widths over the energy continuum.
      // TODO: remove this when you don't need to use it in old testing code anymore
      static constexpr int DEFAULT_CONTINUUM_SUBINTERVALS = 50;

      // Helper functions used internally while computing decay widths using
      // the Hauser-Feshbach statistical model
      static double fragment_continuum_partial_width(int fragment_pid, double Ea,
        int two_s, marley::Parity Pa, int twoJi, marley::Parity Pi,
        marley::SphericalOpticalModel& om, marley::LevelDensityModel& ldm,
        double Exf);

      static double fragment_discrete_partial_width(double Exf_max, double Mconst,
        double Mfgs_ion, double Mi, int twoJi, marley::Parity Pi, int fragment_pid,
        int two_s, marley::Parity Pa, marley::SphericalOpticalModel& om,
        const std::vector<marley::Level*>* sorted_lps);

      static double gamma_continuum_partial_width(int Z, int A, int twoJi,
        marley::Parity Pi, double Exi, marley::LevelDensityModel& ldm,
        double Exf);


      static double gamma_cpw(int Z, int A, int mpol, int twoJf,
        double e_gamma, marley::Parity Pi, marley::LevelDensityModel& ldm,
        double Exf);

      // Helper function used when sampling continuum spin-parities for gamma-ray
      // transitions
      static double store_gamma_pws(double Exf, int twoJf,
        marley::Parity Pi, std::vector<double>& widths, std::vector<int>& twoJfs,
        std::vector<marley::Parity>& Pfs, double tcE, double tcM, int mpol,
        marley::LevelDensityModel& ldm);
  };

}
