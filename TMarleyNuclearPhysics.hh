#pragma once
#include <cmath>
#include <unordered_map>

#include "TMarleyBackshiftedFermiGasModel.hh"
#include "TMarleyFragment.hh"
#include "TMarleyHFTable.hh"
#include "TMarleyLevel.hh"
#include "TMarleyMassTable.hh"
#include "TMarleyParity.hh"
#include "TMarleySphericalOpticalModel.hh"
#include "TMarleyStructureDatabase.hh"
#include "marley_utils.hh"

class TMarleyNuclearPhysics {
  public:
    enum class TransitionType { electric, magnetic };

    static TransitionType determine_gamma_transition_type(int twoJi,
      TMarleyParity Pi, int twoJf, TMarleyParity Pf, int& l);

    static TransitionType determine_gamma_transition_type(int twoJi,
      TMarleyParity Pi, TMarleyLevel* level_f, int& l);

    static TransitionType determine_gamma_transition_type(TMarleyLevel* level_i,
      TMarleyLevel* level_f, int& l);

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
    static inline double recoil_ke(int Zi, int Ai, double Exi, int Zf, int Af,
      double Exf)
    {
      int z = Zi - Zf;
      int a = Ai - Af;
      double Mi = TMarleyMassTable::get_atomic_mass(Zi, Ai) + Exi;
      double Mf = TMarleyMassTable::get_atomic_mass(Zf, Af) + Exf
        + z*TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);
      double m = 0;
      // If there is no change in the nuclear mass number, assume that
      // a photon is being emitted (and therefore the fragment mass
      // m = 0)
      if (a != 0) {
        int fragment_pid;
        if (a == 1 && z == 1) fragment_pid = marley_utils::PROTON;
        else fragment_pid = marley_utils::get_nucleus_pid(z, a);
        m = TMarleyMassTable::get_particle_mass(fragment_pid);
      }
      return recoil_ke(Mi, Mf, m);
    }

    // Computes the partial width (multiplied by 2*pi and the initial nuclear level
    // density) predicted by the Hauser-Feshbach statistical model for decay of a
    // given nuclear level via emission of the fragment f.
    static double hf_fragment_partial_width(int Zi, int Ai, double Ex,
      int twoJi, TMarleyParity Pi, const TMarleyFragment& f,
      /*const*/ TMarleySphericalOpticalModel& om, const TMarleyDecayScheme& ds);

    static double hf_gamma_partial_width(double Ex, int twoJi,
      TMarleyParity Pi, const TMarleyDecayScheme& ds);

    // Table of nuclear fragments that will be considered when computing
    // branching ratios for nuclear de-excitations.
    inline static const std::vector<TMarleyFragment>& get_fragments() {
      return fragments;
    }

    static bool hauser_feshbach_decay(int Zi, int Ai,
      const TMarleyAtom& initial_particle, TMarleyParticle& first_product,
      TMarleyAtom& second_product, double& Ex, int& twoJ, TMarleyParity& Pi,
      TMarleyStructureDatabase& db, TMarleyGenerator& gen);

    static void hf_test(int Zi, int Ai,
      const TMarleyParticle& initial_particle, double Ex,
      int twoJ, TMarleyParity Pi,
      TMarleyStructureDatabase& db, TMarleyGenerator& gen);

    static void hf_test2(int Zi, int Ai,
      const TMarleyParticle& initial_particle, double Ex,
      int twoJ, TMarleyParity Pi,
      TMarleyStructureDatabase& db, TMarleyGenerator& gen,
      size_t num_trials, std::vector<double>& KEs,
      std::vector<int>& f_pids);

    static void hf_test3(int Zi, int Ai, const TMarleyParticle& initial_particle,
        double Ex, int twoJ, TMarleyParity Pi, TMarleyStructureDatabase& db,
        std::unordered_map<const TMarleyFragment*,
          std::function<double(double)> >& funcs,
        std::unordered_map<const TMarleyFragment*, double>& total_c_widths,
        std::unordered_map<const TMarleyFragment*, double>& E_c_mins,
        std::unordered_map<const TMarleyFragment*, double>& Exf_maxes);

    // Based on some initial parameters (including the nuclear 2J and parity)
    // determine the final nuclear spin and parity after a Hauser-Feshbach
    // decay and store them in twoJ and Pi
    static void sample_gamma_spin_parity(int Z, int A, int& twoJ,
      TMarleyParity& Pi, double Exi, double Exf, TMarleyGenerator& gen);
    static void sample_fragment_spin_parity(int& twoJ, TMarleyParity& Pi,
      const TMarleyFragment& f, /*const*/ TMarleySphericalOpticalModel& om,
      TMarleyGenerator& gen, double Exf, double Ea);

    // Create a table that can be used to sample Hauser-Feshbach decay events
    static TMarleyHFTable create_hf_table(int Zi, int Ai,
      const TMarleyParticle& initial_particle,
      double Ex, int twoJ, TMarleyParity Pi, TMarleyStructureDatabase& db,
      TMarleyGenerator& gen);

    static double get_fragment_emission_threshold(const int Zi, const int Ai,
      const TMarleyFragment& f);

  private:
    // Mass of a charged pion
    static constexpr double mpiplus = 139.57018; // MeV
    // Squared pion Compton wavelength
    static constexpr double lambda_piplus2 = std::pow(marley_utils::hbar_c
      / mpiplus, 2); // fm

    // Table of nuclear fragments that will be considered when computing
    // branching ratios for nuclear de-excitations.
    static const std::vector<TMarleyFragment> fragments;

    // Table for looking up optical model objects (used for Hauser-Feshbach
    // fragment evaporation calculations) by PDG particle ID
    static std::unordered_map<int, TMarleySphericalOpticalModel>
      pid_optical_model_table;

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
    static double fragment_continuum_partial_width(double Mconst,
      double Mfgs_ion, double Mi, int twoJi, TMarleyParity Pi, int fragment_pid,
      int two_s, TMarleyParity Pa, /*const*/ TMarleySphericalOpticalModel& om,
      double Exf);

    static double fragment_discrete_partial_width(double Exf_max, double Mconst,
      double Mfgs_ion, double Mi, int twoJi, TMarleyParity Pi, int fragment_pid,
      int two_s, TMarleyParity Pa, /*const*/ TMarleySphericalOpticalModel& om,
      const std::vector<TMarleyLevel*>* sorted_lps);

    static double gamma_continuum_partial_width(int Z, int A, int twoJi,
      double Exi, double Exf);

    static inline double gamma_cpw(int Z, int A, int mpol, int twoJf,
      double e_gamma, double Exf)
    {
      double tcE = gamma_transmission_coefficient(Z, A,
        TransitionType::electric, mpol, e_gamma);
      double tcM = gamma_transmission_coefficient(Z, A,
        TransitionType::magnetic, mpol, e_gamma);
      // Note that since our level density model assumes equipartition of parity,
      // we can multiply both transition types by the same level density function.
      double rho = TMarleyBackshiftedFermiGasModel::level_density(Z, A, Exf,
        twoJf);
      return (tcE + tcM) * rho;
    }

    // Helper function used when sampling continuum spin-parities for gamma-ray
    // transitions
    static inline double store_gamma_pws(int Z, int A, double Exf, int twoJf,
      TMarleyParity Pi, std::vector<double>& widths, std::vector<int>& twoJfs,
      std::vector<TMarleyParity>& Pfs, double tcE, double tcM, int mpol)
    {
      // Note that since our level density model assumes equipartition of parity,
      // we can multiply both transition types by the same level density function.
      double rho = TMarleyBackshiftedFermiGasModel::level_density(Z, A, Exf,
        twoJf);

      TMarleyParity Pf;
      double combined_width = 0.;

      // Compute and store information for the electric transition
      double width = tcE * rho;
      combined_width += width;
      twoJfs.push_back(twoJf);
      // Electric transitions represent a parity flip for odd multipolarities
      if (mpol % 2) Pf = -Pi;
      else Pf = Pi;
      Pfs.push_back(Pf);
      widths.push_back(width);

      // Compute and store information for the magnetic transition
      width = tcM * rho;
      combined_width += width;
      twoJfs.push_back(twoJf);
      // Magnetic transitions have opposite final parity from electric
      // transitions of the same multipolarity
      !Pf;
      Pfs.push_back(Pf);
      widths.push_back(width);

      return combined_width;
    }
};
