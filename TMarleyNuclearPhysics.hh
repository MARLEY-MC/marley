#pragma once
#include <cmath>

#include "TMarleyFragment.hh"
#include "TMarleyLevel.hh"
#include "TMarleyMassTable.hh"
#include "TMarleyParity.hh"
#include "marley_utils.hh"

class TMarleyNuclearPhysics {
  public:
    enum class TransitionType { electric, magnetic };

    static TransitionType determine_gamma_transition_type(int twoJi,
      TMarleyParity Pi, int twoJf, TMarleyParity Pf, int& l);

    static TransitionType determine_gamma_transition_type(TMarleyLevel* level_i,
      TMarleyLevel* level_f, int&l);

    static double weisskopf_partial_decay_width(int A, TransitionType type,
      int l, double e_gamma);

    static double gamma_strength_function(int Z, int A, TransitionType type,
      int l, double e_gamma);

    static inline double gamma_transmission_coefficient(int Z, int A,
      TransitionType type, int l, double e_gamma)
    {
      return 2 * marley_utils::pi * gamma_strength_function(Z, A, type,
        l, e_gamma) * std::pow(e_gamma, 2*l + 1);
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
    // frame) of the recoiling residual nucleus (with mass Mf = Mgsf + Exf,
    // where Mgsf is its ground-state mass and Exf is its excitation energy)
    // from a decay in which the original nucleus (with mass Mi = Mgsi + Exi)
    // emits a fragment with mass m.
    static inline double recoil_ke(double Mi, double Mf, double m)
    {
      return (std::pow(Mi, 2) - std::pow(m, 2) + std::pow(Mf, 2)) / (2 * Mi)
        - Mf;
    }

    // Computes the recoil kinetic energy (in MeV and in the initial nucleus's
    // rest frame) of a residual nucleus with atomic number Zf, mass number Af,
    // and excitation energy Exf that was created when an initial nucleus (Zi,
    // Ai, Exi) emitted a fragment (z = Zi - Zf, a = Ai - Af) which is assumed
    // to have been emitted in its ground state.
    static inline double recoil_ke(int Zi, int Ai, double Exi, int Zf, int Af,
      double Exf)
    {
      double Mi = TMarleyMassTable::get_atomic_mass(Zi, Ai) + Exi;
      double Mf = TMarleyMassTable::get_atomic_mass(Zf, Af) + Exf;
      double m = 0;
      int a = Ai - Af;
      // If there is no change in the nuclear mass number, assume that
      // a photon is being emitted (and therefore the fragment mass
      // m = 0)
      if (a != 0) m = TMarleyMassTable::get_particle_mass(
        marley_utils::get_nucleus_pid(Zi - Zf, a));

      return recoil_ke(Mi, Mf, m);
    }

  private:
    // Mass of a charged pion
    static constexpr double mpiplus = 139.57018; // MeV
    // Squared pion Compton wavelength
    static constexpr double lambda_piplus2 = std::pow(marley_utils::hbar_c
      / mpiplus, 2); // fm

    // Table of nuclear fragments that will be considered when computing
    // branching ratios for nuclear de-excitations.
    static const std::vector<TMarleyFragment> fragments;
};
