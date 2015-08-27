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
