/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see \${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#include <cmath>
#include <string>

#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Logger.hh"
#include "marley/StandardLorentzianModel.hh"

using TrType = marley::GammaStrengthFunctionModel::TransitionType;

marley::StandardLorentzianModel::StandardLorentzianModel(int Z, int A)
  : GammaStrengthFunctionModel(Z, A)
{
  // E1 giant resonance parameterization taken from the empirical fit for
  // spherical nuclei given in the
  // <a href="https://www-nds.iaea.org/RIPL-2/handbook/ripl2.pdf">
  // RIPL-2 handbook</a>, p. 129
  // @todo Consider updating E1 parameters to the new SLO fit from RIPL-3 (see
  // equation 174 in the RIPL-3 Nuclear Data Sheets paper, and compare equation
  // 173)
  e_E1_ = 31.2*std::pow(A_, -1.0/3.0) + 20.6*std::pow(A_, -1.0/6.0); // MeV
  gamma_E1_ = 0.026 * std::pow(e_E1_, 1.91); // MeV
  sigma_E1_ = 1.2 * 120 * (A_ - Z_) * Z_/(A_ * marley_utils::pi * gamma_E1_)
    * marley_utils::mb; // mb

  // E2 giant resonance parameterization taken from the global fit given by
  // Kopecky in the
  // <a href="https://www-nds.iaea.org/ripl/readme/ripl_handbook.ps">RIPL-1
  // handbook</a>, p. 103
  // RIPL-2 and RIPL-3 do not update this parameterization for the SLO
  e_E2_ = 63*std::pow(A_, -1.0/3.0); // MeV
  gamma_E2_ = 6.11 - 0.012*A_; // MeV
  sigma_E2_ = 0.00014 * std::pow(Z_, 2) * e_E2_
    / (std::pow(A_, 1.0/3.0) * gamma_E2_) * marley_utils::mb; // mb

  // M1 parameterization taken from the global SLO model fit given in the
  // <a href="https://www-nds.iaea.org/RIPL-2/handbook/ripl2.pdf">RIPL-2
  // handbook</a>, p. 132
  /// @note Be careful! You <b>must</b> initialize the E1 giant resonance
  /// parameters before the M1 parameters!
  constexpr double e_gamma_ref = 7.0; // MeV
  double factor_m1 = strength_function(TrType::electric, 1, e_gamma_ref)
    / (0.0588 * std::pow(A_, 0.878));
  gamma_M1_ = 4.0; // MeV
  e_M1_ = 41*std::pow(A_, -1.0/3.0); // MeV
  sigma_M1_ = (std::pow(std::pow(e_gamma_ref, 2) - std::pow(e_M1_, 2), 2)
    + std::pow(e_gamma_ref, 2) * std::pow(gamma_M1_, 2))
    * (3 * std::pow(marley_utils::pi, 2) * factor_m1)
    / (e_gamma_ref * std::pow(gamma_M1_, 2)); // mb
}

double marley::StandardLorentzianModel::strength_function_coefficient(
  TrType type, int l, double e_gamma)
{
  check_multipolarity(l);

  // The strength, energy, and width of the giant resonance for a transition of
  // type x (E or M) and multipolarity l
  double e_xl = 0.;
  double gamma_xl = 0.;
  double sigma_xl = 0.;

  if (type == TrType::electric) {

    if (l == 1) {
      e_xl = e_E1_;
      sigma_xl = sigma_E1_;
      gamma_xl = gamma_E1_;
    }

    else if (l > 1) {

      e_xl = e_E2_;
      sigma_xl = sigma_E2_;
      gamma_xl = gamma_E2_;

      // If this is an E2 transition, we're done. Otherwise,
      // compute the giant resonance strength iteratively
      // l > 2 prescription taken from TALYS 1.6 manual
      for (int i = 2; i < l; ++i) sigma_xl *= 8e-4;
    }
  }

  else if (type == TrType::magnetic) {
    e_xl = e_M1_;
    sigma_xl = sigma_M1_;
    gamma_xl = gamma_M1_;
    // If this is an M1 transition, we're done. Otherwise,
    // compute the giant resonance strength iteratively
    // l > 1 prescription taken from TALYS 1.6 manual
    for (int i = 1; i < l; ++i) sigma_xl *= 8e-4;
  }

  else if (type == TrType::unphysical) {
    MARLEY_LOG_WARNING() << "Unphysical EM transition encountered in"
      << " StandardLorentzianModel::strength_function_coefficient()."
      << " The strength function will be set to zero.";
    return 0.;
  }

  /// @todo Improve error message
  else throw marley::Error(std::string("Invalid transition type")
    + " given for gamma ray strength function calculation");

  // Now that we have the appropriate giant resonance parameters,
  // calculate the strength function using the Brink-Axel expression.
  // Note that the strength function has units of MeV^(-3)
  double coeff = (sigma_xl * std::pow(gamma_xl, 2)) / ((2*l + 1)
    * std::pow(marley_utils::pi, 2) * (std::pow(std::pow(e_gamma, 2)
    - std::pow(e_xl, 2), 2) + std::pow(e_gamma, 2) * std::pow(gamma_xl, 2)));

  return coeff;
}

double marley::StandardLorentzianModel::strength_function(TrType type, int l,
  double e_gamma)
{
  return std::pow(e_gamma, 3 - 2*l)
    * strength_function_coefficient(type, l, e_gamma);
}

double marley::StandardLorentzianModel::transmission_coefficient(TrType type,
  int l, double e_gamma)
{
  // Eg^4 = (Eg^[2l + 1] * Eg^[3 - 2l]
  return 2. * marley_utils::pi * strength_function_coefficient(type, l, e_gamma)
    * std::pow(e_gamma, 4);
}
