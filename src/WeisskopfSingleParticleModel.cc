/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/COPYING or
// visit http://opensource.org/licenses/GPL-3.0

#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Logger.hh"
#include "marley/MassTable.hh"
#include "marley/WeisskopfSingleParticleModel.hh"

using TrType = marley::GammaStrengthFunctionModel::TransitionType;

marley::WeisskopfSingleParticleModel::WeisskopfSingleParticleModel(int Z,
  int A, double D0) : marley::GammaStrengthFunctionModel(Z, A), D0_(D0)
{
  if (D0_ <= 0.) throw marley::Error("Invalid"
    " level spacing parameter " + std::to_string(D0_)
    + " MeV passed to the constructor of"
    " marley::WeisskopfGammaStrengthFunctionModel.");
}

// Computes the partial decay width for a gamma transition under the Weisskopf
// single-particle approximation.
double marley::WeisskopfSingleParticleModel::partial_decay_width(TrType type,
  int l, double e_gamma)
{
  return D0_ * std::pow(e_gamma, 2*l + 1)
    * strength_function(type, l, e_gamma);
}

double marley::WeisskopfSingleParticleModel::strength_function(TrType type,
  int l, double /*e_gamma unused*/)
{
  check_multipolarity(l);

  // Compute double factorial of 2l + 1
  int dfact = 1;
  for (int n = 2*l + 1; n > 0; n -= 2) dfact *= n;

  // Multipolarity factor (dimensionless)
  double lambda = (l + 1.) / (l * std::pow(dfact, 2))
    * std::pow(3.0 / (l + 3.), 2);

  // Estimated nuclear radius (fm)
  double R = marley_utils::r0 * std::pow(A_, 1.0/3.0);

  // Electric transition strength function (MeV^[-2l-1])
  double el_sf = 2 * marley_utils::alpha * lambda
    * std::pow(R / marley_utils::hbar_c, 2*l) / D0_;

  if (type == TrType::electric) {
    return el_sf;
  }

  else if (type == TrType::magnetic) {
    static const double mp = marley::MassTable::Instance().get_particle_mass(
      marley_utils::PROTON);
    return 10. * el_sf * std::pow(marley_utils::hbar_c / (mp * R), 2);
  }

  else if (type == TrType::unphysical) {
    MARLEY_LOG_WARNING() << "Unphysical EM transition encountered in"
      << " WeisskopfSingleParticleModel::strength_function()."
      << " The strength function will be set to zero.";
    return 0.;
  }

  // @todo Improve error message
  else throw marley::Error( "Invalid transition type"
    " given for Weisskopf gamma-ray strength function calculation" );
}

double marley::WeisskopfSingleParticleModel::transmission_coefficient(
  TrType type, int l, double e_gamma)
{
  return 2. * marley_utils::pi * strength_function(type, l, e_gamma)
    * std::pow(e_gamma, 2*l + 1);
}
