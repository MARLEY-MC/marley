/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
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

// MARLEY includes
#include "marley/NuclearFormFactor.hh"

double marley::HelmNuclearFormFactor::F( double kappa ) const {
  if ( kappa == 0. ) return 1.;
  double kappa_in_inverse_fm = kappa / marley_utils::hbar_c;
  double x = kappa_in_inverse_fm * R_;
  double j1 = ( std::sin( x ) / x - std::cos( x ) ) / x;
  double F = 3. * j1 / x * std::exp( -kappa_in_inverse_fm
    * kappa_in_inverse_fm * s_ * s_ / 2. );
  return F;
}

double marley::KleinNystrandNuclearFormFactor::F( double kappa ) const {
  if ( kappa == 0. ) return 1.;
  double kappa_in_inverse_fm = kappa / marley_utils::hbar_c;
  double x = kappa_in_inverse_fm * R_;
  double j1 = ( std::sin( x ) / x - std::cos( x ) ) / x;
  double F = 3. * j1 / x * ( 1. / ( 1.
    + kappa_in_inverse_fm*kappa_in_inverse_fm*a_*a_ ) );
  return F;
}
