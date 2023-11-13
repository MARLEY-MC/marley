/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
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

// Functions to calculate the Coulomb wavefunctions using the GNU Scientific
// Library
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_coulomb.h"

#include "marley/coulomb_wavefunctions.hh"
#include "marley/Logger.hh"

std::complex< double > coulomb_H_plus( int l, double eta, double rho ) {

  // Temporarily use the MARLEY gsl error handler (designed to ignore overflows
  // signaled by GSL_EOVRFLW, which are sometimes expected when computing the
  // Coulomb wavefunctions). Keep a pointer to the previous error handler to
  // restore it later in this function.
  auto* previous_handler = gsl_set_error_handler( &marley_gsl_error_handler );

  // H+ = G + i F
  gsl_sf_result F, Fp, G, Gp;
  double exp_F, exp_G;
  int code = gsl_sf_coulomb_wave_FG_e( eta, rho, static_cast<double>(l), 0,
    &F, &Fp, &G, &Gp, &exp_F, &exp_G );

  // Compute the value of H+, correcting for an overflow condition if one was
  // signaled by GSL
  std::complex< double > result;

  if ( code == GSL_EOVRFLW ) {
    double Fl = F.val * std::exp( exp_F );
    double Gl = G.val * std::exp( exp_G );
    result = std::complex<double>( Gl, Fl );
  }
  else {
    result = std::complex<double>( G.val, F.val );
  }

  // Restore the old GSL error handler to play nicely with other packages that
  // may manipulate it
  gsl_set_error_handler( previous_handler );

  return result;
}

/// Custom GSL error handler used by MARLEY when computing the Coulomb
/// wavefunctions
void marley_gsl_error_handler( const char* reason, const char* file, int line,
  int gsl_errno )
{
  // Overflows are expected sometimes, and they are handled explicitly by
  // coulomb_H_plus()
  if ( gsl_errno == GSL_EOVRFLW ) return;

  MARLEY_LOG_ERROR() << "GSL error: " << reason << " in file " << file
    << " at line " << line << " with error code " << gsl_errno;
}
