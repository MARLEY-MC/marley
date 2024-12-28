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
#include "marley/Error.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/NucleonFormFactors.hh"
#include "marley/marley_utils.hh"

double marley::DipoleSachsFormFactors::GEp( double Q2 ) {
  return marley_utils::g_V * this->dipole( Q2 );
}

double marley::DipoleSachsFormFactors::GMp( double Q2 ) {
  return marley_utils::mu_p * this->dipole( Q2 );
}

double marley::DipoleSachsFormFactors::GMn( double Q2 ) {
  return marley_utils::mu_n * this->dipole( Q2 );
}

// Implements Eq. (6) from Nucl. Phys. B Proc. Suppl. 159, 127 (2006)
// https://doi.org/10.1016/j.nuclphysbps.2006.08.028
double marley::BBBA05SachsFormFactors::bbba05_G( double tau,
  const std::vector< double >& a_coeffs,
  const std::vector< double >& b_coeffs )
{
  double numer = 0.;
  double denom = 1.;

  double tau_pow = 1.;
  for ( const auto& a : a_coeffs ) {
    numer += a * tau_pow;
    tau_pow *= tau;
  }

  tau_pow = tau;
  for ( const auto& b : b_coeffs ) {
    denom += b * tau_pow;
    tau_pow *= tau;
  }

  double ff_G = numer / denom;
  return ff_G;
}

double marley::BBBA05SachsFormFactors::GEp( double Q2 ) {
  double tau = this->tau( Q2 );
  // TODO: remove hard-coding here and make the parameters configurable
  return marley_utils::g_V * this->bbba05_G( tau, { 1., -0.0578 },
    { 11.1, 13.6, 33. } );
}

double marley::BBBA05SachsFormFactors::GMp( double Q2 ) {
  double tau = this->tau( Q2 );
  // TODO: remove hard-coding here and make the parameters configurable
  return marley_utils::mu_p * this->bbba05_G( tau, { 1., 0.150 },
    { 11.1, 19.6, 7.54 } );
}

double marley::BBBA05SachsFormFactors::GEn( double Q2 ) {
  double tau = this->tau( Q2 );
  // TODO: remove hard-coding here and make the parameters configurable
  return this->bbba05_G( tau, { 0., 1.25, 1.30 },
    { -9.86, 305., -758., 802. } );
}

double marley::BBBA05SachsFormFactors::GMn( double Q2 ) {
  double tau = this->tau( Q2 );
  // TODO: remove hard-coding here and make the parameters configurable
  return marley_utils::mu_n * this->bbba05_G( tau, { 1., 1.81 },
    { 14.1, 20.7, 68.7 } );
}

double marley::TrivialSachsFormFactors::GEp( double /*Q2*/ )
  { return marley_utils::g_V; }

double marley::TrivialSachsFormFactors::GEn( double /*Q2*/ ) { return 0.; }

double marley::TrivialSachsFormFactors::GMp( double /*Q2*/ )
  { return marley_utils::mu_p; }

double marley::TrivialSachsFormFactors::GMn( double /*Q2*/ )
  { return marley_utils::mu_n; }

double marley::TrivialAxialFormFactors::FA( double /*Q2*/ )
  { return -marley_utils::g_A; }

double marley::TrivialAxialFormFactors::FP( double /*Q2*/ ) {
  return -2. * marley_utils::g_A * marley_utils::m_nucleon
    / marley_utils::m_pion / marley_utils::m_pion;
}

double marley::NucleonFormFactors::F1p( double Q2 ) const {
  double tau = sachs_ff_->tau( Q2 );
  double GEp = sachs_ff_->GEp( Q2 );
  double GMp = sachs_ff_->GMp( Q2 );
  double F1p = ( GEp + tau*GMp ) / ( 1. + tau );
  return F1p;
}

double marley::NucleonFormFactors::F1n( double Q2 ) const {
  double tau = sachs_ff_->tau( Q2 );
  double GEn = sachs_ff_->GEn( Q2 );
  double GMn = sachs_ff_->GMn( Q2 );
  double F1n = ( GEn + tau*GMn ) / ( 1. + tau );
  return F1n;
}

double marley::NucleonFormFactors::F2p( double Q2 ) const {
  double tau = sachs_ff_->tau( Q2 );
  double GEp = sachs_ff_->GEp( Q2 );
  double GMp = sachs_ff_->GMp( Q2 );
  double F2p = ( GMp - GEp ) / ( 1. + tau ) / 2. / marley_utils::m_nucleon;
  return F2p;
}

double marley::NucleonFormFactors::F2n( double Q2 ) const {
  double tau = sachs_ff_->tau( Q2 );
  double GEn = sachs_ff_->GEn( Q2 );
  double GMn = sachs_ff_->GMn( Q2 );
  double F2n = ( GMn - GEn ) / ( 1. + tau ) / 2. / marley_utils::m_nucleon;
  return F2n;
}

marley::NucleonFormFactors::NucleonFormFactors( const marley::JSON& config ) {

  // Flag used to ensure that we send information about the form factor
  // configuration to the logger exactly once when this function is first
  // called
  static bool need_to_log = true;

  // If the user has requested use of the allowed approximation, then
  // configure trivial form factor models and return without logging this
  // choice (it will be handled elsewhere)
  if ( marley::JSONConfig::check_for_allowed_approximation(config) ) {
    sachs_ff_ = std::make_shared< TrivialSachsFormFactors >();
    axial_ff_ = std::make_shared< TrivialAxialFormFactors >();
    need_to_log = false;
    return;
  }

  if ( !config.has_key("sachs_model") ) {
    throw marley::Error( "Missing Sachs form factor model configuration" );
  }
  const auto& sachs_json = config.at( "sachs_model" );
  if ( !sachs_json.is_string() ) {
    throw marley::Error( "Invalid Sachs form factor model "
      + sachs_json.dump_string() );
  }
  std::string sachs_ff_model = sachs_json.to_string();

  if ( !config.has_key("axial_model") ) {
    throw marley::Error( "Missing axial form factor model configuration" );
  }
  const auto& axial_json = config.at( "axial_model" );
  if ( !axial_json.is_string() ) {
    throw marley::Error( "Invalid axial form factor model "
      + axial_json.dump_string() );
  }
  std::string axial_ff_model = axial_json.to_string();

  // Choose the model to use for the Sachs form factors
  if ( sachs_ff_model == "trivial" ) {
    sachs_ff_ = std::make_shared< TrivialSachsFormFactors >();
    if ( need_to_log ) {
      MARLEY_LOG_INFO() << "Using trivial Sachs form factors";
    }
  }
  else if ( sachs_ff_model == "dipole" ) {
    sachs_ff_ = std::make_shared< DipoleSachsFormFactors >( marley_utils::M_V );
    if ( need_to_log ) {
      MARLEY_LOG_INFO() << "Using dipole Sachs form factors";
    }
  }
  else if ( sachs_ff_model == "bbba05" ) {
    sachs_ff_ = std::make_shared< BBBA05SachsFormFactors >();
    if ( need_to_log ) {
      MARLEY_LOG_INFO() << "Using BBBA05 Sachs form factors";
    }
  }
  else throw marley::Error( "Unrecognized Sachs form factor model name \""
    + sachs_ff_model + "\" in constructor of marley::NucleonFormFactors" );

  // Choose the model to use for the axial form factors
  if ( axial_ff_model == "trivial" ) {
    axial_ff_ = std::make_shared< TrivialAxialFormFactors >();
    if ( need_to_log ) {
      MARLEY_LOG_INFO() << "Using trivial axial form factors";
    }
  }
  else if ( axial_ff_model == "dipole" ) {
    axial_ff_ = std::make_shared< DipoleAxialFormFactors >( marley_utils::g_A,
      marley_utils::M_A );
    if ( need_to_log ) {
      MARLEY_LOG_INFO() << "Using dipole axial form factors";
    }
  }
  else throw marley::Error( "Unrecognized axial form factor model name \""
    + axial_ff_model + "\" in constructor of marley::NucleonFormFactors" );

  need_to_log = false;
}

double marley::SachsFormFactors::tau( double Q2 ) {
  return Q2 / 4. / marley_utils::m_nucleon / marley_utils::m_nucleon;
}

double marley::DipoleAxialFormFactors::FP( double Q2 ) {
  return 2. * marley_utils::m_nucleon * this->FA( Q2 )
    / ( marley_utils::m_pion*marley_utils::m_pion + Q2 );
}
