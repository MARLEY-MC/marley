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

// Standard library includes

// MARLEY includes
#include "marley/FormFactors.hh"
#include "marley/marley_utils.hh"
#include "marley/Error.hh"

// The default constructor sets the Q^2 scaling mode
marley::FormFactors::FormFactors( FFScalingMode mode ) : ff_scaling_mode_(mode) {}

// Helper map for converting strings to Q^2 scaling modes
std::map < marley::FormFactors::FFScalingMode, std::string >
  marley::FormFactors::ff_scaling_mode_map_ = {
  { marley::FormFactors::FFScalingMode::DIPOLE, "dipole" },
  { marley::FormFactors::FFScalingMode::FLAT, "flat" }
};

// Convert a string to a Q^2 scaling mode
marley::FormFactors::FFScalingMode marley::FormFactors
  ::ff_scaling_mode_from_string( const std::string& str )
{
  for (const auto& pair : ff_scaling_mode_map_) {
    if (pair.second == str) return pair.first;
  }
  throw marley::Error( "Invalid Q^2 scaling mode string: " + str );
}

// Convert a Q^2 scaling mode to a string
std::string marley::FormFactors::string_from_ff_scaling_mode(
  FFScalingMode mode )
{
  auto it = ff_scaling_mode_map_.find( mode );
  if ( it != ff_scaling_mode_map_.end() ) return it->second;
  else throw marley::Error( "Invalid FFScalingMode value encountered in"
    " marley::FormFactors::string_from_q2_scaling_mode" );
}

// Dipole dependence form factor
double marley::FormFactors::dipole( double Q2, double M ) const {
  return 1.0 / ( ( 1.0 + Q2 / M / M ) * ( 1.0 + Q2 / M / M ) );
}

// F1_n form factor
double marley::FormFactors::F1_n( double Q2, double M ) const {
  if ( ff_scaling_mode_ == FFScalingMode::FLAT ) { Q2 = 0.; }
  else if ( ff_scaling_mode_ != FFScalingMode::DIPOLE ) {
    throw marley::Error("Invalid Q^2 scaling mode.");
  }

  double prefactor = ( marley_utils::mu_n * Q2 )
    / ( 4. * marley_utils::m_nucleon2 + Q2 );

  return prefactor * dipole( Q2, M );
}

// F1_p form factor
double marley::FormFactors::F1_p( double Q2, double M ) const {
  if ( ff_scaling_mode_ == FFScalingMode::FLAT ) { Q2 = 0.; }
  else if ( ff_scaling_mode_ != FFScalingMode::DIPOLE ) {
    throw marley::Error( "Invalid Q^2 scaling mode." );
  }

  double prefactor = ( 4. * marley_utils::m_nucleon2
    + marley_utils::mu_p * Q2 ) / (4. * marley_utils::m_nucleon2 + Q2 );

  return prefactor * dipole( Q2, M );
}

// F2_n form factor
double marley::FormFactors::F2_n( double Q2, double M ) const {
  if ( ff_scaling_mode_ == FFScalingMode::FLAT ) { Q2 = 0.; }
  else if ( ff_scaling_mode_ != FFScalingMode::DIPOLE ) {
    throw marley::Error( "Invalid Q^2 scaling mode." );
  }

  double prefactor = ( 4. * marley_utils::m_nucleon2 * marley_utils::mu_n )
    / ( 4. * marley_utils::m_nucleon2 + Q2 );

  prefactor /= 2. * marley_utils::m_nucleon;

  return prefactor * dipole( Q2, M );
}

// F2_p form factor
double marley::FormFactors::F2_p( double Q2, double M ) const {
  if ( ff_scaling_mode_ == FFScalingMode::FLAT ) { Q2 = 0.; }
  else if ( ff_scaling_mode_ != FFScalingMode::DIPOLE ) {
    throw marley::Error( "Invalid Q^2 scaling mode." );
  }

  double prefactor = ( 4. * marley_utils::m_nucleon2
    * (marley_utils::mu_p - 1) ) / ( 4. * marley_utils::m_nucleon2 + Q2 );

  prefactor /= 2. * marley_utils::m_nucleon;

  return prefactor * dipole( Q2, M );
}

// FS form factor
double marley::FormFactors::FS( double Q2, double M ) const { return 0.0; }

// FA form factor
double marley::FormFactors::FA( double Q2, double M ) const {
  if ( ff_scaling_mode_ == FFScalingMode::FLAT ) { Q2 = 0.; }
  else if ( ff_scaling_mode_ != FFScalingMode::DIPOLE ) {
    throw marley::Error( "Invalid Q^2 scaling mode." );
  }

  return -1. * marley_utils::g_A * dipole( Q2, M );
}

// FP form factor
double marley::FormFactors::FP( double Q2, double M ) const {
  if (ff_scaling_mode_ == FFScalingMode::FLAT) { Q2 = 0.; }
  else if (ff_scaling_mode_ != FFScalingMode::DIPOLE) {
    throw marley::Error( "Invalid Q^2 scaling mode." );
  }

  double prefactor = ( 2. * marley_utils::m_nucleon )
    / ( marley_utils::m_pion * marley_utils::m_pion + Q2 );

  return prefactor * this->FA( Q2, M );
}

// FT form factor
double marley::FormFactors::FT( double Q2, double M ) const { return 0.0; }
