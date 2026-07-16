/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
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

#include "marley/CachedOpticalModel.hh"
#include "marley/LinearInterpolatingFunction.hh"

// Default to using the cache unless it is explicitly disabled
bool marley::CachedOpticalModel::USE_CACHE = true;

double marley::CachedOpticalModel::transmission_coefficient(
  double total_KE_CM, int fragment_pdg, int two_j, int l, int two_s,
  int target_charge )
{
  // If we're in the region of kinetic energies where we use the cache, then
  // check whether we've encountered this combination of quantum numbers before
  if ( USE_CACHE && total_KE_CM >= MIN_TOTAL_KE_CM
    && total_KE_CM <= MAX_TOTAL_KE_CM )
  {
    TCKey temp_key( fragment_pdg, two_j, l, two_s, target_charge );
    const auto cache_iter = tc_cache_.find( temp_key );

    // If we can find a match in the cache, evaluate the transmission
    // coefficient value using the stored interpolating function
    if ( cache_iter != tc_cache_.cend() ) {
      return cache_iter->second->evaluate( total_KE_CM );
    }

    // Otherwise, create an interpolating function object over the caching
    // kinetic energy range
    std::function< double(double) > t_coeff = [ this, fragment_pdg, two_j,
      l, two_s, target_charge ]( double tot_KE_CM ) -> double
    {
      return this->compute_transmission_coefficient( tot_KE_CM, fragment_pdg,
        two_j, l, two_s, target_charge );
    };

    auto t_coeff_interp = std::make_shared<
      marley::LinearInterpolatingFunction >( t_coeff, MIN_TOTAL_KE_CM,
        MAX_TOTAL_KE_CM, marley::DEFAULT_N_LINEAR );

    // Compute the interpolated transmission coefficient for the current
    // call to this function
    double tc_result = t_coeff_interp->evaluate( total_KE_CM );

    // Add the completed interpolating function to the cache before returning
    // the current result
    tc_cache_[ temp_key ] = std::move( t_coeff_interp );
    return tc_result;
  }

  // If we're outside the caching energy range, just compute the transmission
  // coefficient without relying on the cache at all
  return this->compute_transmission_coefficient( total_KE_CM, fragment_pdg,
    two_j, l, two_s, target_charge );
}
