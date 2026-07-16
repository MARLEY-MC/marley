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

#pragma once

#include <complex>
#include <iostream>
#include <map>
#include <memory>

#include "marley/OpticalModel.hh"

namespace marley {

  // Forward-declare the InterpolatingFunction class (used in the cache)
  class InterpolatingFunction;

  /// @brief Nuclear optical model that caches transmission coefficient
  /// calculations for efficiency
  class CachedOpticalModel : public OpticalModel {

    public:

      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      CachedOpticalModel( int Z, int A ) : OpticalModel( Z, A ) {}

      virtual ~CachedOpticalModel() = default;

      virtual double transmission_coefficient( double total_KE_CM,
        int fragment_pdg, int two_j, int l, int two_s,
        int target_charge = 0 ) override;

      struct TCKey {

        TCKey() {}
        TCKey( int fragment_pdg, int two_j, int l, int two_s,
          int target_charge ) : frag_pdg_( fragment_pdg ), two_j_( two_j ),
          l_( l ), two_s_( two_s ), target_charge_( target_charge ) {}

        // Defines an ordering allowing this struct to be used as the key
        // for a std::map
        bool operator<( const TCKey& other ) const {
          if ( frag_pdg_ != other.frag_pdg_ ) {
            return frag_pdg_ < other.frag_pdg_;
          }
          else if ( two_j_ != other.two_j_ ) {
            return two_j_ < other.two_j_;
          }
          else if ( l_ != other.l_ ) {
            return l_ < other.l_;
          }
          else if ( two_s_ != other.two_s_ ) {
            return two_s_ < other.two_s_;
          }
          else if ( target_charge_ != other.target_charge_ ) {
            return target_charge_ < other.target_charge_;
          }
          return false;
        }

        int frag_pdg_ = 0;
        int two_j_ = 0;
        int l_ = 0;
        int two_s_ = 0;
        int target_charge_ = 0;
      };

      inline static void set_use_cache( bool use_it ) {
        USE_CACHE = use_it;
      }

    protected:

      /// @brief Do the actual calculation of the transmission coefficient
      virtual double compute_transmission_coefficient( double total_KE_CM,
        int fragment_pdg, int two_j, int l, int two_s,
        int target_charge = 0 ) = 0;

      /// @brief Saved InterpolatingFunction objects corresponding
      /// to previously-encountered transmission coefficient requests
      std::map< TCKey, std::shared_ptr< InterpolatingFunction > >
        tc_cache_;

      /// @brief Minimum kinetic energy to use when interacting with the cache
      static constexpr double MIN_TOTAL_KE_CM = 1e-10; // MeV
      /// @brief Maximum kinetic energy to use when interacting with the cache
      static constexpr double MAX_TOTAL_KE_CM = 100.; // MeV

      /// @brief Boolean switch that allows global enabling/disabling of the
      /// cache, which is used by default
      static bool USE_CACHE;
  };

}
