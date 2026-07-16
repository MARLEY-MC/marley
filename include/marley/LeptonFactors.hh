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

// MARLEY includes
#include "marley/NuclearResponses.hh"

namespace marley {

  class LeptonFactors {

    public:

      LeptonFactors() {}

      LeptonFactors( double vCC, double vLL, double vCL, double vT,
        double vTprime ) : vCC_( vCC ), vLL_( vLL ), vCL_( vCL ),
        vT_( vT ), vTprime_( vTprime ) {}

      inline double vCC() const { return vCC_; }
      inline double vLL() const { return vLL_; }
      inline double vCL() const { return vCL_; }
      inline double vT() const { return vT_; }
      inline double vTprime() const { return vTprime_; }

      inline void set_vCC( double val ) { vCC_ = val; }
      inline void set_vLL( double val ) { vLL_ = val; }
      inline void set_vCL( double val ) { vCL_ = val; }
      inline void set_vT( double val ) { vT_ = val; }
      inline void set_vTprime( double val ) { vTprime_ = val; }

      // Scalar product with the corresponding nuclear responses.
      // Guard each term so that a zero nuclear response yields a zero
      // contribution even when the corresponding lepton factor is NaN.
      // Note that NaN * 0 == NaN rather than 0, which would silently corrupt
      // the result.
      /// @todo Revisit this hacky fix as you do a full implementation of
      /// inelastic NC scattering at finite momentum transfer
      inline double operator*( const NuclearResponses& nr ) {
        double product = 0.;
        double rcc = nr.RCC();
        double rll = nr.RLL();
        double rcl = nr.RCL();
        double rt  = nr.RT();
        double rtp = nr.RTprime();
        if ( rcc != 0. ) product += this->vCC_ * rcc;
        if ( rll != 0. ) product += this->vLL_ * rll;
        if ( rcl != 0. ) product += this->vCL_ * rcl;
        if ( rt  != 0. ) product += this->vT_ * rt;
        if ( rtp != 0. ) product += this->vTprime_ * rtp;
        return product;
      }

    protected:
      double vCC_ = 0.;
      double vLL_ = 0.;
      double vCL_ = 0.;
      double vT_ = 0.;
      double vTprime_ = 0.;
  };

}
