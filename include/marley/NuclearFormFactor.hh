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

#pragma once

// Standard library includes
#include <cmath>

// MARLEY includes
#include "marley/marley_utils.hh"

namespace marley {

  /// Base class for nuclear form factor models
  class NuclearFormFactor {
    public:

      /// @param Z Proton number
      /// @param A Nucleon number
      NuclearFormFactor( int Z, int A ) : Z_( Z ), A_( A ) {}

      virtual ~NuclearFormFactor() = default;

      /// Evaluate the nuclear form factor at a given 3-momentum transfer
      /// @param kappa Magnitude of the 3-momentum transfer (MeV)
      virtual double F( double kappa ) const = 0;

      inline int Z() const { return Z_; }
      inline int A() const { return A_; }

    protected:

      /// Proton number
      int Z_;

      /// Nucleon number
      int A_;
  };

  // Implements a nuclear form factor that is trivially unity
  class TrivialNuclearFormFactor : public NuclearFormFactor {
    public:

      inline TrivialNuclearFormFactor( int Z, int A )
        : NuclearFormFactor( Z, A ) {}

      inline virtual double F( double kappa ) const override final
        { return 1.; }
  };

  /// Implements a nuclear form factor based on Phys. Rev. 104, 1466 (1956).
  /// Parameter values are taken J. Cosmol. Astropart. Phys. 4, 012 (2007)
  /// which itself took fit results from Astropart. Phys. 6, 87 (1996)
  class HelmNuclearFormFactor : public NuclearFormFactor {
    public:

      inline HelmNuclearFormFactor( int Z, int A ) : NuclearFormFactor( Z, A )
      {
        c_ = 1.23*std::pow( A, marley_utils::ONE_THIRD ) - 0.6; // fm
        this->update_effective_radius();
      }

      virtual double F( double kappa ) const override final;

      inline double s() const { return s_; }
      inline void set_s( double s ) {
        s_ = s;
        this->update_effective_radius();
      }

      inline double a() const { return a_; }

      inline void set_a( double a ) {
        a_ = a;
        this->update_effective_radius();
      }

      inline double c() const { return c_; }
      inline void set_c( double c ) {
        c_ = c;
        this->update_effective_radius();
      }

      inline void update_effective_radius() {
        R_ = marley_utils::real_sqrt( c_*c_ + 7.*marley_utils::pi
          *marley_utils::pi*a_*a_/3. - 5.*s_*s_ );
      }

    protected:

      // Parameters in the fit to muon spectroscopy data
      /// Nuclear skin thickness parameter (fm)
      double s_ = 0.9; // fm
      double a_ = 0.52; // fm
      double c_; // fm

      /// Effective nuclear radius (fm)
      double R_;
  };

  /// Implements a nuclear form factor based on Phys. Rev. C 60, 014903 (1999)
  /// Parameter values for the default and COHERENT's "adapted" version are
  /// described in Universe 9, 207 (2023)
  class KleinNystrandNuclearFormFactor : public NuclearFormFactor {
    public:

      inline KleinNystrandNuclearFormFactor( int Z, int A )
        : NuclearFormFactor( Z, A )
      {
        R_ = 1.23 * std::pow( A_, marley_utils::ONE_THIRD ); // fm
      }

      virtual double F( double kappa ) const override final;

      /// Switches to using a COHERENT-style "adapted" form factor
      /// @param r0 Proton rms radius (fm)
      void use_adapted_version( double r0 ) {
        R_ = marley_utils::real_sqrt( 5.*r0*r0/3. - 10.*a_*a_ );
      }

    protected:

      /// Range (fm) of the assumed Yukawa potential
      double a_ = 0.7;

      /// Effective nuclear radius (fm)
      double R_;
  };

}
