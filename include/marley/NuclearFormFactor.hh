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
#include <limits>

// MARLEY includes
#include "marley/Logger.hh"
#include "marley/marley_utils.hh"

namespace marley {

  // Forward-declare the JSON class
  class JSON;

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

      /// Factory method to create derived class objects
      /// @param Z Proton number
      /// @param A Nucleon number
      /// @param ff_config JSON object containing the model configuration to
      /// use when initializing the new object
      static std::shared_ptr< marley::NuclearFormFactor > create( int Z, int A,
        const JSON& ff_config );

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

      virtual ~TrivialNuclearFormFactor() = default;

      inline virtual double F( double /*kappa*/ ) const override final
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

      virtual ~HelmNuclearFormFactor() = default;

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

      inline KleinNystrandNuclearFormFactor( int Z, int A,
        bool adapted = false, double r0 = DUMMY_r0_VALUE )
        : NuclearFormFactor( Z, A ), adapted_( adapted ), r0_( r0 )
      {
        // COHERENT-style "adapted" treatment uses the rms charge radius r0 to
        // calculate the effective nuclear radius. If the user did not supply
        // one, then look it up from the table of measurements.
        if ( adapted_ ) {
          if ( r0_ == DUMMY_r0_VALUE ) {
            // Ensure that the table has been loaded
            if ( !r0_table_ ) this->initialize_r0_table();

            int pdg = marley_utils::get_nucleus_pid( Z, A );
            auto iter = r0_table_->find( pdg );
            if ( iter != r0_table_->end() ) {
              r0_ = iter->second;
            }
            else {
              throw marley::Error( "Unable to find tabulated rms charge"
                " radius for nucleus with PDG code " + std::to_string( pdg )
                + ". Please specify a value in fm using the \"r0\" JSON key." );
            }
          }

          R_ = marley_utils::real_sqrt( 5.*r0_*r0_/3. - 10.*a_*a_ );
          return;
        }

        // Default treatment assigns the effective nuclear radius based solely
        // on the nucleon number
        R_ = 1.23 * std::pow( A_, marley_utils::ONE_THIRD ); // fm
      }

      virtual ~KleinNystrandNuclearFormFactor() = default;

      virtual double F( double kappa ) const override final;

      /// Returns the value of the rms charge radius used with the adapted
      /// version
      inline double r0() const {
        if ( !adapted_ ) MARLEY_LOG( WARN, "physics.formfactor" )
          << "Requested rms charge radius"
          << " when using default Klein-Nystrand nuclear form factor";
        return r0_;
      }

      /// Dummy value used to signal the need to look up the rms charge radius
      /// from a table of measurements
      static constexpr double DUMMY_r0_VALUE
        = std::numeric_limits< double >::lowest();

    protected:

      static void initialize_r0_table();

      /// @brief Stores measured rms charge radii for many nucleii
      /// @details Keys are nuclear PDG codes, values are rms charge radii (fm)
      inline static std::unique_ptr< std::map< int, double > > r0_table_;

      /// @brief Name of the data file containing the measured rms charge
      /// radii
      inline static const std::string r0_data_file_name_
        = "nuclear_charge_radii.js";

      /// Range (fm) of the assumed Yukawa potential
      double a_ = 0.7;

      /// Effective nuclear radius (fm)
      double R_;

      /// Flag indicating whether we are using the COHERENT-style "adapted"
      /// version
      bool adapted_ = false;

      /// Value of the rms charge radius, used only for the adapted version
      double r0_ = DUMMY_r0_VALUE;
  };

}
