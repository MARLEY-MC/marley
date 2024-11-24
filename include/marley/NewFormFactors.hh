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
#include <memory>
#include <string>
#include <vector>

namespace marley {

  class JSON;

  class SachsFormFactors {
    public:

      SachsFormFactors() {}

      /// Proton electric form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      virtual double GEp( double Q2 ) = 0;

      /// Neutron electric form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      virtual double GEn( double Q2 ) = 0;

      /// Proton magnetic form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      virtual double GMp( double Q2 ) = 0;

      /// Neutron magnetic form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      virtual double GMn( double Q2 ) = 0;

      /// Helper function that computes @f$ \tau = \frac{ Q^{2} }{ 4 m_N^{2} }
      /// @f$
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      virtual double tau( double Q2 );
  };

  class DipoleSachsFormFactors : public SachsFormFactors {

    public:

      /// @param Mv Vector mass parameter (MeV)
      DipoleSachsFormFactors( double Mv ) : Mv_( Mv ) {}

      virtual double GEp( double Q2 ) override final;
      inline virtual double GEn( double /*Q2*/ ) override final { return 0.; }
      virtual double GMp( double Q2 ) override final;
      virtual double GMn( double Q2 ) override final;

      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      inline double dipole( double Q2 )
        { return 1.0 / std::pow( 1.0 + Q2 / Mv_ / Mv_, 2 ); }

    protected:

      /// Vector mass parameter (MeV)
      double Mv_;
  };


  class BBBA05SachsFormFactors : public SachsFormFactors {

    public:

      BBBA05SachsFormFactors() {}

      virtual double GEp( double Q2 ) override final;
      virtual double GEn( double Q2 ) override final;
      virtual double GMp( double Q2 ) override final;
      virtual double GMn( double Q2 ) override final;

    protected:

      /// @brief Helper function for the form factor calculation @details
      /// Implements Eq. (6) from Nucl. Phys. B Proc. Suppl. 159, 127 (2006)
      /// https://doi.org/10.1016/j.nuclphysbps.2006.08.028
      /// @param tau Dimensionless negative square of the four-momentum
      /// transfer @f$ \tau = \frac{ Q^{2} }{ 4 m_N^{2} } @f$
      /// @param a_coeffs Coefficients @f$ a_k @f$ of polynomial involving @f$
      /// \tau^k @f$ in the numerator
      /// @param b_coeffs Coefficients @f$ b_k @f$ of polynomial involving @f$
      /// \tau^k @f$ in the denominator
      double bbba05_G( double tau, const std::vector< double >& a_coeffs,
        const std::vector< double >& b_coeffs );
  };

  class AxialFormFactors {
    public:

      AxialFormFactors() {}

      /// Axial-vector form factor
      virtual double FA( double Q2 ) = 0;

      /// Pseudoscalar form factor
      virtual double FP( double Q2 ) = 0;
  };

  class DipoleAxialFormFactors : public AxialFormFactors {
    public:

      // @param gA Axial-vector coupling constant (dimensionless)
      // @param Ma Axial mass parameter (MeV)
      DipoleAxialFormFactors( double gA, double Ma ) : gA_( gA ), Ma_( Ma ) {}

      inline double FA( double Q2 ) override final
        { return -gA_ / std::pow( 1.0 + Q2 / Ma_ / Ma_, 2 ); }

      double FP( double Q2 ) override final;

    protected:

      /// Axial-vector coupling constant of the nucleon (dimensionless)
      double gA_;

      /// Axial mass parameter (MeV)
      double Ma_;
  };

  class NewFormFactors {
    public:

      NewFormFactors( const JSON& config );

      /// Charged-current vector form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      inline double F1( double Q2 ) const
        { return this->F1p( Q2 ) - this->F1n( Q2 ); }

      /// Charged-current weak magnetic form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      inline double F2( double Q2 ) const
        { return this->F2p( Q2 ) - this->F2n( Q2 ); }

      /// Proton vector form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      double F1p( double Q2 ) const;

      /// Neutron vector form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      double F1n( double Q2 ) const;

      /// Proton weak magnetic form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      double F2p( double Q2 ) const;

      /// Neutron weak magnetic form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      double F2n( double Q2 ) const;

      /// Charged-current axial-vector form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      inline double FA( double Q2 ) const { return axial_ff_->FA( Q2 ); }

      /// Charged-current pseudoscalar form factor
      /// @param Q2 Negative square @f$ Q^{2} @f$ of the four-momentum transfer
      /// (MeV<sup> 2</sup>)
      inline double FP( double Q2 ) const { return axial_ff_->FP( Q2 ); }

      /// Get const access to the owned Sachs form factors
      inline const SachsFormFactors* sachs_ff() const
        { return sachs_ff_.get(); }

      /// Get const access to the owned axial form factors
      inline const AxialFormFactors* axial_ff() const
        { return axial_ff_.get(); }

    protected:

      std::shared_ptr< SachsFormFactors > sachs_ff_;
      std::shared_ptr< AxialFormFactors > axial_ff_;
  };

  /// Computes Sachs form factors in the limit @f$ Q^{2} \to 0 @f$
  class TrivialSachsFormFactors : public SachsFormFactors {

    public:

      TrivialSachsFormFactors() {}

      virtual double GEp( double Q2 ) override final;
      virtual double GEn( double Q2 ) override final;
      virtual double GMp( double Q2 ) override final;
      virtual double GMn( double Q2 ) override final;

      /// Neglect the @f$ Q^{2} @f$ dependence of the form factors entirely
      inline virtual double tau( double /*Q2*/ ) override final { return 0.; }
  };

  /// Computes axial form factors in the limit @f$ Q^{2} \to 0 @f$
  class TrivialAxialFormFactors : public AxialFormFactors {

    public:

      TrivialAxialFormFactors() {}

      double FA( double Q2 ) override final;
      double FP( double Q2 ) override final;

  };

} // marley namespace
