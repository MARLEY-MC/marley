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

// Standard library includes
#include <memory>

// MARLEY includes
#include "marley/NuclearReaction.hh"
#include "marley/TabulatedXSec.hh"

namespace marley {

  /// @brief Generates inclusive scattering events and computes cross
  /// sections using tabulated nuclear responses
  class ContinuumNuclearReaction : public NuclearReaction {

    public:

      /// Defines the approach to handling possible reassignment of
      /// cross-section strength that falls below the excitation energy
      /// threshold for the continuum
      enum class SubContinuumMode { IGNORE, MIRROR, ACCUMULATE };

      ContinuumNuclearReaction( Reaction::ProcessType pt, int pdg_a,
        int pdg_b, int pdg_c, int pdg_d, int q_d,
        const std::shared_ptr<TabulatedXSec>& txsec,
        const std::string& source_file );

      virtual std::shared_ptr< HepMC3::GenEvent > create_event(
        int particle_id_a, double KEa, marley::Generator& gen ) const override;

      virtual double total_xs( int pdg_a, double KEa ) const override;

      inline const TabulatedXSec& get_tabulated_xsec() const
        { return *xsec_; }

      /// Gets the approach to handling sub-continuum cross-section strength
      static inline SubContinuumMode sub_continuum_mode() { return sc_mode_; }

      /// Sets the approach to handling sub-continuum cross-section strength
      static inline void set_sub_continuum_mode( SubContinuumMode scm )
        { sc_mode_ = scm; }

      static SubContinuumMode sub_continuum_mode_from_string(
        const std::string& sc_mode_str );

      static std::string string_from_sub_continuum_mode(
        const SubContinuumMode sc_mode );

    protected:

      virtual void set_description() override;

      /// @brief Helper function for create_event() that potentially reassigns
      /// the value of the energy transfer
      /// @details If reassignment is enabled, moves cross-section strength
      /// falling below the unbound threshold to the continuum by updating the
      /// value of the energy transfer. Returns true if the reassignment was
      /// successful and is kinematically allowed. Returns false otherwise.
      /// @param[in,out] w Energy transfer (MeV, not shifted by
      /// @f$ \Delta_\mathrm{IAS} @f$)
      /// @param[in] ctl Lepton scattering cosine
      /// @param[in] KEa Projectile kinetic energy (MeV)
      /// @return Returns true if the reassignment procedure was successful (or
      /// skipped because it was unnecessary or disabled). Returns false
      /// otherwise.
      bool reassign_sub_continuum( double& w, const double ctl,
        const double KEa ) const;

      /// @brief Helper function for reassign_sub_continuum() that solves
      /// for the outgoing lepton total energy that corresponds to the input
      /// kinematic variables.
      /// @param[in] Ex Nuclear excitation energy (MeV)
      /// @param[in] cos_theta Lepton scattering cosine
      /// @param[in] KEa Projectile kinetic energy
      /// @param[out] jacobian If this argument is not nullptr, then
      /// the target double will be filled with the value of the Jacobian needed
      /// to convert from @f$ d\sigma/dE_\ell @f$ to @f$ d\sigma/dE_x @f$.
      double get_Ec_from_Ex( const double Ex, const double cos_theta,
        const double KEa, double* jacobian = nullptr) const;

      /// @brief Helper object that handles cross section calculations
      std::shared_ptr< TabulatedXSec > xsec_;

      /// @brief Indicates the desired method for handling events with
      /// excitation energies originally sampled below the continuum threshold
      static SubContinuumMode sc_mode_;

      /// @brief Helper map used for conversions between a SubContinuumMode
      /// value and a std::string
      static const std::map< SubContinuumMode, std::string > sc_mode_string_map_;
  };

}
