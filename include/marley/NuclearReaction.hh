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
#include <functional>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "marley/DecayScheme.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/MatrixElement.hh"
#include "marley/Reaction.hh"
#include "marley/StructureDatabase.hh"

namespace marley {

  class Generator;

  /// @brief A neutrino-nucleus reaction
  class NuclearReaction : public Reaction {

    public:

      /// @param pt Type of scattering process represented by this Reaction
      /// @param pdg_a Projectile PDG code
      /// @param pdg_b Target PDG code
      /// @param pdg_c Ejectile PDG code
      /// @param pdg_d Residue PDG code
      /// @param q_d Charge of the residue after the prompt 2->2 scatter
      /// represented by this NuclearReaction object
      NuclearReaction( ProcessType pt, int pdg_a, int pdg_b, int pdg_c,
        int pdg_d, int q_d, const std::string& source_file );

      inline virtual marley::TargetAtom atomic_target() const override final
        { return marley::TargetAtom( pdg_b_ ); }

      /// @brief Get the minimum lab-frame kinetic energy (MeV) of the
      /// projectile that allows this reaction to proceed via a transition to
      /// the residue's ground state
      double threshold_kinetic_energy() const override;

      /// @brief Get the maximum possible excitation energy (MeV) of the
      /// final-state residue that is kinematically allowed
      /// @param KEa Projectile lab-frame kinetic energy (MeV)
      double max_level_energy( double KEa ) const;

      /// Computes the weak nuclear charge @f$ Q_W = N - (1
      /// - 4\sin^2\theta_W)Z @f$ for the target nucleus
      /// @details In the expression above, @f$ N @f$ (@f$Z@f$) is the
      /// neutron (proton) number of the target nucleus and
      /// @f$ \theta_W @f$ is the weak mixing angle.
      double weak_nuclear_charge() const;

    protected:

      /// @brief Creates the description string based on the
      /// PDG code values for the initial and final particles
      virtual void set_description();

      /// @brief Helper function that sets the charges of the target
      /// and residue in an otherwise complete event record
      void set_charge_attributes( std::shared_ptr< HepMC3::GenEvent >& event )
        const;

      /// @brief Helper function that adds the nuclear level attributes
      /// (@f$ E_x @f$, @f$ 2J @f$, and parity) needed to keep track of the
      /// residue's de-excitation state to an otherwise complete event record
      /// @param residue GenParticle object for the residue
      /// @param E_level Residue excitation energy (MeV)
      /// @param twoJ Two times the residue spin
      /// @param P Intrinsic parity of the residue
      void set_nuclear_residue_attributes(
        const std::shared_ptr< HepMC3::GenParticle >& residue,
        double E_level, int twoJ, const marley::Parity& P ) const;

      /// @brief Helper function that makes a complete event object for a
      /// nuclear reaction
      /// @details This function should be called by
      /// marley::NuclearReaction::create_event() after CM frame scattering
      /// angles have been sampled for the ejectile. In addition to creating
      /// the event skeleton (by delegating to
      /// marley::Reaction::make_event_object()), it attaches the charge and
      /// nuclear level attributes needed to keep track of the residue's
      /// de-excitation state.
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param pc_cm Ejectile 3-momentum magnitude (MeV) in the CM frame
      /// @param cos_theta_c_cm Cosine of ejectile's CM frame polar angle
      /// @param phi_c_cm Ejectile's CM frame azimuthal angle (radians)
      /// @param Ec_cm Ejectile total energy (MeV) in the CM frame
      /// @param Ed_cm Residue total energy (MeV) in the CM frame
      /// @param E_level Residue excitation energy (MeV)
      /// @param twoJ Two times the residue spin
      /// @param P Intrinsic parity of the residue
      virtual std::shared_ptr< HepMC3::GenEvent > make_nuclear_event_object(
        double KEa, double pc_cm, double cos_theta_c_cm, double phi_c_cm,
        double Ec_cm, double Ed_cm, double E_level, int twoJ,
        const marley::Parity& P ) const;

      /// @brief Helper function that makes a complete event object for a
      /// nuclear reaction
      /// @details This function expects pre-made HepMC3::GenParticle
      /// objects as input that have four-momenta expressed in the lab frame.
      /// @param KEa Lab-frame kinetic energy (MeV) of the projectile
      /// @param ejectile GenParticle object for the ejectile
      /// @param residue GenParticle object for the residue
      /// @param E_level Residue excitation energy (MeV)
      /// @param twoJ Two times the residue spin
      /// @param P Intrinsic parity of the residue
      virtual std::shared_ptr< HepMC3::GenEvent > make_nuclear_event_object(
        double KEa, const std::shared_ptr< HepMC3::GenParticle >& ejectile,
        const std::shared_ptr< HepMC3::GenParticle >& residue,
        double E_level, int twoJ, const marley::Parity& P ) const;

      double md_gs_; ///< Ground state mass (MeV) of the residue

      int Zi_; ///< Target atomic number
      int Ai_; ///< Target mass number
      int Zf_; ///< Residue atomic number
      int Af_; ///< Residue mass number

      /// @brief Net charge of the residue (in units of the proton charge)
      /// following this reaction
      int q_d_;

      /// @brief Lab-frame kinetic energy of the projectile at threshold for
      /// this reaction (i.e., the residue is produced in its ground state, and
      /// all final-state particles are at rest in the CM frame)
      double KEa_threshold_;
  };

}
