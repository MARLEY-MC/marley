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

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/MatrixElement.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/Parity.hh"

using ME_Type = marley::MatrixElement::TransitionType;

namespace {
  // In cases where no discrete level data are available, a continuum level
  // density is used all the way down to the ground state. To avoid
  // asymptotically approaching Ex = 0 in these cases, the de-excitation cascade
  // will end once the excitation energy of the residual nucleus falls below
  // this (small) value. Excitation energies below this value are considered
  // "close enough" to the ground state for MARLEY not to worry about further
  // de-excitations.
  /// @todo Make this configurable?
  constexpr double CONTINUUM_GS_CUTOFF = 0.001; // MeV

  // The size of a tolerable discrepancy (in MeV) between the excitation energy
  // stored in the event record and other (hopefully consistent) versions of it
  constexpr double EX_TOLERANCE = 1e-5; // MeV
}

void marley::NucleusDecayer::process_event( HepMC3::GenEvent& event,
  marley::Generator& gen )
{
  auto undecayed_residues = marley_hepmc3::get_particles_with_status(
    marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, event );

  MARLEY_LOG( DEBUG, "physics.deexcitation" ) << "NucleusDecayer: processing "
    << undecayed_residues.size() << " undecayed residue(s)";

  // Check the reaction process that created this event. The process types
  // distinguish between discrete and continuum reactions, which is helpful
  // below.
  int proc_id = event.attribute< HepMC3::IntAttribute >(
    "signal_process_id" )->value();
  auto proc_type = marley_hepmc3::from_nuhepmc_proc_id( proc_id );
  bool is_continuum_channel = false;
  if ( proc_type == marley::Reaction::ProcessType::NeutrinoCC_Continuum
    || proc_type == marley::Reaction::ProcessType::AntiNeutrinoCC_Continuum 
    || proc_type == marley::Reaction::ProcessType::NC_Continuum )
  {
    is_continuum_channel = true;
  }

  for ( auto residue : undecayed_residues ) {

    // Get the residue excitation energy from the event. These values represent
    // its state immediately following the initial two-two scattering reaction.
    double Ex = residue->attribute< HepMC3::DoubleAttribute >( "Ex" )->value();
    int twoJ = residue->attribute< HepMC3::IntAttribute >( "twoJ" )->value();
    int p_int = residue->attribute< HepMC3::IntAttribute >( "parity" )->value();
    marley::Parity P( p_int );

    // If the residue is in its ground state, then there's nothing for us to do.
    // Just continue the loop without comment.
    if ( Ex == 0. ) continue;

    MARLEY_LOG( DEBUG, "physics.deexcitation" ) << "De-exciting residue PDG "
      << residue->pid() << ": Ex = " << Ex << " MeV, 2J = " << twoJ
      << ", P = " << P;

    // The excitation energy should be nonnegative. Complain if it's not.
    if ( Ex < 0. ) throw marley::Error("Negative excitation energy Ex = "
      + std::to_string(Ex) + " MeV encountered in marley::NucleusDecayer::"
      "deexcite_residue()");

    // To prevent accidental double application of the de-excitation cascade,
    // check that the residue mass is consistent with the excitation energy
    // stored in the event record (and thus was never decayed).
    const auto& mt = marley::MassTable::Instance();
    int initial_residue_pdg = residue->pid();
    int qIon = marley_hepmc3::get_particle_charge( *residue );

    // Check that the residue PDG code makes sense. If it's not a nucleus,
    // warn the user and refuse to do the cascade.
    if ( !marley_utils::is_ion(initial_residue_pdg) ) {
      MARLEY_LOG( WARN, "physics.deexcitation" )
        << "Unrecognized nuclear PDG code "
        << initial_residue_pdg << " encountered in marley::NucleusDecayer::"
        << "deexcite_residue(). The de-excitation cascade will be skipped";
      continue;
    }

    double residue_mass = residue->generated_mass();

    // Ground-state residue mass
    double gs_residue_mass = mt.get_atomic_mass( initial_residue_pdg )
      - qIon*mt.get_particle_mass( marley_utils::ELECTRON );

    double expected_residue_mass = gs_residue_mass + Ex;

    if ( std::abs(residue_mass - expected_residue_mass) > EX_TOLERANCE ) {

      if ( std::abs(residue_mass - gs_residue_mass) <= EX_TOLERANCE ) {
        MARLEY_LOG( WARN, "physics.deexcitation" )
          << "Encountered ground-state nuclear remnant"
          << " in marley::NucleusDecay::deexcite_residue(). The de-excitation"
          << " cascade has already been applied.";
        continue;
      }

      // If we get here, then the residue is not in its ground state but also
      // not in the initial excited state given in the event record. Something
      // went wrong with a partial application of a de-excitation cascade.
      // Throw an error rather than trying to figure out how to do the right
      // thing.
      /// @todo Revisit this
      throw marley::Error("Partially de-excited nuclear remnant encountered"
        " in marley::NucleusDecay::deexcite_residue().");
    }

    // Decide whether we need to start the de-excitation cascade from a
    // discrete nuclear level or from the continuum. Do this by comparing the
    // excitation energy from the event record to the "unbound threshold" for
    // the residue. If we're above the unbound threshold, do a continuum decay.
    // Also start with a continuum decay if no discrete level data are
    // available for the residue.
    auto* ds = gen.get_structure_db().get_decay_scheme( initial_residue_pdg );
    double unbound_threshold = mt.unbound_threshold( initial_residue_pdg );

    // If Reaction::set_level_ptrs() changes, you'll want to change this too.
    // TODO: find a better way of keeping the two pieces of code in sync
    // TODO: numerical round-off can cause issues with the first test, so you
    // should revisit it again when interfacing MARLEY with other codes
    // that use it solely as a de-excitation model. For now, the third
    // option in the logical OR prevents issues with numerical round-off near
    // the unbound threshold.
    bool continuum = ( Ex > unbound_threshold )
      || ( !ds ) || ( is_continuum_channel );

    // Keep track of whether the cascade was started from the continuum
    // or not. If it was started from a discrete level, we'll double-check that
    // discrete level's excitation energy below.
    bool started_from_continuum = continuum;

    MARLEY_LOG( DEBUG, "physics.deexcitation" ) << "De-excitation path: "
      << ( continuum ? "continuum (Hauser-Feshbach)" : "discrete gamma cascade" );

    if ( continuum ) {

      // Particles used for storage of binary decay products during the
      // de-excitation cascade
      auto first = std::make_shared< HepMC3::GenParticle >();
      auto second = std::make_shared< HepMC3::GenParticle >();

      // The selected level is unbound, so handle its de-excitation using
      // the Hauser-Feshbach statistical model.
      while ( continuum && Ex > CONTINUUM_GS_CUTOFF ) {

        auto& sdb = gen.get_structure_db();

        marley::HauserFeshbachDecay hfd( residue, Ex, twoJ, P, sdb );
        MARLEY_LOG( DEBUG, "physics.deexcitation.hauser" ) << hfd;

        int q_second;
        const auto& exit_channel = hfd.do_decay( Ex, twoJ, P, first, second,
          q_second, gen );

        continuum = exit_channel.is_continuum();

        double width_tot = hfd.total_width();
        double width_ec = exit_channel.width();

        MARLEY_LOG( DEBUG, "physics.deexcitation.hauser" )
          << "Hauser-Feshbach decay to " << first->pid()
          << " and " << second->pid();
        MARLEY_LOG( DEBUG, "physics.deexcitation.hauser" )
          << second->pid() << " is at Ex = " << Ex << " MeV.";

        // Create a new binary decay vertex
        auto decay_vtx = std::make_shared< HepMC3::GenVertex >();
        decay_vtx->set_status( marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX );

        decay_vtx->add_particle_in( residue );
        decay_vtx->add_particle_out( first );
        decay_vtx->add_particle_out( second );

        // Sample a decay time (MeV^{-1}) for emission of the chosen particle
        // and store this timing information in the new binary decay vertex
        marley_hepmc3::store_decay_time( width_ec, gen, decay_vtx, residue );

        // Add the decay vertex to the event record
        event.add_vertex( decay_vtx );

        // We can now set the charge of the daughter ion because it belongs
        // to the parent event (through the decay vertex)
        marley_hepmc3::set_particle_charge( *second, q_second );

        // We can also now set the attributes representing the daughter ion's
        // excitation energy, spin, and parity
        second->add_attribute( "Ex",
          std::make_shared< HepMC3::DoubleAttribute >(Ex) );
        second->add_attribute( "twoJ",
          std::make_shared< HepMC3::IntAttribute >(twoJ) );
        second->add_attribute( "parity",
          std::make_shared< HepMC3::IntAttribute >(static_cast<int>( P )) );

        // The daughter ion now takes the role of the residue for the next loop
        // iteration
        residue.swap( second );

        // Store some information about the total and partial widths of
        // the simulated compound nucleus decay in attributes attached to
        // the decay vertex
        decay_vtx->add_attribute( "TotalWidth",
          std::make_shared< HepMC3::DoubleAttribute >(width_tot) );

        decay_vtx->add_attribute( "ECWidth",
          std::make_shared< HepMC3::DoubleAttribute >(width_ec) );

        // In the case of a transition to the continuum, also store the partial
        // differential width for the chosen spin-parity of the daughter
        // nucleus
        if ( continuum ) {
          const auto& cec = dynamic_cast< const marley::ContinuumExitChannel& >(
            exit_channel );

          const auto* spw_ptr = cec.get_last_sampled_spw();
          double width_sp = spw_ptr->diff_width;

          decay_vtx->add_attribute( "SPWidth",
            std::make_shared< HepMC3::DoubleAttribute >(width_sp) );

          bool is_fragment_emission = exit_channel.emits_fragment();

          if ( is_fragment_emission ) {
            const auto* f_spw = static_cast< const marley
              ::FragmentContinuumExitChannel::FragmentSpinParityWidth* >(
              spw_ptr );
            decay_vtx->add_attribute( "two_j_frag",
              std::make_shared< HepMC3::IntAttribute >(f_spw->two_j_frag) );
            decay_vtx->add_attribute( "orb_l",
              std::make_shared< HepMC3::IntAttribute >(f_spw->orb_l) );
          }

          else {
            // Gamma-ray emission in the continuum
            const auto* g_spw = static_cast< const marley
              ::GammaContinuumExitChannel::GammaSpinParityWidth* >( spw_ptr );

            decay_vtx->add_attribute( "multipolarity",
              std::make_shared< HepMC3::IntAttribute >(g_spw->multipolarity) );
          }
        }
      }
    }

    if ( !continuum ) {
      // Either the selected initial level was bound (so it will only decay via
      // gamma emission) or the Hauser-Feshbach decay process has now accessed
      // a bound level in the residual nucleus. In either case, use gamma-ray
      // decay scheme data to sample the de-excitation gammas and add them to
      // this event's final particle list.
      marley::DecayScheme* dec_scheme = gen.get_structure_db()
        .get_decay_scheme( residue->pid() );

      // Start the gamma cascade from this discrete level
      marley::Level* lev = dec_scheme->get_pointer_to_closest_level( Ex );

      // If we get a null level pointer from the decay scheme, complain
      if ( !lev ) throw marley::Error( "Null nuclear level pointer encountered"
        " in marley::NucleusDecayer::deexcite_residue()" );

      // If we did not simulate any continuum decays before getting to this
      // point, then double-check that the excitation energy from the event
      // record and the initial level are consistent. If they're not, then
      // complain by throwing an error.
      if ( !started_from_continuum ) {
        double Ex_level = lev->energy();
        if ( std::abs(Ex - Ex_level) > EX_TOLERANCE ) {
          throw marley::Error( "Excitation energy mismatch encountered in"
            " marley::NucleusDecayer::deexcite_residue(). Event has Ex = "
            + std::to_string(Ex) + " MeV while the initial discrete level has "
            + std::to_string(Ex_level) + " MeV" );
        }
      }

      dec_scheme->do_cascade( *lev, event, gen, residue );
    }

  } // loop over undecayed residues

}
