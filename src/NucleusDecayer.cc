/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Event.hh"
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

void marley::NucleusDecayer::process_event( marley::Event& event,
  marley::Generator& gen )
{
  // Get the residue excitation energy from the event. These values represent
  // its state immediately following the initial two-two scattering reaction.
  double Ex = event.Ex();
  int twoJ = event.twoJ();
  marley::Parity P = event.parity();

  // If the initial two-two scatter left the residue in its ground state,
  // there's nothing for us to do. Just return without comment.
  if ( Ex == 0. ) return;

  // The excitation energy should be nonnegative. Complain if it's not.
  if ( Ex < 0. ) throw marley::Error("Negative excitation energy Ex = "
    + std::to_string(Ex) + " MeV encountered in marley::NucleusDecayer::"
    "deexcite_residue()");

  // To prevent accidental double application of the de-excitation cascade,
  // check that the residue mass is consistent with the excitation energy
  // stored in the event record (and thus was never decayed).
  const auto& mt = marley::MassTable::Instance();
  int initial_residue_pdg = event.residue().pdg_code();
  int qIon = event.residue().charge();

  // Check that the residue PDG code makes sense. If it's not a nucleus,
  // warn the user and refuse to do the cascade.
  if ( !marley_utils::is_ion(initial_residue_pdg) ) {
    MARLEY_LOG_WARNING() << "Unrecognized nuclear PDG code "
      << initial_residue_pdg << " encountered in marley::NucleusDecayer::"
      << "deexcite_residue(). The de-excitation cascade will be skipped";
    return;
  }

  double residue_mass = event.residue().mass();

  // Ground-state residue mass
  double gs_residue_mass = mt.get_atomic_mass( initial_residue_pdg )
    - qIon*mt.get_particle_mass( marley_utils::ELECTRON );

  double expected_residue_mass = gs_residue_mass + Ex;

  if ( std::abs(residue_mass - expected_residue_mass) > EX_TOLERANCE ) {

    if ( std::abs(residue_mass - gs_residue_mass) <= EX_TOLERANCE ) {
      MARLEY_LOG_WARNING() << "Encountered ground-state nuclear remnant"
        << " in marley::NucleusDecay::deexcite_residue(). The de-excitation"
        << " cascade has already been applied.";
      return;
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

  // Decide whether we need to start the de-excitation cascade from a discrete
  // nuclear level or from the continuum. Do this by comparing the excitation
  // energy from the event record to the "unbound threshold" for the residue.
  // If we're above the unbound threshold, do a continuum decay. Also start
  // with a continuum decay if no discrete level data are available for the
  // residue.
  auto* ds = gen.get_structure_db().get_decay_scheme( initial_residue_pdg );
  double unbound_threshold = mt.unbound_threshold( initial_residue_pdg );

  // If Reaction::set_level_ptrs() changes, you'll want to change this too.
  // TODO: find a better way of keeping the two pieces of code in sync
  bool continuum = ( Ex > unbound_threshold ) || ( !ds );

  // Keep track of whether the cascade was started from the continuum
  // or not. If it was started from a discrete level, we'll double-check that
  // discrete level's excitation energy below.
  bool started_from_continuum = continuum;

  // Get a non-const reference to the nuclear residue. We'll use it to
  // update the event record during each step of the Hauser-Feshbach
  // cascade.
  marley::Particle& residue = event.residue();

  if ( continuum ) {

    // Dummy particles used for temporary storage of binary decay products
    // during the de-excitation cascade
    marley::Particle first, second;

    // The selected level is unbound, so handle its de-excitation using
    // the Hauser-Feshbach statistical model.
    while ( continuum && Ex > CONTINUUM_GS_CUTOFF ) {

      marley::HauserFeshbachDecay hfd( residue, Ex, twoJ, P, gen );
      MARLEY_LOG_DEBUG() << hfd;

      continuum = hfd.do_decay(Ex, twoJ, P, first, second);

      MARLEY_LOG_DEBUG() << "Hauser-Feshbach decay to " << first.pdg_code()
        << " and " << second.pdg_code();
      MARLEY_LOG_DEBUG() << second.pdg_code() << " is at Ex = "
        << Ex << " MeV.";

      residue = second;
      event.add_final_particle( first );
    }
  }

  if ( !continuum ) {
    // Either the selected initial level was bound (so it will only decay via
    // gamma emission) or the Hauser-Feshbach decay process has now accessed a
    // bound level in the residual nucleus. In either case, use gamma-ray decay
    // scheme data to sample the de-excitation gammas and add them to this
    // event's final particle list.
    marley::DecayScheme* dec_scheme = gen.get_structure_db()
      .get_decay_scheme( residue.pdg_code() );

    // Start the gamma cascade from this discrete level
    marley::Level* lev = dec_scheme->get_pointer_to_closest_level( Ex );

    // If we get a null level pointer from the decay scheme, complain
    if ( !lev ) throw marley::Error("Null nuclear level pointer encountered"
      " in marley::NucleusDecayer::deexcite_residue()");

    // If we did not simulate any continuum decays before getting to this
    // point, then double-check that the excitation energy from the event
    // record and the initial level are consistent. If they're not, complain
    // by throwing an error.
    if ( !started_from_continuum ) {
      double Ex_level = lev->energy();
      if ( std::abs(Ex - Ex_level) > EX_TOLERANCE ) {
        throw marley::Error("Excitation energy mismatch encountered in"
          " marley::NucleusDecayer::deexcite_residue(). Event has Ex = "
          + std::to_string(Ex) + " MeV while the initial discrete level has "
          + std::to_string(Ex_level) + " MeV");
      }
    }

    dec_scheme->do_cascade( *lev, event, gen, residue.charge() );
  }

}
