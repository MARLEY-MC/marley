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

using ME_Type = marley::MatrixElement::TransitionType;

namespace {
  // In cases where no discrete level data are available, a continuum level density
  // is used all the way down to the ground state. To avoid asymptotically approaching
  // Ex = 0 in these cases, the de-excitation cascade will end once the excitation energy of
  // the residual nucleus falls below this (small) value. Excitation energies below this
  // value are considered "close enough" to the ground state for MARLEY not to worry about
  // further de-excitations.
  /// @todo Make this configurable?
  constexpr double CONTINUUM_GS_CUTOFF = 0.001; // MeV

  /// @brief The size of a tolerable discrepancy (in MeV) between the excitation energy
  /// stored in the event record and other (hopefully consistent) versions of it
  constexpr double EX_TOLERANCE = 1e-5; // MeV
}

void marley::NucleusDecayer::deexcite_residue( marley::Event& event,
  const marley::Level* plevel, int twoJ, marley::Parity P,
  marley::Generator& gen )
{
  // Get the excitation energy from the event. If the discrete level
  // pointer is not null, check that its excitation energy is consistent
  // with the one stored in the event record.
  double Ex = event.Ex();

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
  int residue_pdg = event.residue().pdg_code();
  int qIon = event.residue().charge();

  // Check that the residue PDG code makes sense. If it's not a nucleus,
  // warn the user and refuse to do the cascade.
  if ( !marley_utils::is_ion(residue_pdg) ) {
    MARLEY_LOG_WARNING() << "Unrecognized nuclear PDG code " << residue_pdg
      << " encountered in marley::NucleusDecayer::deexcite_residue(). The"
      << " de-excitation cascade will be skipped";
    return;
  }

  double residue_mass = event.residue().mass();

  // Ground-state residue mass
  double gs_residue_mass = mt.get_atomic_mass( residue_pdg )
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
    // Throw an error rather than trying to figure out how to do the right thing.
    /// @todo Revisit this
    throw marley::Error("Partially de-excited nuclear remnant encountered"
        " in marley::NucleusDecay::deexcite_residue().");
  }

  // If we've accessed a discrete nuclear level, make sure that its excitation
  // energy matches that in the event record. If it doesn't, complain by
  // throwing an error.
  if ( plevel ) {
    double Ex_level = plevel->energy();
    if ( std::abs(Ex - Ex_level) > EX_TOLERANCE ) throw marley::Error(
      "Excitation energy mismatch encountered in marley::NucleusDecayer::"
      "deexcite_residue(). Event has Ex = " + std::to_string(Ex) + " MeV"
      " while the initial discrete level has " + std::to_string(Ex_level)
      + " MeV");
  }

  // We've passed all the checks. Now decide whether we need to start out
  // with a continuum or discrete decay.
  bool continuum = ( plevel == nullptr );

  int initial_nucleus_pdg = event.residue().pdg_code();
  int Z = marley_utils::get_particle_Z( initial_nucleus_pdg );
  int A = marley_utils::get_particle_A( initial_nucleus_pdg );

  if ( continuum ) {

    // Dummy particles used for temporary storage of binary decay products
    // during the de-excitation cascade
    marley::Particle first, second;

    // Get a non-const reference to the nuclear residue. We'll use it to
    // update the event record during each step of the Hauser-Feshbach
    // cascade.
    marley::Particle& residue = event.residue();

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
      Z = marley_utils::get_particle_Z( residue.pdg_code() );
      A = marley_utils::get_particle_A( residue.pdg_code() );
      event.add_final_particle( first );
    }
  }

  if ( !continuum ) {
    // Either the selected initial level was bound (so it will only decay via
    // gamma emission) or the Hauser-Feshbach decay process has now accessed a
    // bound level in the residual nucleus. In either case, use gamma-ray decay
    // scheme data to sample the de-excitation gammas and add them to this
    // event's final particle list.
    marley::DecayScheme* dec_scheme
      = gen.get_structure_db().get_decay_scheme( Z, A );
    dec_scheme->do_cascade( *dec_scheme->get_pointer_to_closest_level(Ex),
      event, gen, event.residue().charge() );
  }

}

void marley::NucleusDecayer::deexcite_residue( marley::Event& event,
  const marley::MatrixElement& matrix_el, marley::Generator& gen )
{
  // Get a pointer to the selected starting nuclear level
  // (will be nullptr if we're in the continuum)
  const marley::Level* plevel = matrix_el.level();

  bool continuum = ( plevel == nullptr );

  int twoJ = 0;
  marley::Parity P;

  if ( !continuum ) {
    twoJ = plevel->twoJ();
    P = plevel->parity();
  }
  else {
    // The accessed nuclear level is in the unbound continuum

    // Load the initial twoJ and parity values into twoJ and P.  These
    // variables will be changed during every step of the Hauser-Feshbach decay
    // cascade.
    // Right now, matrix element types are represented with 0 <-> Fermi, 1 <->
    // Gamow-Teller.  Since the transitions are from the 40Ar ground state (Jpi
    // = 0+), these also correspond to the 40K* state spins.
    /// @todo Come up with a better way of determining the Jpi values that will
    /// work for forbidden transition operators.
    if ( matrix_el.type() == ME_Type::FERMI) twoJ = 0;
    else if ( matrix_el.type() == ME_Type::GAMOW_TELLER) twoJ = 2;
    else throw marley::Error( "Unrecognized matrix element type encountered"
      " during a continuum decay in marley::NucleusDecayer"
      "::deexcite_residue()" );

    // Fermi transition gives 0+ --> 0+, GT transition gives 0+ --> 1+
    /// @todo handle target nuclei that are not initially 0+
    /// @todo include possibility of negative parity here.
    P = marley::Parity( true ); // positive parity
  }

  // We now have all the information we need to start the cascade, so do it by
  // delegating to an overloaded version of this function.
  this->deexcite_residue( event, plevel, twoJ, P, gen );
}
