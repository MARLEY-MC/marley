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

// Nonstandard but widely-supported (see
// http://en.wikipedia.org/wiki/Pragma_once) preprocessor directive that
// prevents this file from being included multiple times. Another option is an
// include guard (http://en.wikipedia.org/wiki/Include_guard).
#pragma once

// Standard library includes
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// MARLEY includes
#include "marley/Reaction.hh"

namespace marley {
  class Generator;
}

namespace HepMC3 {
  class FourVector;
  class GenCrossSection;
  class GenEvent;
  class GenRunInfo;
  class GenVertex;
  class GenParticle;
}

namespace marley_hepmc3 {

  // G.R.2
  constexpr int NUHEPMC_MAJOR_VERSION = 1;
  constexpr int NUHEPMC_MINOR_VERSION = 0;
  constexpr int NUHEPMC_PATCH_VERSION = 0;

  // G.R.8
  void prepare_process_metadata( HepMC3::GenRunInfo& run_info );

  // E.C.1
  int get_nuhepmc_proc_id( const marley::Reaction::ProcessType pt );
  marley::Reaction::ProcessType from_nuhepmc_proc_id( const int proc_id );

  // G.R.9
  void prepare_vertex_status_metadata( HepMC3::GenRunInfo& run_info );

  // G.R.10
  void prepare_particle_status_metadata( HepMC3::GenRunInfo& run_info );

  // G.R.11
  void prepare_non_standard_pdg_code_metadata( HepMC3::GenRunInfo& run_info );

  // G.R.4, G.R.6, G.C.2, G.C.3
  void apply_nuhepmc_runinfo_conventions( HepMC3::GenRunInfo& run_info,
    const double flux_avg_xsec );

  // V.R.1
  // Vertex status codes for NuHepMC
  constexpr int NUHEPMC_PRIMARY_VERTEX = 1;
  constexpr int NUHEPMC_HF_DECAY_VERTEX = 22;
  constexpr int NUHEPMC_GAMMA_DECAY_VERTEX = 23;

  // P.R.1
  // Particle status codes for NuHepMC
  constexpr int NUHEPMC_FINAL_STATE_STATUS = 1;
  constexpr int NUHEPMC_PROJECTILE_STATUS = 4;
  constexpr int NUHEPMC_TARGET_STATUS = 20;
  constexpr int NUHEPMC_UNDECAYED_RESIDUE_STATUS = 27;
  constexpr int NUHEPMC_INTERMEDIATE_RESIDUE_STATUS = 28;

  constexpr double DUMMY_PARTICLE_MASS
    = std::numeric_limits< double >::lowest();

  void set_particle_charge( HepMC3::GenParticle& particle, int charge );

  int get_particle_charge( HepMC3::GenParticle& particle );

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    const HepMC3::FourVector& mom4, int pdg, int status,
    double mass = DUMMY_PARTICLE_MASS );

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, double px, double py, double pz, double E, int status,
    double mass = DUMMY_PARTICLE_MASS );

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, double px, double py, double pz, int status,
    double mass );

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, int status, double mass = DUMMY_PARTICLE_MASS );


  std::vector< std::shared_ptr< HepMC3::GenVertex > >
    get_vertices_with_status( int status, HepMC3::GenEvent& ev );

  std::vector< std::shared_ptr< HepMC3::GenParticle > >
    get_particles_with_status( int status, HepMC3::GenEvent& ev );

  std::shared_ptr< HepMC3::GenParticle >
    get_first_particle_with_status( int status, HepMC3::GenEvent& ev );

  std::shared_ptr< HepMC3::GenParticle > get_projectile(
    HepMC3::GenEvent& ev );

  std::shared_ptr< HepMC3::GenParticle > get_target(
    HepMC3::GenEvent& ev );

  std::shared_ptr< HepMC3::GenParticle > get_ejectile(
    HepMC3::GenEvent& ev );

  std::shared_ptr< HepMC3::GenParticle > get_residue(
    HepMC3::GenEvent& ev );

  /// Handles sampling and storing a random decay time for a particle decay
  /// vertex
  /// @param partial_width Partial decay width (MeV) for the decay process
  ///   of interest
  void store_decay_time( double partial_width, marley::Generator& gen,
    std::shared_ptr< HepMC3::GenVertex >& decay_vtx,
    const std::shared_ptr< HepMC3::GenParticle >& parent );

  /// @brief Print a human-readable summary of a MARLEY HepMC3 event to
  /// the given output stream in the standard MARLEY display format.
  /// @param ev The event to print (const reference).
  /// @param os The output stream to write to. Defaults to std::cout.
  void print_event( const HepMC3::GenEvent& ev, std::ostream& os = std::cout );

  /// @brief Check whether two GenRunInfo objects are compatible for merging
  /// @param ref The reference run information
  /// @param candidate The candidate run information to compare against
  /// @return An empty string if compatible, or a human-readable description
  /// of the first incompatibility found
  std::string check_run_info_compatibility(
    const HepMC3::GenRunInfo& ref,
    const HepMC3::GenRunInfo& candidate );
};
