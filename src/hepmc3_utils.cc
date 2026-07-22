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

// Standard library includes
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/FourVector.h"
#include "HepMC3/GenCrossSection.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/GenParticle.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/Reaction.hh"

namespace {

  constexpr int DUMMY_NUHEPMC_PROC_ID = 0;

  // G.R.8 and E.C.1
  struct NuHepMCProcess {
    NuHepMCProcess( int procID, std::string name, std::string description )
      : id_( procID ), name_( name ), desc_( description ) {}

    int id_;
    std::string name_;
    std::string desc_;
  };

  /// @todo Vary description for discrete/continuum processes?
  const std::map< marley::Reaction::ProcessType, NuHepMCProcess >
    ptype_to_nuhepmc_proc =
  {
    { marley::Reaction::ProcessType::NeutrinoCC_Discrete,
      { 100, "vCC-discrete", "charged-current neutrino-nucleus"
        " scattering via discrete transitions" } },
    { marley::Reaction::ProcessType::AntiNeutrinoCC_Discrete,
      { 110, "anti-vCC-discrete", "charged-current antineutrino-nucleus"
        " scattering via discrete transitions" } },
    { marley::Reaction::ProcessType::NC_Discrete,
      { 150, "NC-discrete", "neutral-current (anti)neutrino-nucleus"
        " scattering via discrete transitions" } },
    { marley::Reaction::ProcessType::NuElectronElastic,
      { 700, "v-e", "(anti)neutrino-electron elastic scattering" } },
    { marley::Reaction::ProcessType::NeutrinoCC_Continuum,
      { 101, "vCC-continuum", "charged-current neutrino-nucleus"
        " scattering via continuum transitions" } },
    { marley::Reaction::ProcessType::AntiNeutrinoCC_Continuum,
      { 111, "anti-vCC-continuum", "charged-current antineutrino-nucleus"
        " scattering via continuum transitions" } },
    { marley::Reaction::ProcessType::NC_Continuum,
      { 151, "NC-continuum", "neutral-current (anti)neutrino-nucleus"
        " scattering via continuum transitions" } },
    { marley::Reaction::ProcessType::StandaloneDecay,
      { 800, "standalone-decay", "standalone nuclear de-excitation"
        " with no simulated primary reaction" } },
  };

  // G.R.9
  std::map< int, std::pair< std::string, std::string > >
    vertex_status_map =
  {

    { marley_hepmc3::NUHEPMC_PRIMARY_VERTEX, { "Primary",
      "Represents the primary interaction" } },

    { marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX, { "HFDecay",
      "Represents a nuclear de-excitation step simulated using the"
      " Hauser-Feshbach treatment" } },

    { marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX, { "GammaDecay",
      "Represents a nuclear de-excitation step simulated using tabulated"
      " gamma-ray branching ratios" } },

  };

  // G.R.10 and P.R.1
  std::map< int, std::pair< std::string, std::string > >
    particle_status_map =
  {

    { marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, { "Final-state",
      "Undecayed physical particle" } },

    { marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, { "Projectile",
      "Incoming beam particle" } },

    { marley_hepmc3::NUHEPMC_TARGET_STATUS, { "Target",
      "Target particle struck by incoming beam particle in"
      " the primary interaction" } },

    { marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, { "UndecayedRemnant",
      "Nuclear remnant before simulation of nuclear de-excitations" } },

    { marley_hepmc3::NUHEPMC_INTERMEDIATE_RESIDUE_STATUS,
      { "IntermediateRemnant",
        "Nuclear remnant during simulation of nuclear de-excitations" } },

  };

  // G.R.4
  std::vector< std::string > nuhepmc_convention_vec = {
    "G.C.2", // flux-averaged total cross section stored on GenRunInfo
    "G.C.3", // citation metadata for MARLEY publications
    "E.C.1", // process ID categorization
    "E.C.2", // total cross section per event
    "E.C.3", // process-specific cross section per event
  };

}

namespace marley_hepmc3 {

  int get_nuhepmc_proc_id( const marley::Reaction::ProcessType pt ) {
    auto itr = ptype_to_nuhepmc_proc.find( pt );
    if ( itr != ptype_to_nuhepmc_proc.end() ) {
      return itr->second.id_;
    }
    return DUMMY_NUHEPMC_PROC_ID;
  }

  marley::Reaction::ProcessType from_nuhepmc_proc_id( const int proc_id ) {
    marley::Reaction::ProcessType pt = marley::Reaction::ProcessType::Unknown;
    auto itr = std::find_if( ptype_to_nuhepmc_proc.cbegin(),
      ptype_to_nuhepmc_proc.cend(), [ proc_id ]( auto& pair ) -> bool
      { return pair.second.id_ == proc_id; }
    );
    if ( itr != ptype_to_nuhepmc_proc.end() ) {
      pt = itr->first;
    }
    return pt;
  }

  void set_particle_charge( HepMC3::GenParticle& particle, int charge ) {
    bool added_ok = particle.add_attribute( "charge",
      std::make_shared< HepMC3::IntAttribute >(charge) );
    if ( !added_ok ) {
      throw marley::Error( "Failed to set particle charge in marley_hepmc3"
        "::set_particle_charge()" );
    }
  }

  int get_particle_charge( HepMC3::GenParticle& particle ) {
    // Return the charge stored in the particle attributes if it is set
    auto* q_ptr = particle.attribute< HepMC3::IntAttribute >( "charge" ).get();
    if ( q_ptr ) return q_ptr->value();
    // Otherwise, look up the charge based on the PDG code
    return marley_utils::get_particle_charge( particle.pid() );
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    const HepMC3::FourVector& mom4, int pdg, int status,
    double mass )
  {
    auto particle = std::make_shared< HepMC3::GenParticle >(
      mom4, pdg, status );
    if ( mass != DUMMY_PARTICLE_MASS ) {
      particle->set_generated_mass( mass );
    }
    return particle;
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, double px, double py, double pz, double E, int status,
    double mass )
  {
    HepMC3::FourVector mom4( px, py, pz, E );
    return make_particle( mom4, pdg, status, mass );
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, double px, double py, double pz, int status,
    double mass )
  {
    double E = marley_utils::real_sqrt( px*px + py*py + pz*pz + mass*mass );
    HepMC3::FourVector mom4( px, py, pz, E );
    return make_particle( mom4, pdg, status, mass );
  }

  std::shared_ptr< HepMC3::GenParticle > make_particle(
    int pdg, int status, double mass )
  {
    HepMC3::FourVector mom4;
    if ( mass != DUMMY_PARTICLE_MASS ) {
      mom4.set_e( mass );
    }
    return make_particle( mom4, pdg, status, mass );
  }



  std::vector< std::shared_ptr< HepMC3::GenParticle > >
    get_particles_with_status( int status, HepMC3::GenEvent& ev )
  {
    const auto& particles = ev.particles();
    std::vector< std::shared_ptr< HepMC3::GenParticle > > found_particles;
    for ( auto& p : particles ) {
      if ( p->status() == status ) {
        found_particles.push_back( p );
      }
    }
    return found_particles;
  }

  std::vector< std::shared_ptr< HepMC3::GenVertex > >
    get_vertices_with_status( int status, HepMC3::GenEvent& ev )
  {
    const auto& vertices = ev.vertices();
    std::vector< std::shared_ptr< HepMC3::GenVertex > > found_vertices;
    for ( auto& v : vertices ) {
      if ( v->status() == status ) {
        found_vertices.push_back( v );
      }
    }
    return found_vertices;
  }

  std::shared_ptr< HepMC3::GenParticle >
    get_first_particle_with_status( int status, HepMC3::GenEvent& ev )
  {
    auto particle_ptrs = get_particles_with_status( status, ev );
    if ( !particle_ptrs.empty() ) return particle_ptrs.front();
    return nullptr;
  }

  std::shared_ptr< HepMC3::GenParticle > get_projectile(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_PROJECTILE_STATUS, ev );
  }

  std::shared_ptr< HepMC3::GenParticle > get_target(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_TARGET_STATUS, ev );
  }

  std::shared_ptr< HepMC3::GenParticle > get_ejectile(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS, ev );
  }

  std::shared_ptr< HepMC3::GenParticle > get_residue(
    HepMC3::GenEvent& ev )
  {
    return get_first_particle_with_status(
      marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, ev );
  }

  // G.R.8
  void prepare_process_metadata( HepMC3::GenRunInfo& run_info ) {
    std::vector< int > proc_id_vec;
    for ( const auto& pair : ptype_to_nuhepmc_proc ) {
      int proc_id = pair.second.id_;
      proc_id_vec.push_back( proc_id );

      std::string attr_prefix = "NuHepMC.ProcessInfo["
        + std::to_string( proc_id ) + "].";

      run_info.add_attribute( attr_prefix + "Name",
        std::make_shared< HepMC3::StringAttribute >(pair.second.name_) );

      run_info.add_attribute( attr_prefix + "Description",
        std::make_shared< HepMC3::StringAttribute >(pair.second.desc_) );
    }

    run_info.add_attribute( "NuHepMC.ProcessIDs",
      std::make_shared< HepMC3::VectorIntAttribute >(proc_id_vec) );
  }

  // G.R.9
  void prepare_vertex_status_metadata( HepMC3::GenRunInfo& run_info ) {
    std::vector< int > status_vec;
    for ( const auto& pair : vertex_status_map ) {
      int status = pair.first;
      status_vec.push_back( status );

      std::string attr_prefix = "NuHepMC.VertexStatusInfo["
        + std::to_string( status ) + "].";

      run_info.add_attribute( attr_prefix + "Name",
        std::make_shared< HepMC3::StringAttribute >(pair.second.first) );

      run_info.add_attribute( attr_prefix + "Description",
        std::make_shared< HepMC3::StringAttribute >(pair.second.second) );
    }

    run_info.add_attribute( "NuHepMC.VertexStatusIDs",
      std::make_shared< HepMC3::VectorIntAttribute >(status_vec) );
  }

  // G.R.10
  void prepare_particle_status_metadata( HepMC3::GenRunInfo& run_info ) {
    std::vector< int > status_vec;
    for ( const auto& pair : particle_status_map ) {
      int status = pair.first;
      status_vec.push_back( status );

      std::string attr_prefix = "NuHepMC.ParticleStatusInfo["
        + std::to_string( status ) + "].";

      run_info.add_attribute( attr_prefix + "Name",
        std::make_shared< HepMC3::StringAttribute >(pair.second.first) );

      run_info.add_attribute( attr_prefix + "Description",
        std::make_shared< HepMC3::StringAttribute >(pair.second.second) );
    }

    run_info.add_attribute( "NuHepMC.ParticleStatusIDs",
      std::make_shared< HepMC3::VectorIntAttribute >(status_vec) );
  }

  // G.R.11
  void prepare_non_standard_pdg_code_metadata( HepMC3::GenRunInfo& run_info )
  {
    // PDG code 0 is used as a dummy/absent projectile in standalone
    // nuclear de-excitation events generated by "marley decay"
    const std::vector< int > non_standard_PDGs = { 0 };
    run_info.add_attribute( "NuHepMC.AdditionalParticleNumbers",
      std::make_shared< HepMC3::VectorIntAttribute >(non_standard_PDGs) );

    run_info.add_attribute( "NuHepMC.AdditionalParticleNumbers[0].Name",
      std::make_shared< HepMC3::StringAttribute >( "Absent" ) );

    run_info.add_attribute( "NuHepMC.AdditionalParticleNumbers[0].Description",
      std::make_shared< HepMC3::StringAttribute >( "Dummy particle"
        " representing an absent projectile in a standalone nuclear"
        " de-excitation event (used by \"marley decay\")" ) );
  }

  // NOTE: the input flux-averaged total cross section is assumed to be in
  // natural units (MeV^{-2})
  void apply_nuhepmc_runinfo_conventions( HepMC3::GenRunInfo& run_info,
    const double flux_avg_xsec )
  {

    // G.R.4
    run_info.add_attribute( "NuHepMC.Conventions",
      std::make_shared< HepMC3::VectorStringAttribute >(
        nuhepmc_convention_vec )
    );

    // G.R.6
    run_info.add_attribute( "NuHepMC.Units.CrossSection.Unit",
      std::make_shared< HepMC3::StringAttribute >( "pb" )
    );

    run_info.add_attribute( "NuHepMC.Units.CrossSection.TargetScale",
      std::make_shared< HepMC3::StringAttribute >( "PerAtom" )
    );

    // G.C.2
    double xsec_picobarn = flux_avg_xsec * marley_utils::hbar_c2
      * marley_utils::fm2_to_picobarn;
    run_info.add_attribute( "NuHepMC.FluxAveragedTotalCrossSection",
      std::make_shared< HepMC3::DoubleAttribute >( xsec_picobarn )
    );

    // G.C.3
    std::vector< std::string > marley_DOIs = {
      "10.1103/PhysRevC.103.044604",
      "10.1016/j.cpc.2021.108123"
    };

    run_info.add_attribute( "NuHepMC.Citations.Generator.DOI",
      std::make_shared< HepMC3::VectorStringAttribute >( marley_DOIs )
    );

    std::vector< std::string > marley_arXivs = {
      "2010.02393",
      "2101.11867",
      "2604.26801"
    };

    run_info.add_attribute( "NuHepMC.Citations.Generator.arXiv",
      std::make_shared< HepMC3::VectorStringAttribute >( marley_arXivs )
    );

    std::vector< std::string > marley_INSPIREs = {
      "Gardiner:2020ulp",
      "Gardiner:2021qfr"
    };

    run_info.add_attribute( "NuHepMC.Citations.Generator.InspireHEP",
      std::make_shared< HepMC3::VectorStringAttribute >( marley_INSPIREs )
    );

  }

  std::string check_run_info_compatibility(
    const HepMC3::GenRunInfo& ref,
    const HepMC3::GenRunInfo& candidate )
  {
    if ( ref.weight_names() != candidate.weight_names() ) {
      return "weight names differ between files";
    }

    auto ref_attr = ref.attribute< HepMC3::StringAttribute >(
      "MARLEY.JSONconfig" );
    auto cand_attr = candidate.attribute< HepMC3::StringAttribute >(
      "MARLEY.JSONconfig" );

    if ( ref_attr && cand_attr ) {
      try {
        marley::JSON ref_json = marley::JSON::load( ref_attr->value() );
        marley::JSON cand_json = marley::JSON::load( cand_attr->value() );

        auto strip_run_keys = []( marley::JSON& j ) -> marley::JSON {
          marley::JSON result = marley::JSON::object();
          for ( const auto& [key, value] : j.object_range() ) {
            if ( key != "seed" && key != "generate" ) {
              result[key] = value;
            }
          }
          return result;
        };

        marley::JSON ref_stripped = strip_run_keys( ref_json );
        marley::JSON cand_stripped = strip_run_keys( cand_json );

        if ( ref_stripped.dump_string() != cand_stripped.dump_string() ) {
          return "MARLEY JSON configuration differs between files";
        }
      } catch ( const std::exception& e ) {
        return "failed to parse MARLEY JSON configuration: "
          + std::string( e.what() );
      }
    } else if ( static_cast< bool >( ref_attr )
      != static_cast< bool >( cand_attr ) )
    {
      return "one file has MARLEY JSON configuration and the other does not";
    }

    return {};
  }

}

// Handles sampling and storing a random decay time for a binary decay vertex
void marley_hepmc3::store_decay_time( double partial_width,
  marley::Generator& gen, std::shared_ptr< HepMC3::GenVertex >& decay_vtx,
  const std::shared_ptr< HepMC3::GenParticle >& parent )
{
  // Sample a decay time (MeV^{-1}) to assign to the decay vertex
  double decay_time = gen.sample_decay_time( partial_width );
  MARLEY_LOG( TRACE, "physics.deexcitation.gamma" ) << "decay_time = "
    << marley_utils::hbar * decay_time << " s";

  // Convert to the appropriate time units (cm) for a NuHepMC 4-position.
  // See marley::Reaction::make_event_object() where the units are defined.
  decay_time *= marley_utils::hbar_c * marley_utils::fm_to_cm;

  // The decay width treatment above assumes that the parent particle is
  // at rest. Apply a (typically very small) time dilation correction
  // since it may be moving in the laboratory frame.
  const HepMC3::FourVector& mom4_parent = parent->momentum();
  double E2_parent = std::pow( mom4_parent.e(), 2 );
  double beta2_parent = mom4_parent.length2() / E2_parent;
  double gamma_parent = 1. / marley_utils::real_sqrt( 1. - beta2_parent );
  decay_time *= gamma_parent;

  // Get the creation time of the parent particle from its starting vertex
  const auto parent_prod_vtx = parent->production_vertex();
  if ( !parent_prod_vtx ) throw marley::Error( "Could not access parent"
    " particle production vertex in marley_hepmc3::store_decay_time()" );
  double old_time = parent_prod_vtx->position().t(); // cm

  // Set and store the absolute time in the decay vertex
  // TODO: add spatial information as needed
  double new_time = old_time + decay_time; // cm
  HepMC3::FourVector decay_pos4;
  decay_pos4.set_t( new_time );

  decay_vtx->set_position( decay_pos4 );
}

namespace {

  // ──────────────────────────────────────────────────────────────────
  // Helpers for marley_hepmc3::print_event()
  // ──────────────────────────────────────────────────────────────────

  // Named UTF-8 constants for the Unicode characters used in the display.
  const std::string BOX_HEAVY   = "━";  // U+2501 BOX DRAWINGS HEAVY HORIZONTAL
  const std::string THIN_CHAR   = "─";  // U+2500 BOX DRAWINGS LIGHT HORIZONTAL
  const std::string ARROW_RIGHT = "►";  // U+25BA BLACK RIGHT-POINTING POINTER
  const std::string BOX_VERT    = "│";  // U+2502 BOX DRAWINGS LIGHT VERTICAL

  // Thick separator: 71 × BOX_HEAVY (━), programmatically generated.
  const std::string THICK_SEP = []() {
    std::string s;
    s.reserve( 71 * 3 );
    for ( int k = 0; k < 71; ++k ) s += BOX_HEAVY;
    return s;
  }();

  // Compute display width of a UTF-8 string.
  // Counts Unicode code points and skips combining characters (U+0300–U+036F)
  // since they do not advance the terminal cursor.
  int utf8_display_width( const std::string& s ) {
    int w = 0;
    size_t i = 0;
    while ( i < s.size() ) {
      auto c = static_cast< unsigned char >( s[i] );
      uint32_t cp = 0;
      size_t len = 1;
      if ( c < 0x80u ) {
        cp = c; len = 1;
      } else if ( c < 0xE0u ) {
        cp = static_cast< uint32_t >( c & 0x1Fu ) << 6;
        if ( i + 1 < s.size() )
          cp |= static_cast< unsigned char >( s[i+1] ) & 0x3Fu;
        len = 2;
      } else if ( c < 0xF0u ) {
        cp = static_cast< uint32_t >( c & 0x0Fu ) << 12;
        if ( i + 1 < s.size() )
          cp |= static_cast< uint32_t >(
            static_cast< unsigned char >( s[i+1] ) & 0x3Fu ) << 6;
        if ( i + 2 < s.size() )
          cp |= static_cast< unsigned char >( s[i+2] ) & 0x3Fu;
        len = 3;
      } else {
        // 4-byte sequence: decode all three continuation bytes
        cp = static_cast< uint32_t >( c & 0x07u ) << 18;
        if ( i + 1 < s.size() )
          cp |= static_cast< uint32_t >(
            static_cast< unsigned char >( s[i+1] ) & 0x3Fu ) << 12;
        if ( i + 2 < s.size() )
          cp |= static_cast< uint32_t >(
            static_cast< unsigned char >( s[i+2] ) & 0x3Fu ) << 6;
        if ( i + 3 < s.size() )
          cp |= static_cast< unsigned char >( s[i+3] ) & 0x3Fu;
        len = 4;
      }
      i += len;
      // Combining characters (U+0300–U+036F) contribute zero display width
      if ( cp < 0x0300u || cp > 0x036Fu ) ++w;
    }
    return w;
  }

  // Build a thin separator line of total display width `total` characters.
  // The prefix may contain multi-byte UTF-8 characters; fill uses THIN_CHAR.
  std::string make_thin_sep( const std::string& prefix, int total = 71 ) {
    std::string s = prefix;
    int fill = total - utf8_display_width( prefix );
    for ( int k = 0; k < fill; ++k ) s += THIN_CHAR;
    return s;
  }

  // Left-justify s in a field of target display-width chars, padding with
  // spaces
  std::string left_pad( const std::string& s, int target ) {
    int dw = utf8_display_width( s );
    int pad = ( target > dw ) ? ( target - dw ) : 0;
    return s + std::string( static_cast< size_t >( pad ), ' ' );
  }

  // Format the spin-parity string from twoJ and parity integer
  std::string format_jp( int twoJ, int par ) {
    std::string j;
    if ( twoJ % 2 == 0 ) j = std::to_string( twoJ / 2 );
    else j = std::to_string( twoJ ) + "/2";
    return j + ( par >= 0 ? "+" : "-" );
  }

  // Return the 6-display-char tag for the interaction-chain primary line.
  // All returned strings are exactly 6 display characters wide.
  std::string chain_tag_primary( marley::Reaction::ProcessType pt,
    int proj_pdg )
  {
    std::string tag = marley_utils::get_particle_symbol( proj_pdg );
    switch ( pt ) {
      case marley::Reaction::ProcessType::NeutrinoCC_Discrete:
      case marley::Reaction::ProcessType::NeutrinoCC_Continuum:
      case marley::Reaction::ProcessType::AntiNeutrinoCC_Discrete:
      case marley::Reaction::ProcessType::AntiNeutrinoCC_Continuum:
        tag += " CC ";
        break;
      case marley::Reaction::ProcessType::NC_Discrete:
      case marley::Reaction::ProcessType::NC_Continuum:
        tag += " NC ";
        break;
      case marley::Reaction::ProcessType::NuElectronElastic:
        tag += "+e⁻ ";
        break;
      case marley::Reaction::ProcessType::StandaloneDecay:
        return " dcay ";
      default:
        return " ???? ";
    }
    return tag;
  }

  // Helper: get double attribute from particle, returns false if absent
  bool get_double_attr( const HepMC3::ConstGenParticlePtr& p,
    const std::string& name, double& val )
  {
    auto attr = p->attribute< HepMC3::DoubleAttribute >( name );
    if ( !attr ) return false;
    val = attr->value();
    return true;
  }

  // Helper: get double attribute from vertex, returns false if absent
  bool get_double_attr_vtx( const HepMC3::ConstGenVertexPtr& v,
    const std::string& name, double& val )
  {
    auto attr = v->attribute< HepMC3::DoubleAttribute >( name );
    if ( !attr ) return false;
    val = attr->value();
    return true;
  }

  // Helper: get int attribute from particle, returns false if absent
  bool get_int_attr( const HepMC3::ConstGenParticlePtr& p,
    const std::string& name, int& val )
  {
    auto attr = p->attribute< HepMC3::IntAttribute >( name );
    if ( !attr ) return false;
    val = attr->value();
    return true;
  }

  // Helper: get int attribute from vertex, returns false if absent
  bool get_int_attr_vtx( const HepMC3::ConstGenVertexPtr& v,
    const std::string& name, int& val )
  {
    auto attr = v->attribute< HepMC3::IntAttribute >( name );
    if ( !attr ) return false;
    val = attr->value();
    return true;
  }

  // Format a momentum component with explicit sign and 3 decimal places
  std::string fmtpm( double v ) {
    std::ostringstream ss;
    ss << std::showpos << std::fixed << std::setprecision( 3 ) << v;
    return ss.str();
  }

  // Print one particle line in a vertex IN or OUT block.
  //
  // Rules (applied in order):
  //  1. If the particle has an "Ex" attribute and we are in the IN block
  //     (is_in_block == true): print simplified nuclear excitation line only.
  //  2. Otherwise print the full PDG/E/momentum line, then optionally the
  //     excitation annotation below.
  //
  // Full-line momentum display (for OUT and primary-vertex IN):
  //  a. If |p3| < 1e-6 MeV/c:                   show "(at rest)"
  //  b. If px≈0 and py≈0 (projectile on z-axis): show "pz = <+val>"
  //  c. If photon (PDG 22) in gamma-cascade vertex with multipolarity > 0:
  //     show "M<L> (l = <L>)" replacing momentum
  //  d. If nuclear (PDG > 1e9) and Ex attribute absent or zero,
  //     and vertex is a de-excitation vertex (is_deex_vtx == true):
  //     show "[ground state]" replacing momentum
  //  e. Otherwise: show "p = (px, py, pz) MeV/c"
  //
  // Name field: 11 display chars (left-justified).
  // PDG field: 10 chars (right-justified).
  // Excitation annotation (31-space indent): follows OUT particle if Ex > 0.
  void print_particle_line( std::ostream& os,
    const HepMC3::ConstGenParticlePtr& p,
    bool is_in_block,
    bool is_deex_vtx,
    bool is_gamma_vtx,
    int  vtx_multipolarity )
  {
    int pdg = p->pid();

    double Ex  = 0.;
    int twoJ   = -1;
    int par    = 1;
    bool has_ex = get_double_attr( p, "Ex", Ex );
    get_int_attr( p, "twoJ",   twoJ );
    get_int_attr( p, "parity", par  );
    bool excited = has_ex && ( Ex > 0. );

    std::string name = marley_utils::get_particle_symbol( pdg, excited );

    // ── Rule 1: simplified IN-block line for nuclear excited states ───────
    if ( is_in_block && has_ex ) {
        os << "    " << left_pad( name, 11 )
           << "Ex = " << std::fixed << std::setprecision( 2 ) << Ex << " MeV";
        if ( twoJ >= 0 )
          os << "   Jπ = " << format_jp( twoJ, par ); // Jπ = ...
      os << "\n";
      return;
    }

    // ── Rules 2+: full particle line ─────────────────────────────────────
    const HepMC3::FourVector& mom = p->momentum();
    double E  = mom.e();
    double px = mom.px();
    double py = mom.py();
    double pz = mom.pz();
    double p3mag = std::sqrt( px*px + py*py + pz*pz );

    std::ostringstream line;
    line << "    " << left_pad( name, 11 )
         << "E = " << std::fixed << std::setprecision( 3 )
         << std::setw( 9 ) << E << " MeV  ";

    // Determine the momentum display.
    // "ground state nucleus" in the OUT block of a de-excitation vertex:
    //  - Must be a nuclear PDG code (> 1e9)
    //  - Must be in the OUT block (not IN)
    //  - Must be in a de-excitation vertex (HF or gamma)
    //  - Must not have an excitation energy set
    //  - Excludes final-state particles in HF vertices: those are emitted
    //    fragments (e.g. α, d, t) with genuine momenta, not daughter residues.
    //    In gamma-cascade vertices all nuclear OUT particles are daughters.
    bool is_nuclear = ( pdg > 1000000000 );
    bool is_final_state_particle =
      ( p->status() == marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS );
    bool ground_state_nucleus = is_nuclear && !is_in_block && is_deex_vtx
      && ( !has_ex || Ex == 0. )
      && ( is_gamma_vtx || !is_final_state_particle );

    if ( ground_state_nucleus ) {
      // Rule 2d: nuclear ground state in de-excitation vertex
      line << "[ground state]";
    } else if ( p3mag < 1e-6 ) {
      // Rule 2a: at rest
      line << "(at rest)";
    } else if ( std::abs(px) < 1e-6 && std::abs(py) < 1e-6 ) {
      // Rule 2b: projectile aligned along z-axis
      line << "pz = " << std::showpos << std::fixed << std::setprecision( 3 )
           << pz << std::noshowpos;
    } else if ( pdg == 22 && is_gamma_vtx && vtx_multipolarity > 0 ) {
      // Rule 2c: photon with known multipolarity in gamma-cascade vertex
      line << "M" << vtx_multipolarity
           << " (l = " << vtx_multipolarity << ")";
    } else {
      // Rule 2e: general 3-momentum
      line << "p = (" << fmtpm(px) << ", " << fmtpm(py) << ", "
           << fmtpm(pz) << ") MeV/c";
    }

    os << line.str() << "\n";

    // ── Excitation annotation (31-space indent) ───────────────────────────
    // Applies to OUT particles with Ex > 0, but not to photons in gamma vertices.
    if ( !is_in_block && excited && !ground_state_nucleus
      && !( pdg == 22 && is_gamma_vtx ) && twoJ >= 0 )
    {
        os << "                               Ex = "
           << std::fixed << std::setprecision( 2 ) << Ex
           << " MeV   Jπ = " << format_jp( twoJ, par ) << "\n";
    }
  }

  // Get the display name for a particle in the interaction-chain summary.
  // Appends "*" to nuclei when Ex > 0.
  std::string chain_particle_name( const HepMC3::ConstGenParticlePtr& p ) {
    int pdg = p->pid();
    double Ex = 0.;
    bool excited = get_double_attr( p, "Ex", Ex ) && ( Ex > 0. );
    return marley_utils::get_particle_symbol( pdg, excited );
  }

  // Print the INTERACTION CHAIN section.
  void print_interaction_chain( std::ostream& os, const HepMC3::GenEvent& ev,
    marley::Reaction::ProcessType proc_type )
  {
    const auto& verts = ev.vertices();
    if ( verts.empty() ) return;

    // Primary vertex
    const auto& pv     = verts.front();
    const auto& pv_in  = pv->particles_in();
    const auto& pv_out = pv->particles_out();

    // Identify projectile and target from IN particles
    HepMC3::ConstGenParticlePtr proj = nullptr;
    HepMC3::ConstGenParticlePtr target = nullptr;
    if ( pv_in.size() == 2u ) {
      const auto& pv_in1 = pv_in.front();
      const auto& pv_in2 = pv_in.back();
      if ( pv_in1->status() == marley_hepmc3::NUHEPMC_PROJECTILE_STATUS ) {
        proj = pv_in1;
        target = pv_in2;
      }
      else {
        target = pv_in1;
        proj = pv_in2;
      }
    }
    else throw marley::Error( "Primary vertex without two incoming particles"
      " encountered in marley_hepmc3::print_event()" );

    // Identify ejectile and residue from OUT particles.
    HepMC3::ConstGenParticlePtr ejectile = nullptr;
    HepMC3::ConstGenParticlePtr residue  = nullptr;
    if ( pv_out.size() == 2u ) {
      const auto& pv_out1 = pv_out.front();
      const auto& pv_out2 = pv_out.back();
      // TODO: double-check that ES events look reasonable with this recipe
      if ( pv_out1->status()
        == marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS )
      {
        residue = pv_out1;
        ejectile = pv_out2;
      }
      else {
        ejectile = pv_out1;
        residue = pv_out2;
      }
    }
    else throw marley::Error( "Primary vertex without two outgoing particles"
      " encountered in marley_hepmc3::print_event()" );

    // 6-character chain tag
    std::string tag = chain_tag_primary( proc_type, proj->pid() );

    // Build primary chain line
    std::ostringstream pline;
    pline << "  " << chain_particle_name( proj )
      << " (" << std::fixed << std::setprecision( 2 )
      << proj->momentum().e() << " MeV)";
    pline << " + " << chain_particle_name( target );
    // ──[<tag>]──►
    pline << "  " << THIN_CHAR << THIN_CHAR << "[" << tag
          << "]" << THIN_CHAR << THIN_CHAR << ARROW_RIGHT << "  ";
    pline << chain_particle_name( ejectile );
    pline << "  +  " << chain_particle_name( residue );

    // Residue excitation bracket
    {
      double Ex = 0.; int twoJ = -1; int par = 1;
      bool has_ex = get_double_attr( residue, "Ex", Ex );
      get_int_attr( residue, "twoJ",   twoJ );
      get_int_attr( residue, "parity", par  );
      if ( has_ex && Ex > 0. ) {
        pline << "  [Ex=" << std::fixed << std::setprecision( 2 ) << Ex
              << " MeV, " << format_jp( twoJ, par ) << "]";
      }
    }
    os << pline.str() << "\n";

    // De-excitation chain lines
    for ( size_t vi = 1; vi < verts.size(); ++vi ) {
      const auto& vtx   = verts[vi];
      int vstatus = vtx->status();
      const auto& vin   = vtx->particles_in();
      const auto& vout  = vtx->particles_out();

      // Choose the decay tag (6 display chars)
      // "  γ   " uses γ = U+03B3 = 0xCE 0xB3 (1 display char)
      std::string dtag;
      if ( vstatus == marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX ) {
        dtag = "  HF  ";
      }
      else if ( vstatus == marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX ) {
        dtag = "  γ   ";
      }
      else {
        dtag = "  ??  ";
      }

      // IN: first particle is the decaying residue
      HepMC3::ConstGenParticlePtr dec_in = vin.front();

      // OUT: light fragment + heavy daughter nucleus
      HepMC3::ConstGenParticlePtr frag = nullptr;
      HepMC3::ConstGenParticlePtr daughter = nullptr;
      if ( vout.size() != 2u ) {
        throw marley::Error( "Binary decay vertex does not have exactly"
          " two outgoing particles" );
      }

      auto prod1 = vout.front();
      auto prod2 = vout.back();

      if ( prod1->generated_mass() <= prod2->generated_mass() ) {
        frag = prod1;
        daughter = prod2;
      }
      else {
        frag = prod2;
        daughter = prod1;
      }

      std::ostringstream dline;
      dline << "  ";
      if ( dec_in  ) dline << chain_particle_name( dec_in );
      dline << "  " << THIN_CHAR << THIN_CHAR << "[" << dtag
            << "]" << THIN_CHAR << THIN_CHAR << ARROW_RIGHT << "  ";
      if ( frag    ) dline << chain_particle_name( frag );
      dline << "  +  ";
      if ( daughter ) dline << chain_particle_name( daughter );

      // Daughter excitation bracket or ground-state label
      if ( daughter ) {
        double Ex = 0.; int twoJ = -1; int par = 1;
        bool has_ex = get_double_attr( daughter, "Ex", Ex );
        get_int_attr( daughter, "twoJ",   twoJ );
        get_int_attr( daughter, "parity", par  );
        if ( has_ex && Ex > 0. ) {
          dline << "  [Ex=" << std::fixed << std::setprecision( 2 ) << Ex
                << " MeV, " << format_jp( twoJ, par ) << "]";
        } else if ( daughter->pid() > 1000000000 ) {
          dline << "   [ground state]";
        }
      }
      os << dline.str() << "\n";
    }
  }

  // Compute the gamma-ray multipolarity string (e.g., "E1", "M2") by
  // inspecting the spin-parity attributes of the IN and OUT particles
  // of a gamma-cascade vertex. Returns an empty string when the needed
  // attributes are missing.
  //
  // The multipolarity ℓ is approximated as the lowest value allowed
  // by angular momentum conservation: ℓ = max(1, |J_i - J_f|).
  // The type (E or M) is determined from the parity rule:
  //   Electric if π_i = (-1)^ℓ × π_f, Magnetic otherwise.
  std::string compute_gamma_multipolarity_string(
    const HepMC3::ConstGenVertexPtr& vtx )
  {
    const auto& in_particles = vtx->particles_in();
    if ( in_particles.empty() ) return {};

    int twoJ_i = 0, Pi = 0;
    if ( !get_int_attr( in_particles.front(), "twoJ",   twoJ_i ) ) return {};
    if ( !get_int_attr( in_particles.front(), "parity", Pi     ) ) return {};

    // Find the daughter nucleus (non-photon) among OUT particles
    const auto& out_particles = vtx->particles_out();
    HepMC3::ConstGenParticlePtr daughter = nullptr;
    for ( const auto& p : out_particles ) {
      if ( p->pid() != 22 ) { daughter = p; break; }
    }
    if ( !daughter ) return {};

    int twoJ_f = 0, Pf = 0;
    if ( !get_int_attr( daughter, "twoJ",   twoJ_f ) ) return {};
    if ( !get_int_attr( daughter, "parity", Pf     ) ) return {};

    // ℓ = max(1, |J_i - J_f|).  The factor 1/2 converts twoJ to J.
    int ell = std::abs( twoJ_i - twoJ_f ) / 2;
    if ( ell < 1 ) ell = 1;

    int phase = ( ell % 2 == 0 ) ? 1 : -1; // (-1)^ℓ
    char type = ( Pi == phase * Pf ) ? 'E' : 'M';

    return type + std::to_string( ell );
  }

  // Print the full detail block for one vertex.
  void print_vertex_block( std::ostream& os,
    const HepMC3::ConstGenVertexPtr& vtx, int vtx_index )
  {
    int vstatus = vtx->status();
    bool is_deex  = ( vstatus == marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX
                   || vstatus == marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX );
    bool is_gamma = ( vstatus == marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX );

    // ── Vertex header ─────────────────────────────────────────────────────
    std::string type_name;
    if ( vstatus == marley_hepmc3::NUHEPMC_PRIMARY_VERTEX )
      type_name = "PRIMARY INTERACTION";
    else if ( vstatus == marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX )
      type_name = "HAUSER-FESHBACH";
    else if ( vstatus == marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX )
      type_name = "GAMMA CASCADE";
    else
      type_name = "VERTEX (status " + std::to_string( vstatus ) + ")";

    std::ostringstream hdr;
    hdr << "  [V" << vtx_index << "]  " << type_name;

    // Width attributes for Hauser-Feshbach decay vertices
    if ( vstatus == marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX ) {
      double width_tot = 0., width_ec = 0.;
      bool has_tot = get_double_attr_vtx( vtx, "TotalWidth", width_tot );
      bool has_ec  = get_double_attr_vtx( vtx, "ECWidth",    width_ec  );
      if ( has_tot && has_ec ) {
        hdr << "  Γ_ec = "
            << std::scientific << std::setprecision( 2 ) << width_ec
            << " MeV  Γ_tot = "
            << std::scientific << std::setprecision( 2 ) << width_tot
            << " MeV";
      }
    }
    // Additional info for gamma-cascade vertices
    else if ( is_gamma ) {
      {
        std::string xl = compute_gamma_multipolarity_string( vtx );
        if ( !xl.empty() ) hdr << "  " << xl;
      }
      double br = 0.;
      if ( get_double_attr_vtx( vtx, "GammaBranchingRatio", br ) ) {
        hdr << "  BR = " << std::scientific << std::setprecision( 2 ) << br;
      }
      double width_tot = 0.;
      if ( get_double_attr_vtx( vtx, "TotalWidth", width_tot ) ) {
        hdr << "  Γ_tot = "
            << std::scientific << std::setprecision( 2 ) << width_tot
            << " MeV";
      }
    }
    os << hdr.str() << "\n";

    // Retrieve multipolarity if this is a gamma-cascade vertex (used by
    // print_particle_line for the photon momentum display)
    int multi = -1;
    if ( is_gamma ) get_int_attr_vtx( vtx, "multipolarity", multi );

    // ── IN particles ──────────────────────────────────────────────────────
    os << make_thin_sep( THIN_CHAR + THIN_CHAR + " IN " ) << "\n";
    for ( const auto& p : vtx->particles_in() ) {
      print_particle_line( os, p,
        /*is_in_block=*/true, is_deex, is_gamma, multi );
    }

    // ── OUT particles ─────────────────────────────────────────────────────
    os << make_thin_sep( THIN_CHAR + THIN_CHAR + " OUT " ) << "\n";
    for ( const auto& p : vtx->particles_out() ) {
      print_particle_line( os, p,
        /*is_in_block=*/false, is_deex, is_gamma, multi );
    }
  }

} // end anonymous namespace (print_event helpers)

void marley_hepmc3::print_event( const HepMC3::GenEvent& ev,
  std::ostream& os )
{
  // Retrieve the process type for the input event
  marley::Reaction::ProcessType proc_type =
    marley::Reaction::ProcessType::Unknown;

  auto attr = ev.attribute< HepMC3::IntAttribute >( "signal_process_id" );
  if ( attr ) proc_type = marley_hepmc3::from_nuhepmc_proc_id( attr->value() );
  std::string proc_name = marley::Reaction::proc_type_to_string( proc_type );

  // Top thick separator
  os << THICK_SEP << "\n";

  // Header line: "  MARLEY  │  Event #N  │  <process>"
  os << "  MARLEY  " << BOX_VERT << "  Event #" << ev.event_number()
     << "  " << BOX_VERT << "  " << proc_name << "\n";

  // Thick separator
  os << THICK_SEP << "\n";

  // INTERACTION CHAIN
  os << "  INTERACTION CHAIN\n";
  print_interaction_chain( os, ev, proc_type );

  // Vertex detail blocks
  int vtx_idx = 1;
  for ( const auto& vtx : ev.vertices() ) {
    os << THICK_SEP << "\n";
    print_vertex_block( os, vtx, vtx_idx++ );
  }

  // ── FINAL STATE ───────────────────────────────────────────────────────
  os << THICK_SEP << "\n";
  os << "  FINAL STATE\n";
  for ( const auto& p : ev.particles() ) {
    if ( p->status() == marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS ) {
      // Final-state particles are ground-state; show full momentum.
      // Re-use print_particle_line with is_in_block=false, is_deex=false.
      print_particle_line( os, p,
        /*is_in_block=*/false,
        /*is_deex_vtx=*/false,
        /*is_gamma_vtx=*/false,
        /*vtx_multipolarity=*/-1 );
    }
  }

  // ── Cross-section footer ──────────────────────────────────────────────
  // σ = U+03C3 = 0xCF 0x83
  os << THICK_SEP << "\n";
  {
    auto proc_attr = ev.attribute< HepMC3::DoubleAttribute >( "proc_xs" );
    auto tot_attr  = ev.attribute< HepMC3::DoubleAttribute >( "tot_xs"  );
    double proc_xs = proc_attr ? proc_attr->value() : 0.;
    double tot_xs  = tot_attr  ? tot_attr->value()  : 0.;
    os << "  σ[" << proc_name << "] = "
       << std::scientific << std::setprecision( 3 ) << proc_xs
       << " pb     σ[tot] = "
       << std::scientific << std::setprecision( 3 ) << tot_xs << " pb\n";
  }

  // ── Bottom thick separator ────────────────────────────────────────────
  os << THICK_SEP << "\n";
}
