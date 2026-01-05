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

// Standard library includes
#include <algorithm>

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/hepmc3_utils.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/JSON.hh"
#include "marley/Logger.hh"
#include "marley/OMPWeightCalculator.hh"
#include "marley/StructureDatabase.hh"

namespace {
  constexpr double PRETTY_SMALL = 1e-6;
}

marley::OMPWeightCalculator::OMPWeightCalculator( const marley::JSON& config )
  : marley::WeightCalculator( config )
{
  if ( !config.is_object() ) {
    throw marley::Error( "Non-object JSON configuration passed to constructor"
      " of marley::OMPWeightCalculator" );
  }

  if ( !config.has_key("opt_mod") ) {
    throw marley::Error( "Missing \"opt_mod\" key in marley::OMPWeight"
      "Calculator JSON configuration" );
  }

  // Load the JSON configuration for the optical model parameters
  const auto& om_config = config.at( "opt_mod" );
  if ( !om_config.is_object() ) {
    throw marley::Error( "Non-object associated with \"opt_mod\" key in"
      " marley::OMPWeightCalculator JSON configuration" );
  }

  // Create a structure database with a custom set of optical model parameters
  sdb_ = std::make_shared< marley::StructureDatabase >();
  sdb_->load_optical_model_params( &om_config );
}

double marley::OMPWeightCalculator::weight( HepMC3::GenEvent& ev,
  marley::Generator& /*gen*/ ) const
{
  // Start with a result of unity
  double weight = 1.;

  // Get a vector of pointers to all decay vertices in the event that
  // were handled using the Hauser-Feshbach model
  auto hf_vtx_vec = marley_hepmc3::get_vertices_with_status(
    marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX, ev );

  // Include each vertex that was found in the reweighting calculation
  for ( const auto& vtx : hf_vtx_vec ) {

    // Start with the default assumption that the decay was to a discrete
    // nuclear level
    bool decayed_to_continuum = false;

    const auto& in_vec = vtx->particles_in();
    const auto& out_vec = vtx->particles_out();

    // Double-check that the expected numbers of particles are present
    if ( in_vec.size() != 1u || out_vec.size() != 2u ) {
    throw marley::Error( "Hauser-Feshbach vertex encountered that"
      " is not a binary decay" );
    }

    // Retrieve the decay widths computed in the original simulation
    auto width_tot_attr = vtx->attribute< HepMC3::DoubleAttribute >(
      "TotalWidth" );
    double width_tot = width_tot_attr->value();

    auto width_ec_attr = vtx->attribute< HepMC3::DoubleAttribute >(
      "ECWidth" );
    double width_ec = width_ec_attr->value();

    // This attribute is only present for decays to the continuum
    auto width_sp_attr = vtx->attribute< HepMC3::DoubleAttribute >(
      "SPWidth" );
    double width_sp = 0.;
    if ( width_sp_attr ) {
      decayed_to_continuum = true;
      width_sp = width_sp_attr->value();
    }

    // Get access to the mother nucleus for the current binary decay
    const auto& mother = in_vec.front();

    // Retrieve the starting nuclear excitation energy, spin, and parity
    // TODO: add error handling here for missing attributes
    auto Exi_attr = mother->attribute< HepMC3::DoubleAttribute >( "Ex" );
    double Exi = Exi_attr->value();

    auto twoJi_attr = mother->attribute< HepMC3::IntAttribute >( "twoJ" );
    int twoJi = twoJi_attr->value();

    auto Pi_attr = mother->attribute< HepMC3::IntAttribute >( "parity" );
    marley::Parity Pi( Pi_attr->value() );

    // NucleusDecayer populates the outgoing particles in the decay
    // vertex in a specific order. The first is the emitted nuclear
    // fragment or gamma-ray, while the second is the daughter nucleus.
    const auto& emitted_particle = out_vec.front();
    const auto& daughter = out_vec.back();

    int emitted_pdg = emitted_particle->pid();

    // Retrieve the final nuclear excitation energy, spin, and parity
    // TODO: add error handling here for missing attributes
    auto Exf_attr = daughter->attribute< HepMC3::DoubleAttribute >( "Ex" );
    double Exf = Exf_attr->value();

    auto twoJf_attr = daughter->attribute< HepMC3::IntAttribute >( "twoJ" );
    int twoJf = twoJf_attr->value();

    auto Pf_attr = daughter->attribute< HepMC3::IntAttribute >( "parity" );
    marley::Parity Pf( Pf_attr->value() );

    // Construct a set of decay widths calculated with alternative settings
    marley::HauserFeshbachDecay hf_alt( mother, Exi, twoJi, Pi, *sdb_ );

    // Get the total decay width under the alternative calculation
    double width_tot_alt = hf_alt.total_width();

    // Find the ExitChannel in the alternative calculation corresponding
    // to the original decay that was sampled
    const auto& ec_vec = hf_alt.exit_channels();
    auto ec_iter = std::find_if( ec_vec.cbegin(), ec_vec.cend(),
    [ emitted_pdg, decayed_to_continuum, Exf ](
      const std::unique_ptr< marley::ExitChannel >& test_ec ) -> bool
    {
      if ( emitted_pdg != test_ec->emitted_particle_pdg() ) return false;

      if ( decayed_to_continuum ) {
      if ( !test_ec->is_continuum() ) return false;
      else return true;
      }

      // If we get to here, then we're dealing with a discrete transition
      if ( test_ec->is_continuum() ) return false;

      const auto* dec = dynamic_cast<
      const marley::DiscreteExitChannel* >( test_ec.get() );
      if ( !dec ) return false;

      // Check the excitation energy of the final discrete level as
      // a last confirmation that we've found the correct ExitChannel
      const auto& lev = dec->get_final_level();
      double Ex_level = lev.energy();

      // To deal with possible numerical precision issues, this check of
      // the final excitation energy allows for small differences
      if ( std::abs(Exf - Ex_level) > PRETTY_SMALL ) return false;
      return true;
    }
    );

    if ( ec_iter == ec_vec.cend() ) {
      MARLEY_LOG_WARNING() << "Could not find ExitChannel"
        << " during reweighting";
      weight = 0.;
      continue;
    }

    // Partial width for the chosen ExitChannel under the alternative
    // calculation
    double width_ec_alt = ( *ec_iter )->width();

    // Partial differential width for the chosen spin-parity under the
    // alternative calculation (applies only to decays to the continuum)
    double width_sp_alt = 0.;
    if ( decayed_to_continuum ) {
      // Evaluate the differential width at the sampled excitation energy.
      // This will populate the table of SpinParityWidth objects with the
      // correct values.
      const auto& cec = dynamic_cast<
        const marley::ContinuumExitChannel& >( *(*ec_iter) );
      cec.differential_width( Exf, true );

      // Retrieve extra information based on the kind of emitted particle
      int mpol, two_j_frag, orb_l;
      bool emitted_gamma = false;
      auto mpol_attr = vtx->attribute< HepMC3::IntAttribute >(
        "multipolarity" );
      if ( mpol_attr ) {
        emitted_gamma = true;
        mpol = mpol_attr->value();
      }
      else {
        auto two_j_frag_attr = vtx->attribute< HepMC3::IntAttribute >(
          "two_j_frag" );
        auto orb_l_attr = vtx->attribute< HepMC3::IntAttribute >( "orb_l" );

        // TODO: add error handling for missing attributes here

        two_j_frag = two_j_frag_attr->value();
        orb_l = orb_l_attr->value();
      }

      // Find the SpinParityWidth object corresponding to the spin-parity
      // value that was actually sampled
      const auto& spw_vec = cec.get_spw_table();
      auto spw_iter = std::find_if( spw_vec.cbegin(), spw_vec.cend(),
        [ twoJf, Pf, emitted_gamma, mpol, two_j_frag, orb_l ](
          const std::unique_ptr< marley::ContinuumExitChannel
          ::SpinParityWidth >& spw ) -> bool
        {
          if ( Pf != spw->Pf ) return false;
          if ( twoJf != spw->twoJf ) return false;
          if ( emitted_gamma ) {
            const auto* g_spw = static_cast< const marley
              ::GammaContinuumExitChannel::GammaSpinParityWidth* >(
              spw.get() );
            if ( !g_spw ) return false;
            if ( mpol != g_spw->multipolarity ) return false;
          }
          else {
            // emitted fragment
            const auto* f_spw = static_cast< const marley
              ::FragmentContinuumExitChannel::FragmentSpinParityWidth* >(
              spw.get() );
            if ( !f_spw ) return false;
            if ( two_j_frag != f_spw->two_j_frag ) return false;
            if ( orb_l != f_spw->orb_l ) return false;
          }
          return true;
        }
      );

      if ( spw_iter == spw_vec.cend() ) {
        MARLEY_LOG_WARNING() << "Could not find SpinParityWidth during"
          << " reweighting";
        weight = 0.;
        continue;
      }

      // Store the partial differential width to this spin-parity state
      // under the alternative calculation
      width_sp_alt = spw_iter->get()->diff_width;
    }

    // We've accumulated all the width values we need to compute an
    // event weight. Do some sanity checks beforehand.
    std::string bad_width_name;
    if ( width_tot_alt <= 0. ) {
      bad_width_name = "Alternate total_width";
    }
    else if ( width_tot <= 0. ) {
      bad_width_name = "Original total width";
    }
    else if ( width_ec <= 0. ) {
      bad_width_name = "Original exit channel";
    }
    else if ( width_ec_alt <= 0. ) {
      bad_width_name = "Alternate exit channel";
    }
    else if ( decayed_to_continuum ) {
      if ( width_sp_alt <= 0. ) bad_width_name = "Alternate spin-parity";
      else if ( width_sp <= 0. ) bad_width_name = "Original spin-parity";
    }

    if ( !bad_width_name.empty() ) {
      MARLEY_LOG_WARNING() << bad_width_name << " is non-positive";
      weight = 0.;
      continue;
    }

    // Everything looks okay, so do the weight calculation using the
    // appropriate expression for a decay to the continuum or a discrete
    // nuclear level
    double  w = width_tot / width_tot_alt;
    if ( decayed_to_continuum ) {
      w *= width_sp_alt / width_sp;
    }
    else {
      w *= width_ec_alt / width_ec;
    }

    // Multiply the weight for the current decay vertex into the overall
    // event weight
    weight *= w;

  } // decay vertex loop

  return weight;
}
