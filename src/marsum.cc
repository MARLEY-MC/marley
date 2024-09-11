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

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef USE_ROOT
  #error Building the marsum executable requires ROOT
#endif

// HepMC3 includes
#include "HepMC3/FourVector.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// MARLEY includes
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"
#include "marley/EventFileReader.hh"

int main( int argc, char* argv[] ) {

  // If the user has not supplied enough command-line arguments, display the
  // standard help message and exit
  if ( argc <= 2 ) {
    std::cout << "Usage: " << argv[0] << " OUTPUT_FILE INPUT_FILE...\n";
    return 0;
  }

  // Temporary storage for output TTree branch variables
  double Ex; // nuclear excitation energy
  int twoJ; // two times the residue spin immediately after the two-two
            // scattering reaction
  int par; // integer representation of the intrinsic parity of the
           // residue immediately following the two-two reaction
  double flux_avg_tot_xsec; // flux-averaged total cross section
  double Ev, KEv, pxv, pyv, pzv; // projectile
  double Mt; // target mass
  double El, KEl, pxl, pyl, pzl; // ejectile
  double Er, KEr, pxr, pyr, pzr; // residue (after de-excitations)
  int pdgv, pdgt, pdgl, pdgr; // PDG codes
  int np; // number of de-excitation products (final-state particles other
          // than the ejectile and residue)

  // Information about each of the other final-state particles
  std::vector<int> PDGs;
  std::vector<double> Es, KEs, pXs, pYs, pZs;

  // Event weights
  double cv_weight; // Central-value weight (normally unity, required by
                    // NuHepMC standard)
  std::vector< double > other_weights; // Other user-defined weights

  // Check whether the output file exists and warn the user before
  // overwriting it if it does
  std::ifstream temp_stream( argv[1] );
  if ( temp_stream ) {
    bool overwrite = marley_utils::prompt_yes_no(
      "Really overwrite " + std::string(argv[1]) + '?');
    if ( !overwrite ) {
      std::cout << "Action aborted.\n";
      return 0;
    }
  }

  TFile out_tfile( argv[1], "recreate" );
  TTree* out_tree = new TTree( "mst", "MARLEY summary tree" );

  // projectile branches
  out_tree->Branch( "pdgv", &pdgv, "pdgv/I" );
  out_tree->Branch( "Ev", &Ev, "Ev/D" );
  out_tree->Branch( "KEv", &KEv, "KEv/D" );
  out_tree->Branch( "pxv", &pxv, "pxv/D" );
  out_tree->Branch( "pyv", &pyv, "pyv/D" );
  out_tree->Branch( "pzv", &pzv, "pzv/D" );

  // target branches
  out_tree->Branch( "pdgt", &pdgt, "pdgt/I" );
  out_tree->Branch( "Mt", &Mt, "Mt/D" );

  // ejectile branches
  out_tree->Branch( "pdgl", &pdgl, "pdgl/I" );
  out_tree->Branch( "El", &El, "El/D" );
  out_tree->Branch( "KEl", &KEl, "KEl/D" );
  out_tree->Branch( "pxl", &pxl, "pxl/D" );
  out_tree->Branch( "pyl", &pyl, "pyl/D" );
  out_tree->Branch( "pzl", &pzl, "pzl/D" );

  // residue branches
  out_tree->Branch( "pdgr", &pdgr, "pdgr/I" );
  out_tree->Branch( "Er", &Er, "Er/D" );
  out_tree->Branch( "KEr", &KEr, "KEr/D" );
  out_tree->Branch( "pxr", &pxr, "pxr/D" );
  out_tree->Branch( "pyr", &pyr, "pyr/D" );
  out_tree->Branch( "pzr", &pzr, "pzr/D" );

  // Nuclear excitation energy branch
  out_tree->Branch( "Ex", &Ex, "Ex/D" );

  // Spin and parity branches
  out_tree->Branch( "twoJ", &twoJ, "twoJ/I" );
  out_tree->Branch( "parity", &par, "parity/I" );

  // De-excitation products (final-state particles other than the
  // ejectile and ground-state residue)
  out_tree->Branch( "np", &np, "np/I" );
  out_tree->Branch( "pdgp", &PDGs );
  out_tree->Branch( "Ep",  &Es );
  out_tree->Branch( "KEp", &KEs );
  out_tree->Branch( "pxp", &pXs );
  out_tree->Branch( "pyp", &pYs );
  out_tree->Branch( "pzp", &pZs );

  // Flux-averaged total cross section
  out_tree->Branch( "xsec", &flux_avg_tot_xsec, "xsec/D" );

  // Event weights
  out_tree->Branch( "cv_weight", &cv_weight, "cv_weight/D" );
  out_tree->Branch( "other_weights", &other_weights );

  // Prepare to read the input file(s)
  std::vector<std::string> input_file_names;
  for ( int i = 2; i < argc; ++i ) input_file_names.push_back( argv[i] );

  // File loop
  for ( const auto& file_name : input_file_names ) {

    // Open the current file for reading
    marley::EventFileReader efr( file_name );
    std::cout << "Opened file \"" << file_name << "\"\n";

    // Temporary object to use for reading in saved events
    HepMC3::GenEvent ev;

    // Stores the number of elements expected in the other_weights vector
    // (determined on the first iteration of the event loop below)
    int num_other_weights = 0;

    // Event loop
    int event_num = 0;
    while ( efr >> ev ) {

      // Write the vector of custom weight names to the output file on the
      // first iteration of the event loop. We only need to do this once since
      // the list is stored in the run information and is common to all events.
      if ( event_num == 0 ) {
        // TODO: add error handling for when the number of weights changes,
        // indicating that the events have inconsistent run information
        auto run_info = ev.run_info();
        auto wgt_names = run_info->weight_names();

        // Drop the first weight name since it will always be the central-value
        // weight, which is stored separately in the output TTree
        wgt_names.erase( wgt_names.begin() );

        // Store the names in the output TFile
        out_tfile.WriteObject( &wgt_names,
          "MARLEY_other_weight_names", "WriteDelete" );

        // Record the expected number of other (i.e., non-CV) weights
        num_other_weights = wgt_names.size();
      }

      if ( event_num % 1000 == 0 ) std::cout << "Event " << event_num << '\n';

      PDGs.clear();
      Es.clear();
      KEs.clear();
      pXs.clear();
      pYs.clear();
      pZs.clear();

      auto projectile = marley_hepmc3::get_projectile( ev );
      pdgv = projectile->pid();

      double mv = projectile->generated_mass();
      const HepMC3::FourVector& p4v = projectile->momentum();

      Ev = p4v.e();
      KEv = std::max( 0., Ev - mv );
      pxv = p4v.px();
      pyv = p4v.py();
      pzv = p4v.pz();

      auto target = marley_hepmc3::get_target( ev );
      pdgt = target->pid();
      Mt = target->generated_mass();

      auto ejectile = marley_hepmc3::get_ejectile( ev );
      pdgl = ejectile->pid();

      // Object ID number for the ejectile in the event (*not* the same as the
      // PDG code)
      int ej_id = ejectile->id();

      double ml = ejectile->generated_mass();
      const HepMC3::FourVector& p4l = ejectile->momentum();

      El = p4l.e();
      KEl = std::max( 0., El - ml );
      pxl = p4l.px();
      pyl = p4l.py();
      pzl = p4l.pz();

      auto residue = marley_hepmc3::get_residue( ev );
      pdgr = residue->pid();

      double mr = residue->generated_mass();
      const HepMC3::FourVector& p4r = residue->momentum();

      Er = p4r.e();
      KEr = std::max( 0., Er - mr );
      pxr = p4r.px();
      pyr = p4r.py();
      pzr = p4r.pz();

      // TODO: add error handling here for missing attributes
      auto Ex_attr = residue->attribute< HepMC3::DoubleAttribute >( "Ex" );
      Ex = Ex_attr->value();

      auto twoJ_attr = residue->attribute< HepMC3::IntAttribute >( "twoJ" );
      twoJ = twoJ_attr->value();

      auto parity_attr = residue->attribute< HepMC3::IntAttribute >( "parity" );
      par = parity_attr->value();

      flux_avg_tot_xsec = efr.flux_averaged_xsec();

      np = 0;
      const auto& particles = ev.particles();
      for ( const auto& p : particles ) {
        // Only save information about final-state particles in the array of
        // nuclear de-excitation products
        if ( p->status() != marley_hepmc3::NUHEPMC_FINAL_STATE_STATUS ) {
          continue;
        }

        // Skip the ejectile since its information is already saved elsewhere in
        // the output branches. We use the unique object ID number in the event
        // (not the same as the PDG code) to check for this
        if ( p->id() == ej_id ) continue;

        ++np;

        PDGs.push_back( p->pid() );

        double mp = p->generated_mass();
        const HepMC3::FourVector& p4p = p->momentum();

        double Ep = p4p.e();
        Es.push_back( Ep );

        double KEp = std::max( 0., Ep - mp );
        KEs.push_back( KEp );

        pXs.push_back( p4p.px() );
        pYs.push_back( p4p.py() );
        pZs.push_back( p4p.pz() );
      }

      // Make a copy of the vector of weights for the current event
      other_weights = ev.weights();

      // Drop the first element and store it in the central-value weight
      // instead
      cv_weight = other_weights.front();
      other_weights.erase( other_weights.begin() );

      // All variables are ready. Fill the output TTree and advance to the next
      // input event
      out_tree->Fill();

      ++event_num;

    } // event loop
  } // file loop

  out_tfile.cd();
  out_tree->Write();
  out_tfile.Close();
  return 0;
}
