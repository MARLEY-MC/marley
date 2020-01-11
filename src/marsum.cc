/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see \${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef USE_ROOT
  #error Building the marsum executable requires ROOT
#endif

// ROOT includes
#include "TFile.h"
#include "TTree.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Event.hh"
#include "marley/RootEventFileReader.hh"
#include "marley/Particle.hh"

namespace {
  // Index of the first final-state particle that is not the
  // ejectile or residue (first de-excitation product).
  // TODO: Make this easier to maintain. If you change the layout of
  // marley::Event, this will break.
  constexpr size_t FIRST_PROD_IDX = 2u;
}

int main(int argc, char* argv[]) {

  // If the user has not supplied any command-line arguments, display the
  // standard help message and exit
  if (argc <= 2) {
    std::cout << "Usage: " << argv[0] << " OUTPUT_FILE INPUT_FILE...\n";
    return 0;
  }

  // Temporary storage for output TTree branch variables
  double Ex; // nuclear excitation energy
  double flux_avg_tot_xsec; // flux-averaged total cross section
  double Ev, KEv, pxv, pyv, pzv; // projectile
  double Mt; // target mass
  double El, KEl, pxl, pyl, pzl; // ejectile
  double Er, KEr, pxr, pyr, pzr; // residue (after de-excitations)
  int pdgv, pdgt, pdgl, pdgr; // PDG codes
  int np; // number of de-excitation products (final-state particles other
          // than the ejectile and residue)

  // Inpormation about each of the other final-state particles
  std::vector<int> PDGs;
  std::vector<double> Es, KEs, pXs, pYs, pZs;

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
  TTree* out_tree = new TTree("mst", "MARLEY summary tree");

  // projectile branches
  out_tree->Branch("pdgv", &pdgv, "pdgv/I");
  out_tree->Branch("Ev", &Ev, "Ev/D");
  out_tree->Branch("KEv", &KEv, "KEv/D");
  out_tree->Branch("pxv", &pxv, "pxv/D");
  out_tree->Branch("pyv", &pyv, "pyv/D");
  out_tree->Branch("pzv", &pzv, "pzv/D");

  // target branches
  out_tree->Branch("pdgt", &pdgt, "pdgt/I");
  out_tree->Branch("Mt", &Mt, "Mt/D");

  // ejectile branches
  out_tree->Branch("pdgl", &pdgl, "pdgl/I");
  out_tree->Branch("El", &El, "El/D");
  out_tree->Branch("KEl", &KEl, "KEl/D");
  out_tree->Branch("pxl", &pxl, "pxl/D");
  out_tree->Branch("pyl", &pyl, "pyl/D");
  out_tree->Branch("pzl", &pzl, "pzl/D");

  // residue branches
  out_tree->Branch("pdgr", &pdgr, "pdgr/I");
  out_tree->Branch("Er", &Er, "Er/D");
  out_tree->Branch("KEr", &KEr, "KEr/D");
  out_tree->Branch("pxr", &pxr, "pxr/D");
  out_tree->Branch("pyr", &pyr, "pyr/D");
  out_tree->Branch("pzr", &pzr, "pzr/D");

  // Nuclear excitation energy branch
  out_tree->Branch("Ex", &Ex, "Ex/D");

  // De-excitation products (final-state particles other than the
  // ejectile and ground-state residue)
  out_tree->Branch("np", &np, "np/I");
  out_tree->Branch("pdgp", PDGs.data(), "pdgp[np]/I");
  out_tree->Branch("Ep",  Es.data(), "Ep[np]/D");
  out_tree->Branch("KEp", KEs.data(), "KEp[np]/D");
  out_tree->Branch("pxp", pXs.data(), "pxp[np]/D");
  out_tree->Branch("pyp", pYs.data(), "pyp[np]/D");
  out_tree->Branch("pzp", pZs.data(), "pzp[np]/D");

  // Flux-averaged total cross section
  out_tree->Branch("xsec", &flux_avg_tot_xsec, "xsec/D");

  // Prepare to read the input file(s)
  std::vector<std::string> input_file_names;
  for ( int i = 2; i < argc; ++i ) input_file_names.push_back( argv[i] );

  // File loop
  for ( const auto& file_name : input_file_names ) {

    // Open the current file for reading
    marley::RootEventFileReader refr( file_name );
    std::cout << "Opened file \"" << file_name << "\"\n";

    // Temporary object to use for reading in saved events
    marley::Event ev;

    // Event loop
    int event_num = 0;
    while ( refr >> ev ) {

      if ( event_num % 1000 == 0 ) std::cout << "Event " << event_num << '\n';

      PDGs.clear();
      Es.clear();
      KEs.clear();
      pXs.clear();
      pYs.clear();
      pZs.clear();

      pdgv = ev.projectile().pdg_code();
      Ev = ev.projectile().total_energy();
      KEv = ev.projectile().kinetic_energy();
      pxv = ev.projectile().px();
      pyv = ev.projectile().py();
      pzv = ev.projectile().pz();

      pdgt = ev.target().pdg_code();
      Mt = ev.target().mass();

      pdgl = ev.ejectile().pdg_code();
      El = ev.ejectile().total_energy();
      KEl = ev.ejectile().kinetic_energy();
      pxl = ev.ejectile().px();
      pyl = ev.ejectile().py();
      pzl = ev.ejectile().pz();

      pdgr = ev.residue().pdg_code();
      Er = ev.residue().total_energy();
      KEr = ev.residue().kinetic_energy();
      pxr = ev.residue().px();
      pyr = ev.residue().py();
      pzr = ev.residue().pz();

      Ex = ev.Ex();
      flux_avg_tot_xsec = refr.flux_averaged_xsec();

      const auto& fparts = ev.get_final_particles();
      np = fparts.size() - FIRST_PROD_IDX;
      for ( size_t j = FIRST_PROD_IDX; j < fparts.size(); ++j ) {
        const auto* fp = fparts.at( j );
        PDGs.push_back( fp->pdg_code() );
        Es.push_back( fp->total_energy() );
        KEs.push_back( fp->kinetic_energy() );
        pXs.push_back( fp->px() );
        pYs.push_back( fp->py() );
        pZs.push_back( fp->pz() );
      }

      // Update the branch addresses (manipulating the vectors may have
      // invalidated them)
      out_tree->SetBranchAddress("pdgp", PDGs.data());
      out_tree->SetBranchAddress("Ep",  Es.data());
      out_tree->SetBranchAddress("KEp", KEs.data());
      out_tree->SetBranchAddress("pxp", pXs.data());
      out_tree->SetBranchAddress("pyp", pYs.data());
      out_tree->SetBranchAddress("pzp", pZs.data());

      out_tree->Fill();

      ++event_num;
    } // event loop
  } // file loop

  out_tfile.cd();
  out_tree->Write();
  out_tfile.Close();
  return 0;
}
