#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifndef USE_ROOT
  #error Building the marsum executable requires ROOT
#endif

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Event.hh"
#include "marley/Particle.hh"

int main(int argc, char* argv[]) {

  // If the user has not supplied any command-line
  // arguments, display the standard help message
  // and exit
  if (argc <= 2) {
    std::cout << "Usage: " << argv[0] << " OUTPUT_FILE INPUT_FILE...\n";
    return 0;
  }

  // Temporary storage for output TTree branch variables
  double nu_E, nu_KE, nu_px, nu_py, nu_pz;
  int nu_pdg, target_pdg, num_final;
  std::vector<int> PDGs;
  std::vector<double> Es, KEs, pXs, pYs, pZs;
  double Ex;

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

  TFile out_tfile(argv[1], "recreate");
  TTree* out_tree = new TTree("summary_tree", "MARLEY summary tree");
  out_tree->Branch("nu_pdg", &nu_pdg, "nu_pdg/I");
  out_tree->Branch("target_pdg", &target_pdg, "target_pdg/I");
  out_tree->Branch("nu_E", &nu_E, "nu_E/D");
  out_tree->Branch("nu_KE", &nu_KE, "nu_KE/D");
  out_tree->Branch("nu_px", &nu_px, "nu_px/D");
  out_tree->Branch("nu_py", &nu_py, "nu_py/D");
  out_tree->Branch("nu_pz", &nu_pz, "nu_pz/D");
  out_tree->Branch("Ex", &Ex, "Ex/D");
  out_tree->Branch("num_final", &num_final, "num_final/I");
  out_tree->Branch("pdg", PDGs.data(), "pdg[num_final]/I");
  out_tree->Branch("E",  Es.data(), "E[num_final]/D");
  out_tree->Branch("KE", KEs.data(), "KE[num_final]/D");
  out_tree->Branch("px", pXs.data(), "px[num_final]/D");
  out_tree->Branch("py", pYs.data(), "py[num_final]/D");
  out_tree->Branch("pz", pZs.data(), "pz[num_final]/D");

  // Prepare to read the input file(s)
  marley::Event* ev = nullptr;
  TChain in_chain("MARLEY_event_tree");
  std::vector<std::string> input_file_names;
  for (int i = 2; i < argc; ++i) {
    in_chain.Add( argv[i] );
  }
  in_chain.SetBranchAddress("event", &ev);

  // Event loop
  // Don't use TChain::GetEntries() because it can be slow for large chains
  int entry = 0;
  while ( true ) {
    int local_entry = in_chain.LoadTree( entry );
    if ( local_entry < 0 ) break;
    in_chain.GetEntry( entry );

    if ( entry % 1000 == 0 ) std::cout << "Entry " << entry << '\n';

    PDGs.clear();
    Es.clear();
    KEs.clear();
    pXs.clear();
    pYs.clear();
    pZs.clear();

    nu_pdg = ev->projectile().pdg_code();
    nu_E = ev->projectile().total_energy();
    nu_KE = ev->projectile().kinetic_energy();
    nu_px = ev->projectile().px();
    nu_py = ev->projectile().py();
    nu_pz = ev->projectile().pz();

    target_pdg = ev->target().pdg_code();
    Ex = ev->Ex();

    num_final = ev->get_final_particles().size();
    for ( const auto& fp : ev->get_final_particles() ) {
      PDGs.push_back( fp->pdg_code() );
      Es.push_back( fp->total_energy() );
      KEs.push_back( fp->kinetic_energy() );
      pXs.push_back( fp->px() );
      pYs.push_back( fp->py() );
      pZs.push_back( fp->pz() );
    }

    // Update the branch addresses (manipulating the vectors may have
    // invalidated them)
    out_tree->SetBranchAddress("pdg", PDGs.data());
    out_tree->SetBranchAddress("E",  Es.data());
    out_tree->SetBranchAddress("KE", KEs.data());
    out_tree->SetBranchAddress("px", pXs.data());
    out_tree->SetBranchAddress("py", pYs.data());
    out_tree->SetBranchAddress("pz", pZs.data());

    out_tree->Fill();

    ++entry;
  }

  out_tree->Write();
  out_tfile.Close();
  return 0;
}
