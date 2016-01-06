#ifndef USE_ROOT
#error "This MARLEY code must be compiled with ROOT support."
#endif

#include <fstream>
#include <iomanip>
#include <iostream>

#include "TMarleyEvent.hh"

#include "TFile.h"
#include "TTree.h"

// Function that dumps a TMarleyParticle to an output stream in HEPEVT format
void dump_particle(const TMarleyParticle& p, std::ostream& os,
  bool track = true)
{
  if (track) os << "1 ";
  else os << "0 ";

  os << p.get_id() << " 0 0 0 0 " << p.get_px() << " " << p.get_py() << " "
    << p.get_pz() << " " << p.get_total_energy() << " " << p.get_mass()
    << " " << " 0.0 0.0 0.0 0.0 " << std::endl;
}

// Script to dump events from a MARLEY event tree file (in ROOT format) to
// stdout in HEPEVT format
int main(){

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a TMarleyEvent pointer.
  TFile* file = new TFile("event_tree.root", "READ");
  TTree* t = nullptr;
  file->GetObject("MARLEY Event Tree", t);
  TMarleyEvent* e = new TMarleyEvent;
  t->GetBranch("events")->SetAddress(&e);

  // Use 16-digit precision for floating point numbers
  std::cout << std::setprecision(16) << std::scientific;

  // Write each event in the tree to stdout
  for (int i = 0, n = t->GetEntries(); i < n; ++i) {
    t->GetEntry(i);
    const auto& fp_list = e->get_final_particles();
    std::cout << i << " " << fp_list.size() + 1 << std::endl;
    dump_particle(*(e->get_projectile()), std::cout, false);
    for (const auto& fp : fp_list) dump_particle(fp, std::cout);
    //std::cout << *e << std::endl;
  }

  // Close the event tree ROOT file
  file->Close();
  delete file;

  return 0;
}
