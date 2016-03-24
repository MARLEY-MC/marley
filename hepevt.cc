#ifndef USE_ROOT
#error "This MARLEY code must be compiled with ROOT support."
#endif

#include <fstream>
#include <iomanip>
#include <iostream>

#include "Event.hh"

#include "TFile.h"
#include "TTree.h"

// Script to dump events from a MARLEY event tree file (in ROOT format) to
// stdout in HEPEvt format
int main(){

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a marley::Event pointer.
  TFile* file = new TFile("events.root", "READ");
  TTree* t = nullptr;
  file->GetObject("MARLEY Event Tree", t);
  marley::Event* e = new marley::Event;
  t->GetBranch("events")->SetAddress(&e);

  // Use 16-digit precision for floating point numbers
  std::cout << std::setprecision(16) << std::scientific;

  // Write each event in the tree to stdout
  for (int i = 0, n = t->GetEntries(); i < n; ++i) {
    t->GetEntry(i);
    e->write_hepevt(i, std::cout);
  }

  // Close the event tree ROOT file
  file->Close();
  delete file;

  return 0;
}
