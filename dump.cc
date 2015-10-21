#ifndef USE_ROOT
#error "This MARLEY code must be compiled with ROOT support."
#endif

#include <iomanip>
#include <iostream>

#include "TMarleyEvent.hh"

#include "TFile.h"
#include "TTree.h"

// Script to dump events from a MARLEY event tree file (in ROOT format) to stdout
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
    std::cout << *e << std::endl;
  }

  // Close the event tree ROOT file
  file->Close();
  delete file;

  return 0;
}
