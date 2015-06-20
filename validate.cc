#ifndef USE_ROOT
#error "This MARLEY validation code currently must be compiled with ROOT support."
#endif

#include <iostream>

#include "TMarleyEvent.hh"

#include "TFile.h"
#include "TTree.h"

int main(){

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a TMarleyEvent pointer.
  TFile file("event_tree.root");
  TTree* t = static_cast<TTree*>(file.Get("event_tree"));
  TMarleyEvent* e = new TMarleyEvent;
  t->GetBranch("events")->SetAddress(&e);

  // Cycle through each event in the tree
  for (int i = 0, n = t->GetEntries(); i < n; ++i) {

    // Load the current event
    t->GetEntry(i);

    // Put event analysis code here
    std::cout << e->get_ejectile()->get_total_energy() << std::endl;
  }

  return 0;
}
