#ifndef USE_ROOT
#error "This MARLEY validation code currently must be compiled with ROOT support."
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyEvent.hh"

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

int main(){

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a TMarleyEvent pointer.
  TFile* file = new TFile("event_tree.root", "READ");
  TTree* t = nullptr;
  file->GetObject("MARLEY Event Tree", t);
  TMarleyEvent* e = new TMarleyEvent;
  t->GetBranch("events")->SetAddress(&e);

  // Change in energy between the initial and final states for each event
  std::vector<double> delta_Es;

  // Cycle through each event in the tree
  for (int i = 0, n = t->GetEntries(); i < n; ++i) {

    // Load the current event
    t->GetEntry(i);
    //std::cout << std::endl << "Event " << i << std::endl;

    // Compute the sum of all of the initial particle total energies
    // for the current event
    double Ei = 0.0;
    std::list<TMarleyParticle>* iparts = e->get_initial_particles();
    for(std::list<TMarleyParticle>::iterator it = iparts->begin();
      it != iparts->end(); ++it)
    {
      Ei += it->get_total_energy();
    }

    // Do the same thing for the final particles
    double Ef = 0.0;
    std::list<TMarleyParticle>* fparts = e->get_final_particles();
    for(std::list<TMarleyParticle>::iterator it = fparts->begin();
      it != fparts->end(); ++it)
    {
      Ef += it->get_total_energy();
      //std::cout << "Final particle: ID = " << it->get_id() << std::endl;
    }

    double change_in_E = Ef - Ei;

    delta_Es.push_back(change_in_E);
  }

  // Close the event tree ROOT file
  file->Close();
  delete file;

  // Find the minimum and maximum changes in the total energy of an event
  double max_change_in_E = *(std::max_element(delta_Es.begin(), delta_Es.end()));
  double min_change_in_E = *(std::min_element(delta_Es.begin(), delta_Es.end()));

  std::cout << "Change in energy ranges from " << min_change_in_E
    << " to " << max_change_in_E << std::endl;

  // Create a histogram to store the change in total energy for each event.
  // Use 100 bins, and adjust the upper limit of the last bin so that
  // the maximum observed value just barely falls within it.
  TH1D delta_E_histogram("delta_E_histogram",
    "Conservation of Energy in MARLEY Events; #Delta E [MeV]; Counts",
    100, min_change_in_E, 1.2*std::nextafter(max_change_in_E, DBL_MAX));

  // Fill the energy change histogram with values
  for(double dE : delta_Es) delta_E_histogram.Fill(dE);

  // Write the energy change histogram to new a ROOT file
  file = new TFile("validation.root","RECREATE");
  delta_E_histogram.Write();
  file->Close();
  delete file;

  return 0;
}
