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
  // Changes in 3-momentum components for each event
  std::vector<double> delta_pxs;
  std::vector<double> delta_pys;
  std::vector<double> delta_pzs;

  // Cycle through each event in the tree
  for (int i = 0, n = t->GetEntries(); i < n; ++i) {

    // Load the current event
    t->GetEntry(i);
    //std::cout << std::endl << "Event " << i << std::endl;

    // Compute the sum of all of the initial particle total energies
    // for the current event
    double Ei = 0.0;
    // Compute the total 3-momentum components as well
    double pxi = 0., pyi = 0., pzi = 0.;
    std::list<TMarleyParticle>* iparts = e->get_initial_particles();
    for(std::list<TMarleyParticle>::iterator it = iparts->begin();
      it != iparts->end(); ++it)
    {
      Ei += it->get_total_energy();
      pxi += it->get_px();
      pyi += it->get_py();
      pzi += it->get_pz();
    }

    // Do the same thing for the final particles
    double Ef = 0.0;
    double pxf = 0., pyf = 0., pzf = 0.;
    std::list<TMarleyParticle>* fparts = e->get_final_particles();
    for(std::list<TMarleyParticle>::iterator it = fparts->begin();
      it != fparts->end(); ++it)
    {
      //std::cout << "Final particle: ID = " << it->get_id() << std::endl;
      Ef += it->get_total_energy();
      pxf += it->get_px();
      pyf += it->get_py();
      pzf += it->get_pz();

    }

    double change_in_E = Ef - Ei;
    double change_in_px = pxf - pxi;
    double change_in_py = pyf - pyi;
    double change_in_pz = pzf - pzi;

    delta_Es.push_back(change_in_E);
    delta_pxs.push_back(change_in_px);
    delta_pys.push_back(change_in_py);
    delta_pzs.push_back(change_in_pz);
  }

  // Close the event tree ROOT file
  file->Close();
  delete file;

  // Find the minimum and maximum changes in the total energy of an event
  double max_change_in_E = *(std::max_element(delta_Es.begin(), delta_Es.end()));
  double min_change_in_E = *(std::min_element(delta_Es.begin(), delta_Es.end()));
  // Do the same for the 3-momentum components
  double max_change_in_px = *(std::max_element(delta_pxs.begin(), delta_pxs.end()));
  double min_change_in_px = *(std::min_element(delta_pxs.begin(), delta_pxs.end()));
  double max_change_in_py = *(std::max_element(delta_pys.begin(), delta_pys.end()));
  double min_change_in_py = *(std::min_element(delta_pys.begin(), delta_pys.end()));
  double max_change_in_pz = *(std::max_element(delta_pzs.begin(), delta_pzs.end()));
  double min_change_in_pz = *(std::min_element(delta_pzs.begin(), delta_pzs.end()));

  std::cout << "Change in energy ranges from " << min_change_in_E
    << " to " << max_change_in_E << std::endl;
  std::cout << "Change in px ranges from " << min_change_in_px
    << " to " << max_change_in_px << std::endl;
  std::cout << "Change in py ranges from " << min_change_in_py
    << " to " << max_change_in_py << std::endl;
  std::cout << "Change in pz ranges from " << min_change_in_pz
    << " to " << max_change_in_pz << std::endl;

  // Create a histogram to store the change in total energy for each event.
  // Use 100 bins, and adjust the upper limit of the last bin so that
  // the maximum observed value just barely falls within it.
  TH1D delta_E_histogram("delta_E_histogram",
    "Conservation of Energy in MARLEY Events; #Delta E [MeV]; Counts",
    100, min_change_in_E, 1.2*std::nextafter(max_change_in_E, DBL_MAX));
  // Do the same for the 3-momentum components
  TH1D delta_px_histogram("delta_px_histogram",
    "Conservation of Momentum in MARLEY Events; #Delta p_x [MeV]; Counts",
    100, min_change_in_px, 1.2*std::nextafter(max_change_in_px, DBL_MAX));
  TH1D delta_py_histogram("delta_py_histogram",
    "Conservation of Momentum in MARLEY Events; #Delta p_y [MeV]; Counts",
    100, min_change_in_py, 1.2*std::nextafter(max_change_in_py, DBL_MAX));
  TH1D delta_pz_histogram("delta_pz_histogram",
    "Conservation of Momentum in MARLEY Events; #Delta p_z [MeV]; Counts",
    100, min_change_in_pz, 1.2*std::nextafter(max_change_in_pz, DBL_MAX));

  // Fill the histograms with values
  for(double dE : delta_Es) delta_E_histogram.Fill(dE);
  for(double px : delta_pxs) delta_px_histogram.Fill(px);
  for(double py : delta_pys) delta_py_histogram.Fill(py);
  for(double pz : delta_pzs) delta_pz_histogram.Fill(pz);

  // Write the energy change histogram to new a ROOT file
  file = new TFile("validation.root","RECREATE");
  delta_E_histogram.Write();
  delta_px_histogram.Write();
  delta_py_histogram.Write();
  delta_pz_histogram.Write();
  file->Close();
  delete file;

  return 0;
}
