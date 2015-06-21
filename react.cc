#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleyEvent.hh"
#include "TMarleyMassTable.hh"
#include "TMarleyReaction.hh"

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

double fermi_dirac_distribution(double C, bool e_flavor, bool anti, double nu_energy){
  double eta = 0;
  double T = 0;
  double N_nu;

  if(e_flavor && !anti){
    T = 3.5; // temperature in MeV
    N_nu = 2.8;  // total number of electron neutrinos expected (x10^57)
  }
  else if(e_flavor && anti) {
    T = 5.0; // temperature in MeV
    N_nu = 1.9;  // total number of electron anti-neutrinos expected (x10^57)
  }
  else { // !e_flavor
    T = 8.0; // temperature in MeV
    N_nu = 5.0; // total number of mu+tau neutrinos + anti-neutrinos expected (x10^57)
  }
  return (C/std::pow(T,3))*(std::pow(nu_energy,2)/(1+std::exp(nu_energy/(T-eta))))*N_nu;
}

int main(){

  //for (int i = 1; i < 8; i++)
  //  TMarleyMassTable::print_separation_energies(19, 40, i);
  //std::cout << std::endl << std::endl;

  // Sample from electron flavor supernova neutrinos for C = 0.55
  std::function<double(double)> f = std::bind(fermi_dirac_distribution, 0.55,
    true, false, std::placeholders::_1);

  #ifdef USE_ROOT
  // Create a pointer to an event object. This will
  // be used to fill the event tree with the
  // events generated in the loop below
  TMarleyEvent* p_event = nullptr;

  // Create a ROOT file to store the event tree
  TFile treeFile("event_tree.root","RECREATE");

  // Create a ROOT tree to store the events
  TTree event_tree("event_tree", "A tree of TMarleyEvent objects");

  // Create a branch in this ROOT tree, and associate
  // it with the event pointer we made before
  event_tree.Branch("events", &p_event, 32000, 99);

  // Number of MB written to the ROOT file
  double MB_written = 0;
  #endif

  // Select the isotope and ENSDF file to use for the simulation
  std::string nuc_id = marley_utils::nuc_id(19, 40); // 40K
  std::string filename = "ensdf.040";

  // Create a decay scheme object to store data
  // imported from the ENSDF file
  TMarleyDecayScheme ds(nuc_id, filename);

  TMarleyReaction r("ve40ArCC.react", &ds);

  // TODO: debug numerical errors that arise when
  // Ea = E_threshold

  // Simulate a charged current reaction
  int n_events = 1000;
  double Ea; // Incident neutrino energy

  // Set the precision for outputting floating-point numbers
  std::cout << std::fixed;
  std::cout.precision(3);

  for (int i = 1; i <= n_events; ++i) {

    // Sample a supernova neutrino energy
    Ea = marley_utils::rejection_sample(f, 4.36, 50, 1e-8);

    // Create an event using the charged current reaction
    TMarleyEvent e = r.create_event(Ea);

    std::cout << "Event Count = " << i << "/" << n_events << std::endl;

    #ifdef USE_ROOT
    // Get the address of this event object
    p_event = new TMarleyEvent;
    *p_event = e;

    // Store this event in the ROOT tree
    MB_written += event_tree.Fill()/1e6;
    std::cout << "MB written = " << MB_written << std::endl;
    //std::cout << "Elapsed time = "
    // Move up one line in std::cout
    std::cout << "\033[F";
    #endif

    // Move up one line in std::cout
    std::cout << "\033[F";
  }

  #ifdef USE_ROOT
  event_tree.Write();
  treeFile.Close();
  #endif

  return 0;
}
