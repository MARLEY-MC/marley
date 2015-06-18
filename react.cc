#include <iostream>
#include <functional>
#include <string>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleyEvent.hh"
#include "TMarleyReaction.hh"
#include "TMarleyMassTable.hh"

#ifdef USE_ROOT
// ROOT Stuff
#include "TH2D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLatex.h"
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
  else if(e_flavor && anti){
    T = 5.0; // temperature in MeV
    N_nu = 1.9;  // total number of electron anti-neutrinos expected (x10^57)
  }
  else if(!e_flavor && (!anti || anti)){
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
  TH2D hist = TH2D("hist", "E_{e^-} vs E_level",
    100, 0, 50, 25, 0, 6); // MeV

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

  // Incident neutrino energy
  double Ea;

  // Simulate a charged current reaction
  double E_level, e_eminus;
  for (int i = 1; i <= 100; i++) {
    std::cout << "Event Count = " << i << std::endl;

    // Sample a supernova neutrino energy
    Ea = marley_utils::rejection_sample(f, 4.36, 50, 1e-8);

    // Create an event using the charged current reaction
    TMarleyEvent e = r.create_event(Ea);

    E_level = e.get_E_level();
    e_eminus = e.get_ejectile()->get_total_energy();

    e.print_event();
    //std::cout << "neutrino energy = " << Ea << std::endl;
    std::cout << std::endl;
    std::cout << "E_level = " << E_level << std::endl;
    std::cout << "e- energy = " << e_eminus << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << std::endl << std::endl;

    #ifdef USE_ROOT
    hist.Fill(e_eminus, E_level);

    // Get the address of this event object
    p_event = new TMarleyEvent;
    *p_event = e;
    // Store this event in the ROOT tree
    int bytes_written = event_tree.Fill();
    std::cout << "bytes written = " << bytes_written << std::endl;
    #endif
  }

  #ifdef USE_ROOT
  hist.GetXaxis()->SetTitle("E_{e^-} [MeV]");
  hist.GetYaxis()->SetTitle("E_{level} [MeV]");

  event_tree.Write();
  treeFile.Close();

  TFile histFile("hist.root","RECREATE");
  hist.Write();
  histFile.Close();

  // Read the tree back from the file and
  // display each event
  TMarleyEvent* pe = new TMarleyEvent;
  TFile file("event_tree.root");
  TTree* t = static_cast<TTree*>(file.Get("event_tree"));
  t->GetBranch("events")->SetAddress(&pe);
  for (int i = 0, n = t->GetEntries(); i < n; ++i) {
    t->GetEntry(i);
    std::cout << "n_children = " << pe->get_residue()->get_children()->size() << std::endl;
  }
  #endif

  return 0;
}
