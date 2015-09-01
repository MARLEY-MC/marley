#ifndef USE_ROOT
#error "This MARLEY analysis code currently must be compiled with ROOT support."
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyEvent.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleyReaction.hh"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"

int main(){

  // Create plot of total reaction cross section
  std::vector<double> Eas;
  std::vector<double> tot_xs;
  std::cout << std::setprecision(16) << std::scientific;
  TMarleyDecayScheme ds(19, 40, "ensdf.040");
  TMarleyReaction r("ve40ArCC.react", &ds);
  for (double Ea = 4.4; Ea < 60; Ea += 0.1) {
    Eas.push_back(Ea);
    tot_xs.push_back(1e12 * r.total_xs(Ea) / marley_utils::mb); // fb
  }

  TGraph xs(tot_xs.size(), &Eas.front(), &tot_xs.front());
  xs.SetLineColor(4);
  xs.SetLineWidth(3);

  TCanvas canvas;

  xs.Draw("AL");

  xs.SetTitle("Total Cross Section for CC   #nu_{e} on  ^{40}Ar");
  xs.GetXaxis()->SetTitle("Neutrino Energy (MeV)");
  xs.GetXaxis()->CenterTitle();
  xs.GetXaxis()->SetTitleOffset(1.3);
  xs.GetYaxis()->SetTitle("Cross Section (fb)");
  xs.GetYaxis()->CenterTitle();

  canvas.SetLogy();

  canvas.SaveAs("xs.pdf");

  canvas.Clear();
  //return 0;

  // Create plot of B(GT) strength
  double rescale_factor = 1.4;
  std::vector<double> cheoun2012_Exs_40K;
  std::vector<double> cheoun2012_integrated_BGTs;
  std::vector<double> rescaled_cheoun2012_integrated_BGTs;
  std::vector<double> bhattacharya_1998_Exs_40K;
  std::vector<double> bhattacharya_1998_integrated_BGTs;
  std::vector<double> bhattacharya_2009_Exs_40K;
  std::vector<double> bhattacharya_2009_integrated_BGTs;
  std::ifstream file1("plot_prep/cheoun2012_qrpa_BGTs");
  double Ex, integrated_BGT;
  while (file1.good()) {
    file1 >> Ex >> integrated_BGT;
    cheoun2012_Exs_40K.push_back(Ex);
    cheoun2012_integrated_BGTs.push_back(integrated_BGT);
    rescaled_cheoun2012_integrated_BGTs.push_back(rescale_factor
      * integrated_BGT);
  }
  ifstream file2("plot_prep/bhattacharya1998_BGTs");
  while (file2.good()) {
    file2 >> Ex >> integrated_BGT;
    bhattacharya_1998_Exs_40K.push_back(Ex);
    bhattacharya_1998_integrated_BGTs.push_back(integrated_BGT);
  }
  ifstream file3("plot_prep/bhattacharya2009_BGTs");
  while (file3.good()) {
    file3 >> Ex >> integrated_BGT;
    bhattacharya_2009_Exs_40K.push_back(Ex);
    bhattacharya_2009_integrated_BGTs.push_back(integrated_BGT);
  }

  int num = cheoun2012_integrated_BGTs.size();
  TGraph cheoun2012(num, &cheoun2012_Exs_40K.front(), &cheoun2012_integrated_BGTs.front());
  cheoun2012.SetLineColor(3);
  cheoun2012.SetLineWidth(3);

  num = rescaled_cheoun2012_integrated_BGTs.size();
  TGraph cheoun2012_rescaled(num, &cheoun2012_Exs_40K.front(), &rescaled_cheoun2012_integrated_BGTs.front());
  cheoun2012_rescaled.SetLineColor(1);
  cheoun2012_rescaled.SetLineWidth(3);

  num = bhattacharya_1998_integrated_BGTs.size();
  TGraph bhattacharya1998(num, &bhattacharya_1998_Exs_40K.front(), &bhattacharya_1998_integrated_BGTs.front());
  bhattacharya1998.SetLineColor(4);
  bhattacharya1998.SetLineWidth(3);

  num = bhattacharya_2009_integrated_BGTs.size();
  TGraph bhattacharya2009(num, &bhattacharya_2009_Exs_40K.front(), &bhattacharya_2009_integrated_BGTs.front());
  bhattacharya2009.SetLineColor(2);
  bhattacharya2009.SetLineWidth(3);

  cheoun2012_rescaled.Draw("AL");
  cheoun2012_rescaled.SetTitle("Integrated Gamow-Teller Strength for CC   #nu_{e} on  ^{40}Ar");
  cheoun2012_rescaled.GetXaxis()->SetTitle("^{40}K* Excitation Energy (MeV)");
  cheoun2012_rescaled.GetXaxis()->CenterTitle();
  cheoun2012_rescaled.GetXaxis()->SetTitleOffset(1.3);
  cheoun2012_rescaled.GetYaxis()->SetTitle("Integrated B(GT)");
  cheoun2012_rescaled.GetYaxis()->CenterTitle();

  cheoun2012.Draw("L");
  bhattacharya1998.Draw("L");
  bhattacharya2009.Draw("L");

  TLegend legend(0.5, 0.1, 0.9, 0.4);
  legend.AddEntry(&cheoun2012, "QRPA from Cheoun, et al. (2012)", "l");
  legend.AddEntry(&bhattacharya2009, "#splitline{(p,n) Data from Bhattacharya,}{et al. (2009)}", "l");
  legend.AddEntry(&cheoun2012_rescaled, "Rescaled QRPA", "l");
  legend.AddEntry(&bhattacharya1998, "#splitline{Analog Decay Experiment Data}{from Bhattacharya, et al. (1998)}", "l");
  legend.Draw();

  canvas.SaveAs("bgt.pdf");
  canvas.Clear();

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a TMarleyEvent pointer.
  TFile* file = new TFile("event_tree.root", "READ");
  TTree* t = nullptr;
  file->GetObject("event_tree", t);
  TMarleyEvent* e = new TMarleyEvent;
  t->GetBranch("events")->SetAddress(&e);

  // Particle energies from all events
  std::vector<double> gamma_Es, electron_Es, neutrino_Es;
  // Total number of gamma rays emitted
  int num_gammas = 0;
  // Total number of events
  int n = t->GetEntries();

  // Cycle through each event in the tree
  for (int i = 0; i < n; ++i) {
    if (i % 5000 == 0) std::cout << "Event " << i << std::endl;

    // Load the current event
    t->GetEntry(i);

    // Loop over the initial particles. For each neutrino, store
    // its energy.
    std::list<TMarleyParticle>* iparts = e->get_initial_particles();
    for(const auto& particle : *iparts) {
      int pid = particle.get_id();
      if (pid == marley_utils::ELECTRON_NEUTRINO) {
        neutrino_Es.push_back(particle.get_total_energy());
      }
    }

    // Loop over the final particles. For each gamma ray, store
    // its energy and increment the gamma counter. Also store the
    // energy of each electron.
    std::list<TMarleyParticle>* fparts = e->get_final_particles();
    for(const auto& particle : *fparts) {
      int pid = particle.get_id();
      if (pid == marley_utils::PHOTON) {
        gamma_Es.push_back(particle.get_total_energy());
        ++num_gammas;
      }
      else if (pid == marley_utils::ELECTRON) {
        electron_Es.push_back(particle.get_total_energy());
      }
    }
  }

  std::cout << "There were " << num_gammas << " produced in " << n
    << " events." << std::endl;

  // Close the event tree ROOT file
  file->Close();
  delete file;

  // Find the maximum gamma-ray energy
  double E_gmax = 0, E_emax = 0, E_nmax = 0;
  for (const double& eg : gamma_Es) {
    if (eg > E_gmax) E_gmax = eg;
  }
  for (const double& ee : electron_Es) {
    if (ee > E_emax) E_emax = ee;
  }
  for (const double& en : neutrino_Es) {
    if (en > E_nmax) E_nmax = en;
  }

  // Create a histogram to store the gamma ray energies.
  // Use 100 bins, and adjust the upper limit of the last bin so that
  // the maximum observed value just barely falls within it.
  TH1D gamma_E_histogram("gamma_E_histogram",
    "De-excitation  #gamma-ray spectrum for CC SN   #nu_{e} on  ^{40}Ar; Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(E_gmax, DBL_MAX));
  // Do the same thing for the other two particle types
  TH1D electron_E_histogram("electron_E_histogram",
    "Electron spectrum for CC SN   #nu_{e} on  ^{40}Ar; Total Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(E_emax, DBL_MAX));
  TH1D neutrino_E_histogram("neutrino_E_histogram",
    "Supernova  #nu_{e}  spectrum; Total Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(E_nmax, DBL_MAX));

  // Fill the histograms with values
  for(const double& eg : gamma_Es) gamma_E_histogram.Fill(eg);
  for(const double& ee : electron_Es) electron_E_histogram.Fill(ee);
  for(const double& en : neutrino_Es) neutrino_E_histogram.Fill(en);

  //// Write the histograms to new a ROOT file
  //file = new TFile("plot_histograms.root","RECREATE");
  //gamma_E_histogram.Write();
  //electron_E_histogram.Write();
  //neutrino_E_histogram.Write();
  //file->Close();
  //delete file;

  // Create PDFs plotting the spectra histograms
  gamma_E_histogram.Draw();
  gamma_E_histogram.GetYaxis()->SetTitleOffset(1.3);
  gamma_E_histogram.GetXaxis()->CenterTitle();
  gamma_E_histogram.GetYaxis()->CenterTitle();
  gamma_E_histogram.GetXaxis()->SetTitleOffset(1.3);
  gamma_E_histogram.SetLineWidth(1.8);
  canvas.SaveAs("E_gamma.pdf");

  canvas.Clear();
  electron_E_histogram.Draw();
  electron_E_histogram.GetXaxis()->SetTitleOffset(1.3);
  electron_E_histogram.GetYaxis()->SetTitleOffset(1.5);
  electron_E_histogram.GetXaxis()->CenterTitle();
  electron_E_histogram.GetYaxis()->CenterTitle();
  electron_E_histogram.SetLineWidth(1.8);
  canvas.SaveAs("E_electron.pdf");

  canvas.Clear();
  neutrino_E_histogram.Draw();
  neutrino_E_histogram.SetLineWidth(1.8);
  neutrino_E_histogram.GetXaxis()->SetTitleOffset(1.2);
  neutrino_E_histogram.GetYaxis()->SetTitleOffset(1.5);
  neutrino_E_histogram.GetXaxis()->CenterTitle();
  neutrino_E_histogram.GetYaxis()->CenterTitle();
  canvas.SaveAs("E_neutrino.pdf");

  return 0;
}
