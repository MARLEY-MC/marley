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
#include "TMarleyDecayScheme.hh"
#include "TMarleyEvent.hh"
#include "TMarleyNeutrinoSource.hh"
#include "TMarleyReaction.hh"

#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"

void prepare_integrated_bgts(const std::vector<double>& Es,
  const std::vector<double>& bgts, std::vector<double>& plot_Es,
  std::vector<double>& plot_ibgts, const double max_E)
{
  plot_Es.clear();
  plot_ibgts.clear();
  for (double e = 0.; e <= max_E; e += 0.1) {
    double ibgt = 0.;
    plot_Es.push_back(e);
    for (int j = 0; j < Es.size(); ++j) {
      double Enew = Es.at(j);
      if (Enew >= e) break;
      ibgt += bgts.at(j);
    }
    plot_ibgts.push_back(ibgt);
  }
}

int main(){

  //// Create plot of total reaction cross section
  //std::vector<double> Eas;
  //std::vector<double> tot_xs;
  //std::cout << std::setprecision(16) << std::scientific;
  //TMarleyDecayScheme ds(19, 40, "ensdf.040");
  //TMarleyReaction r("ve40ArCC_1998.react", &ds);
  //for (double Ea = 4.4; Ea < 60; Ea += 0.1) {
  //  Eas.push_back(Ea);
  //  tot_xs.push_back(1e12 * r.total_xs(Ea) / marley_utils::mb); // fb
  //}

  //TGraph xs(tot_xs.size(), &Eas.front(), &tot_xs.front());
  //xs.SetLineColor(4);
  //xs.SetLineWidth(3);

  TCanvas canvas;

  //xs.Draw("AL");

  //xs.SetTitle("Total Cross Section for CC   #nu_{e} on  ^{40}Ar");
  //xs.GetXaxis()->SetTitle("Neutrino Energy (MeV)");
  //xs.GetXaxis()->CenterTitle();
  //xs.GetXaxis()->SetTitleOffset(1.3);
  //xs.GetYaxis()->SetTitle("Cross Section (fb)");
  //xs.GetYaxis()->CenterTitle();

  ////canvas.SetLogy();

  //canvas.SaveAs("xs.pdf");

  //canvas.Clear();
  ////return 0;

  //// Create plot of B(GT) strength
  //std::vector<double> cheoun2012_Exs_40K;
  //std::vector<double> cheoun2012_bgts;
  //std::vector<double> bhattacharya1998_Exs_40K;
  //std::vector<double> bhattacharya1998_bgts;
  //std::vector<double> bhattacharya2009_Exs_40K;
  //std::vector<double> bhattacharya2009_bgts;
  //std::vector<double> accepted1998_Exs_40K;
  //std::vector<double> accepted1998_bgts;
  //std::vector<double> accepted2009_Exs_40K;
  //std::vector<double> accepted2009_bgts;
  //double Ex, bgt;

  //std::ifstream file1("cheoun2012/bgt_qrpa_calc_cheoun2012");
  //while (file1 >> Ex >> bgt) {
  //  cheoun2012_Exs_40K.push_back(Ex);
  //  cheoun2012_bgts.push_back(bgt);
  //}

  //ifstream file2("cheoun2012/bgt_data_bhattacharya1998");
  //while (file2 >> Ex >> bgt) {
  //  bhattacharya1998_Exs_40K.push_back(Ex);
  //  bhattacharya1998_bgts.push_back(bgt);
  //}

  //ifstream file3("cheoun2012/bgt_data_bhattacharya2009");
  //while (file3 >> Ex >> bgt) {
  //  bhattacharya2009_Exs_40K.push_back(Ex);
  //  bhattacharya2009_bgts.push_back(bgt);
  //}

  //ifstream file4("cheoun2012/accepted_bgt_data_exp1998");
  //while (file4 >> Ex >> bgt) {
  //  accepted1998_Exs_40K.push_back(Ex);
  //  accepted1998_bgts.push_back(bgt);
  //}

  //ifstream file5("cheoun2012/accepted_bgt_data_exp2009");
  //while (file5 >> Ex >> bgt) {
  //  accepted2009_Exs_40K.push_back(Ex);
  //  accepted2009_bgts.push_back(bgt);
  //}

  //std::vector<double> cheoun2012_plot_Es;
  //std::vector<double> cheoun2012_plot_ibgts;
  //prepare_integrated_bgts(cheoun2012_Exs_40K, cheoun2012_bgts,
  //  cheoun2012_plot_Es, cheoun2012_plot_ibgts, 60.1);
  //TGraph cheoun2012(cheoun2012_plot_Es.size(),
  //  &cheoun2012_plot_Es.front(), &cheoun2012_plot_ibgts.front());
  //cheoun2012.SetLineColor(3);
  //cheoun2012.SetLineWidth(3);

  //std::vector<double> bhattacharya1998_plot_Es;
  //std::vector<double> bhattacharya1998_plot_ibgts;
  //prepare_integrated_bgts(bhattacharya1998_Exs_40K, bhattacharya1998_bgts,
  //  bhattacharya1998_plot_Es, bhattacharya1998_plot_ibgts,
  //  bhattacharya1998_Exs_40K.back() + 0.2);
  //TGraph bhattacharya1998(bhattacharya1998_plot_Es.size(),
  //  &bhattacharya1998_plot_Es.front(), &bhattacharya1998_plot_ibgts.front());
  //bhattacharya1998.SetLineColor(4);
  //bhattacharya1998.SetLineWidth(3);

  //std::vector<double> bhattacharya2009_plot_Es;
  //std::vector<double> bhattacharya2009_plot_ibgts;
  //prepare_integrated_bgts(bhattacharya2009_Exs_40K, bhattacharya2009_bgts,
  //  bhattacharya2009_plot_Es, bhattacharya2009_plot_ibgts,
  //  bhattacharya2009_Exs_40K.back() + 0.2);
  //TGraph bhattacharya2009(bhattacharya2009_plot_Es.size(),
  //  &bhattacharya2009_plot_Es.front(), &bhattacharya2009_plot_ibgts.front());
  //bhattacharya2009.SetLineColor(2);
  //bhattacharya2009.SetLineWidth(3);

  //TColor saddle_brown(1756, 0.545098, 0.270588, 0.0745098);
  //std::vector<double> accepted1998_plot_Es;
  //std::vector<double> accepted1998_plot_ibgts;
  //prepare_integrated_bgts(accepted1998_Exs_40K, accepted1998_bgts,
  //  accepted1998_plot_Es, accepted1998_plot_ibgts, 60.1);
  //TGraph accepted1998(accepted1998_plot_Es.size(),
  //  &accepted1998_plot_Es.front(), &accepted1998_plot_ibgts.front());
  //accepted1998.SetLineColor(1756);
  //accepted1998.SetLineWidth(3);
  //accepted1998.SetLineStyle(2);

  //std::vector<double> accepted2009_plot_Es;
  //std::vector<double> accepted2009_plot_ibgts;
  //prepare_integrated_bgts(accepted2009_Exs_40K, accepted2009_bgts,
  //  accepted2009_plot_Es, accepted2009_plot_ibgts, 60.1);
  //TGraph accepted2009(accepted2009_plot_Es.size(),
  //  &accepted2009_plot_Es.front(), &accepted2009_plot_ibgts.front());
  //accepted2009.SetLineColor(1);
  //accepted2009.SetLineWidth(3);
  //accepted2009.SetLineStyle(2);

  //cheoun2012.Draw("AL");
  //cheoun2012.SetTitle("Integrated Gamow-Teller Strength for CC   #nu_{e} on  ^{40}Ar");
  //cheoun2012.GetXaxis()->SetTitle("^{40}K* Excitation Energy (MeV)");
  //cheoun2012.GetXaxis()->CenterTitle();
  //cheoun2012.GetXaxis()->SetTitleOffset(1.3);
  //cheoun2012.GetYaxis()->SetTitle("Integrated B(GT)");
  //cheoun2012.GetYaxis()->CenterTitle();

  //bhattacharya1998.Draw("L");
  //bhattacharya2009.Draw("L");
  //accepted1998.Draw("L");
  //accepted2009.Draw("L");

  //TLegend legend(0.37, 0.1, 0.9, 0.5);
  //legend.SetMargin(0.2);
  //legend.AddEntry(&bhattacharya1998, "#splitline{Analog Decay Experiment Data}{from Bhattacharya, et al. (1998)}", "l");
  //legend.AddEntry(&cheoun2012, "QRPA from Cheoun, et al. (2012)", "l");
  //legend.AddEntry(&bhattacharya2009, "#splitline{(p,n) Data from Bhattacharya,}{et al. (2009)}", "l");
  //legend.AddEntry(&accepted1998, "MARLEY B(GT) based on 1998 data", "l");
  //legend.AddEntry(&accepted2009, "MARLEY B(GT) based on 2009 data", "l");
  //legend.Draw();

  //canvas.SaveAs("bgt.pdf");
  //canvas.Clear();

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a TMarleyEvent pointer.
  TFile* file = new TFile("event_tree.root", "READ");
  TTree* t = nullptr;
  file->GetObject("MARLEY Event Tree", t);
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
    "De-excitation  #gamma-ray spectrum for CC FD   #nu_{e} on  ^{40}Ar (2009 data); Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(E_gmax, DBL_MAX));
  // Do the same thing for the other two particle types
  TH1D electron_E_histogram("electron_E_histogram",
    "Electron spectrum for CC FD   #nu_{e} on  ^{40}Ar (2009 data); Total Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(E_emax, DBL_MAX));
  TH1D neutrino_E_histogram("neutrino_E_histogram",
    "Fermi-Dirac  #nu_{e}  spectrum (2009 data); Total Energy [MeV]; Counts",
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

  //std::function<double(double)> f = std::bind(
  //  &TMarleyNeutrinoSource::fermi_dirac_distribution,
  //  3.5, 2.8, std::placeholders::_1 /*nu_energy*/);

  //double integral = marley_utils::num_integrate(f, 0, 100, 1e5);
  //std::vector<double> nu_E_dist;
  //std::vector<double> nu_Es;
  //for (double e = 0.; e <= 60; e += 0.1) {
  //  nu_Es.push_back(e);
  //  nu_E_dist.push_back(1e6 * f(e) / integral);
  //}

  //TGraph fd(nu_E_dist.size(), &nu_Es.front(), &nu_E_dist.front());
  //fd.Draw("L");

  canvas.SaveAs("E_neutrino.pdf");

  return 0;
}
