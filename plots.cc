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
#include "TMarleyGenerator.hh"
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
    for (size_t j = 0; j < Es.size(); ++j) {
      double Enew = Es.at(j);
      if (Enew >= e) break;
      ibgt += bgts.at(j);
    }
    plot_ibgts.push_back(ibgt);
  }
}

int main(){

  TMarleyGenerator gen("config.txt");

  // Create plot of total reaction cross section
  std::vector<double> Eas;
  std::vector<double> tot_xs_1998;
  std::vector<double> tot_xs_2009;
  //std::vector<double> tot_xs_diffs;
  std::cout << std::setprecision(16) << std::scientific;
  TMarleyDecayScheme ds(19, 40, "ensdf.040");
  TMarleyReaction r1998("ve40ArCC_1998.react", gen.get_structure_db());
  TMarleyReaction r2009("ve40ArCC_2009.react", gen.get_structure_db());
  for (double Ea = 0.; Ea < 100; Ea += 0.1) {
    Eas.push_back(Ea);
    tot_xs_1998.push_back(1e15 * r1998.total_xs_cm(Ea)
      / marley_utils::mb); // 1e-42 cm^2
    tot_xs_2009.push_back(1e15 * r2009.total_xs_cm(Ea)
      / marley_utils::mb); // 1e-42 cm^2
  }

  //for (int i = 0; i < tot_xs.size(); ++i)
  //  tot_xs_diffs.push_back(tot_xs[i] - tot_xs_cm[i]);

  TGraph xs1998(tot_xs_1998.size(), &Eas.front(), &tot_xs_1998.front());
  xs1998.SetLineColor(4);
  xs1998.SetLineWidth(3);

  TGraph xs2009(tot_xs_2009.size(), &Eas.front(), &tot_xs_2009.front());
  xs2009.SetLineColor(2);
  xs2009.SetLineWidth(3);

  //TGraph xs_diffs(tot_xs_diffs.size(), &Eas.front(), &tot_xs_diffs.front());
  //xs_diffs.SetLineColor(4);
  //xs_diffs.SetLineWidth(3);

  TCanvas canvas;

  canvas.SetLogy(1);

  //xs_diffs.Draw("AL");
  xs1998.Draw("AL");
  xs1998.SetTitle("Total Cross Section for CC   #nu_{e} on  ^{40}Ar");
  xs1998.GetXaxis()->SetRangeUser(0.,100.);
  xs1998.GetYaxis()->SetRangeUser(1e-3,10000);
  xs1998.GetXaxis()->SetTitle("Neutrino Energy (MeV)");
  xs1998.GetXaxis()->CenterTitle();
  xs1998.GetXaxis()->SetTitleOffset(1.3);
  xs1998.GetYaxis()->SetTitle("Cross Section (10^{ -42} cm^{2})");
  xs1998.GetYaxis()->CenterTitle();

  xs2009.Draw("L");

  TLegend xs_legend(0.2, 0.15, 0.9, 0.3);
  xs_legend.SetMargin(0.2);
  xs_legend.SetTextSize(0.03);
  xs_legend.AddEntry(&xs1998, "MARLEY #nu_{e}ArCC total cross section based on   ^{40}Ti decay data", "l");
  xs_legend.AddEntry(&xs2009, "MARLEY #nu_{e}ArCC total cross section based on (p,n) data", "l");
  xs_legend.Draw();

  canvas.SaveAs("xs.pdf");

  canvas.Clear();
  canvas.SetLogy(0);

  // Create plot of B(GT) strength
  std::vector<double> cheoun2012_Exs_40K;
  std::vector<double> cheoun2012_bgts;
  std::vector<double> bhattacharya1998_Exs_40K;
  std::vector<double> bhattacharya1998_bgts;
  std::vector<double> bhattacharya2009_Exs_40K;
  std::vector<double> bhattacharya2009_bgts;
  std::vector<double> accepted1998_Exs_40K;
  std::vector<double> accepted1998_bgts;
  std::vector<double> accepted2009_Exs_40K;
  std::vector<double> accepted2009_bgts;
  double Ex, bgt;

  std::ifstream file1("cheoun2012/bgt_qrpa_calc_cheoun2012");
  while (file1 >> Ex >> bgt) {
    cheoun2012_Exs_40K.push_back(Ex);
    cheoun2012_bgts.push_back(bgt);
  }

  ifstream file2("cheoun2012/bgt_data_bhattacharya1998");
  while (file2 >> Ex >> bgt) {
    bhattacharya1998_Exs_40K.push_back(Ex);
    bhattacharya1998_bgts.push_back(bgt);
  }

  ifstream file3("cheoun2012/bgt_data_bhattacharya2009");
  while (file3 >> Ex >> bgt) {
    bhattacharya2009_Exs_40K.push_back(Ex);
    bhattacharya2009_bgts.push_back(bgt);
  }

  ifstream file4("cheoun2012/accepted_bgt_data_exp1998");
  while (file4 >> Ex >> bgt) {
    accepted1998_Exs_40K.push_back(Ex);
    accepted1998_bgts.push_back(bgt);
  }

  ifstream file5("cheoun2012/accepted_bgt_data_exp2009");
  while (file5 >> Ex >> bgt) {
    accepted2009_Exs_40K.push_back(Ex);
    accepted2009_bgts.push_back(bgt);
  }

  std::vector<double> cheoun2012_plot_Es;
  std::vector<double> cheoun2012_plot_ibgts;
  prepare_integrated_bgts(cheoun2012_Exs_40K, cheoun2012_bgts,
    cheoun2012_plot_Es, cheoun2012_plot_ibgts, 60.1);
  TGraph cheoun2012(cheoun2012_plot_Es.size(),
    &cheoun2012_plot_Es.front(), &cheoun2012_plot_ibgts.front());
  cheoun2012.SetLineColor(3);
  cheoun2012.SetLineWidth(3);

  std::vector<double> bhattacharya1998_plot_Es;
  std::vector<double> bhattacharya1998_plot_ibgts;
  prepare_integrated_bgts(bhattacharya1998_Exs_40K, bhattacharya1998_bgts,
    bhattacharya1998_plot_Es, bhattacharya1998_plot_ibgts,
    bhattacharya1998_Exs_40K.back() + 0.2);
  TGraph bhattacharya1998(bhattacharya1998_plot_Es.size(),
    &bhattacharya1998_plot_Es.front(), &bhattacharya1998_plot_ibgts.front());
  bhattacharya1998.SetLineColor(4);
  bhattacharya1998.SetLineWidth(3);

  std::vector<double> bhattacharya2009_plot_Es;
  std::vector<double> bhattacharya2009_plot_ibgts;
  prepare_integrated_bgts(bhattacharya2009_Exs_40K, bhattacharya2009_bgts,
    bhattacharya2009_plot_Es, bhattacharya2009_plot_ibgts,
    bhattacharya2009_Exs_40K.back() + 0.2);
  TGraph bhattacharya2009(bhattacharya2009_plot_Es.size(),
    &bhattacharya2009_plot_Es.front(), &bhattacharya2009_plot_ibgts.front());
  bhattacharya2009.SetLineColor(2);
  bhattacharya2009.SetLineWidth(3);

  TColor saddle_brown(1756, 0.545098, 0.270588, 0.0745098);
  std::vector<double> accepted1998_plot_Es;
  std::vector<double> accepted1998_plot_ibgts;
  prepare_integrated_bgts(accepted1998_Exs_40K, accepted1998_bgts,
    accepted1998_plot_Es, accepted1998_plot_ibgts, 60.1);
  TGraph accepted1998(accepted1998_plot_Es.size(),
    &accepted1998_plot_Es.front(), &accepted1998_plot_ibgts.front());
  accepted1998.SetLineColor(1756);
  accepted1998.SetLineWidth(3);
  accepted1998.SetLineStyle(2);

  std::vector<double> accepted2009_plot_Es;
  std::vector<double> accepted2009_plot_ibgts;
  prepare_integrated_bgts(accepted2009_Exs_40K, accepted2009_bgts,
    accepted2009_plot_Es, accepted2009_plot_ibgts, 60.1);
  TGraph accepted2009(accepted2009_plot_Es.size(),
    &accepted2009_plot_Es.front(), &accepted2009_plot_ibgts.front());
  accepted2009.SetLineColor(1);
  accepted2009.SetLineWidth(3);
  accepted2009.SetLineStyle(2);

  cheoun2012.Draw("AL");
  cheoun2012.SetTitle("Integrated Gamow-Teller Strength for CC   #nu_{e} on  ^{40}Ar");
  cheoun2012.GetXaxis()->SetTitle("^{40}K* Excitation Energy (MeV)");
  cheoun2012.GetXaxis()->CenterTitle();
  cheoun2012.GetXaxis()->SetTitleOffset(1.3);
  cheoun2012.GetYaxis()->SetTitle("Integrated B(GT)");
  cheoun2012.GetYaxis()->CenterTitle();

  bhattacharya1998.Draw("L");
  bhattacharya2009.Draw("L");
  accepted1998.Draw("L");
  accepted2009.Draw("L");

  TLegend legend(0.37, 0.1, 0.9, 0.5);
  legend.SetMargin(0.2);
  legend.AddEntry(&bhattacharya1998, "#splitline{Analog Decay Experiment Data}{from Bhattacharya, et al. (1998)}", "l");
  legend.AddEntry(&cheoun2012, "QRPA from Cheoun, et al. (2012)", "l");
  legend.AddEntry(&bhattacharya2009, "#splitline{(p,n) Data from Bhattacharya,}{et al. (2009)}", "l");
  legend.AddEntry(&accepted1998, "MARLEY B(GT) based on 1998 data", "l");
  legend.AddEntry(&accepted2009, "MARLEY B(GT) based on 2009 data", "l");
  legend.Draw();

  canvas.SaveAs("bgt.pdf");
  canvas.Clear();

  // Open a file containing a ROOT tree filled with MARLEY events.
  // Associate the tree with a TMarleyEvent pointer.
  TFile* file = new TFile("event_tree.root", "READ");
  TTree* t = nullptr;
  file->GetObject("MARLEY Event Tree", t);
  TMarleyEvent* e = new TMarleyEvent;
  t->GetBranch("events")->SetAddress(&e);

  // Particle energies from all events
  std::vector<double> gamma_Es, electron_Es, neutrino_Es;
  // Excitation energies of residual nucleus
  std::vector<double> nuc_Exs;
  // Total number of events
  int n = t->GetEntries();

  // Final particle kinetic energies (lookup table by PID)
  std::unordered_map<int, std::vector<double> > fparticle_KEs;
  // Final particle counts (lookup table by PID)
  std::unordered_map<int, int> fparticle_counts;

  // Storage for 3-momentum components used for momentum transfer plot
  double nu_px = 0, nu_py = 0, nu_pz = 0, l_px = 0, l_py = 0, l_pz = 0;
  // Vector of 3-momentum transfer values
  std::vector<double> momentum_transfers;

  // Cycle through each event in the tree
  for (int i = 0; i < n; ++i) {
    if (i % 5000 == 0) std::cout << "Event " << i << std::endl;

    // Load the current event
    t->GetEntry(i);

    bool found_nu = false, found_l = false;

    nuc_Exs.push_back(e->get_E_level());

    // Loop over the initial particles. For each neutrino, store
    // its energy.
    std::list<TMarleyParticle>* iparts = e->get_initial_particles();
    for(const auto& particle : *iparts) {
      int pid = particle.get_id();
      if (pid == marley_utils::ELECTRON_NEUTRINO) {
        neutrino_Es.push_back(particle.get_total_energy());
        nu_px = particle.get_px();
        nu_py = particle.get_py();
        nu_pz = particle.get_pz();
        found_nu = true;
      }
    }

    // Loop over the final particles.
    std::list<TMarleyParticle>* fparts = e->get_final_particles();
    for(const auto& particle : *fparts) {

      int pid = particle.get_id();

      if (pid == marley_utils::ELECTRON) {
        l_px = particle.get_px();
        l_py = particle.get_py();
        l_pz = particle.get_pz();
        found_l = true;
      }

      std::unordered_map<int, std::vector<double> >::iterator it
        = fparticle_KEs.find(pid);

      if (it == fparticle_KEs.end()) {
        // This is a particle that doesn't have a vector of kinetic energies yet
        fparticle_KEs[pid] = std::vector<double>(1, particle.get_kinetic_energy());
        fparticle_counts[pid] = 1;
      }
      // This particle has a vector of kinetic energies already, so add it
      else {
        (it->second).push_back(particle.get_kinetic_energy());
        ++fparticle_counts[pid];
      }
    }

    if (found_nu && found_l) {
      double qx = l_px - nu_px;
      double qy = l_py - nu_py;
      double qz = l_pz - nu_pz;
      double q = std::sqrt(std::pow(qx, 2) + std::pow(qy, 2) + std::pow(qz, 2));
      momentum_transfers.push_back(q);
    }
  }

  // Close the event tree ROOT file
  file->Close();
  delete file;

  // Find the maximum neutrino energy
  double E_nmax = 0;
  for (const double& en : neutrino_Es) {
    if (en > E_nmax) E_nmax = en;
  }

  // Find the maximum momentum transfer
  double qmax = 0;
  for (const double& q : momentum_transfers) {
    if (q > qmax) qmax = q;
  }

  // Create a histogram to store the neutrino energies.
  // Use 100 bins, and adjust the upper limit of the last bin so that
  // the maximum observed value just barely falls within it.
  TH1D neutrino_E_histogram("neutrino_E_histogram",
    "Interacting  #nu_{e}  from Fermi-Dirac source; Total Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(E_nmax, DBL_MAX));
  // Do the same thing for the residual nucleus excitation energies
  TH1D ex_residue_histogram("Ex_40K",
    "Strength to ^{40}K* states; Excitation Energy [MeV]; Counts",
    100, 0., 50.);
  // Do the same thing for the momentum transfers
  TH1D q_histogram("q_hist",
    "Momentum Transfer; q [MeV]; Counts",
    100, 0., 1.1*std::nextafter(qmax, DBL_MAX));

  // Fill the histograms with values
  for(const double& en : neutrino_Es) neutrino_E_histogram.Fill(en);
  for(const double& ex : nuc_Exs) ex_residue_histogram.Fill(ex);
  for(const double& q : momentum_transfers) q_histogram.Fill(q);

  //// Write the histograms to new a ROOT file
  //file = new TFile("plot_histograms.root","RECREATE");
  //gamma_E_histogram.Write();
  //electron_E_histogram.Write();
  //neutrino_E_histogram.Write();
  //ex_residue_histogram.Write();
  //file->Close();
  //delete file;

  // Create PDFs plotting the spectra histograms
  canvas.Clear();
  neutrino_E_histogram.Draw();
  neutrino_E_histogram.SetLineWidth(1.8);
  neutrino_E_histogram.GetXaxis()->SetTitleOffset(1.2);
  neutrino_E_histogram.GetYaxis()->SetTitleOffset(1.5);
  neutrino_E_histogram.GetXaxis()->CenterTitle();
  neutrino_E_histogram.GetYaxis()->CenterTitle();

  //std::vector<double> nu_E_dist;
  //std::vector<double> nu_Es;
  //double Emin = gen.get_nu_source().get_Emin();
  //double Emax = gen.get_nu_source().get_Emax();
  //double step = (Emax - Emin) / 1e4;
  //for (double e = Emin; e <= Emax; e += step) {
  //  nu_Es.push_back(e);
  //  nu_E_dist.push_back(n * gen.normalized_Ea_pdf(e));
  //}

  //TGraph fd(nu_E_dist.size(), &nu_Es.front(), &nu_E_dist.front());
  //fd.Draw("L");

  canvas.SaveAs("E_neutrino.pdf");

  canvas.Clear();
  ex_residue_histogram.Draw();
  ex_residue_histogram.SetLineWidth(1.8);
  ex_residue_histogram.GetXaxis()->SetTitleOffset(1.2);
  ex_residue_histogram.GetYaxis()->SetTitleOffset(1.5);
  ex_residue_histogram.GetXaxis()->CenterTitle();
  ex_residue_histogram.GetYaxis()->CenterTitle();

  canvas.SaveAs("E_residue.pdf");

  canvas.Clear();
  q_histogram.Draw();
  q_histogram.SetLineWidth(1.8);
  q_histogram.GetXaxis()->SetTitleOffset(1.2);
  q_histogram.GetYaxis()->SetTitleOffset(1.5);
  q_histogram.GetXaxis()->CenterTitle();
  q_histogram.GetYaxis()->CenterTitle();

  canvas.SaveAs("momentum_transfer.pdf");

  // Loop over all of the final particles, making kinetic energy spectrum
  // plots for each one
  for (const auto& pair : fparticle_KEs) {
    int pid = pair.first;
    std::string symbol;
    if (pid > marley_utils::ALPHA) {
      int Z = TMarleyMassTable::get_particle_Z(pid);
      int A = TMarleyMassTable::get_particle_A(pid);
      symbol = std::to_string(A) + marley_utils::element_symbols.at(Z);
    }
    else symbol = marley_utils::particle_symbols.at(pid);

    std::cout << symbol << " multiplicity per event = "
      << static_cast<double>(fparticle_counts.at(pid)) / n
      << std::endl;

    // Get the maximum kinetic energy for this particle
    double Emax = 0.;
    for (const double& e : pair.second)
      if (e > Emax) Emax = e;
    // Create a histogram to store the gamma ray energies.
    // Use 100 bins, and adjust the upper limit of the last bin so that
    // the maximum observed value just barely falls within it.
    std::string plot_title;
    if (pid == marley_utils::ELECTRON) plot_title = std::string(
      "Electron spectrum for CC FD   #nu_{e} on ")
      + "  ^{40}Ar; Kinetic Energy [MeV]; Counts";
    else if (pid == marley_utils::PHOTON)
       plot_title = std::string("De-excitation  #gamma-ray")
       + " spectrum for CC FD   #nu_{e}"
       + " on  ^{40}Ar; Energy [MeV]; Counts";
    else plot_title = "De-excitation " + symbol
      + " spectrum for CC FD   #nu_{e}"
      + " on  ^{40}Ar; Kinetic Energy [MeV]; Counts";
    TH1D KE_histogram("KE_histogram", plot_title.c_str(),
      100, 0, 1.1*std::nextafter(Emax, DBL_MAX));
    for(const double& e : pair.second) KE_histogram.Fill(e);
    canvas.Clear();
    KE_histogram.Draw();
    KE_histogram.GetXaxis()->SetTitleOffset(1.3);
    KE_histogram.GetYaxis()->SetTitleOffset(1.5);
    KE_histogram.GetXaxis()->CenterTitle();
    KE_histogram.GetYaxis()->CenterTitle();
    KE_histogram.SetLineWidth(1.8);

    canvas.SaveAs(("KE_" + symbol + ".pdf").c_str());
  }

  return 0;
}
