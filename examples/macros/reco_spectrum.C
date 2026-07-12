#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

void reco_spectrum(const std::string& filename) {

  const int NUM_BINS = 110;
  const double E_MIN = 0.;
  const double E_MAX = 55.;
  const double BIN_WIDTH = ( E_MAX - E_MIN ) / NUM_BINS;

  const int NEUTRON = 2112;
  const int ALPHA = 1000020040;
  const double Q_ground_state = 1.5044;

  TH1D* true_Es = new TH1D("true_Es", "E_{#nu,true}",
    NUM_BINS, E_MIN, E_MAX);
  TH1D* eq_Es = new TH1D("eq_Es", "E_{#nu,reco 2}",
    NUM_BINS, E_MIN, E_MAX);
  TH1D* reco_Es = new TH1D("reco_Es", "energy distributions;"
    " neutrino energy E_{#nu} (MeV); #left[d#sigma / dE_{#nu}"
    "#right]_{flux} (10^{-42} cm^{2} / MeV)", NUM_BINS, E_MIN, E_MAX);

  true_Es->SetStats(false);
  eq_Es->SetStats(false);
  reco_Es->SetStats(false);

  const int all_gamma_nuc_pdg = 1000190400; // 40K
  const int one_n_nuc_pdg = 1000190390; // 39K
  const int one_p_nuc_pdg = 1000180390; // 39Ar

  size_t all_gamma_count = 0;
  size_t one_n_count = 0;
  size_t one_p_count = 0;
  size_t other_count = 0;

  TFile* f = TFile::Open(filename.c_str());
  TTree* t = static_cast<TTree*>(f->Get("mst"));

  double Ev, KEl, xsec;
  int np;
  std::vector<int>* pdgp = nullptr;
  std::vector<double>* KEp = nullptr;

  t->SetBranchAddress("Ev", &Ev);
  t->SetBranchAddress("KEl", &KEl);
  t->SetBranchAddress("xsec", &xsec);
  t->SetBranchAddress("np", &np);
  t->SetBranchAddress("pdgp", &pdgp);
  t->SetBranchAddress("KEp", &KEp);

  Long64_t num_events = t->GetEntries();

  for (Long64_t i = 0; i < num_events; ++i) {
    t->GetEntry(i);
    if (i % 1000 == 0) std::cout << "Event " << i << '\n';

    double KE = 0.;

    KE += KEl;

    double KE_e_plus_Qgs = KEl + Q_ground_state;
    eq_Es->Fill(KE_e_plus_Qgs);

    for (int j = 0; j < np; ++j) {
      int pdg = (*pdgp)[j];
      if (pdg != NEUTRON) KE += (*KEp)[j];
      if (pdg > ALPHA) {
        if (pdg == all_gamma_nuc_pdg) ++all_gamma_count;
        else if (pdg == one_n_nuc_pdg) ++one_n_count;
        else if (pdg == one_p_nuc_pdg) ++one_p_count;
        else ++other_count;
      }
    }

    true_Es->Fill(Ev);
    double E_reco = KE + Q_ground_state;
    reco_Es->Fill(E_reco);
  }

  double event_count = static_cast<double>(num_events);
  double scale_factor = xsec / event_count / BIN_WIDTH;
  true_Es->Scale(scale_factor);
  reco_Es->Scale(scale_factor);
  eq_Es->Scale(scale_factor);

  TCanvas* c = new TCanvas;
  c->cd();

  TLegend* legend = new TLegend(0.65, 0.65, 0.9, 0.9);
  legend->AddEntry(true_Es, true_Es->GetTitle(), "l");
  legend->AddEntry(reco_Es, "E_{#nu,reco 1}", "l");
  legend->AddEntry(eq_Es, eq_Es->GetTitle(), "l");

  double max_reco = reco_Es->GetMaximum();
  double max_true = true_Es->GetMaximum();
  double max_eq = eq_Es->GetMaximum();
  double ymax = 1.05 * std::max(max_reco, std::max(max_true, max_eq));

  reco_Es->GetXaxis()->SetTitleOffset(1.3);
  reco_Es->GetYaxis()->SetTitleOffset(1.3);

  reco_Es->SetLineColor(kBlue);
  reco_Es->SetLineWidth(2);
  reco_Es->Draw("hist");
  reco_Es->GetYaxis()->SetRangeUser(0., ymax);
  true_Es->SetLineColor(kBlack);
  true_Es->SetLineWidth(2);
  true_Es->Draw("hist same");
  eq_Es->SetLineColor(kRed);
  eq_Es->SetLineWidth(2);
  eq_Es->Draw("hist same");

  legend->Draw();

  std::cout << "** Summary **" << '\n';
  std::cout << "e- + gammas only: " << all_gamma_count << " events ("
    << all_gamma_count / event_count * 100 << "%)" << '\n';
  std::cout << "single n: " << one_n_count << " events ("
    << one_n_count / event_count * 100 << "%)" << '\n';
  std::cout << "single p: " << one_p_count << " events ("
    << one_p_count / event_count * 100 << "%)" << '\n';
  std::cout << "other: " << other_count << " events ("
    << other_count / event_count * 100 << "%)" << '\n';
}
