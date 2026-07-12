#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

void nu_spect(const std::string& file_name) {

  const int NUM_BINS = 100;
  const double E_MIN = 0.;
  const double E_MAX = 60.;
  const double BIN_WIDTH = ( E_MAX - E_MIN ) / NUM_BINS;

  TH1D* Ev_hist = new TH1D("Ev_hist", "reacting neutrino spectrum;"
    "neutrino energy E_{#nu} (MeV); #left[ d#sigma/dE_{#nu} #right]_{flux}"
    " (10^{-42} cm^{2} / MeV)", NUM_BINS, E_MIN, E_MAX);

  TFile* f = TFile::Open(file_name.c_str());
  TTree* t = static_cast<TTree*>(f->Get("mst"));

  double Ev, xsec;
  t->SetBranchAddress("Ev", &Ev);
  t->SetBranchAddress("xsec", &xsec);

  Long64_t num_events = t->GetEntries();
  for (Long64_t i = 0; i < num_events; ++i) {
    t->GetEntry(i);
    if (i % 1000 == 0) std::cout << "Event " << i << '\n';
    Ev_hist->Fill(Ev);
  }

  double scale_factor = xsec / num_events / BIN_WIDTH;
  Ev_hist->Scale(scale_factor);

  TCanvas* c = new TCanvas;
  c->cd();

  gStyle->SetOptStat();
  Ev_hist->SetStats(false);
  Ev_hist->GetXaxis()->SetTitleOffset(1.2);
  Ev_hist->GetYaxis()->SetTitleOffset(1.2);
  Ev_hist->SetLineColor(kBlue);
  Ev_hist->SetLineWidth(2);
  Ev_hist->Draw("hist");
}
