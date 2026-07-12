#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

void Ex(const std::string& filename) {

  TFile* f = TFile::Open(filename.c_str());
  TTree* t = static_cast<TTree*>(f->Get("mst"));

  double Ex;
  t->SetBranchAddress("Ex", &Ex);

  Long64_t n = t->GetEntries();

  double Ex_max = -1;
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (Ex > Ex_max) Ex_max = Ex;
  }

  TH1D* Ex_hist = new TH1D("Ex_hist", "Nuclear excitation energies;"
    " E_{x} (MeV); events", 100, Ex_max * 1.1, 0.);
  Ex_hist->SetDirectory(NULL);

  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (i % 1000 == 0) std::cout << "Event " << i << '\n';
    Ex_hist->Fill(Ex);
  }

  TCanvas* c = new TCanvas;
  c->cd();

  Ex_hist->SetStats(false);
  Ex_hist->SetLineColor(kBlack);
  Ex_hist->SetLineWidth(2);
  Ex_hist->GetXaxis()->SetTitleOffset(1.2);
  Ex_hist->GetYaxis()->SetTitleOffset(1.2);
  Ex_hist->Draw("hist");
}
