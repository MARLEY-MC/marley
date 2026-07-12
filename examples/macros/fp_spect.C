#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"

void fp_spect(const std::string& file_name, int pdg) {

  TFile* f = TFile::Open(file_name.c_str());
  TTree* t = static_cast<TTree*>(f->Get("mst"));

  int pdgl, np;
  double KEl;
  std::vector<int>* pdgp = nullptr;
  std::vector<double>* KEp = nullptr;

  t->SetBranchAddress("pdgl", &pdgl);
  t->SetBranchAddress("KEl", &KEl);
  t->SetBranchAddress("np", &np);
  t->SetBranchAddress("pdgp", &pdgp);
  t->SetBranchAddress("KEp", &KEp);

  Long64_t n = t->GetEntries();
  std::vector<double> KE_vec;

  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    if (i % 1000 == 0) std::cout << "Event " << i << '\n';

    if (pdgl == pdg) KE_vec.push_back(KEl);

    for (int j = 0; j < np; ++j) {
      if ((*pdgp)[j] == pdg) KE_vec.push_back((*KEp)[j]);
    }
  }

  double KE_max = -1e30;
  double KE_min = 1e30;
  for (auto ke : KE_vec) {
    if (ke > KE_max) KE_max = ke;
    if (ke < KE_min) KE_min = ke;
  }

  TString title_str;
  title_str.Form("kinetic energies for pdg = %d; kinetic energy (MeV);"
    "events", pdg);

  TH1D* KEs = new TH1D("KEs", title_str.Data(), 100, KE_max, KE_min);

  KEs->SetDirectory(NULL);

  for (auto ke : KE_vec) {
    KEs->Fill(ke);
  }

  TCanvas* c = new TCanvas;
  c->cd();

  gStyle->SetOptStat();

  KEs->SetStats(false);
  KEs->SetLineColor(kBlue);
  KEs->SetLineWidth(2);
  KEs->Draw("hist");

  std::cout << "Found " << KE_vec.size() << " particles with"
    << " pdg = " << pdg << " in " << n << " events\n";
}
