#include <iostream>
#include <string>

void nu_spect(const std::string& filename) {

  TFile* file = new TFile(filename.c_str(), "read");
  TTree* tree = NULL;
  file->GetObject("MARLEY Event Tree", tree);
  if (!tree) {
    std::cout << "MARLEY event tree not found" << '\n';
    return;
  }

  marley::Event* ev = new marley::Event;
  tree->SetBranchAddress("events", &ev);

  size_t num_events = tree->GetEntries();

  std::vector<double> E_vec;

  for (size_t i = 0; i < num_events; ++i) {

    tree->GetEntry(i);

    E_vec.push_back(ev->projectile().total_energy());

    if (i % 1000 == 0) std::cout << "Event " << i << '\n';
  }

  double E_max = -1e30;
  double E_min = 1e30;
  for (size_t k = 0; k < E_vec.size(); ++k) {
    double e = E_vec.at(k);
    if (e > E_max) E_max = e;
    else if (e < E_min) E_min = e;
  }

  TString title_str;

  TH1D* Es = new TH1D("nu_Es", "reacting neutrino spectrum", 100,
    E_max, E_min);

  for (size_t j = 0; j < E_vec.size(); ++j) {
    Es->Fill(E_vec.at(j));
  }

  TCanvas* c = new TCanvas;
  c->cd();

  gStyle->SetOptStat();

  Es->SetStats(true);
  Es->SetLineColor(kBlue);
  Es->SetLineWidth(2);
  Es->Draw();

  //c->SaveAs("nu_Es.pdf");
}
