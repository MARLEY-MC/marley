#include <iostream>
#include <string>

void Ex(const std::string& filename) {

  TFile* file = new TFile(filename.c_str(), "read");
  TTree* tree = NULL;
  file->GetObject("MARLEY_event_tree", tree);
  if (!tree) {
    std::cout << "MARLEY event tree not found" << '\n';
    return;
  }

  marley::Event* ev = new marley::Event;
  tree->SetBranchAddress("event", &ev);

  size_t num_events = tree->GetEntries();

  TVectorD* Ex_vec = new TVectorD(num_events);

  for (size_t i = 0; i < num_events; ++i) {

    tree->GetEntry(i);

    double Ex = ev->Ex();
    (*Ex_vec)[i] = Ex;

    if (i % 1000 == 0) std::cout << "Event " << i << '\n';
  }

  TH1D* Ex_hist = new TH1D("Ex_hist", "residual nucleus excitation energies",
    100, Ex_vec->Max()*1.1, Ex_vec->Min()*0.9);

  for (size_t j = 0; j < Ex_vec->GetNoElements(); ++j) {
    Ex_hist->Fill((*Ex_vec)[j]);
  }

  TCanvas* c = new TCanvas;
  c->cd();

  gStyle->SetOptStat();

  Ex_hist->SetStats(false);
  Ex_hist->SetLineColor(kBlue);
  Ex_hist->SetLineWidth(2);
  Ex_hist->Draw();

  //c->SaveAs("Ex_hist.pdf");

  //TFile* outfile = new TFile("Ex_hist.root", "recreate");
  //outfile->cd();

  //Ex_hist->Write();
  //outfile->WriteTObject(Ex_vec, "Ex_vec");
  //outfile->Close();
}
