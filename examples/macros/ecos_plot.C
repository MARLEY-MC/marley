#include <iostream>
#include <string>

const double Ex_IAS = 4.384;

void ecos_plot(const std::string& filename) {

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

  std::vector<double> e_cos_theta_vec;

  for (size_t i = 0; i < num_events; ++i) {

    tree->GetEntry(i);

    const marley::Particle& nu = ev->projectile();
    const marley::Particle& e = ev->ejectile();

    double pnu_dot_pe = nu.px() * e.px() + nu.py() * e.py()
      + nu.pz() * e.pz();

    double norm_pnu = nu.momentum_magnitude();
    double norm_pe = e.momentum_magnitude();

    double cos_theta = pnu_dot_pe / (norm_pnu * norm_pe);
    /*if (ev.Ex() == Ex_IAS)*/ e_cos_theta_vec.push_back(cos_theta);

    if (i % 1000 == 0) std::cout << "Event " << i << '\n';
  }

  TH1D* e_cos_theta_hist = new TH1D("e_cos_theta_hist",
    "electron cos(#theta)", 100, -1., 1.);

  for (size_t j = 0; j < e_cos_theta_vec.size(); ++j) {
    e_cos_theta_hist->Fill(e_cos_theta_vec.at(j));
  }

  TCanvas* c = new TCanvas;
  c->cd();

  e_cos_theta_hist->SetStats(true);
  e_cos_theta_hist->SetLineColor(kBlue);
  e_cos_theta_hist->SetLineWidth(2);
  e_cos_theta_hist->Draw();

  //c->SaveAs("e_cos_theta_hist.pdf");

  //TFile* outfile = new TFile("e_cos_theta_hist.root", "recreate");
  //outfile->cd();
}
