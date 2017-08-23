#include <iostream>
#include <string>

void fp_spect(const std::string& filename, int pdg) {

  TFile file(filename.c_str(), "read");
  TTree* tree = NULL;
  file.GetObject("MARLEY_event_tree", tree);
  if (!tree) {
    std::cout << "MARLEY event tree not found" << '\n';
    return;
  }

  marley::Event* ev = NULL;
  tree->SetBranchAddress("event", &ev);

  size_t num_events = tree->GetEntries();

  std::vector<double> KE_vec;

  for (size_t i = 0; i < num_events; ++i) {

    tree->GetEntry(i);

    const std::vector<marley::Particle*>& finals = ev->get_final_particles();

    for (size_t j = 0; j < finals.size(); ++j) {
      const marley::Particle* fp = finals.at(j);
      if (fp->pdg_code() == pdg) {
        KE_vec.push_back(fp->kinetic_energy());
      }
    }

    if (i % 1000 == 0) std::cout << "Event " << i << '\n';
  }

  double KE_max = -1e30;
  double KE_min = 1e30;
  for (size_t k = 0; k < KE_vec.size(); ++k) {
    double ke = KE_vec.at(k);
    if (ke > KE_max) KE_max = ke;
    else if (ke < KE_min) KE_min = ke;
  }

  TString title_str;
  title_str.Form("kinetic energies for pdg = %d", pdg);

  TH1D* KEs = new TH1D("KEs", title_str.Data(), 100, KE_max, KE_min);
  // Prevent this histogram from being automatically deleted when the
  // current TFile is closed
  KEs->SetDirectory(NULL);

  for (size_t j = 0; j < KE_vec.size(); ++j) {
    KEs->Fill(KE_vec.at(j));
  }

  TCanvas* c = new TCanvas;
  c->cd();

  gStyle->SetOptStat();

  KEs->SetStats(true);
  KEs->SetLineColor(kBlue);
  KEs->SetLineWidth(2);
  KEs->Draw();

  //c->SaveAs("KEs.pdf");

  std::cout << "Found " << KE_vec.size() << " particles with"
    << " pdg = " << pdg << " in " << num_events << " events\n";
}
