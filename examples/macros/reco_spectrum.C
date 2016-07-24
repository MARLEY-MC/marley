#include <algorithm>
#include <iostream>
#include <string>

void reco_spectrum(const std::string& filename) {

  // Difference between ground-state masses of 40K and 40Ar
  const double Q_ground_state = 1.5044; // MeV

  TFile* file = new TFile(filename.c_str(), "read");
  TTree* tree = NULL;
  file->GetObject("MARLEY Event Tree", tree);
  if (!tree) {
    std::cout << "MARLEY event tree not found" << '\n';
    return;
  }

  marley::Event* ev = new marley::Event;
  tree->SetBranchAddress("events", &ev);

  TH1D* true_Es = new TH1D("true_Es", "true neutrino energies", 110, 0., 60.);
  TH1D* eq_Es = new TH1D("eq_Es", "e- KE + g.s.->g.s. Q value", 110, 0., 60.);
  TH1D* reco_Es = new TH1D("reco_Es", "reco neutrino energies", 110, 0., 60.);

  true_Es->SetStats(false);
  eq_Es->SetStats(false);
  reco_Es->SetStats(false);

  size_t num_events = tree->GetEntries();

  TVectorD* true_E_vec = new TVectorD(num_events);
  TVectorD* eq_E_vec = new TVectorD(num_events);
  TVectorD* reco_E_vec = new TVectorD(num_events);

  const int all_gamma_nuc_pid = 1000190400; // 40K
  const int one_n_nuc_pid = 1000190390; // 39K
  const int one_p_nuc_pid = 1000180390; // 39Ar

  size_t all_gamma_count = 0;
  size_t one_n_count = 0;
  size_t one_p_count = 0;
  size_t other_count = 0;

  for (size_t i = 0; i < num_events; ++i) {

    tree->GetEntry(i);

    double KE = 0.;

    marley::Particle* fp = NULL;

    const std::vector<marley::Particle*>& finals = ev->get_final_particles();
    for (size_t j = 0; j < finals.size(); ++j) {
      const marley::Particle* fp = finals.at(j);
      int pid = fp->pdg_code();
      if (pid != 2112 /* neutron */) KE += fp->kinetic_energy();
      if (pid == 11 /* electron */) {
        double KE_e_plus_Qgs = fp->kinetic_energy() + Q_ground_state; // e- KE + Q_gs
        eq_Es->Fill(KE_e_plus_Qgs);
        //eq_E_vec->operator[](i) = KE_e_plus_Qgs;
        (*eq_E_vec)[i] = KE_e_plus_Qgs;
      }
      if (pid > 1000020040 /* alpha */) {
        // nucleus
        if (pid == all_gamma_nuc_pid) ++all_gamma_count;
        else if (pid == one_n_nuc_pid) ++one_n_count;
        else if (pid == one_p_nuc_pid) ++one_p_count;
        else ++other_count;
      }
    }

    double E_true = ev->projectile().kinetic_energy();
    double E_reco = KE + Q_ground_state;

    true_Es->Fill(E_true);
    reco_Es->Fill(E_reco);
    if (i % 1000 == 0) std::cout << "Event " << i << '\n';

    (*true_E_vec)[i] = E_true;
    (*reco_E_vec)[i] = E_reco;
  }

  TCanvas* c = new TCanvas;
  c->cd();

  TLegend* legend = new TLegend(0.8, 0.8, 1., 1.);
  legend->AddEntry(true_Es, true_Es->GetTitle(), "l");
  legend->AddEntry(reco_Es, "#splitline{perfect reconstruction"
    " missing all neutrons}{and nuclear fragment Q-values}", "l");
  legend->AddEntry(eq_Es, eq_Es->GetTitle(), "l");

  double max_reco = reco_Es->GetMaximum();
  double max_true = true_Es->GetMaximum();
  double max_eq = eq_Es->GetMaximum();

  double ymax = 1.05*std::max(max_reco, std::max(max_true, max_eq));

  reco_Es->SetLineColor(kBlue);
  reco_Es->SetLineWidth(2);
  reco_Es->Draw();
  reco_Es->GetYaxis()->SetRangeUser(0., ymax);
  true_Es->SetLineColor(kBlack);
  true_Es->SetLineWidth(2);
  true_Es->Draw("same");
  eq_Es->SetLineColor(kRed);
  eq_Es->SetLineWidth(2);
  eq_Es->Draw("same");

  legend->Draw();

  //c->SaveAs("results.pdf");

  //TFile* outfile = new TFile("results.root", "recreate");
  //outfile->cd();

  //true_Es->Write();
  //reco_Es->Write();
  //eq_Es->Write();
  //outfile->WriteTObject(true_E_vec, "true_E_vec");
  //outfile->WriteTObject(reco_E_vec, "reco_E_vec");
  //outfile->WriteTObject(eq_E_vec, "eq_E_vec");
  //outfile->Close();

  std::cout << "** Summary **" << '\n';
  // One operand must be a floating-point number to get a float
  // when doing division
  double event_num = static_cast<double>(num_events);
  std::cout << "e- + gammas only: " << all_gamma_count << " events ("
    << all_gamma_count / event_num * 100
    << "%)" << '\n';
  std::cout << "single n: " << one_n_count << " events ("
    << one_n_count / event_num * 100
    << "%)" << '\n';
  std::cout << "single p: " << one_p_count << " events ("
    << one_p_count / event_num * 100
    << "%)" << '\n';
  std::cout << "other: " << other_count << " events ("
    << other_count / event_num * 100
    << "%)" << '\n';
}
