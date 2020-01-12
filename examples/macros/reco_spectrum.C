#include <algorithm>
#include <iostream>
#include <string>

// Histogram settings
const int NUM_BINS = 110;
const double E_MIN = 0.; // MeV
const double E_MAX = 55.; // MeV
const double BIN_WIDTH = ( E_MAX - E_MIN ) / NUM_BINS;

// PDG codes
const int NEUTRON = 2112;
const int ELECTRON = 11;
const int ALPHA = 1000020040;

void reco_spectrum(const std::string& filename) {

  // Difference between ground-state masses of 40K and 40Ar
  const double Q_ground_state = 1.5044; // MeV

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

  marley::R5EFR efr( filename );
  marley::Event ev;

  long num_events = 0;

  while ( efr >> ev ) {

    if (num_events % 1000 == 0) std::cout << "Event " << num_events << '\n';

    double KE = 0.;

    marley::Particle* fp = NULL;

    const std::vector<marley::Particle*>& finals = ev.get_final_particles();
    for (size_t j = 0; j < finals.size(); ++j) {
      const marley::Particle* fp = finals.at(j);
      int pdg = fp->pdg_code();
      if ( pdg != NEUTRON ) KE += fp->kinetic_energy();
      if ( pdg == ELECTRON ) {
        double KE_e_plus_Qgs = fp->kinetic_energy() + Q_ground_state; // e- KE + Q_gs
        eq_Es->Fill( KE_e_plus_Qgs );
      }
      if ( pdg > ALPHA ) {
        // nucleus
        if (pdg == all_gamma_nuc_pdg) ++all_gamma_count;
        else if (pdg == one_n_nuc_pdg) ++one_n_count;
        else if (pdg == one_p_nuc_pdg) ++one_p_count;
        else ++other_count;
      }
    }

    double E_true = ev.projectile().total_energy();
    double E_reco = KE + Q_ground_state;

    true_Es->Fill( E_true );
    reco_Es->Fill( E_reco );

    ++num_events;
  }

  // One operand must be a floating-point number to get a float
  // when doing division
  double event_count = static_cast<double>(num_events);

  // Normalize the event histograms so that they represent
  // differential cross sections
  double xsec = efr.flux_averaged_xsec(); // 10^{-42} cm^2
  double scale_factor = xsec / event_count / BIN_WIDTH;
  true_Es->Scale( scale_factor );
  reco_Es->Scale( scale_factor );
  eq_Es->Scale( scale_factor );

  TCanvas* c = new TCanvas;
  c->cd();

  TLegend* legend = new TLegend(0.65, 0.65, 0.9, 0.9);
  legend->AddEntry(true_Es, true_Es->GetTitle(), "l");
  legend->AddEntry(reco_Es, "E_{#nu,reco 1}", "l");
  legend->AddEntry(eq_Es, eq_Es->GetTitle(), "l");

  double max_reco = reco_Es->GetMaximum();
  double max_true = true_Es->GetMaximum();
  double max_eq = eq_Es->GetMaximum();

  double ymax = 1.05*std::max(max_reco, std::max(max_true, max_eq));

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
    << all_gamma_count / event_count * 100
    << "%)" << '\n';
  std::cout << "single n: " << one_n_count << " events ("
    << one_n_count / event_count * 100
    << "%)" << '\n';
  std::cout << "single p: " << one_p_count << " events ("
    << one_p_count / event_count * 100
    << "%)" << '\n';
  std::cout << "other: " << other_count << " events ("
    << other_count / event_count * 100
    << "%)" << '\n';
}
