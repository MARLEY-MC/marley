#include <iostream>
#include <string>

const int NUM_BINS = 100;
const double E_MIN = 0.;
const double E_MAX = 60.;
const double BIN_WIDTH = ( E_MAX - E_MIN ) / NUM_BINS;

void nu_spect(const std::string& file_name) {

  TH1D* Ev_hist = new TH1D("Ev_hist", "reacting neutrino spectrum;"
    "neutrino energy E_{#nu} (MeV); #left[ d#sigma/dE_{#nu} #right]_{flux}"
    " (10^{-42} cm^{2} / MeV)", NUM_BINS, E_MIN, E_MAX);

  marley::MacroEventFileReader reader( file_name );
  marley::Event ev;

  long num_events = 0;

  while ( reader >> ev ) {

    if (num_events % 1000 == 0) std::cout << "Event " << num_events << '\n';

    Ev_hist->Fill( ev.projectile().total_energy() );

    ++num_events;

  }

  double xsec = reader.flux_averaged_xsec();
  double scale_factor = xsec / num_events / BIN_WIDTH;

  Ev_hist->Scale( scale_factor );

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
