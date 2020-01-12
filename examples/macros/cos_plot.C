#include <iostream>
#include <string>

const int NUM_BINS = 10;
const double COS_MIN = -1.;
const double COS_MAX = 1.;
const double BIN_WIDTH = ( COS_MAX - COS_MIN ) / NUM_BINS;

void cos_plot(const std::string& file_name) {

  TH1D* cos_theta_hist = new TH1D("cos_theta_hist",
    "scattering cosine distribution; cos#theta; #left[ d#sigma"
    "/dcos#theta #right]_{flux} (10^{-42} cm^{2})", NUM_BINS,
    COS_MIN, COS_MAX);

  marley::MacroEventFileReader reader( file_name );
  marley::Event ev;

  int num_events = 0;
  while ( reader >> ev ) {

    if ( num_events % 1000 == 0 ) std::cout << "Event " << num_events << '\n';

    const marley::Particle& nu = ev.projectile();
    const marley::Particle& e = ev.ejectile();

    double pnu_dot_pe = nu.px() * e.px() + nu.py() * e.py()
      + nu.pz() * e.pz();

    double norm_pnu = nu.momentum_magnitude();
    double norm_pe = e.momentum_magnitude();

    double cos_theta = pnu_dot_pe / ( norm_pnu * norm_pe );

    cos_theta_hist->Fill( cos_theta );

    ++num_events;
  }

  double xsec = reader.flux_averaged_xsec();
  double scale_factor = xsec / BIN_WIDTH / num_events;

  cos_theta_hist->Scale( scale_factor );

  TCanvas* c = new TCanvas;
  c->cd();

  cos_theta_hist->GetXaxis()->SetTitleOffset(1.2);
  cos_theta_hist->GetYaxis()->SetTitleOffset(1.2);
  cos_theta_hist->SetStats( false );
  cos_theta_hist->SetLineColor( kBlue );
  cos_theta_hist->SetLineWidth( 2 );
  cos_theta_hist->Draw("hist");
}
