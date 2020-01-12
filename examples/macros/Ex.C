#include <cfloat>
#include <iostream>
#include <string>

void Ex(const std::string& filename) {

  double Ex_max = DBL_MIN;

  std::vector<double> Ex_values;

  marley::R5EFR efr( filename );
  marley::Event ev;

  int counter = 0;

  while ( efr >> ev ) {

    if (counter % 1000 == 0) std::cout << "Event " << counter << '\n';

    double Ex = ev.Ex();
    if ( Ex < Ex_max ) Ex_max = Ex;

    Ex_values.push_back( Ex );

    ++counter;
  }

  TH1D* Ex_hist = new TH1D("Ex_hist", "Nuclear excitation energies;"
    " E_{x} (MeV); events", 100, Ex_max*1.1, 0.);
  Ex_hist->SetDirectory(NULL);

  for ( size_t j = 0u; j < Ex_values.size(); ++j ) {
    Ex_hist->Fill( Ex_values.at(j) );
  }

  TCanvas* c = new TCanvas;
  c->cd();

  Ex_hist->SetStats( false );
  Ex_hist->SetLineColor( kBlack );
  Ex_hist->SetLineWidth( 2 );
  Ex_hist->GetXaxis()->SetTitleOffset(1.2);
  Ex_hist->GetYaxis()->SetTitleOffset(1.2);
  Ex_hist->Draw("hist");

}
