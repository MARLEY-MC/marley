#include <functional>
#include <iostream>

#include "Generator.hh"
#include "NeutrinoSource.hh"

#include "TH1D.h"
#include "TCanvas.h"

int main() {

  marley::Generator gen("config.txt");


  double m_mu = marley::MassTable::get_particle_mass(marley_utils::MUON);

  TH1D hist("neutrino_E_histogram", "Stopped #pi beam #nu_{e} spectrum; Total Energy [MeV]; Counts",
    100, 0., m_mu/2 * 1.1);

  marley::NeutrinoSource* source = gen.get_nu_source();

  std::function<double(double)> pdf = [source](double E)
    -> double { return source->pdf(E); };

  //std::vector<double> energies;
  for (size_t j = 0; j < 1e6; ++j) {
    double e = gen.rejection_sample(pdf, source->get_Emin(),
      source->get_Emax());
    hist.Fill(e);
  }

  TCanvas canvas;
  hist.Draw();
  hist.SetLineWidth(1.8);
  hist.GetXaxis()->SetTitleOffset(1.2);
  hist.GetYaxis()->SetTitleOffset(1.5);
  hist.GetXaxis()->CenterTitle();
  hist.GetYaxis()->CenterTitle();
  canvas.SaveAs("hist.pdf");

  canvas.Clear();
  return 0;
}
