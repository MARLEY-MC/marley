#ifndef USE_ROOT
#error "This MARLEY analysis code currently must be compiled with ROOT support."
#endif

#include <vector>

#include "TMarleyGenerator.hh"
#include "TMarleyNeutrinoSource.hh"

#include "TH1D.h"
#include "TCanvas.h"
#include "TColor.h"
//#include "TLegend.h"

int main() {

  size_t num_trials = 1e6;

  TMarleyGenerator gen("config.txt");

  TMarleyNeutrinoSource& source = gen.get_nu_source();

  // Sample a bunch of neutrino energies for the histogram plot
  std::vector<double> nu_Es;
  double nu_E_max = 0.;
  for (size_t j = 0; j < num_trials; ++j) {
    double nu_E = source.sample_neutrino_energy(gen);
    if (nu_E_max < nu_E) nu_E_max = nu_E;
    nu_Es.push_back(nu_E);
  }

  // Create a histogram to store the neutrino energies.
  // Use 100 bins, and adjust the upper limit of the last bin so that
  // the maximum observed value just barely falls within it.
  TH1D neutrino_E_histogram("neutrino_E_histogram",
    "Fermi-Dirac  #nu_{e}  spectrum; Total Energy [MeV]; Counts",
    100, 0, 1.1*std::nextafter(nu_E_max, DBL_MAX));

  // Fill the histogram with the sampled neutrino energies
  for(const double e : nu_Es) neutrino_E_histogram.Fill(e);

  TCanvas canvas;
  canvas.Clear();
  neutrino_E_histogram.Draw();
  neutrino_E_histogram.SetLineWidth(1.8);
  neutrino_E_histogram.GetXaxis()->SetTitleOffset(1.2);
  neutrino_E_histogram.GetYaxis()->SetTitleOffset(1.5);
  neutrino_E_histogram.GetXaxis()->CenterTitle();
  neutrino_E_histogram.GetYaxis()->CenterTitle();
  canvas.SaveAs("E_neutrino_source.pdf");

  return 0;
}
