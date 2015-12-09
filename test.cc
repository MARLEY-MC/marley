#include <functional>
#include <iostream>

#include "TMarleyGenerator.hh"
#include "TMarleyNeutrinoSource.hh"

#include "TH1D.h"
#include "TCanvas.h"

int main() {

  TMarleyGenerator gen("config.txt");

  std::vector<double> bins = { 0., 10., 20. };
  std::vector<double> weights = { 1., 5., 1. };

  TMarleyGridNeutrinoSource source(bins, weights,
    marley_utils::ELECTRON_NEUTRINO,
    InterpolationGrid<double>::InterpolationMethod::LinearLinear);

  TH1D hist("hist", "hist; bin; counts",
    100, -10., 30.);

  std::function<double(double)> pdf = std::bind(&TMarleyGridNeutrinoSource::pdf, &source,
    std::placeholders::_1);

  //std::vector<double> energies;
  for (size_t j = 0; j < 1e6; ++j) {
    double e = gen.rejection_sample(pdf, source.get_Emin(),
      source.get_Emax());
    hist.Fill(e);
  }

  TCanvas canvas;
  hist.Draw();
  canvas.SaveAs("hist.pdf");

  canvas.Clear();
  return 0;
}
