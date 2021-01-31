/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#include <memory>

#include "marley/marley_root.hh"

std::unique_ptr<marley::GridNeutrinoSource>
  marley_root::make_root_neutrino_source(int pdg_code, const TH1* th1)
{
  if (!th1) {
    throw marley::Error(std::string("Null TH1* passed")
      + " to marley_root::make_root_neutrino_source");
    return std::unique_ptr<marley::GridNeutrinoSource>(nullptr);
  }

  // Read in the (energy low edge, bin weight) ordered pairs. Keep the
  // overflow bin so that we can use its left edge as the maximum energy
  // value. Skip the unneeded underflow bin.
  size_t n_bins = th1->GetNbinsX() + 1;
  std::vector<double> Es(n_bins);
  std::vector<double> PDs(n_bins);

  // underflow bin is bin 0
  for (size_t b = 1; b <= n_bins; ++b) {
    Es.at(b - 1) = th1->GetBinLowEdge(b);
    // the content of the overflow bin is reset to zero below
    PDs.at(b - 1) = th1->GetBinContent(b);
  }

  // Convert the bin weights to probability densities (used by our
  // GridNeutrinoSource object) by dividing each bin weight by the
  // corresponding bin width.
  for (size_t c = 0; c < n_bins - 1; ++c) {
    double width = Es.at(c + 1) - Es.at(c);
    if (width <= 0) throw marley::Error(std::string("Invalid bin width")
      + std::to_string(width) + " encountered when creating a TH1 neutrino"
      + " source");
    PDs.at(c) /= width;
  }

  // Assign zero probability density to the overflow bin's
  // left edge. This ensures that neutrino energies will be
  // sampled on the half-open interval [Elow, Ehigh), where
  // Elow is the left edge of the first bin and Ehigh is the
  // left edge of the overflow bin.
  PDs.back() = 0.;

  // Now that we've processed grid points, create the grid neutrino
  // source
  auto source = std::make_unique<marley::GridNeutrinoSource>(
    Es, PDs, pdg_code, marley::InterpolationGrid<double>
    ::InterpolationMethod::Constant);
  return source;
}

std::unique_ptr<marley::GridNeutrinoSource>
  marley_root::make_root_neutrino_source(int pdg_code, const TGraph* tg)
{
  if (!tg) {
    throw marley::Error(std::string("Null TGraph* passed")
      + " to marley_root::make_root_neutrino_source");
    return std::unique_ptr<marley::GridNeutrinoSource>(nullptr);
  }

  size_t num_points = tg->GetN();
  std::vector<double> Es(num_points);
  std::vector<double> PDs(num_points);

  // Load the energies and PDF values into our vectors
  for(size_t p = 0; p < num_points; ++p) {
    tg->GetPoint(p, Es.at(p), PDs.at(p));
  }

  // Create a neutrino source based on the grid
  auto source = std::make_unique<marley::GridNeutrinoSource>(Es, PDs,
    pdg_code, marley::InterpolationGrid<double>::InterpolationMethod
    ::LinearLinear);
  return source;
}
