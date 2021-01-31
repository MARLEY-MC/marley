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

#include <string>
#include <vector>
#include <regex>

#include "marley/marley_utils.hh"
#include "marley/Generator.hh"
#include "marley/Level.hh"
#include "marley/Parity.hh"

marley::Level::Level(double E, int twoJ, marley::Parity pi) : energy_(E),
  twoJ_(twoJ), parity_(pi) {}

const marley::Gamma* marley::Level::sample_gamma(marley::Generator& gen)
{
  if (gammas_.empty()) return nullptr;
  else {
    // Get the index of the gamma to return by randomly sampling from the
    // discrete distribution gamma_dist using the standard marley_utils random
    // number generator.
    size_t g_index = gen.sample_from_distribution(gamma_dist_);
    // Return a pointer to the corresponding gamma
    return &(gammas_[g_index]);
  }
}

marley::Gamma& marley::Level::add_gamma(const marley::Gamma& gamma) {
  // Update the vector of gamma objects
  gammas_.push_back(gamma);

  // Update the distribution for sampling gammas
  update_gamma_distribution();

  // Return a reference to the newly added gamma
  return gammas_.back();
}

marley::Gamma& marley::Level::add_gamma(double energy, double branching_ratio,
  marley::Level* end_lev)
{
  // Update the vector of gamma objects
  gammas_.emplace_back(energy, branching_ratio, this, end_lev);

  // Update the distribution for sampling gammas
  update_gamma_distribution();

  // Return a reference to the newly added gamma
  return gammas_.back();
}

void marley::Level::clear_gammas() {
  gammas_.clear();

  // The discrete distribution will be cleared by this command because the
  // vector of gammas is now empty.
  update_gamma_distribution();
}

std::string marley::Level::spin_parity_string() const {
  std::string str = std::to_string(twoJ_ / 2);
  // If 2*J is odd, then the level has half-integer spin
  if (twoJ_ % 2) str += "/2";
  return str + parity_.to_char();
}

void marley::Level::update_gamma_distribution() {
  // Get iterators to the relative intensities of the gammas owned by this
  // level.
  auto ri_begin = marley::Gamma::make_intensity_iterator(gammas_.begin());
  auto ri_end = marley::Gamma::make_intensity_iterator(gammas_.end());

  // Update the discrete distribution used to sample gammas
  std::discrete_distribution<size_t>::param_type params(ri_begin, ri_end);
  gamma_dist_.param(params);
}
