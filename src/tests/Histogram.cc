/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#include "marley/tests/Histogram.hh"

// Operator for writing a Histogram to a std::ostream
std::ostream& operator<<(std::ostream& out, const marley::tests::Histogram& h) {
  h.print( out );
  return out;
}

// Operator for reading in a Histogram from a std::istream
std::istream& operator>>(std::istream& in, marley::tests::Histogram& h) {
  in >> h.N_bins_;
  in >> h.x_min_;
  in >> h.x_max_;
  in >> h.num_entries_;
  in >> h.x_step_;
  h.bin_counts_.clear();

  int count;
  for (size_t b = 0; b < h.N_bins_; ++b) {
    in >> count;
    h.bin_counts_.push_back( count );
  }

  return in;
}
