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
