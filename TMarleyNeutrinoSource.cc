#include "TMarleyGenerator.hh"
#include "TMarleyNeutrinoSource.hh"

// Particle IDs of each neutrino that could possibly be produced by a
// source object
const std::set<int> TMarleyNeutrinoSource::pids = {
  marley_utils::ELECTRON_NEUTRINO,
  marley_utils::ELECTRON_ANTINEUTRINO,
  marley_utils::MUON_NEUTRINO,
  marley_utils::MUON_ANTINEUTRINO,
  marley_utils::TAU_NEUTRINO,
  marley_utils::TAU_ANTINEUTRINO
};

//double TMarleyFermiDiracNeutrinoSource::sample_energy(TMarleyGenerator& gen)
//{
//  return gen.rejection_sample(fd_dist, E_min, E_max);
//}
//
//double TMarleyHistogramNeutrinoSource::sample_energy(TMarleyGenerator& gen) {
//  return gen.piecewise_constant_sample(energy_dist);
//}
//
//double TMarleyGridNeutrinoSource::sample_energy(TMarleyGenerator& gen) {
//  return gen.piecewise_linear_sample(energy_dist);
//}

//double TMarleyHistogram::sample_value(TMarleyGenerator& gen) {
//  size_t bin_index = gen.discrete_sample(bin_dist);
//  double xmin = lower_bounds.at(bin_index);
//  double xmax;
//  if (bin_index < last_bin_index) xmax = lower_bounds.at(bin_index + 1);
//  else if (bin_index == last_bin_index) xmax = x_max;
//  else throw std::runtime_error(std::string("Bin")
//    + " indexing error in TMarleyHistogram::sample_value()");
//  return gen.uniform_random_double(xmin, xmax, true);
//}
