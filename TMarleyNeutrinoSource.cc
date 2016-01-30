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
