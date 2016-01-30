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

// Constants used for a decay at rest source
const double TMarleyDecayAtRestNeutrinoSource::m_mu =
  TMarleyMassTable::get_particle_mass(marley_utils::MUON);
const double TMarleyDecayAtRestNeutrinoSource::m_mu_to_the_minus_four
  = std::pow(m_mu, -4);
const double TMarleyDecayAtRestNeutrinoSource::E_max = m_mu / 2.;
