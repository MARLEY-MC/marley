#include "Generator.hh"
#include "NeutrinoSource.hh"

// Particle IDs of each neutrino that could possibly be produced by a
// source object
const std::set<int> marley::NeutrinoSource::pids = {
  marley_utils::ELECTRON_NEUTRINO,
  marley_utils::ELECTRON_ANTINEUTRINO,
  marley_utils::MUON_NEUTRINO,
  marley_utils::MUON_ANTINEUTRINO,
  marley_utils::TAU_NEUTRINO,
  marley_utils::TAU_ANTINEUTRINO
};

// Constants used for a decay at rest source
const double marley::DecayAtRestNeutrinoSource::m_mu =
  marley::MassTable::get_particle_mass(marley_utils::MUON);
const double marley::DecayAtRestNeutrinoSource::m_mu_to_the_minus_four
  = std::pow(m_mu, -4);
const double marley::DecayAtRestNeutrinoSource::E_max = m_mu / 2.;
