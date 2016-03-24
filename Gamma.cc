#include "Gamma.hh"

marley::Gamma::Gamma(double energy, double ri, marley::Level* start_level) {
  fEnergy = energy;
  fRI = ri;
  pStartLevel = start_level;
  pEndLevel = nullptr;
}

void marley::Gamma::set_start_level(marley::Level* start_level) {
  pStartLevel = start_level;
}

void marley::Gamma::set_end_level(marley::Level* end_level) {
  pEndLevel = end_level;
}

marley::Level* marley::Gamma::get_start_level() const {
  return pStartLevel;
}

marley::Level* marley::Gamma::get_end_level() const {
  return pEndLevel;
}

double marley::Gamma::get_energy() const {
  return fEnergy;
}

double marley::Gamma::get_ri() const {
  return fRI;
}
