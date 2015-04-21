#include "TMarleyGamma.hh"

TMarleyGamma::TMarleyGamma(double energy, double ri, TMarleyLevel* start_level) {
  fEnergy = energy;
  fRI = ri;
  fCC = 0;
  fTI = 0;
  pStartLevel = start_level;
  pEndLevel = nullptr;
}

void TMarleyGamma::set_start_level(TMarleyLevel* start_level) {
  pStartLevel = start_level;
}

void TMarleyGamma::set_end_level(TMarleyLevel* end_level) {
  pEndLevel = end_level;
}

TMarleyLevel* TMarleyGamma::get_start_level() const {
  return pStartLevel;
}

TMarleyLevel* TMarleyGamma::get_end_level() const {
  return pEndLevel;
}

double TMarleyGamma::get_energy() const {
  return fEnergy;
}

double TMarleyGamma::get_ri() const {
  return fRI;
}
