#include "TMarleyEvaporationThreshold.hh"

TMarleyEvaporationThreshold::TMarleyEvaporationThreshold(double E_separation,
  int pid)
{
  separation_energy = E_separation;
  fragment_pid = pid;
}

double TMarleyEvaporationThreshold::get_separation_energy() const {
  return separation_energy;
}

int TMarleyEvaporationThreshold::get_fragment_pid() const {
  return fragment_pid;
}

void TMarleyEvaporationThreshold::set_separation_energy(double E_separation) {
  separation_energy = E_separation;
}

void TMarleyEvaporationThreshold::set_fragment_pid(int pid) {
  fragment_pid = pid;
}
