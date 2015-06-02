#pragma once
#include <unordered_map>

class TMarleyMassTable {

  public:
    double get_particle_mass(int particle_id);
    double get_atomic_mass(int nucleus_pid);
    double get_atomic_mass(int Z, int A);

  private:
  // Factor to use when converting from micro-amu to MeV
  static constexpr double micro_amu = 0.000931494061;

  // Lookup table for particle masses. Keys are PDG particle
  // ID numbers, values are masses in micro-amu.
  static const std::unordered_map<int, double> particle_masses;

  // Lookup table for atomic masses. Keys are PDG particle
  // ID numbers for the nuclei, values are masses in micro-amu.
  static const std::unordered_map<int, double> atomic_masses;
};
