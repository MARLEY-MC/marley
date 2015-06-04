#pragma once
#include <unordered_map>

class TMarleyMassTable {

  public:
    static double get_particle_mass(int particle_id);
    static double get_atomic_mass(int nucleus_pid);
    static double get_atomic_mass(int Z, int A);
    static double get_proton_separation_energy(int Z, int A);
    static double get_neutron_separation_energy(int Z, int A);

  private:
  // Factor to use when converting from micro-amu to MeV
  static constexpr double micro_amu = 0.000931494061;

  // Lookup table for particle masses. Keys are PDG particle
  // ID numbers, values are masses in micro-amu.
  static const std::unordered_map<int, double> particle_masses;

  // Lookup table for atomic masses. Keys are PDG particle
  // ID numbers for the nuclei, values are masses in micro-amu.
  static const std::unordered_map<int, double> atomic_masses;

  // Convert the atomic number Z and the mass number A into
  // a PDG-compliant particle ID number
  inline static int get_nucleus_pid(int Z, int A) {
    return 10000*Z + 10*A + 1000000000;
  }
};
