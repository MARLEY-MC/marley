#pragma once
#include <unordered_map>

class TMarleyMassTable {

  public:
    static double get_particle_mass(int particle_id);
    static double get_atomic_mass(int nucleus_pid);
    static double get_atomic_mass(int Z, int A);
    static double get_proton_separation_energy(int Z, int A);
    static double get_neutron_separation_energy(int Z, int A);
    static double get_particle_separation_energy(int Z, int A, int pid);
    static void print_separation_energies(int Z, int A, unsigned n);

    // Convert the atomic number Z and the mass number A into
    // a PDG-compliant particle ID number
    static inline int get_nucleus_pid(int Z, int A) {
      return 10000*Z + 10*A + 1000000000;
    }

    static inline int get_particle_Z(int pid) {
      // proton
      if (pid == 2212) return 1;
      // neutron
      else if (pid == 2112) return 0;
      // nuclear fragment
      else if (pid > 1000000000) return (pid % 10000000)/10000;
      // other particle
      else return 0;
    }

    static inline int get_particle_A(int pid) {
      // proton
      if (pid == 2212) return 1;
      // neutron
      else if (pid == 2112) return 1;
      // nuclear fragment
      else if (pid > 1000000000) return (pid % 10000)/10;
      // other particle
      else return 0;
    }

  private:
    // Factor to use when converting from micro-amu to MeV
    static constexpr double micro_amu = 0.000931494061;

    // Lookup table for particle masses. Keys are PDG particle
    // ID numbers, values are masses in micro-amu.
    static const std::unordered_map<int, double> particle_masses;

    // Lookup table for atomic masses. Keys are PDG particle
    // ID numbers for the nuclei, values are masses in micro-amu.
    static const std::unordered_map<int, double> atomic_masses;



    // Particle IDs for all of the nuclear fragments that will
    // be considered when calculating separation energies
    static const std::vector<int> fragment_pids;

    static void iterate_multicombination_counts(unsigned n_objects,
      std::vector<unsigned>& vec, std::function<void(std::vector<unsigned>&)> f,
      unsigned n_boxes);

    static void iterate_multicombination_counts(unsigned n_boxes,
      unsigned n_objects, std::function<void(std::vector<unsigned>&)> f);

    static void print_separation_energy(int Z, int A, std::vector<unsigned>& fragment_counts);
};
