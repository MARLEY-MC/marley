#pragma once
#include <unordered_map>
#include <functional>
#include <vector>

#include "marley_utils.hh"

namespace marley {

  // TODO: Consider using a namespace (possibly with an anonymous namespace
  // in the source file to encapsulate the current private members)
  // rather than a class for the mass table. For more thoughts on this
  // possibility, see http://www.cplusplus.com/forum/general/121429/
  class MassTable {
  
    public:
      static double get_particle_mass(int particle_id);
      static double get_atomic_mass(int nucleus_pid, bool theory_ok = true);
      static double get_atomic_mass(int Z, int A, bool theory_ok = true);
      static double get_fragment_separation_energy(int Z, int A, int pid,
        bool theory_ok = true);
      static double get_binding_energy(int Z, int A, bool theory_ok = true);
      static double get_mass_excess(int Z, int A, bool theory_ok = true);
      static double liquid_drop_model_mass_excess(int Z, int A);
  
      static inline double liquid_drop_model_atomic_mass(int Z, int A) {
        return liquid_drop_model_mass_excess(Z, A) + micro_amu*1e6*A;
      }
  
      static inline int get_particle_Z(int pid) {
        if (pid == marley_utils::PROTON) return 1;
        else if (pid == marley_utils::NEUTRON) return 0;
        // nuclear fragment
        else if (pid > 1000000000) return (pid % 10000000)/10000;
        // other particle
        else return 0;
      }
  
      static inline int get_particle_A(int pid) {
        if (pid == marley_utils::PROTON) return 1;
        else if (pid == marley_utils::NEUTRON) return 1;
        // nuclear fragment
        else if (pid > 1000000000) return (pid % 10000)/10;
        // other particle
        else return 0;
      }
  
      static inline const std::vector<int>& get_fragment_pids() {
        return fragment_pids;
      }
  
      // Factor to use when converting from micro-amu to MeV
      static constexpr double micro_amu = 0.000931494061;
  
    private:
      // Function used internally by the mass table. Returns an atomic mass and
      // loads the boolean value exp with true if it is an experimental value from
      // the lookup table (and therefore has units of micro-amu) or false if it is
      // a theoretical value computed using the liquid drop model (and has units
      // of MeV).  This method of looking up masses is used to avoid unnecessary
      // unit conversions that can lead to losses of precision.
      static double lookup_atomic_mass(int nucleus_pid, bool& exp,
        bool theory_ok = true);
      static double lookup_atomic_mass(int Z, int A, bool& exp,
        bool theory_ok = true);
  
      // Lookup table for particle masses. Keys are PDG particle
      // ID numbers, values are masses in micro-amu.
      static const std::unordered_map<int, double> particle_masses;
  
      // Lookup table for atomic masses. Keys are PDG particle
      // ID numbers for the nuclei, values are masses in micro-amu.
      static const std::unordered_map<int, double> atomic_masses;
  
      // Particle IDs for all of the nuclear fragments that will
      // be considered when calculating separation energies
      static const std::vector<int> fragment_pids;
  
      // Liquid drop model parameters (taken from A. J. Koning, et al., Nucl.
      // Phys. A 810 (2008) pp. 13-76)
      static constexpr double Mn = 8.07144; // MeV
      static constexpr double MH = 7.28899; // MeV
      static constexpr double a1 = 15.677; // MeV
      static constexpr double a2 = 18.56; // MeV
      static constexpr double kappa = 1.79;
      static constexpr double c3 = 0.717; // MeV
      static constexpr double c4 = 1.21129; // MeV
  };

}
