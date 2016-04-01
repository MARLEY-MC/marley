#pragma once
#include <unordered_map>
#include <functional>
#include <vector>
#include <utility>

#include "marley_utils.hh"
#include "Error.hh"

namespace marley {

  // TODO: Consider using a namespace (possibly with an anonymous namespace
  // in the source file to encapsulate the current private members)
  // rather than a class for the mass table. For more thoughts on this
  // possibility, see http://www.cplusplus.com/forum/general/121429/
  class MassTable {
  
    public:
      static constexpr double get_particle_mass(int particle_id) {
        int id = particle_id;
        // The lookup table only includes entries for particles (as opposed to
        // antiparticles), so flip the sign of the input particle id for the
        // lookup if it represents an antiparticle.
        if (id < 0) id *= -1;
        // Find the particle's mass in the lookup table, and convert its
        // value from micro-amu to MeV
        return micro_amu * find_particle_mass(id);
      }

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
      static constexpr std::pair<int, double> particle_masses[] = {
        {11, 548.57990946}, // e-
        {12, 0.0}, // nu_e
        {13, 113428.9267}, // mu-
        {14, 0.0}, // nu_mu
        {15, 1907490}, // tau-
        {16, 0.0}, // nu_tau
        {22, 0.0}, // photon
        {2112, 1008664.91585}, // neutron
        {2212, 1007276.466812}, // proton
        {1000010020, 2013553.212712}, // deuteron
        {1000010030, 3015500.7134}, // triton
        {1000020030, 3014932.2468}, // helion
        {1000020040, 4001506.179125}, // alpha
      };

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

      // Helper code for get_particle_mass(). This approach to accessing a
      // constexpr lookup table is based on the nifty trick described here:
      // http://stackoverflow.com/a/26079954/4081973
      static constexpr size_t particle_masses_size = sizeof(marley::MassTable::particle_masses)
        / sizeof(marley::MassTable::particle_masses[0]);
      
      static constexpr double find_particle_mass(int pid,
        size_t range = particle_masses_size)
      {
        return (range == 0) ? throw marley::Error(std::string("Unknown particle in marley::")
          + "MassTable::find_particle_mass()") :
          (particle_masses[range - 1].first == pid) ? particle_masses[range - 1].second :
          find_particle_mass(pid, range - 1);
      }
      
  };

}
