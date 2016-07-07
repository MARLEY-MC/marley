#pragma once
#include <functional>
#include <unordered_map>
#include <vector>

#include "marley_utils.hh"
#include "Error.hh"

namespace marley {

  /// @brief Lookup table for particle and atomic masses
  class MassTable {

    public:

      // delete copy and move constructors and assign operators
      // Copy constructor
      MassTable(const MassTable&) = delete;
      // Move constructor
      MassTable(MassTable&&) = delete;
      // Copy assignment operator
      MassTable& operator=(const MassTable&) = delete;
      // Move assignment operator
      MassTable& operator=(MassTable&&) = delete;

      static const MassTable& Instance();

      double get_particle_mass(int particle_id) const;
      double get_atomic_mass(int nucleus_pid, bool theory_ok = true) const;
      double get_atomic_mass(int Z, int A, bool theory_ok = true) const;
      double get_fragment_separation_energy(int Z, int A, int pid,
        bool theory_ok = true) const;
      double get_binding_energy(int Z, int A, bool theory_ok = true) const;
      double get_mass_excess(int Z, int A, bool theory_ok = true) const;
      double liquid_drop_model_mass_excess(int Z, int A) const;
      double liquid_drop_model_atomic_mass(int Z, int A) const;

    protected:

      MassTable();

    private:

      // Function used internally by the mass table. Returns an atomic mass and
      // loads the boolean value exp with true if it is an experimental value
      // from the lookup table (and therefore has units of micro-amu) or false
      // if it is a theoretical value computed using the liquid drop model (and
      // has units of MeV).  This method of looking up masses is used to avoid
      // unnecessary unit conversions that can lead to losses of precision.
      double lookup_atomic_mass(int nucleus_pid, bool& exp,
        bool theory_ok = true) const;
      double lookup_atomic_mass(int Z, int A, bool& exp,
        bool theory_ok = true) const;

      // Lookup table for particle masses. Keys are PDG particle
      // ID numbers, values are masses in micro-amu.
      const std::unordered_map<int, double> particle_masses_;

      // Lookup table for atomic masses. Keys are PDG particle
      // ID numbers for the nuclei, values are masses in micro-amu.
      const std::unordered_map<int, double> atomic_masses_;

      // Factor to use when converting from micro-amu to MeV
      static constexpr double micro_amu_ = 0.000931494061;
  };

}
