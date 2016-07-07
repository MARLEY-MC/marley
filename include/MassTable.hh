#pragma once
#include <functional>
#include <unordered_map>
#include <vector>

#include "marley_utils.hh"
#include "Error.hh"

namespace marley {

  /// @brief Singleton lookup table for particle and atomic masses
  class MassTable {

    public:

      /// @brief Deleted copy constructor
      MassTable(const MassTable&) = delete;

      /// @brief Deleted move constructor
      MassTable(MassTable&&) = delete;

      /// @brief Deleted copy assignment operator
      MassTable& operator=(const MassTable&) = delete;

      /// @brief Deleted move assignment operator
      MassTable& operator=(MassTable&&) = delete;

      /// @brief Get a const reference to the singleton instance of the
      /// MassTable
      static const MassTable& Instance();

      /// @brief Get the mass of a particle
      /// @param pdg_code PDG code identifying the type of particle
      /// @return %Particle mass (MeV)
      double get_particle_mass(int pdg_code) const;

      /// @brief Get the mass of an atom
      /// @param pdg_code PDG code identifying the nucleus of the atom
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Atomic mass (MeV)
      double get_atomic_mass(int pdg_code, bool theory_ok = true) const;

      /// @brief Get the mass of an atom
      /// @param Z Atomic number of the atom
      /// @param A Mass number of the atom
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Atomic mass (MeV)
      double get_atomic_mass(int Z, int A, bool theory_ok = true) const;

      /// @brief Get the separation energy for emission of a nuclear fragment
      /// from a nucleus
      /// @param Z Atomic number for the mother nucleus
      /// @param A Mass number for the mother nucleus
      /// @param pdg PDG code for the emitted nuclear fragment
      /// @param theory_ok Whether to calculate a theoretical mass using
      /// the liquid drop model if an experimental mass cannot be found.
      /// If theory_ok is false, then a marley::Error will be thrown
      /// if an experimental mass cannot be found.
      /// @return Separation energy (MeV)
      double get_fragment_separation_energy(int Z, int A, int pdg,
        bool theory_ok = true) const;

      /// @brief Get the binding energy of a nucleus
      /// @param Z Atomic number
      /// @param A Mass number
      /// @param theory_ok Whether to use a theoretical (liquid drop model)
      /// atomic mass during the calculation if an experimental atomic mass
      /// cannot be found. If theory_ok is false, then a marley::Error will be
      /// thrown if an experimental atomic mass cannot be found.
      /// @return Binding energy (MeV)
      double get_binding_energy(int Z, int A, bool theory_ok = true) const;

      /// @brief Get the <a href="https://en.wikipedia.org/wiki/Mass_excess">
      /// mass excess</a> of a nucleus
      /// @param Z Atomic number
      /// @param A Mass number
      /// @param theory_ok Whether to use a theoretical (liquid drop model)
      /// atomic mass during the calculation if an experimental atomic mass
      /// cannot be found. If theory_ok is false, then a marley::Error will be
      /// thrown if an experimental atomic mass cannot be found.
      /// @return Mass excess (MeV)
      double get_mass_excess(int Z, int A, bool theory_ok = true) const;

      /// @brief Calculate a theoretical
      /// <a href="https://en.wikipedia.org/wiki/Mass_excess">mass excess</a>
      /// for a nucleus using the liquid drop model
      /// @param Z Atomic number
      /// @param A Mass number
      /// @return Mass excess (MeV)
      double liquid_drop_model_mass_excess(int Z, int A) const;

      /// @brief Calculate a theoretical atomic mass using the liquid drop model
      /// @param Z Atomic number
      /// @param A Mass number
      /// @return Atomic mass (MeV)
      double liquid_drop_model_atomic_mass(int Z, int A) const;

    protected:

      /// @brief Create the singleton MassTable object
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
      static constexpr double micro_amu_ = 0.000931494061; // MeV/uAMU
  };

}
