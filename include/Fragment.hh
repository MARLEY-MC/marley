#pragma once
#include "marley_utils.hh"
#include "MassTable.hh"
#include "Parity.hh"

namespace marley {

  /// @brief Simple container for storing reference data about
  /// each of the nuclear fragments considered by MARLEY's
  /// compound nucleus de-excitation models
  class Fragment {

    public:

      /// @param pid PDG particle ID
      /// @param twoS Two times the fragment's spin
      /// @param pi The fragment's parity
      inline Fragment(int pid, int twoS, marley::Parity pi);

      /// Get the PDG particle ID for this fragment
      inline int get_pid() const;

      /// Get two times the spin of this fragment
      inline int get_two_s() const;

      /// Get the parity of this fragment
      inline marley::Parity get_parity() const;

      /// Get the atomic number of this fragment
      inline int get_Z() const;

      /// Get the mass number of this fragment
      inline int get_A() const;

      /// Get the mass (in MeV) of this fragment
      inline double get_mass() const;

    private:
      int pid_; ///< PDG particle ID

      /// Two times the fragment spin (allows for half-integer spins)
      int two_s_;

      marley::Parity parity_; ///< intrinsic parity
      int Z_; ///< atomic number
      int A_; ///< mass number
  };

// Inline function definitions
inline Fragment::Fragment(int pid, int twoS, marley::Parity pi)
  : pid_(pid), two_s_(twoS), parity_(pi)
{
  Z_ = marley_utils::get_particle_Z(pid);
  A_ = marley_utils::get_particle_A(pid);
}

inline int Fragment::get_pid() const
  { return pid_; }

inline int Fragment::get_two_s() const
  { return two_s_; }

inline marley::Parity Fragment::get_parity() const
  { return parity_; }

inline int Fragment::get_Z() const
  { return Z_; }

inline int Fragment::get_A() const
  { return A_; }

inline double Fragment::get_mass() const
  { return marley::MassTable::Instance().get_particle_mass(pid_); }

}
