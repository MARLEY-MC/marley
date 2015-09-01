#pragma once

#include "TMarleyMassTable.hh"
#include "TMarleyParity.hh"

// Stores data about a nuclear fragment
class TMarleyFragment {
  public:
    inline TMarleyFragment(int p_id, int twoS, TMarleyParity pi) {
      pid = p_id;
      two_s = twoS;
      parity = pi;
      Z = TMarleyMassTable::get_particle_Z(p_id);
      A = TMarleyMassTable::get_particle_A(p_id);
    }

    inline int get_pid() const {
      return pid;
    }

    inline int get_two_s() const {
      return two_s;
    }

    inline TMarleyParity get_parity() const {
      return parity;
    }

    inline int get_Z() const {
      return Z;
    }

    inline int get_A() const {
      return A;
    }

    inline double get_mass() const {
      return TMarleyMassTable::get_particle_mass(pid);
    }

  private:
    int pid; // Particle ID for this fragment
    int two_s; // Two times the spin of the fragment (allows for
               // half-integer spins)
    TMarleyParity parity; // Intrinsic parity of the fragment
    int Z, A; // Atomic and mass numbers (extracted from pid upon construction)
};
