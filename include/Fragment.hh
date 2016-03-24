#pragma once

#include "MassTable.hh"
#include "Parity.hh"

namespace marley {

  // Stores data about a nuclear fragment
  class Fragment {
    public:
      inline Fragment(int p_id, int twoS, marley::Parity pi) {
        pid = p_id;
        two_s = twoS;
        parity = pi;
        Z = marley::MassTable::get_particle_Z(p_id);
        A = marley::MassTable::get_particle_A(p_id);
      }
  
      inline int get_pid() const {
        return pid;
      }
  
      inline int get_two_s() const {
        return two_s;
      }
  
      inline marley::Parity get_parity() const {
        return parity;
      }
  
      inline int get_Z() const {
        return Z;
      }
  
      inline int get_A() const {
        return A;
      }
  
      inline double get_mass() const {
        return marley::MassTable::get_particle_mass(pid);
      }
  
    private:
      int pid; // Particle ID for this fragment
      int two_s; // Two times the spin of the fragment (allows for
                 // half-integer spins)
      marley::Parity parity; // Intrinsic parity of the fragment
      int Z, A; // Atomic and mass numbers (extracted from pid upon construction)
  };

}
