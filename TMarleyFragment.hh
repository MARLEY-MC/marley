#pragma once

#include "TMarleyParity.hh"

// Stores data about a nuclear fragment
class TMarleyFragment {
  public:
    inline TMarleyFragment(int p_id, int twoS, TMarleyParity pi) {
      pid = p_id;
      two_s = twoS;
      parity = pi;
    }

    inline int get_pid() {
      return pid;
    }

    inline int get_two_s() {
      return two_s;
    }

    inline TMarleyParity get_parity() {
      return parity;
    }

  private:
    int pid; // Particle ID for this fragment
    int two_s; // Two times the spin of the fragment (allows for
               // half-integer spins)
    TMarleyParity parity; // Intrinsic parity of the fragment
};
