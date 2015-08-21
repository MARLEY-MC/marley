#pragma once
#include <random>
#include "TMarleyGamma.hh"
#include "TMarleyParity.hh"

/// A class that represents an ENSDF nuclear energy level

class TMarleyLevel {
  public:
    /// @param jpi a string containing the spin
    /// and parity of the level (e.g., 0+)
    TMarleyLevel(double E, int twoJ, TMarleyParity pi);
    TMarleyGamma* add_gamma(const TMarleyGamma& gamma);
    void clear_gammas();
    std::vector<TMarleyGamma>* get_gammas();
    bool get_gamma_status() const;
    double get_energy() const;
    int get_two_J() const;
    TMarleyParity get_parity() const;
    void set_energy(double E);
    void set_two_J(int twoJ);
    void set_parity(TMarleyParity pi);
    TMarleyGamma* sample_gamma();
    std::string get_spin_parity_string() const;

  private:
    double energy; // MeV
    int two_J; // Two times the level spin (allows us to represent half-integer
               // spins as integers)
    TMarleyParity parity;
  
    bool gammas_known; // Determining whether or not the gammas are known

    std::vector<TMarleyGamma> gammas;
    std::vector<double> gamma_intensities;
    std::discrete_distribution<int> gamma_dist;
};
