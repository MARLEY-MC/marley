#pragma once
#include <random>
#include "TMarleyGamma.hh"

/// A class that represents an ENSDF nuclear energy level

class TMarleyLevel {
  public:
    /// @param energy a string containing the
    /// level energy in keV
    /// @param jpi a string containing the spin
    /// and parity of the level (e.g., 0+)
    TMarleyLevel(std::string energy = "0", std::string jpi = "");
    void add_gamma(const TMarleyGamma& gamma);
    void clear_gammas();
    std::vector<TMarleyGamma>* get_gammas();
    double get_numerical_energy() const;
    std::string get_string_energy() const;
    std::string get_spin_parity() const;
    void set_energy(std::string energy);
    void set_spin_parity(std::string jpi);
    TMarleyGamma* sample_gamma();

  private:
    std::string sEnergy;
    double fEnergy;
    std::string spin_parity;
    std::vector<TMarleyGamma> gammas;
    std::vector<double> gamma_intensities;
    std::discrete_distribution<int> gamma_dist;
};
