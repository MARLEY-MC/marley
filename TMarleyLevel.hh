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
    bool get_gamma_status() const;
    double get_numerical_energy() const;
    std::string get_string_energy() const;
    std::string get_spin_parity() const;
    int get_ispin() const;
    int get_iparity() const;
    void set_energy(std::string energy);
    void set_spin_parity(std::string jpi);
    TMarleyGamma* sample_gamma();

  private:
    std::string sEnergy;
    double fEnergy;
  
    std::string spin_parity;
    int ispin; // Spin as an integer, not isospin
    int iparity;

    bool gammas_known; // Determining whether or not the gammas are known
    std::vector<TMarleyGamma> gammas;
    std::vector<double> gamma_intensities;
    std::discrete_distribution<int> gamma_dist;
};
