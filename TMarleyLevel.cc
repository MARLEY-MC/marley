#include <string>
#include <vector>
#include "marley_utils.hh"
#include "TMarleyLevel.hh"

TMarleyLevel::TMarleyLevel(std::string energy, std::string jpi) {
  sEnergy = energy;
  fEnergy = std::stod(energy); // Converts the energy string into a double
  spin_parity = jpi;
}

/// Choose a gamma owned by this level randomly based on the relative
/// intensities of all of the gammas.  Return a pointer to the gamma that was
/// chosen.  If this level doesn't have any gammas, return a null pointer.
TMarleyGamma* TMarleyLevel::sample_gamma() {
  if (gammas.empty()) {
    return nullptr;
  }
  else {
    // Get the index of the gamma to return by randomly sampling from the
    // discrete distribution gamma_dist using the standard marley_utils random
    // number generator.
    int g_index = gamma_dist(marley_utils::rand_gen);
    // Return a pointer to the corresponding gamma
    return &(gammas[g_index]);
  }
}

void TMarleyLevel::add_gamma(const TMarleyGamma& gamma) {
  // Update the vector of gamma objects
  gammas.push_back(gamma);

  // Update the vector of gamma intensities
  gamma_intensities.push_back(gamma.get_ri());

  // Update the discrete distribution used to sample gammas
  // when simulating a gamma cascade
  std::discrete_distribution<int>::param_type
    params(gamma_intensities.begin(), gamma_intensities.end());
  gamma_dist.param(params);

}

double TMarleyLevel::get_numerical_energy() const {
  return fEnergy; 
}

std::string TMarleyLevel::get_string_energy() const {
  return sEnergy; 
}

void TMarleyLevel::set_energy(std::string energy) {
  sEnergy = energy;
  fEnergy = std::stod(energy);
}

std::string TMarleyLevel::get_spin_parity() const {
  return spin_parity;
}

void TMarleyLevel::set_spin_parity(std::string jpi) {
  spin_parity = jpi;
}

void TMarleyLevel::clear_gammas() {
  gammas.clear();
  gamma_intensities.clear();
  // The discrete distribution will be cleared by these
  // commands because gamma_intensities is now empty
  std::discrete_distribution<int>::param_type
    params(gamma_intensities.begin(), gamma_intensities.end());
  gamma_dist.param(params);
}

std::vector<TMarleyGamma>* TMarleyLevel::get_gammas() {
  return &gammas;
}
