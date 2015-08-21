#include <string>
#include <vector>
#include <regex>
#include "marley_utils.hh"
#include "TMarleyLevel.hh"
#include "TMarleyParity.hh"

TMarleyLevel::TMarleyLevel(double E, int twoJ, TMarleyParity pi) {
  energy = E;
  two_J = twoJ;
  parity = pi;

  gammas_known = false; // The default is that gammas are not known
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

// Add a new gamma object to the level and return a pointer to it
TMarleyGamma* TMarleyLevel::add_gamma(const TMarleyGamma& gamma) {
  // Update the vector of gamma objects
  gammas.push_back(gamma);

  // Update the vector of gamma intensities
  gamma_intensities.push_back(gamma.get_ri());

  // Update the discrete distribution used to sample gammas
  // when simulating a gamma cascade
  std::discrete_distribution<int>::param_type
    params(gamma_intensities.begin(), gamma_intensities.end());
  gamma_dist.param(params);

  // Update the gamma status to known
  gammas_known = true;

  // Return a pointer to the newly added gamma
  return &gammas.back();
}

double TMarleyLevel::get_energy() const {
  return energy; 
}

void TMarleyLevel::set_energy(double E) {
  energy = E;
}

int TMarleyLevel::get_two_J() const {
  return two_J;
}

void TMarleyLevel::set_two_J(int twoJ) {
  two_J = twoJ;
}

TMarleyParity TMarleyLevel::get_parity() const
{
  return parity;
}

void TMarleyLevel::set_parity(TMarleyParity pi) {
  parity = pi;
}

void TMarleyLevel::clear_gammas() {
  gammas.clear();
  gamma_intensities.clear();
  // The discrete distribution will be cleared by these
  // commands because gamma_intensities is now empty
  std::discrete_distribution<int>::param_type
    params(gamma_intensities.begin(), gamma_intensities.end());
  gamma_dist.param(params);
  // Update the gamma status to unknown
  gammas_known = false;
}

std::vector<TMarleyGamma>* TMarleyLevel::get_gammas() {
  return &gammas;
}

bool TMarleyLevel::get_gamma_status() const
{
  return gammas_known;
}

std::string TMarleyLevel::get_spin_parity_string() const {
  std::string str = std::to_string(two_J / 2);
  // If 2*J is odd, then the level has half-integer spin
  if (two_J % 2) str += "/2";
  return str + parity.str();
}
