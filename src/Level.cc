#include <string>
#include <vector>
#include <regex>

#include "marley_utils.hh"
#include "Generator.hh"
#include "Level.hh"
#include "Parity.hh"

marley::Level::Level(double E, int twoJ, marley::Parity pi) {
  energy = E;
  two_J = twoJ;
  parity = pi;

  gammas_known = false; // The default is that gammas are not known
}

/// Choose a gamma owned by this level randomly based on the relative
/// intensities of all of the gammas.  Return a pointer to the gamma that was
/// chosen.  If this level doesn't have any gammas, return a null pointer.
const marley::Gamma* marley::Level::sample_gamma(marley::Generator& gen) {
  if (gammas.empty()) {
    return nullptr;
  }
  else {
    // Get the index of the gamma to return by randomly sampling from the
    // discrete distribution gamma_dist using the standard marley_utils random
    // number generator.
    int g_index = gen.discrete_sample(gamma_dist);
    // Return a pointer to the corresponding gamma
    return &(gammas[g_index]);
  }
}

// Add a new gamma object to the level and return a pointer to it
marley::Gamma* marley::Level::add_gamma(const marley::Gamma& gamma) {
  // Update the vector of gamma objects
  gammas.push_back(gamma);

  // Update the distribution for sampling gammas
  update_gamma_distribution();

  // Update the gamma status to known
  gammas_known = true;

  // Return a pointer to the newly added gamma
  return &gammas.back();
}

double marley::Level::get_energy() const {
  return energy;
}

void marley::Level::set_energy(double E) {
  energy = E;
}

int marley::Level::get_two_J() const {
  return two_J;
}

void marley::Level::set_two_J(int twoJ) {
  two_J = twoJ;
}

marley::Parity marley::Level::get_parity() const
{
  return parity;
}

void marley::Level::set_parity(marley::Parity pi) {
  parity = pi;
}

void marley::Level::clear_gammas() {
  gammas.clear();

  // The discrete distribution will be cleared by this command because the
  // vector of gammas is now empty.
  update_gamma_distribution();

  // Update the gamma status to unknown
  gammas_known = false;
}

std::vector<marley::Gamma>* marley::Level::get_gammas() {
  return &gammas;
}

bool marley::Level::get_gamma_status() const
{
  return gammas_known;
}

std::string marley::Level::get_spin_parity_string() const {
  std::string str = std::to_string(two_J / 2);
  // If 2*J is odd, then the level has half-integer spin
  if (two_J % 2) str += "/2";
  return str + parity.str();
}
