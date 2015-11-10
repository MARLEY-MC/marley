#pragma once
#include <random>
#include <sstream>
#include <vector>

#include "TMarleyConfigFile.hh"
#include "TMarleyNeutrinoSource.hh"
#include "TMarleyParity.hh"
#include "TMarleyNuclearReaction.hh"
#include "TMarleyStructureDatabase.hh"

class TMarleyGenerator {
  public:
    inline TMarleyGenerator(const TMarleyConfigFile& cf) {
      init(cf);
    }

    inline TMarleyGenerator(const std::string& filename) {
      TMarleyConfigFile cf(filename);
      init(cf);
    }

    inline TMarleyEvent create_event() {
      double E_nu;
      size_t r_index;
      sample_reaction(E_nu, r_index);
      return reactions.at(r_index).create_event(marley_utils::ELECTRON_NEUTRINO,
        E_nu, *this);
    }

    inline uint_fast64_t get_seed() const {
      return seed;
    }

    // TODO: add error handling here (state_string may be nullptr)
    // Use a state string to set the MARLEY RNG
    inline void seed_using_state_string(std::string* state_string) {
      std::stringstream strstr(*state_string);
      strstr >> rand_gen;
    }

    // Get a string that represents the current internal state of the random
    // number generator
    inline std::string get_state_string() {
      std::stringstream ss;
      ss << rand_gen;
      return ss.str();
    }

    // Sample a random number uniformly on either [min, max) or [min, max]
    double uniform_random_double(double min, double max, bool inclusive);

    // Sample a random number x from the pdf f(x) on the interval [xmin, xmax]
    double rejection_sample(std::function<double(double)> f, double xmin,
      double xmax, double max_search_tolerance
      = DEFAULT_REJECTION_SAMPLING_TOLERANCE);

    // Sample from a discrete distribution object
    template <typename numType> inline numType
      discrete_sample(std::discrete_distribution<numType>& disc_dist)
    {
      return disc_dist(rand_gen);
    }

    // Sample from a discrete distribution object using the parameters params
    template <typename numType> inline numType
      discrete_sample(std::discrete_distribution<numType>& disc_dist,
      const typename std::discrete_distribution<numType>::param_type& params)
    {
      return disc_dist(rand_gen, params);
    }

    inline TMarleyStructureDatabase& get_structure_db() {
      return structure_db;
    }

    inline const std::vector<TMarleyNuclearReaction>& get_reactions() const {
      return reactions;
    }

    // Loads Ea with the energy of a reacting neutrino, and loads
    // r_index with the index of the reaction it undergoes. Input
    // values for Ea and r_index are ignored.
    void sample_reaction(double& Ea, size_t& r_index);

    // Returns the PDF that reacting neutrino energies obey. It is normalized
    // between nu_source.E_min and nu_source.E_max.
    inline double normalized_Ea_pdf(double Ea) {
      return Ea_pdf(Ea);
    }

    inline TMarleyNeutrinoSource& get_nu_source() {
      return nu_source;
    }

  private:
    // Seed for the random number generator
    uint_fast64_t seed;
    // 64-bit Mersenne Twister random number generator
    std::mt19937_64 rand_gen;
    // Default tolerance for rejection sampling
    static constexpr double DEFAULT_REJECTION_SAMPLING_TOLERANCE = 1e-8;
    // Initialization code shared by multiple constructors
    void init(const TMarleyConfigFile& cf);

    TMarleyNeutrinoSource nu_source;
    TMarleyStructureDatabase structure_db;
    std::vector<TMarleyNuclearReaction> reactions;

    // Use total cross sections for each reaction as weights for sampling a
    // reaction type
    std::vector<double> total_xs_values;
    // Discrete distribution for sampling reaction types
    std::discrete_distribution<size_t> r_index_dist;
    // PDF used for sampling reacting neutrino energies
    std::function<double(double)> Ea_pdf;

    // Helper function for sampling reacting neutrino energies
    double unnormalized_Ea_pdf(double Ea);
};
