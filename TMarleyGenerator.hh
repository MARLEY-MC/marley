#pragma once
#include <random>
#include <sstream>
#include <vector>

#include "TMarleyConfigFile.hh"
#include "TMarleyNeutrinoSource.hh"
#include "TMarleyParity.hh"
#include "TMarleyReaction.hh"
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
      double E_nu = nu_source.sample_neutrino_energy(*this);
      return reactions.front().create_event(E_nu, *this);
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

    // Assume equipartition of parity and use a Bernoulli distribution
    // with p = 0.5 to generate random parities
    inline TMarleyParity sample_parity() {
      return TMarleyParity(bernoulli_dist(rand_gen));
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

    // Bernoulli distribution used to generate random boolean values.
    // Default constructor (used here) generates true 50% of the time.
    std::bernoulli_distribution bernoulli_dist;

    TMarleyNeutrinoSource nu_source;
    TMarleyStructureDatabase structure_db;
    std::vector<TMarleyReaction> reactions;
};
