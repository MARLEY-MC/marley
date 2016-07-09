#pragma once
#include <memory>
#include <random>
#include <sstream>
#include <vector>

#include "ConfigFile.hh"
#include "NeutrinoSource.hh"
#include "Parity.hh"
#include "NuclearReaction.hh"
#include "StructureDatabase.hh"

namespace marley {

  class Generator {

    public:

      Generator(marley::ConfigFile& cf);

      Generator(const std::string& filename);

      marley::Event create_event();

      inline uint_fast64_t get_seed() const;

      // Use a state string to set the MARLEY RNG
      void seed_using_state_string(std::string& state_string);

      // Get a string that represents the current internal state of the random
      // number generator
      std::string get_state_string();

      // Sample a random number uniformly on either [min, max) or [min, max]
      double uniform_random_double(double min, double max, bool inclusive);

      // Sample a random number x from the pdf f(x) on the interval [xmin, xmax]
      double rejection_sample(std::function<double(double)> f, double xmin,
        double xmax, double max_search_tolerance
        = DEFAULT_REJECTION_SAMPLING_TOLERANCE);

      inline marley::StructureDatabase& get_structure_db();

      inline const std::vector<std::unique_ptr<marley::Reaction> >&
        get_reactions() const;

      // Loads Ea with the energy of a reacting neutrino, and loads
      // r_index with the index of the reaction it undergoes. Input
      // values for Ea and r_index are ignored.
      void sample_reaction(double& Ea, size_t& r_index);

      // Returns the PDF that reacting neutrino energies obey. It is normalized
      // between nu_source->E_min and nu_source->E_max.
      inline double normalized_Ea_pdf(double Ea);

      inline marley::NeutrinoSource* get_nu_source();

      // Sample from a discrete distribution object
      template <typename numType> inline numType
        discrete_sample(std::discrete_distribution<numType>& disc_dist);

      // Sample from a discrete distribution object using the parameters params
      template <typename numType> inline numType
        discrete_sample(std::discrete_distribution<numType>& disc_dist,
        const typename std::discrete_distribution<numType>::param_type& params);

    private:

      // Seed for the random number generator
      uint_fast64_t seed;
      // 64-bit Mersenne Twister random number generator
      std::mt19937_64 rand_gen;
      // Default tolerance for rejection sampling
      static constexpr double DEFAULT_REJECTION_SAMPLING_TOLERANCE = 1e-8;
      // Initialization code shared by multiple constructors
      void init(marley::ConfigFile& cf);

      std::unique_ptr<marley::NeutrinoSource> nu_source;
      marley::StructureDatabase structure_db;
      std::vector<std::unique_ptr<marley::Reaction> > reactions;

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

  // Inline function definitions
  inline uint_fast64_t Generator::get_seed() const { return seed; }

  inline marley::StructureDatabase& Generator::get_structure_db()
    { return structure_db; }

  inline const std::vector<std::unique_ptr<marley::Reaction> >&
    Generator::get_reactions() const { return reactions; }

  inline double Generator::normalized_Ea_pdf(double Ea) { return Ea_pdf(Ea); }

  inline marley::NeutrinoSource* Generator::get_nu_source()
    { return nu_source.get(); }

  template <typename numType> inline numType
    Generator::discrete_sample(std::discrete_distribution<numType>& disc_dist)
    { return disc_dist(rand_gen); }

  template <typename numType> inline numType
    Generator::discrete_sample(std::discrete_distribution<numType>& disc_dist,
    const typename std::discrete_distribution<numType>::param_type& params)
    { return disc_dist(rand_gen, params); }
}
