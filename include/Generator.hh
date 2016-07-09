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

  /// @brief The MARLEY Event generator
  class Generator {

    public:

      /// @brief Create a Generator using a ConfigFile object
      Generator(marley::ConfigFile& cf);

      /// @brief Create a Generator using a configuration file
      /// @param filename Name of the configuration file
      Generator(const std::string& filename);

      /// @brief Create an Event using the NeutrinoSource, Reaction, and
      /// StructureDatabase objects owned by this Generator
      marley::Event create_event();

      /// @brief Get the seed used to initialize this Generator
      inline uint_fast64_t get_seed() const;

      /// @brief Use a string to set this Generator's internal state
      /// @details This function is typically used to restore a Generator
      /// to a state saved using get_state_string().
      void seed_using_state_string(std::string& state_string);

      /// @brief Get a string that represents the current internal state of
      /// this Generator
      std::string get_state_string() const;

      /// @brief Sample a random number uniformly on either [min, max) or
      /// [min, max]
      /// @param min Lower bound of the sampling interval
      /// @param max Upper bound of the sampling interval
      /// @param inclusive Whether the upper bound should be included (true)
      /// or excluded (false) from the possible sampling outcomes
      double uniform_random_double(double min, double max, bool inclusive);

      /// @brief Sample from a given 1D probability density function f(x) on
      /// the interval [xmin, xmax] using a simple rejection method
      /// @param f Probability density function to use for sampling
      /// @param xmin Lower bound of the sampling interval
      /// @param xmax Upper bound of the sampling interval
      /// @param max_search_tolerance Tolerance to use when finding the
      /// maximum of f(x) using <a href="http://tinyurl.com/ntqkfck">Brent's
      /// method</a>
      /// @return Sampled value of x
      double rejection_sample(std::function<double(double)> f, double xmin,
        double xmax, double max_search_tolerance
        = DEFAULT_REJECTION_SAMPLING_TOLERANCE_);

      /// @brief Get a reference to the StructureDatabase owned by this
      /// Generator
      inline marley::StructureDatabase& get_structure_db();

      /// @brief Get a const reference to the vector of Reaction objects
      /// owned by this Generator
      inline const std::vector<std::unique_ptr<marley::Reaction> >&
        get_reactions() const;

      /// @brief Sample a Reaction and an energy for the reacting neutrino
      /// @param[out] E Total energy of the neutrino undergoing the reaction
      /// @return Reference to the sampled Reaction owned by this Generator
      marley::Reaction& sample_reaction(double& E);

      /// @brief Probability density function that describes the distribution
      /// of reacting neutrino energies
      /// @details This function computes the cross-section weighted neutrino
      /// flux (normalized to unity between source_->E_min and
      /// source_->E_max) including cross-section contributions from all
      /// Reactions owned by this Generator. For the distribution of
      /// <i>incident</i> neutrino energies, use marley::NeutrinoSource::pdf()
      /// @param E Total energy of the reacting neutrino
      /// @return Probability density (MeV<sup> -1</sup>)
      double E_pdf(double E);

      /// @brief Get a reference to the NeutrinoSource owned by this Generator
      /// @details Throws a marley::Error if this Generator does not own a
      /// NeutrinoSource object.
      /// @note The sampling interval [Emin, Emax] for the incident neutrino
      /// energy distribution described by marley::NeutrinoSource::pdf() should
      /// not be changed using the reference returned by this function. Doing
      /// so could silently invalidate the normalization of E_pdf(). This is
      /// currently impossible to do using the public interface of
      /// marley::NeutrinoSource, but future changes to the NeutrinoSource or
      /// Generator classes should continue to protect against this
      /// possibility.
      marley::NeutrinoSource& get_source();

      /// @brief Sample from a std::discrete_distribution
      template <typename numType> inline numType
        discrete_sample(std::discrete_distribution<numType>& disc_dist);

      /// @brief Sample from a std::discrete_distribution using the parameters
      /// params
      template <typename numType> inline numType
        discrete_sample(std::discrete_distribution<numType>& disc_dist,
        const typename std::discrete_distribution<numType>::param_type& params);

    private:

      /// @brief Helper function that contains initialization code shared by
      /// multiple constructors
      void init(marley::ConfigFile& cf);

      /// @brief Helper function that updates the normalization factor to
      /// use in E_pdf()
      void normalize_E_pdf();

      /// @brief Seed for the random number generator
      uint_fast64_t seed_;

      /// @brief 64-bit Mersenne Twister random number generator
      std::mt19937_64 rand_gen_;

      /// @brief Default stopping tolerance for rejection sampling
      static constexpr double DEFAULT_REJECTION_SAMPLING_TOLERANCE_ = 1e-8;

      /// @brief Normalizaton factor for E_pdf()
      double norm_ = 1.;

      /// @brief NeutrinoSource used to sample reacting neutrino energies
      std::unique_ptr<marley::NeutrinoSource> source_;

      /// @brief StructureDatabase used to simulate nuclear de-excitations
      /// when creating Event objects
      marley::StructureDatabase structure_db_;

      /// @brief Reaction(s) used to sample reacting neutrino energies
      std::vector<std::unique_ptr<marley::Reaction> > reactions_;

      /// @brief Total cross section values to use as weights for sampling a
      /// Reaction in sample_reaction()
      /// @details These total cross section values are energy dependent. They
      /// are therefore updated with every call to E_pdf().
      std::vector<double> total_xs_values_;

      /// @brief Discrete distribution used for Reaction sampling
      std::discrete_distribution<size_t> r_index_dist_;
  };

  // Inline function definitions
  inline uint_fast64_t Generator::get_seed() const { return seed_; }

  inline marley::StructureDatabase& Generator::get_structure_db()
    { return structure_db_; }

  inline const std::vector<std::unique_ptr<marley::Reaction> >&
    Generator::get_reactions() const { return reactions_; }

  template <typename numType> inline numType
    Generator::discrete_sample(std::discrete_distribution<numType>& disc_dist)
    { return disc_dist(rand_gen_); }

  template <typename numType> inline numType
    Generator::discrete_sample(std::discrete_distribution<numType>& disc_dist,
    const typename std::discrete_distribution<numType>::param_type& params)
    { return disc_dist(rand_gen_, params); }
}
