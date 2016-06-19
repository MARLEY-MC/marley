#pragma once
#include <random>

#include "Gamma.hh"
#include "Parity.hh"

namespace marley {

  class Generator;

  /// @brief A discrete nuclear energy level
  class Level {
    public:

      /// @param E excitation energy of this level (MeV)
      /// @param twoJ two times the level spin
      /// @param pi level parity
      Level(double E, int twoJ, marley::Parity pi);

      /// @brief Retrieve a const reference to the vector of gamma rays owned
      /// by this level
      inline const std::vector<marley::Gamma>& get_gammas() const;

      /// Retrieve a reference to the vector of gamma rays owned by this level
      inline std::vector<marley::Gamma>& get_gammas();

      /// Get the excitation energy of this level (MeV)
      inline double get_energy() const;
      /// Set the excitation energy of this level (MeV)
      inline void set_energy(double E);

      /// Get two times the level spin
      inline int get_twoJ() const;
      /// Set two times the level spin
      inline void set_twoJ(int twoJ);

      /// Get the level parity
      inline marley::Parity get_parity() const;
      /// Set the level parity
      inline void set_parity(marley::Parity pi);

      /// @return true if the level owns at least one Gamma object
      inline bool has_gammas() const;

      /// @brief Add a new gamma-ray transition to this level
      /// @return a reference to the newly added gamma
      marley::Gamma& add_gamma(const marley::Gamma& gamma);

      /// Remove all gamma ray information from this level
      void clear_gammas();

      /// @brief Choose a gamma owned by this level randomly based on the
      /// relative intensities of all of the gammas.
      /// @return a pointer to the selected Gamma object, or nullptr
      /// if the level doesn't have any gammas
      const marley::Gamma* sample_gamma(marley::Generator& gen);

      /// Returns the level spin-parity as a string
      std::string get_spin_parity_string() const;

    private:

      double energy_; ///< excitation energy (MeV)

      /// @brief two times the level spin
      /// @note the factor of two allows us to represent half-integer spins as
      /// integers
      int twoJ_;

      marley::Parity parity_; ///< level parity

      /// @brief gamma-ray transitions owned by this level
      std::vector<marley::Gamma> gammas_;

      /// @brief discrete distribution object used to sample gamma-ray
      /// de-excitations
      std::discrete_distribution<size_t> gamma_dist_;

      /// @brief helper function that updates gamma-ray distribution when Gamma
      /// objects are added or removed from the level
      void update_gamma_distribution();
  };

  // Inline function definitions
  inline double Level::get_energy() const { return energy_; }
  inline void Level::set_energy(double E) { energy_ = E; }

  inline int Level::get_twoJ() const { return twoJ_; }
  inline void Level::set_twoJ(int twoJ) { twoJ_ = twoJ; }

  inline marley::Parity Level::get_parity() const { return parity_; }
  inline void Level::set_parity(marley::Parity pi) { parity_ = pi; }

  inline const std::vector<marley::Gamma>& Level::get_gammas() const
    { return gammas_; }
  inline std::vector<marley::Gamma>& Level::get_gammas() { return gammas_; }

  inline bool Level::has_gammas() const { return gammas_.empty(); }
}
