/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#pragma once

// Standard library includes
#include <memory>
#include <random>
#include <string>
#include <vector>

// MARLEY includes
#include "marley/Reaction.hh"
#include "marley/WeightCalculator.hh"

namespace HepMC3 {
  class GenEvent;
}

namespace marley {

  class DiscreteNuclearReaction;
  class Generator;
  class JSON;

  /// @brief WeightCalculator that varies nuclear matrix element strengths
  /// according to their experimental uncertainties. Three modes are
  /// available: "multisim" (dimidiated Gaussian random draws),
  /// "sigma_shift" (deterministic systematic shifts of all matrix
  /// elements), and "unisim" (deterministic systematic shifts of a
  /// single matrix element at a time).
  class StrengthVariationWeightCalculator : public WeightCalculator {

    public:

      /// @brief Supported variation modes
      enum class VariationMode { multisim, sigma_shift, unisim };

      /// @brief Static factory: validates JSON configuration and creates
      /// all variation instances. Called by the Weighter for each
      /// strength_variation entry in the weights config array.
      /// @param config JSON configuration for this calculator
      /// @param gen Generator whose loaded reactions are used by the
      /// unisim mode to determine the number of matrix elements
      /// @return Vector of shared pointers to the created instances
      static std::vector< std::shared_ptr<
        StrengthVariationWeightCalculator > >
        create_instances( const marley::JSON& config,
          marley::Generator& gen );

      virtual ~StrengthVariationWeightCalculator() = default;

      virtual double weight( HepMC3::GenEvent& event,
        marley::Generator& gen ) const override;

    private:

      /// @brief Constructor for multisim mode
      StrengthVariationWeightCalculator( const std::string& name,
        std::shared_ptr< std::mt19937_64 > rng,
        const std::string& resolved_reaction_file );

      /// @brief Constructor for sigma_shift mode
      StrengthVariationWeightCalculator( const std::string& name,
        double sigma_factor,
        const std::string& resolved_reaction_file );

      /// @brief Constructor for unisim mode
      StrengthVariationWeightCalculator( const std::string& name,
        double sigma_factor,
        size_t matrix_element_index,
        const std::string& resolved_reaction_file );

    protected:

      /// @brief Lazy initialization of per-instance varied strengths
      void ensure_initialized( marley::Generator& gen ) const;

      /// @brief Shared RNG (seeded once, shared across N instances)
      std::shared_ptr< std::mt19937_64 > rng_;

      /// @brief Resolved path of the reaction input file of interest
      std::string resolved_reaction_file_;

      /// @brief Variation mode for this instance
      VariationMode mode_;

      /// @brief Signed sigma factor for systematic shifts (sigma_shift
      /// and unisim modes)
      double sigma_factor_ = 0.;

      /// @brief Index of the single matrix element to vary (unisim mode)
      size_t me_idx_ = 0;

      /// @brief Whether lazy initialization has been completed
      mutable bool initialized_ = false;

      /// @brief Pointer to the matched DiscreteNuclearReaction
      mutable const DiscreteNuclearReaction* dnr_ = nullptr;

      /// @brief Process type of the matched Reaction
      mutable Reaction::ProcessType process_type_
        = Reaction::ProcessType::Unknown;

      /// @brief Target nucleus PDG code for the matched Reaction
      mutable int target_pdg_ = 0;

      /// @brief Pre-generated varied strengths (one per matrix element
      /// in the matched Reaction)
      mutable std::vector< double > varied_;
  };

}
