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

  /// @brief WeightCalculator that propagates experimental uncertainties
  /// on allowed nuclear matrix elements using a dimidiated (bifurcated)
  /// Gaussian probability density function
  class StrengthVariationWeightCalculator : public WeightCalculator {

    public:

      /// @param config JSON configuration for this calculator
      /// @param instance_index Zero-based index of this variation instance
      /// @param rng Shared random number generator (seeded externally,
      /// shared by all instances from the same config entry)
      /// @param resolved_reaction_file Resolved path to the .react file
      /// for the reaction whose matrix elements should be varied
      StrengthVariationWeightCalculator( const marley::JSON& config,
        long instance_index,
        std::shared_ptr< std::mt19937_64 > rng,
        const std::string& resolved_reaction_file );

      virtual ~StrengthVariationWeightCalculator() = default;

      virtual double weight( HepMC3::GenEvent& event,
        marley::Generator& gen ) const override;

    protected:

      /// @brief Lazy initialization of per-instance varied strengths
      void ensure_initialized( marley::Generator& gen ) const;

      /// @brief Shared RNG (seeded once, shared across N instances)
      std::shared_ptr< std::mt19937_64 > rng_;

      /// @brief Resolved path of the reaction input file of interest
      std::string resolved_reaction_file_;

      /// @brief Standard normal distribution for dimidiated Gaussian
      mutable std::normal_distribution< double > normal_dist_;

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
