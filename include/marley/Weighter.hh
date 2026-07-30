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
#include <string>
#include <vector>

// MARLEY includes
#include "marley/EventProcessor.hh"

namespace marley {

  // Forward-declare the WeightCalculator class
  class WeightCalculator;

  /// @brief EventProcessor that assigns event weights
  class Weighter : public EventProcessor {

    public:

      Weighter( const marley::JSON& config,
        marley::Generator& gen );

      inline virtual ~Weighter() = default;

      virtual void process_event( HepMC3::GenEvent& event,
        marley::Generator& gen ) override;

      std::vector< double > compute_weights( HepMC3::GenEvent& event,
        marley::Generator& gen ) const;

      /// @brief Returns the names of the weights associated with all
      /// configured weight calculators
      std::vector< std::string > get_weight_names() const;

      /// @brief Toggles inclusion of the central-value weight
      void set_use_cv_weight( bool use_it );

      /// @brief Provides non-const access to the owned vector
      /// of weight calculators
      inline std::vector< std::shared_ptr< WeightCalculator > >&
        get_weight_calculators() { return calc_vec_; }

    protected:

      /// @brief Reserved name for the central-value weights (required for
      /// compliance with the NuHepMC standard)
      static const std::string CV_WEIGHT_NAME;

      /// @brief Owned vector of weight calculators used to actually compute
      /// weights
      std::vector< std::shared_ptr< WeightCalculator > > calc_vec_;

      /// @brief Helper function to create the weight calculator that
      /// handles the central-value weights during event generation
      std::shared_ptr< marley::WeightCalculator > make_cv_weight_calc();
  };

}
