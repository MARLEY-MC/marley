/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
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

#include "marley/WeightCalculator.hh"

namespace marley {

  /// @brief Computes event weights based on variations to the optical model
  /// potential
  class OMPWeightCalculator : public WeightCalculator {

    public:

      OMPWeightCalculator( const marley::JSON& config );
      virtual ~OMPWeightCalculator() = default;

      virtual double weight( HepMC3::GenEvent& event,
        marley::Generator& gen ) const override;

    protected:

      /// @brief Owned StructureDatabase with a custom set of optical model
      /// parameters
      std::shared_ptr< marley::StructureDatabase > sdb_;
  };

}
