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
#include <string>

namespace HepMC3 {
  class GenEvent;
}

namespace marley {

  class Generator;
  class JSON;

  /// @brief Abstract base class for objects that compute numerical weights for
  /// an input event
  class WeightCalculator {

    public:

      WeightCalculator( const marley::JSON& config );
      WeightCalculator( const std::string& name );
      virtual ~WeightCalculator() = default;

      virtual double weight( HepMC3::GenEvent& event,
        marley::Generator& gen ) const = 0;

      const std::string& name() const { return name_; }

    protected:

      /// @brief Name used to label the output event weight
      std::string name_;
  };

  /// @brief WeightCalculator object that unconditionally returns unit weights
  class TrivialWeightCalculator : public WeightCalculator {

    public:

      TrivialWeightCalculator( const marley::JSON& config )
        : WeightCalculator( config ) {}
      virtual ~TrivialWeightCalculator() = default;

      inline virtual double weight( HepMC3::GenEvent& /*event*/,
        marley::Generator& /*gen*/ ) const override { return 1.; }
  };

}
