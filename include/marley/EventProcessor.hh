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

namespace HepMC3 {
  class GenEvent;
}

namespace marley {

  // Forward-declare some needed classes
  class Generator;

  /// @brief Abstract base class for entities that take a HepMC3::GenEvent as
  /// input and possibly modify it
  class EventProcessor {

    public:

      inline EventProcessor() {}

      inline virtual ~EventProcessor() = default;

      /// @brief Processes an input GenEvent object
      /// @param[in,out] ev The GenEvent object to be processed
      /// @param[in] gen If needed, the Generator object to use during
      /// processing (e.g., for obtaining random numbers)
      virtual void process_event( HepMC3::GenEvent& ev,
        marley::Generator& gen ) = 0;
  };

}
