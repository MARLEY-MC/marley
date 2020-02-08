/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/COPYING or
// visit http://opensource.org/licenses/GPL-3.0

#pragma once

namespace marley {

  // Forward-declare some needed classes
  class Event;
  class Generator;

  /// @brief Abstract base class for entities that take a marley::Event as
  /// input and possibly modify it
  class EventProcessor {

    public:

      inline EventProcessor() {}

      inline virtual ~EventProcessor() = default;

      /// @brief Processes an input MARLEY Event object
      /// @param[in,out] ev The Event object to be processed
      /// @param[in] gen If needed, the Generator object to use during
      /// processing (e.g., for obtaining random numbers)
      virtual void process_event( marley::Event& ev,
        marley::Generator& gen ) = 0;
  };

}
