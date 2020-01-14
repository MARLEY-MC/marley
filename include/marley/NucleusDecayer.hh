/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

// MARLEY includes
#include "marley/EventProcessor.hh"

namespace marley {

  /// @brief EventProcessor that handles nuclear de-excitations
  class NucleusDecayer : public EventProcessor {

    public:

      inline NucleusDecayer() {}

      inline virtual ~NucleusDecayer() = default;

      virtual void process_event( marley::Event& event,
        marley::Generator& gen ) override;
  };

}
