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

// MARLEY includes
#include "marley/JSONConfig.hh"

namespace marley {

  class RootJSONConfig : public JSONConfig {

    public:

      explicit RootJSONConfig(const marley::JSON& object)
        : JSONConfig(object) {}

      explicit RootJSONConfig(const std::string& json_filename)
        : JSONConfig(json_filename) {}

      virtual bool process_extra_source_types(const std::string& type,
        const marley::JSON& source_spec, int pdg_code,
        std::unique_ptr<marley::NeutrinoSource>& source) const override;
  };

}
