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

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
}

void for_each_event(
  const std::vector< std::string >& input_files,
  std::function< void( HepMC3::GenEvent&, bool,
    double, const std::shared_ptr< HepMC3::GenRunInfo >& ) > callback );
