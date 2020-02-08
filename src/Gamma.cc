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

#include "marley/Gamma.hh"

marley::Gamma::Gamma(double energy, double rel_intensity,
  marley::Level* start_lev, marley::Level* end_lev) : energy_(energy),
  relative_intensity_(rel_intensity), start_level_(start_lev),
  end_level_(end_lev)
{
}
