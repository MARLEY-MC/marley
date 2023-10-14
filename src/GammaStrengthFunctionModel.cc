/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
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

#include "marley/Error.hh"
#include "marley/GammaStrengthFunctionModel.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/Parity.hh"

using TrType = marley::GammaStrengthFunctionModel::TransitionType;

marley::GammaStrengthFunctionModel::GammaStrengthFunctionModel(int Z, int A)
  : Z_(Z), A_(A) {}

void marley::GammaStrengthFunctionModel::check_multipolarity(int l) {
  /// @todo Improve error message
  if (l < 1) throw marley::Error( "Invalid multipolarity "
    + std::to_string(l) + " given for gamma ray strength"
    " function calculation" );
}
