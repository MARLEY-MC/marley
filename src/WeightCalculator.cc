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

// MARLEY includes
#include "marley/Error.hh"
#include "marley/JSON.hh"
#include "marley/WeightCalculator.hh"

marley::WeightCalculator::WeightCalculator( const marley::JSON& config )
{
  if ( !config.is_object() ) {
    throw marley::Error( "Non-object JSON configuration passed to constructor"
      " of marley::WeightCalculator" );
  }

  if ( !config.has_key("name") ) {
    throw marley::Error( "Missing \"name\" key in weight calculator JSON"
      " configuration" );
  }

  name_ = config.at( "name" ).to_string();
}

marley::WeightCalculator::WeightCalculator( const std::string& name )
  : name_( name ) {}
