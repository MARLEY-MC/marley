/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see \${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

// MARLEY includes
#include "marley/Event.hh"
#include "marley/R5EFR.hh"
#include "marley/RootEventFileReader.hh"

marley::R5EFR::R5EFR( const std::string& file_name ) {
  event_file_reader_ = new marley::RootEventFileReader( file_name );
}

marley::R5EFR::~R5EFR() {
  auto* efr = static_cast<marley::RootEventFileReader*>( event_file_reader_ );
  if ( efr ) delete efr;
}

bool marley::R5EFR::next_event(marley::Event& ev) {
  auto* efr = static_cast<marley::RootEventFileReader*>( event_file_reader_ );
  return efr->next_event( ev );
}

marley::R5EFR::operator bool() const {
  auto* efr = static_cast<marley::RootEventFileReader*>( event_file_reader_ );
  return static_cast<bool>( *efr );
}
