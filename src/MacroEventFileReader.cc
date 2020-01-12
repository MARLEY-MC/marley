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
#include "marley/MacroEventFileReader.hh"
#include "marley/RootEventFileReader.hh"

namespace {

  marley::RootEventFileReader* get_refr_pointer(void* ptr) {
    return static_cast<marley::RootEventFileReader*>( ptr );
  }

}


marley::MacroEventFileReader::MacroEventFileReader( const std::string& file_name ) {
  event_file_reader_ = new marley::RootEventFileReader( file_name );
}

marley::MacroEventFileReader::~MacroEventFileReader() {
  auto* efr = get_refr_pointer( event_file_reader_ );
  if ( efr ) delete efr;
}

bool marley::MacroEventFileReader::next_event(marley::Event& ev) {
  auto* efr = get_refr_pointer( event_file_reader_ );
  return efr->next_event( ev );
}

marley::MacroEventFileReader::operator bool() const {
  auto* efr = get_refr_pointer( event_file_reader_ );
  return static_cast<bool>( *efr );
}

double marley::MacroEventFileReader::flux_averaged_xsec( bool natural_units ) const {
  auto* efr = get_refr_pointer( event_file_reader_ );
  return efr->flux_averaged_xsec( natural_units );
}
