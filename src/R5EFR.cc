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
