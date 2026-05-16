/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
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

#include <iostream>
#include <exception>

#include "marley/Error.hh"
#include "marley/Logger.hh"

const char* marley::Logger::loglevel_to_str(LogLevel lev)
{
  switch (lev) {
    case LogLevel::ERROR:
      return "[ERROR]: ";
      break;
    case LogLevel::WARNING:
      return "[WARNING]: ";
      break;
    case LogLevel::DEBUG:
      return "[DEBUG]: ";
      break;
    default:
      return "";
      break;
  }
}

marley::Logger::OutStream::OutStream( std::shared_ptr< std::ostream > os,
  LogLevel lev ) : stream_( os ), level_( lev )
{
  constexpr char error_1[] = "std::shared_ptr to ";
  constexpr char error_2[] = " passed to constructor of marley::Logger"
    "::OutStream. Please use the alternate constructor (which takes a"
    " std::ostream&) instead.";
  if ( stream_.get() == &std::cout ) throw marley::Error( error_1
    + std::string( "std::cout" ) + error_2 );
  else if ( stream_.get() == &std::cerr ) throw marley::Error( error_1
    + std::string( "std::cerr" ) + error_2 );
}

marley::Logger::OutStream::OutStream( std::ostream& os,
  LogLevel lev ) : level_( lev )
{
  // Avoid any deletion problems when the std::shared_ptr goes out of scope by
  // providing a custom deleter that does nothing.
  stream_ = std::shared_ptr< std::ostream >( &os,
    [](std::ostream*) -> void {} );
}

void marley::Logger::clear_streams() {
  // Empty the vector of output streams for the Logger
  streams_.clear();

  // Unless the user changes the stderr logging level via a call to
  // add_stream(), the Logger will always write warning and error messages to
  // stderr.
  add_stream( std::cerr, LogLevel::WARNING );
}

marley::Logger::Logger( LogLevel lev ) : level_( lev )
{
  // Unless the user changes the stderr logging level via a call to
  // add_stream(), the Logger will always write warning and error messages to
  // stderr.
  add_stream( std::cerr, LogLevel::WARNING );

  // By default, also add std::cout at the info level
  add_stream( std::cout, LogLevel::INFO );
}

marley::Logger& marley::Logger::Instance() {
  static Logger instance( LogLevel::INFO );
  return instance;
}

bool marley::Logger::has_stream( const std::ostream& os ) const {
  auto stream = get_stream( &os );
  if ( stream ) return true;
  // A nullptr was returned, so the stream couldn't be found
  else return false;
}

const marley::Logger::OutStream* marley::Logger::get_stream(
  const std::ostream* os ) const
{
  auto end = streams_.end();
  auto iter = std::find_if( streams_.begin(), end,
    [ os ]( const OutStream& s ) -> bool { return s.stream_.get() == os; }
  );
  if ( iter == end ) return nullptr;
  else return &( *iter );
}

marley::Logger::OutStream* marley::Logger::get_stream(
  const std::ostream* os )
{
  auto end = streams_.end();
  auto iter = std::find_if( streams_.begin(), end,
    [ os ]( const OutStream& s ) -> bool { return s.stream_.get() == os; }
  );
  if ( iter == end ) return nullptr;
  else return &( *iter );
}

void marley::Logger::add_stream( std::shared_ptr< std::ostream > stream,
  LogLevel level )
{
  // Check to see whether we have already added this stream to the logger
  marley::Logger::OutStream* os = get_stream( stream.get() );

  // If we don't have it yet, then add it with the requested logging level
  if ( !os ) streams_.emplace_back( stream, level );

  // Otherwise, just update the level of the existing one
  else os->level_ = level;
}

void marley::Logger::add_stream( std::ostream& stream, LogLevel level )
{
  // Check to see whether we have already added this stream to the logger
  marley::Logger::OutStream* os = get_stream( &stream );

  // If we don't have it yet, then add it with the requested logging level
  if ( !os ) streams_.emplace_back( stream, level );

  // Otherwise, just update the level of the existing one
  else os->level_ = level;
}

marley::Logger::Message marley::Logger::log( LogLevel lev ) {

  std::vector< OutStream* > active_streams_;
  if ( lev <= level_ ) {
    // Include each stream in the output Message if it needs to receive
    // output based on its configured logging level
    for( auto& s : streams_ ) {
      bool active = ( lev <= s.level_ );
      // Because writing error and warning messages to stdout will cause
      // duplication in a terminal when we're also writing to stderr, prevent
      // stdout from receiving any logger messages that are at the WARNING
      // log level or below. The user may suppress warnings/errors entirely
      // by adjusting the log level for std::cerr.
      if ( s.stream_.get() == &std::cout ) {
        active = active && ( lev > LogLevel::WARNING );
      }

      if ( active ) active_streams_.push_back( &s );
    }
  }

  // Returns a Message object to receive the output and relay it to the active
  // OutStreams. Preprends a prefix based on the relevant logging level before
  // accepting other output.
  marley::Logger::Message msg( active_streams_, &mutex_ );
  msg << loglevel_to_str( lev );
  return msg;
}

marley::Logger::Message& marley::Logger::Message::operator<<( std::ostream&
  (*manip)(std::ostream&) )
{
  buffer_ << manip;
  return *this;
}

marley::Logger::Message& marley::Logger::Message::operator<<( std::ios_base&
  (*manip)(std::ios_base&) )
{
  buffer_ << manip;
  return *this;
}
