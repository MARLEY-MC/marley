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

#include <iostream>

#include "marley/Error.hh"
#include "marley/JSON.hh"
#include "marley/Logger.hh"
#include "marley/marley_utils.hh"

const char* marley::Logger::loglevel_to_str(LogLevel lev)
{
  switch (lev) {
    case LogLevel::FATAL:
      return "[FATAL]: ";
      break;
    case LogLevel::ERROR:
      return "[ERROR]: ";
      break;
    case LogLevel::WARN:
      return "[WARNING]: ";
      break;
    case LogLevel::DEBUG:
      return "[DEBUG]: ";
      break;
    case LogLevel::TRACE:
      return "[TRACE]: ";
      break;
    default:
      return "";
      break;
  }
}
marley::Logger::LogLevel marley::Logger::string_to_loglevel(
  const std::string& str )
{
  LogLevel result = LogLevel::INFO;

  // Convert the string to all-lowercase to achieve case-insensitive matching
  std::string low = marley_utils::to_lowercase( str );
  static const std::unordered_map< std::string, LogLevel > conversion_map = {
    { "trace", LogLevel::TRACE },
    { "debug", LogLevel::DEBUG },
    { "info", LogLevel::INFO },
    { "notice", LogLevel::NOTICE },
    { "warn", LogLevel::WARN },
    { "error", LogLevel::ERROR },
    { "fatal", LogLevel::FATAL }
  };
  auto iter = conversion_map.find( low );
  if ( iter != conversion_map.end() ) result = iter->second;
  else throw marley::Error( "Unrecognized marley::Logger logging level \""
    + str + '\"' );
  return result;
}

marley::Logger::OutStream::OutStream( std::shared_ptr< std::ostream > os,
  LogLevel min, LogLevel max ) : stream_( os )
{
  if ( min > max ) throw marley::Error( "Minimum severity more than maximum"
    " severity in constructor of marley::Logger::OutStream" );
  min_level_ = min;
  max_level_ = max;
}

// Alternative constructor that takes a reference to the std::ostream. Problems
// when the std::shared_ptr goes out of scope are avoided by providing a custom
// deleter that does nothing.
marley::Logger::OutStream::OutStream( std::ostream& os, LogLevel min,
  LogLevel max ) : OutStream( std::shared_ptr< std::ostream >( &os,
    [](std::ostream*) -> void {} ), min, max )
{
}

void marley::Logger::configure( const marley::JSON& config ) {

  if ( config.has_key("out") ) {
    auto out_array = config.at( "out" );
    if ( !out_array.is_array() ) {
      throw marley::Error( "JSON array expected for \"out\" key in"
        " marley::Logger JSON configuration. Read \""
        + out_array.dump_string() + '\"' );
    }
    auto streams = out_array.array_range();
    for ( const auto& s : streams ) {

      // Set defaults here that may be overwritten by the configuration
      LogLevel min = LogLevel::TRACE;
      LogLevel max = LogLevel::FATAL;

      if ( s.has_key("min") ) {
        auto jmin = s.at( "min" );
        if ( !jmin.is_string() ) throw marley::Error("Invalid marley::Logger"
          " minimum logging level specification " + jmin.dump_string() );
        min = this->string_to_loglevel( jmin.to_string() );
      }

      if ( s.has_key("max") ) {
        auto jmax = s.at( "max" );
        if ( !jmax.is_string() ) throw marley::Error("Invalid marley::Logger"
          " maximum logging level specification " + jmax.dump_string() );
        max = this->string_to_loglevel( jmax.to_string() );
      }

      bool is_file = s.has_key( "file" );
      bool is_stream = s.has_key( "stream" );
      if ( is_file ) {
        if ( is_stream ) {
          throw marley::Error( "marley::Logger output stream definition uses"
            "both the \"stream\" and \"file\" keys" );
        }
        auto fname = s.at( "file" );
        if ( !fname.is_string() ) throw marley::Error("Invalid marley::Logger"
          " output file specification " + fname.dump_string() );

        // If the user specified a value for the "overwrite" key, use it
        // to determine whether we should append to the file (false) or
        // overwrite it (true). Otherwise, assume we want to append to it.
        auto file_mode = std::ios::out;
        if ( s.has_key("overwrite") ) {

          auto ow = s.at( "overwrite" );

          bool ok;
          bool overwrite = ow.to_bool( ok );
          if ( !ok ) throw marley::Error( "Invalid log file overwrite"
            " setting \"" + ow.dump_string() + '\"' );

          if ( overwrite ) file_mode |= std::ios::trunc;
          else file_mode |= std::ios::app;
        }
        else file_mode |= std::ios::app;

        auto outfile = std::make_shared< std::ofstream >( fname.to_string(),
          file_mode );
        if ( !outfile || (!outfile->good()) ) throw marley::Error( "Unable"
          " to open the log file \"" + fname.to_string() + "\"" );

        // Create the output stream for the file
        this->add_stream( outfile, min, max );
      }
      else {
        auto st = s.at( "stream" );
        if ( !st.is_string() ) throw marley::Error("Invalid marley::Logger"
          " output stream specification " + st.dump_string() );
        auto stream_name = st.to_string();
        if ( stream_name == "stdout" ) {
          this->add_stream( std::cout, min, max );
        }
        else if ( stream_name == "stderr" ) {
          this->add_stream( std::cerr, min, max );
        }
        else throw marley::Error( "Unrecognized stream name \"" + stream_name
          + "\" encountered in marley::Logger::configure()" );
      }

    } // stream definitions
  } // handling of "out" key

  bool set_default_categ = false;
  if ( !config.has_key("categories") ) {
    // If the user has not configured any categories then globally set the
    // logging level to INFO and move on
    default_level_ = LogLevel::INFO;
    return;
  }

  auto categ_spec = config.at( "categories" );
  if ( !categ_spec.is_object() ) {
    throw marley::Error( "marley::Logger categories must be specified as"
      " a JSON object." );
  }
  auto categs = categ_spec.object_range();
  for ( const auto& [ cat, lev ] : categs ) {
    if ( !lev.is_string() ) {
      throw marley::Error( "Invalid marley::Logger level specification \""
        + lev.dump_string() + '\"' );
    }
    LogLevel ll = this->string_to_loglevel( lev.to_string() );
    // The "default" category doesn't appear in the internal map. Instead,
    // it has a dedicated class member to use as the ultimate fallback.
    if ( cat == "default" ) {
      set_default_categ = true;
      default_level_ = ll;
    }
    else if ( cat.empty() ) {
      throw marley::Error( "Empty name encountered in the category"
        " configuration for marley::Logger" );
    }
    else category_map_[ cat ] = ll;
  }

  if ( !set_default_categ ) {
    throw marley::Error( "Missing \"default\" logging level in the category"
      " configuration for marley::Logger" );
  }

}

marley::Logger::Logger() {
  // This is usually done with the FileManager, but we need to avoid logging
  // messages in this constructor (to avoid recursive initialization). So we
  // do it without the FileManager here.
  char* mar = std::getenv( "MARLEY" );
  if ( !mar ) throw marley::Error( "The MARLEY environment variable is not"
    " set. Please set it (e.g., by sourcing the setup_marley.sh script) and"
    " try again." );

  // This works OK because the marley::JSON class does not use the Logger
  std::string config_file_name = std::string( mar ) + "/data/config/logger.js";
  auto json_config = marley::JSON::load_file( config_file_name );
  this->configure( json_config );
}

marley::Logger& marley::Logger::Instance() {
  static Logger instance;
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
  LogLevel min, LogLevel max )
{
  // Check to see whether we have already added this stream to the logger
  marley::Logger::OutStream* os = get_stream( stream.get() );

  // If we don't have it yet, then add it with the requested logging level
  if ( !os ) streams_.emplace_back( stream, min, max );

  // Otherwise, just update the level of the existing one
  else {
    os->min_level_ = min;
    os->max_level_ = max;
  }
}

void marley::Logger::add_stream( std::ostream& stream,
  LogLevel min, LogLevel max )
{
  // Check to see whether we have already added this stream to the logger
  marley::Logger::OutStream* os = get_stream( &stream );

  // If we don't have it yet, then add it with the requested logging level
  if ( !os ) streams_.emplace_back( stream, min, max );

  // Otherwise, just update the level of the existing one
  else {
    os->min_level_ = min;
    os->max_level_ = max;
  }
}

marley::Logger::Message marley::Logger::log( LogLevel lev,
  const std::string& category )
{
  std::vector< std::ostream* > active_streams;
  // The current severity level determines whether the output Message will
  // accept streamed content at all
  if ( this->should_emit( category, lev) ) {
    for( auto& s : streams_ ) {
      // The Message will send the streamed content only to OutStreams whose
      // severity levels are configured to accept it
      if ( lev <= s.max_level_ && lev >= s.min_level_ ) {
        // Store a pointer to the std::ostream object that will receive output
        active_streams.push_back( s.stream_.get() );
      }
    }
  }

  // Returns a Message object to receive the output and relay it to the active
  // OutStreams. Preprends a prefix based on the relevant logging level before
  // accepting other output.
  marley::Logger::Message msg( active_streams, &mutex_ );
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

marley::Logger::LogLevel marley::Logger::category_level(
  const std::string& category )
{
  // First check the cache. If we've already resolved the level for this
  // category, then just use the result. This avoids unnecessary string
  // splitting to check parent category settings.
  auto iter = resolved_level_cache_.find( category );
  if ( iter != resolved_level_cache_.end() ) return iter->second;

  // Position of the delimiter used to mark category hierarchy separations
  size_t delim_pos = std::string::npos;

  // We need to resolve the severity level for a new category. Copy the
  // input so that we can iteratively trim it to scan up the hierarchy.
  std::string categ( category );
  do {

    // A value of delim_pos other than std::string::npos signals that we
    // need to erase the rightmost category name so that we can look up
    // the severity setting for the immediate parent category below.
    if ( delim_pos != std::string::npos ) categ.erase( delim_pos );

    // Check for a setting for the current category. The most specific setting
    // wins, so cache the result and return immediately if one is found
    const auto cit = category_map_.find( categ );
    if ( cit != category_map_.cend() ) {
      LogLevel resolved = cit->second;
      resolved_level_cache_[ categ ] = resolved;
      return resolved;
    }

    // Search for the last category delimiter in the current category string
    delim_pos = categ.rfind( CATEG_DELIM_ );

    // If one was not found, then exit the loop so we can fall back to the
    // default severity level
  } while( delim_pos != std::string::npos );

  // A specific category setting was not found, so fall back to the default
  // severity level
  return default_level_;
}
