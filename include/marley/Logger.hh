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

#pragma once
#include <algorithm>
#include <fstream>
#include <memory>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <vector>

// Define numerical values for the severity levels in a form that can be
// understood by the preprocessor (enums aren't available yet)
#define MARLEY_LOGGER_LEVEL_TRACE  0
#define MARLEY_LOGGER_LEVEL_DEBUG  1
#define MARLEY_LOGGER_LEVEL_INFO   2
#define MARLEY_LOGGER_LEVEL_NOTICE 3
#define MARLEY_LOGGER_LEVEL_WARN   4
#define MARLEY_LOGGER_LEVEL_ERROR  5
#define MARLEY_LOGGER_LEVEL_FATAL  6

// Double-expansion helper that allows a severity level token (e.g., INFO) to
// be assigned to MARLEY_COMPILED_LOG_LEVEL and matched to one of the numerical
// values above. This is used below to enable numerical comparison in
// MARLEY_LOG_IMPL without needing to reference the numerical values above
// when configuring CMake or GNU Make.
#define MARLEY_LOG_LEVEL_NUM(level) MARLEY_LOG_LEVEL_NUM_IMPL(level)
#define MARLEY_LOG_LEVEL_NUM_IMPL(level) MARLEY_LOGGER_LEVEL_##level

// Default to maximum severity of INFO (this can be overriden via injection
// of a different definition for this macro at compile time). Anything below
// this is blocked from execution at compile time, preventing
// debugging messages from impairing runtime performance when they are not
// enabled.
#ifndef MARLEY_COMPILED_LOG_LEVEL
  #define MARLEY_COMPILED_LOG_LEVEL INFO
#endif

// Main user-facing macro (accepts one or two arguments depending on
// whether a category is specified)
#define MARLEY_LOG(...) \
    MARLEY_LOG_SELECT(__VA_ARGS__, \
        MARLEY_LOG_2, \
        MARLEY_LOG_1) \
    (__VA_ARGS__)

#define MARLEY_LOG_SELECT(_1,_2,NAME,...) NAME

#define MARLEY_LOG_1(level) MARLEY_LOG_IMPL(level,"")
#define MARLEY_LOG_2(level, category) MARLEY_LOG_IMPL(level,category)

#define MARLEY_LOG_IMPL(level, category) \
  ( MARLEY_LOG_LEVEL_NUM(MARLEY_COMPILED_LOG_LEVEL) \
    <= MARLEY_LOGGER_LEVEL_##level ) && \
  marley::Logger::Instance().log( marley::Logger::LogLevel::level, category )

// Forward declare some MARLEY classes and their operator<< functions so that
// we can stream them to the Logger
namespace marley {
  class HauserFeshbachDecay;
  class JSON;
  class Parity;
  class Target;
  class TargetAtom;
}

std::ostream& operator<<( std::ostream& out,
  const marley::HauserFeshbachDecay& hfd );

std::ostream& operator<<( std::ostream& os, const marley::JSON& json );
std::ostream& operator<<( std::ostream& out, const marley::Parity& p );
std::ostream& operator<<( std::ostream& out, const marley::Target& t );
std::ostream& operator<<( std::ostream& out, const marley::TargetAtom& ta );

namespace marley {

  /// @brief Simple singleton logging class
  /// @details This class is based on code found <a
  /// href="http://tinyurl.com/h5oevxo">here</a>.
  class Logger {

    public:

      /// @brief Defines the logging levels recognized by the marley::Logger.
      /// @details Note that we tie the numerical values of the levels
      /// to the previously-defined pre-processor macros to ensure consistency.
      /// The levels appear in increasing numerical order (increasing severity).
      enum class LogLevel {
        TRACE  = MARLEY_LOGGER_LEVEL_TRACE,
        DEBUG  = MARLEY_LOGGER_LEVEL_DEBUG,
        INFO   = MARLEY_LOGGER_LEVEL_INFO,
        NOTICE = MARLEY_LOGGER_LEVEL_NOTICE,
        WARN   = MARLEY_LOGGER_LEVEL_WARN,
        ERROR  = MARLEY_LOGGER_LEVEL_ERROR,
        FATAL  = MARLEY_LOGGER_LEVEL_FATAL,
      };

    private:

      /// @brief Returns the string to prepend to the logger message
      /// based on the given logging level
      const char* loglevel_to_str( LogLevel lev );

      /// @brief Wrapped std::ostream together with its logging configuration
      class OutStream {

        friend class Logger;

        public:

          /// @param os std::shared_ptr that points to a std::ostream object
          /// that will receive logging messages
          /// @param lev Logging level for this stream
          /// @note Because std::shared_ptr will auto-delete the pointed-to
          /// object when the use_count falls to zero (unless a custom deleter
          /// is used), this version of the OutStream constructor <b>should
          /// not</b> be used with the std::cout or std::cerr streams. For
          /// streams that should never be deleted when the OutStream is
          /// destroyed, use @see OutStream::(std::ostream& os)
          /// @exception marley::Error std::cout or std::cerr was passed to
          /// this function. This could still be permissible if a suitable
          /// custom deleter was used to define the shared_ptr, but for safety,
          /// one should use @see OutStream::(std::ostream& os)
          OutStream( std::shared_ptr< std::ostream > os,
            LogLevel min, LogLevel max );

          /// @param [in] os Reference to std::ostream object that will receive
          /// logging messages
          /// @param lev Logging level for this stream
          /// @note Because this version of the constructor uses a bare
          /// reference instead of a std::shared_ptr, the std::ostream is not
          /// guaranteed to exist throughout the life of the OutStream. If the
          /// std::ostream goes out of scope before the OutStream does, logging
          /// to the std::ostream will silently stop. The alternative
          /// constructor @see OutStream::(std::shared_ptr<std::ostream> os)
          /// is recommended for use with all streams except for std::cout and
          /// std::cerr.
          OutStream( std::ostream& os, LogLevel min, LogLevel max );

        private:

          /// @brief Pointer to a std::ostream that will receive logging
          /// messages
          std::shared_ptr< std::ostream > stream_;

          /// @brief Minimum severity to accept for this OutStream
          LogLevel min_level_;

          /// @brief Maximum severity to accept for this OutStream
          LogLevel max_level_;
      };

    public:

      /// @brief Temporary object used for forming logger messages
      /// @details Upon destruction, a newline is appended to the message. This
      /// is based on a trick from https://stackoverflow.com/a/57553824
      class Message {

        public:

          Message( std::vector< std::ostream* >& vec, std::mutex* mtx )
            : osvec_( vec ), mtx_( mtx ) {}

          // Copy constructors cannot be defaulted because the class owns
          // a std::ostringstream, which has a deleted copy constructor.
          // However, the stream can be moved, so we take advantage of this
          // here to allow returning a Message by value.
          Message( Message&& other ) = default;
          Message& operator=( Message&& other ) = default;

          // Allows conversion to bool to get the types to work in the MARLEY_LOG_IMPL macro.
          // The return value is not intended to be used anywhere.
          inline explicit operator bool() { return true; }

          inline ~Message() {
            // If we have no active OutStreams, then just destroy the Message
            // without any output
            if ( osvec_.empty() ) return;
            // Otherwise, append a new line, apply a lock to protect against
            // activity from other threads, and write to each active OutStream
            buffer_ << '\n';
            std::string msg = buffer_.str();
            std::lock_guard< std::mutex > lock( *mtx_ );
            for ( auto& s : osvec_ ) *s << msg;
          }

          template< typename OutputType > Message&
            operator<<( const OutputType& out )
          {
            // If we have no active OutStreams, then don't both to store the
            // streaming output (since it will not be sent anywhere by the
            // destructor
            if ( !osvec_.empty() ) buffer_ << out;
            return *this;
          }

          /// @brief Allows the Message to accept standard stream manipulators
          /// (like std::endl) via operator<<
          /// @note Code for this function is based on a trick discussed here:
          /// http://www.cplusplus.com/forum/general/54588/#msg294798
          Message& operator<<( std::ostream& (*manip)(std::ostream&) );

          /// @brief Backup overload for extra output manipulators
          Message& operator<<( std::ios_base& (*manip)(std::ios_base&) );

        protected:

          /// Used to cache output from multiple calls to operator<< until we
          /// have the full message. Then we just emit it in a single-step
          /// for each active OutStream, all at once. This together with a
          /// std::mutex prevents splitting logging messages due to race
          /// conditions between parallel threads.
          std::ostringstream buffer_;

          /// Active output streams that should receive the Message
          std::vector< std::ostream* > osvec_;

          /// std::mutex to use to lock streaming output when emitting the
          /// final message
          std::mutex* mtx_;
      };

      /// @brief Create the singleton Logger
      Logger();

      /// @brief Initialize the Logger using settings expressed as a JSON
      /// object
      void configure( const marley::JSON& json );

      static LogLevel string_to_loglevel( const std::string& str );

      /// @brief Get the singleton instance of the Logger class
      static Logger& Instance();

      /// @brief Returns true if stream is already registered with the
      /// Logger, or false otherwise
      bool has_stream( const std::ostream& stream ) const;

      /// @brief Prepare the Logger to receive a log message via
      /// the << stream operator
      /// @param lev marley::Logger::LogLevel of the incoming message
      /// @param category Named category of the incoming message
      Message log( LogLevel lev, const std::string& category = "" );

      /// @brief Returns whether the logger should emit a message for the
      /// given category and level
      inline bool should_emit( const std::string& category, LogLevel lev );

      /// @brief Looks up the severity setting for the input category,
      /// including inheritance from the hierarchy
      LogLevel category_level( const std::string& category );

      // Make the singleton Logger uncopyable and unmovable
      /// @brief Deleted copy constructor
      Logger( const Logger& ) = delete;
      /// @brief Deleted copy assignment operator
      Logger& operator=( const Logger& ) = delete;
      /// @brief Deleted move constructor
      Logger( Logger&& ) = delete;
      /// @brief Deleted move assignment operator
      Logger& operator=( Logger&& ) = delete;

    private:

      /// @brief Add a std::ostream to the vector of streams that will receive
      /// Logger output
      /// @note The stream owned by the std::shared_ptr should have been
      /// dynamically allocated since std::shared_ptr will auto-delete it when
      /// use_count falls to zero (unless a suitable custom deleter was used).
      /// For adding std::cout or std::cerr to the Logger, please use @see
      /// add_stream(std::ostream& stream, LogLevel level) instead. See also
      /// the documentation for @see marley::Logger::OutStream::OutStream.
      void add_stream( std::shared_ptr< std::ostream > stream,
        LogLevel min, LogLevel max );

      /// @brief Add a std::ostream to the vector of streams that will receive
      /// Logger output
      /// @note For streams other than std::cout and std::cerr, using @see
      /// add_stream( std::shared_ptr< std::ostream > stream, LogLevel level )
      /// instead of this function is recommended. See also the documentation
      /// for @see marley::Logger::OutStream::OutStream.
      void add_stream( std::ostream& stream, LogLevel min, LogLevel max );

      // @brief Returns a pointer to the given stream's OutStream object if
      // it has been added to the Logger, or nullptr otherwise.
      const OutStream* get_stream( const std::ostream* os ) const;

      // @brief Returns a pointer to the given stream's OutStream object if
      // it has been added to the Logger, or nullptr otherwise.
      OutStream* get_stream( const std::ostream* os );

      /// @brief Vector of wrapped std::ostream objects that will
      /// receive the log messages
      std::vector< OutStream > streams_;

      /// @brief Default logging level for emitting messages
      /// @details This value is the ultimate fallback for missing
      /// category-specific levels
      LogLevel default_level_ = LogLevel::INFO;

      /// @brief Stores configuration of severity levels for each configured
      /// category
      std::unordered_map< std::string, LogLevel > category_map_;

      /// @brief Caches resolved fully-qualified severity levels (which can
      /// differ when using hierarchical categories). This avoids
      /// string splitting beyond the first lookup.
      /// @note The cache becomes invalidated and must be cleared in response to
      /// a change of either category_map_ or default_level_
      std::unordered_map< std::string, LogLevel > resolved_level_cache_;

      /// @brief Used to avoid race conditions when emitting logging messages
      std::mutex mutex_;

      /// @brief Delimiter used to separate sub-categories in labels
      static constexpr char CATEG_DELIM_ = '.';
  };

}

// Inline function definitions
inline bool marley::Logger::should_emit( const std::string& category,
  LogLevel lev )
{
  LogLevel cl = this->category_level( category );
  return ( lev >= cl );
}

// Convenient shortcut functions for recording log messages
inline auto MARLEY_LOG_ERROR() {
  return marley::Logger::Instance().log( marley::Logger::LogLevel::ERROR );
}

inline auto MARLEY_LOG_WARNING() {
  return marley::Logger::Instance().log( marley::Logger::LogLevel::WARNING );
}

inline auto MARLEY_LOG_INFO() {
  return marley::Logger::Instance().log( marley::Logger::LogLevel::INFO );
}

inline auto MARLEY_LOG_DEBUG() {
  return marley::Logger::Instance().log( marley::Logger::LogLevel::DEBUG );
}
