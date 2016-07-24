#pragma once
#include <algorithm>
#include <fstream>
#include <vector>

#include "marley/HauserFeshbachDecay.hh"

namespace marley {

  /// @brief Simple singleton logging class
  /// @details This class is based on code found <a
  /// href="http://tinyurl.com/h5oevxo">here</a>.
  class Logger {

    public:

      /// @brief Defines the logging levels recognized by the marley::Logger.
      /// @details Note that C++ automatically assigns ascending values to the
      /// enum class members (which have an underlying integral type) in the
      /// order that they are written in the definition, so ERROR < WARNING <
      /// INFO < DEBUG.
      enum class LogLevel { ERROR, WARNING, INFO, DEBUG };

    private:

      /// @brief Returns the string to prepend to the logger message
      /// based on the given logging level
      const char* loglevel_to_str(LogLevel lev);

      /// @brief Wrapped std::ostream together with its logging configuration
      class OutStream {

        friend class Logger;
        friend class OutStreamVector;

        public:

          /// @param os std::ostream object that will receive logging messages
          /// @param lev marley::Logger::LogLevel specifier that indicates
          /// what the logging level should be for this stream
          /// @param enable Whether to enable (true) or disable (false)
          /// logging to this stream upon construction
          OutStream(std::ostream& os, LogLevel lev, bool enable = true);

        private:

          /// @brief Pointer to a std::ostream that will receive logging
          /// messages
          std::ostream* stream_;

          /// @brief Logging level
          LogLevel level_;

          /// @brief Whether logging to the stream is enabled (true) or
          /// disabled (false)
          bool enabled_;

          /// @brief Whether the stream has previously received at least one
          /// logging message
          bool previously_used_;
      };

      /// @brief A std::vector of OutStream objects that provides a stream
      /// output << operator for sending logging messages to the streams
      class OutStreamVector : public std::vector<OutStream> {
        public:

          OutStreamVector() : std::vector<OutStream>() {}

          template<typename OutputType> OutStreamVector&
            operator<<(const OutputType& ot);

          /// @brief Allows the logger to use the << operator on output
          /// manipulators like std::endl.
          /// @note Code for this function is based on a trick discussed here:
          /// http://www.cplusplus.com/forum/general/54588/#msg294798
          OutStreamVector& operator<<(std::ostream& (*manip)(std::ostream&));

          /// @brief Backup overload for extra output manipulators
          OutStreamVector& operator<<(std::ios_base& (*manip)(std::ios_base&));
      };

      /// @brief Create the singleton Logger
      /// @param log_enabled Create the Logger so that it is
      /// initially enabled (true) or disabled (false)
      Logger(bool log_enabled = true);

    public:

      /// @brief Get the singleton instance of the Logger class
      static Logger& Instance();

      /// @brief Add a std::ostream to the vector of streams that will receive
      /// Logger output
      void add_stream(std::ostream& stream,
        LogLevel level = LogLevel::WARNING);

      /// @brief Clear the vector of streams that receive Logger output
      inline void clear_streams();

      /// @brief Enable or disable the Logger
      /// @param log_enabled Enable (true) or disable (false) the Logger
      void enable(bool log_enabled = true);

      /// @brief Prepare the Logger to receive a log message via
      /// the << stream operator
      /// @param lev marley::Logger::LogLevel of the incoming message
      OutStreamVector& log(LogLevel lev = LogLevel::WARNING);

      // Make the singleton Logger uncopyable and unmovable
      /// @brief Deleted copy constructor
      Logger(const Logger&) = delete;
      /// @brief Deleted copy assignment operator
      Logger& operator=(const Logger&) = delete;
      /// @brief Deleted move constructor
      Logger(Logger&&) = delete;
      /// @brief Deleted move assignment operator
      Logger& operator=(Logger&&) = delete;

    private:

      /// @brief Vector of wrapped std::ostream objects that will
      /// receive the log messages
      OutStreamVector streams_;

      /// @brief Whether the Logger is currently enabled
      bool enabled_;

      /// @brief LogLevel of the last log message
      LogLevel old_level_;
  };

}

// Inline function definitions
inline void marley::Logger::clear_streams() { streams_.clear(); }

template<typename OutputType> marley::Logger::OutStreamVector&
  marley::Logger::OutStreamVector::operator<<(const OutputType& ot)
{
  for (auto s : *this) if (s.enabled_ && s.stream_) *(s.stream_) << ot;
  return *this;
}

// Convenient shortcut functions for recording log messages
inline auto MARLEY_LOG_ERROR() {
  return marley::Logger::Instance().log(marley::Logger::LogLevel::ERROR);
}

inline auto MARLEY_LOG_WARNING() {
  return marley::Logger::Instance().log(marley::Logger::LogLevel::WARNING);
}

inline auto MARLEY_LOG_INFO() {
  return marley::Logger::Instance().log(marley::Logger::LogLevel::INFO);
}

inline auto MARLEY_LOG_DEBUG() {
  return marley::Logger::Instance().log(marley::Logger::LogLevel::DEBUG);
}
