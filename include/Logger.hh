#pragma once
#include <algorithm>
#include <fstream>
#include <vector>

#include "HauserFeshbachDecay.hh"

namespace marley {

  // Simple logging class based on code found here: http://tinyurl.com/h5oevxo
  // If you need something more sophisticated in the future, a number of
  // logging packages are available. Try spdlog
  // (https://github.com/gabime/spdlog), for instance.
  class Logger {

    public:

      // The different logging levels are defined here. Note that C++
      // automatically assigns ascending values to the enum class members
      // (which have an underlying integral type) in the order that they are
      // written here, so ERROR < WARNING < INFO < DEBUG.
      enum class LogLevel { ERROR, WARNING, INFO, DEBUG };

    private:

      const char* loglevel_to_str(LogLevel lev);

      class OutStream {
        friend class Logger;
        friend class OutStreamVector;
        public:
          OutStream(std::ostream& os, LogLevel lev, bool enable = true);
        private:
          std::ostream* stream_;
          LogLevel level_;
          bool enabled_;
          bool previously_used_;
      };

      class OutStreamVector : public std::vector<OutStream> {
        public:
          OutStreamVector() {}

          template<typename OutputType> OutStreamVector&
            operator<<(const OutputType& ot);

          // Allows the logger to use the << operator on output manipulators
          // like std::endl. Based on a trick discussed here:
          // http://www.cplusplus.com/forum/general/54588/#msg294798
          OutStreamVector& operator<<(std::ostream& (*manip)(std::ostream&));

	  // Backup overload for extra output manipulators
          OutStreamVector& operator<<(std::ios_base& (*manip)(std::ios_base&));
      };

      Logger(bool log_enabled = true);

    public:

      static Logger& Instance();

      void add_stream(std::ostream& stream,
        LogLevel level = LogLevel::WARNING);

      inline void clear_streams();

      void enable(bool log_enabled = true);

      OutStreamVector& log(LogLevel lev = LogLevel::WARNING);

      // Make the logger uncopyable and unmovable
      Logger(const Logger&) = delete;
      Logger& operator=(const Logger&) = delete;
      Logger(Logger&&) = delete;
      Logger& operator=(Logger&&) = delete;

    private:

      OutStreamVector streams_;
      bool enabled_;
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
