#pragma once
#include <algorithm>
#include <fstream>
#include <vector>

namespace marley {

  // Simple logging class based on code found here:
  // https://github.com/SFML/SFML/wiki/Source:-Simple-File-Logger-%28by-mateandmetal%29
  // If you need something more sophisticated in the future, a number of logging
  // packages are available. Try spdlog (https://github.com/gabime/spdlog), for
  // instance.
  class Logger {

    public:

      // The different logging levels are defined here. Note that C++
      // automatically assigns ascending values to the enum class members (which
      // have an underlying integral type) in the order that they are written
      // here, so ERROR < WARNING < INFO < DEBUG.
      enum class LogLevel { ERROR, WARNING, INFO, DEBUG };

    private:

      inline const char* loglevel_to_str(LogLevel lev) {
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

      class OutStream {
        friend class Logger;
        friend class OutStreamVector;
        public:
          inline OutStream(std::ostream& os, LogLevel lev, bool enable = true) :
            stream(&os), level(lev), enabled(enable), previously_used(false) {}

        private:
          std::ostream* stream;
          LogLevel level;
          bool enabled;
          bool previously_used;
      };

      class OutStreamVector : public std::vector<OutStream> {
        public:
          inline OutStreamVector() {}

          template<typename OutputType> OutStreamVector&
            operator<<(const OutputType& ot)
          {
            for (auto s : *this) if (s.enabled && s.stream) *(s.stream) << ot;
            return *this;
          }

          // Allows the logger to use the << operator on output manipulators like
          // std::endl. Based on a trick discussed here:
          // http://www.cplusplus.com/forum/general/54588/#msg294798 
          inline OutStreamVector& operator<<(std::ostream& (*manip)(std::ostream&))
          {
            for (auto s : *this) if (s.enabled && s.stream) manip(*s.stream);
            return *this;
          }

	  // Backup overload for extra output manipulators
          inline OutStreamVector& operator<<(std::ios_base& (*manip)(std::ios_base&))
          {
            for (auto s : *this) if (s.enabled && s.stream) manip(*s.stream);
            return *this;
          }
      };

      inline Logger(bool log_enabled = true) : enabled(log_enabled),
        old_level(LogLevel::INFO) {}

      inline virtual ~Logger() {}

    public:

      static Logger& Instance() {
        static Logger instance;
        return instance;
      }

      inline void add_stream(std::ostream& stream,
        LogLevel level = LogLevel::WARNING)
      {
        // Add the new stream if we don't have an OutStream that contains
        // a pointer to the std::ostream object already.
        std::ostream* os = &stream;
        auto end = streams.end();
        if (std::find_if(streams.begin(), end, [os](const OutStream& s)
          -> bool { return s.stream == os; }) == end)
        {
          bool stream_enabled = (enabled && (level >= old_level));
          streams.emplace_back(stream, level, stream_enabled);
        }
      }

      inline void clear_streams() { streams.clear(); }

      inline void enable(bool log_enabled = true) {
        enabled = log_enabled;
        for (auto& s : streams) s.enabled = enabled && (s.level >= old_level);
      }

      inline OutStreamVector& log(LogLevel lev = LogLevel::WARNING) {
        bool level_changed = (lev != old_level);
        if (enabled) {
          // Update the state of each stream based on the logging level (if needed)
          // and start a new line on each enabled stream for the new message.
          for(auto& s : streams) {
            if (level_changed) s.enabled = (s.level >= lev);
            if (s.enabled) {
	      // Insert a new line at the beginning of each logging message
	      // unless this is the first message.
              if (s.previously_used) *s.stream << std::endl;
              else s.previously_used = true;
              *s.stream << loglevel_to_str(lev);
            }
          }
          // Update the old logging level
          old_level = lev;
        }
        return streams;
      }

      // Make the logger uncopyable
      Logger(const Logger&) = delete;
      Logger& operator=(const Logger&) = delete;
  
    private:
  
      OutStreamVector streams;

      bool enabled;

      LogLevel old_level;
  };

}

#define LOG_ERROR marley::Logger::Instance().log(marley::Logger::LogLevel::ERROR)
#define LOG_WARNING marley::Logger::Instance().log(marley::Logger::LogLevel::WARNING)
#define LOG_INFO marley::Logger::Instance().log(marley::Logger::LogLevel::INFO)
#define LOG_DEBUG marley::Logger::Instance().log(marley::Logger::LogLevel::DEBUG)
