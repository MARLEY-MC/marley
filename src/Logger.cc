#include "Logger.hh"

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

marley::Logger::OutStream::OutStream(std::ostream& os, LogLevel lev,
  bool enable) : stream_(&os), level_(lev), enabled_(enable),
  previously_used_(false) {}

marley::Logger::OutStreamVector&
  marley::Logger::OutStreamVector::operator<<(std::ostream&
  (*manip)(std::ostream&))
{
  for (auto s : *this) if (s.enabled_ && s.stream_) manip(*s.stream_);
  return *this;
}

marley::Logger::OutStreamVector&
  marley::Logger::OutStreamVector::operator<<(std::ios_base&
  (*manip)(std::ios_base&))
{
  for (auto s : *this) if (s.enabled_ && s.stream_) manip(*s.stream_);
  return *this;
}

marley::Logger::Logger(bool log_enabled) : enabled_(log_enabled),
  old_level_(LogLevel::INFO) {}

marley::Logger& marley::Logger::Instance() {
  static Logger instance;
  return instance;
}

void marley::Logger::add_stream(std::ostream& stream, LogLevel level)
{
  // Add the new stream if we don't have an OutStream that contains
  // a pointer to the std::ostream object already.
  std::ostream* os = &stream;
  auto end = streams_.end();
  if (std::find_if(streams_.begin(), end, [os](const OutStream& s)
    -> bool { return s.stream_ == os; }) == end)
  {
    bool stream_enabled = (enabled_ && (level >= old_level_));
    streams_.emplace_back(stream, level, stream_enabled);
  }
}

void marley::Logger::enable(bool log_enabled) {
  enabled_ = log_enabled;
  for (auto& s : streams_) s.enabled_ = enabled_ && (s.level_ >= old_level_);
}

marley::Logger::OutStreamVector& marley::Logger::log(LogLevel lev)
{
  bool level_changed = (lev != old_level_);

  if (enabled_) {
    // Update the state of each stream based on the logging level (if
    // needed) and start a new line on each enabled stream for the new
    // message.
    for(auto& s : streams_) {
      if (level_changed) s.enabled_ = (s.level_ >= lev);
      if (s.enabled_) {
        // Insert a new line at the beginning of each logging message
        // unless this is the first message.
        if (s.previously_used_) *s.stream_ << std::endl;
        else s.previously_used_ = true;
        *s.stream_ << loglevel_to_str(lev);
      }
    }
    // Update the old logging level
    old_level_ = lev;
  }
  return streams_;
}
