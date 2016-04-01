#pragma once
#include <exception>
#include <string>

namespace marley {

  // Base class for all exceptions thrown by MARLEY functions.
  // TODO: create an exception class hierarchy as needed.
  class Error : public std::exception {
    public:
      inline explicit Error(const char* message) : msg(message) {}
      inline explicit Error(const std::string& message) : msg(message) {}
      inline virtual ~Error() {}
      inline virtual const char* what() const noexcept { return msg.c_str(); }
    protected:
      std::string msg;
  };
}
