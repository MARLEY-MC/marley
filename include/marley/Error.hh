#pragma once
#include <exception>
#include <string>

namespace marley {

  /// @brief Base class for all exceptions thrown by MARLEY functions
  /// @details A distinct exception class allows the user to differentiate
  /// between exceptions thrown by MARLEY and those thrown by other libraries.
  /// @todo Create an MARLEY exception class hierarchy as needed.
  class Error : public std::exception {
    public:
      /// @param message
      /// An error message to display if the exception is not caught
      inline explicit Error(const char* message) : msg_(message) {}
      /// @param message
      /// An error message to display if the exception is not caught
      inline explicit Error(const std::string& message) : msg_(message) {}
      inline virtual ~Error() {}
      /// Method called by the C++ standard library to display the error message
      inline virtual const char* what() const noexcept { return msg_.c_str(); }
    protected:
      std::string msg_; //< error message
  };
}
