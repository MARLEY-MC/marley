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

#pragma once
#include <exception>
#include <string>

namespace marley {

  /// @brief Base class for all exceptions thrown by MARLEY functions
  /// @details A distinct exception class allows the user to differentiate
  /// between exceptions thrown by MARLEY and those thrown by other libraries.
  class Error : public std::exception {

    public:

      /// @param message
      /// An error message to display if the exception is not caught
      inline explicit Error( const char* message ) : msg_(message) {}

      /// @param message
      /// An error message to display if the exception is not caught
      inline explicit Error( const std::string& message ) : msg_(message) {}

      inline virtual ~Error() {}

      /// Method called by the C++ standard library to display the error message
      inline virtual const char* what() const noexcept { return msg_.c_str(); }

    protected:

      std::string msg_; //< error message
  };
}
