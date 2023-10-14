/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
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

// Standard library includes
#include <ostream>

namespace marley {

  /// @brief Type-safe representation of a parity value (either +1 or -1)
  class Parity {

    public:

      /// @brief Default constructor chooses positive parity
      constexpr Parity() : is_positive_(true) {}

      /// @brief Create a Parity object from a boolean value
      /// @param is_positive Boolean parameter that is true for +1,
      /// or false for -1.
      constexpr explicit Parity(bool is_positive)
        : is_positive_(is_positive) {}

      /// @brief Create a Parity object from an integer
      /// @details If an int value other than +1 or -1 is passed to this
      /// constructor, a marley::Error will be thrown.
      explicit Parity(int i);

      /// @brief Convert the Parity object to a bool (true for +1, false for -1)
      constexpr explicit operator bool() const { return is_positive_; }

      /// @brief Creates a copy of the Parity object with a flipped value
      constexpr marley::Parity operator-() const {
        return marley::Parity(!is_positive_);
      }

      /// @brief Unary operator that flips a Parity value in place
      inline Parity& operator!() {
        is_positive_ = !is_positive_;
        return *this;
      }

      /// @brief Assigns a parity value using an integer
      /// @details If an int value other than +1 or -1 is used, a marley::Error
      /// will be thrown.
      marley::Parity& operator=(const int& i);

      inline bool operator==(const marley::Parity& p) const
      {
        return is_positive_ == p.is_positive_;
      }

      inline bool operator!=(const marley::Parity& p) const
      {
        return is_positive_ != p.is_positive_;
      }

      inline marley::Parity operator*(const marley::Parity& p) const
      {
        return marley::Parity(is_positive_ == p.is_positive_);
      }

      inline int operator*(const int& i) const
      {
        if (is_positive_) return i;
        else return -i;
      }

      // Allows explicit casts of marley::Parity to int
      /// @todo See comments for marley::TargetAtom::operator int().
      /// A similar situation occurs for this function.
      inline explicit operator int() const {
        if (is_positive_) return 1;
        else return -1;
      }

      /// @brief Convert the Parity object to a char
      inline char to_char() const {
        if (is_positive_) return '+';
        else return '-';
      }

      /// @brief Assign a value to this Parity object from a char
      void from_char(const char c);

    protected:

      /// @brief Boolean representation of the parity value that is true when
      /// the parity is +1 and false when it is -1.
      bool is_positive_;
  };

}

inline int operator*(const int& i, const marley::Parity& p)
{
  if (static_cast<bool>(p)) return i;
  else return -i;
}

inline std::ostream& operator<<(std::ostream& out, const marley::Parity& p)
{
  out << p.to_char();
  return out;
}

std::istream& operator>> (std::istream& in, marley::Parity& p);
