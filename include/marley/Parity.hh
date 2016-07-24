#pragma once
#include <string>

#include "Error.hh"

namespace marley {

  /// @brief Type-safe representation of a parity value (either +1 or -1)
  class Parity {

    public:

      /// @brief Default constructor chooses positive parity
      inline Parity() : is_positive_(true) {}

      /// @brief Create a Parity object from a boolean value
      /// @param is_positive Boolean parameter that is true for +1,
      /// or false for -1.
      inline explicit Parity(bool is_positive)
        : is_positive_(is_positive) {}

      /// @brief Create a Parity object from an integer
      /// @details If an int value other than +1 or -1 is passed to this
      /// constructor, a marley::Error will be thrown.
      inline explicit Parity(int i) {
        if (i == 1) is_positive_ = true;
        else if (i == -1) is_positive_ = false;
        else throw marley::Error(std::string("Invalid parity ")
          + std::to_string(i) + " passed to constructor of marley::Parity");
      }

      /// @brief Convert the Parity object to a bool (true for +1, false for -1)
      inline explicit operator bool() const { return is_positive_; }

      /// @brief Creates a copy of the Parity object with a flipped value
      inline marley::Parity operator-() const {
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
      inline marley::Parity& operator=(const int& i) {
        // Do the assignment
        if (i == 1) is_positive_ = true;
        else if (i == -1) is_positive_ = false;
        else throw marley::Error(std::string("Invalid parity ")
          + std::to_string(i) + " assigned to variable of type marley::Parity");

        // Return the existing object
        return *this;
      }

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
      inline explicit operator int() const {
        if (is_positive_) return 1;
        else return -1;
      }

      /// @brief Convert the Parity object to a char
      inline char to_char() const {
        if (is_positive_) return '+';
        else return '-';
      }

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

inline std::istream& operator>> (std::istream& in,
  marley::Parity& p)
{
  char c;
  in >> c;

  if (c == '+') p = marley::Parity(true);
  else if (c == '-') p = marley::Parity(false);
  else throw marley::Error(std::string("Invalid parity ")
    + c + " assigned via the >> operator to"
    + " a variable of type marley::Parity");

  return in;
}
