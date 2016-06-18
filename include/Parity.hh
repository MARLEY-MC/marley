#pragma once
#include <iostream>
#include <string>

#include "Error.hh"

namespace marley {

  // Class that represents a parity value (either +1 or -1) in a type-safe way
  class Parity {
    public:
      // Default constructor chooses positive parity
      inline Parity() : is_positive(true) {}

      inline explicit Parity(bool is_pos) { is_positive = is_pos; }

      inline explicit Parity(int i) {
        if (i == 1) is_positive = true;
        else if (i == -1) is_positive = false;
        else throw marley::Error(std::string("Invalid parity ")
          + std::to_string(i) + " passed to constructor of marley::Parity");
      }

      // Copy constructor
      inline Parity(const marley::Parity& p) { is_positive = p.is_positive; }

      inline explicit operator bool() const { return is_positive; }

      // The minus unary operator creates a copy of a parity
      // object with a flipped value
      inline marley::Parity operator-() const {
        return marley::Parity(!is_positive);
      }

      // The not unary operator flips a parity value in place
      inline Parity& operator!() {
        is_positive = !is_positive;
        return *this;
      }

      inline marley::Parity& operator=(const marley::Parity& p) {
        // Do the copy
        is_positive = p.is_positive;

        // Return the existing object
        return *this;
      }

      inline marley::Parity& operator=(const int& i) {
        // Do the assignment
        if (i == 1) is_positive = true;
        else if (i == -1) is_positive = false;
        else throw marley::Error(std::string("Invalid parity ")
          + std::to_string(i) + " assigned to variable of type marley::Parity");

        // Return the existing object
        return *this;
      }

      inline bool operator==(const marley::Parity& p) const
      {
        return is_positive == p.is_positive;
      }

      inline bool operator!=(const marley::Parity& p) const
      {
        return is_positive != p.is_positive;
      }

      inline marley::Parity operator*(const marley::Parity& p) const
      {
        return marley::Parity(is_positive == p.is_positive);
      }

      inline int operator*(const int& i) const
      {
        if (is_positive) return i;
        else return -i;
      }

      // Allows explicit casts of marley::Parity to int
      inline explicit operator int() const {
        if (is_positive) return 1;
        else return -1;
      }

      // Converts the parity to a char
      inline char to_char() const {
        if (is_positive) return '+';
        else return '-';
      }

    private:
      // Boolean value that is true when the parity is +1 and false when it is -1
      bool is_positive;
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
