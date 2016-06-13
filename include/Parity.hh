#pragma once
#include <iostream>
#include <string>

#include "Error.hh"

// Forward-declare the class and friend operators so that we can use the
// operators in the global scope
namespace marley { class Parity; }

inline marley::Parity& operator!(marley::Parity& parity);
inline bool operator== (const marley::Parity& p1, const marley::Parity& p2);
inline bool operator!= (const marley::Parity& p1, const marley::Parity& p2);
inline int operator* (const marley::Parity& p, const int& i);
inline int operator* (const int& i, const marley::Parity& p);
inline std::ostream& operator<< (std::ostream& out, const marley::Parity& p);
inline std::istream& operator>> (std::istream& in, marley::Parity& p);

namespace marley {

  // Class that represents a parity value (either +1 or -1) in a type-safe way
  class Parity {
    public:
      // Default constructor chooses positive parity
      inline Parity() {
        is_positive = true;
      }

      inline Parity(bool is_pos) {
        is_positive = is_pos;
      }

      inline Parity(int i) {
        if (i == 1) is_positive = true;
        else if (i == -1) is_positive = false;
        else throw marley::Error(std::string("Invalid parity ")
          + std::to_string(i) + " passed to constructor of marley::Parity");
      }

      // Copy constructor
      inline Parity(const marley::Parity& parity) {
        is_positive = parity.is_positive;
      }

      // The minus unary operator creates a copy of a parity
      // object with a flipped value
      inline marley::Parity operator-(const marley::Parity& parity) {
        return marley::Parity(!parity.is_positive);
      }

      // The not unary operator flips a parity value in place
      friend inline Parity& ::operator!(marley::Parity& parity) {
        parity.is_positive = !parity.is_positive;
        return parity;
      }

      inline marley::Parity& operator= (const marley::Parity& parity) {
        // Do the copy
        is_positive = parity.is_positive;

        // Return the existing object
        return *this;
      }

      inline marley::Parity& operator= (const int& i) {
        // Do the assignment
        if (i == 1) is_positive = true;
        else if (i == -1) is_positive = false;
        else throw marley::Error(std::string("Invalid parity ")
          + std::to_string(i) + " assigned to variable of type marley::Parity");

        // Return the existing object
        return *this;
      }

      friend inline bool ::operator== (const marley::Parity& p1,
        const marley::Parity& p2)
      {
        return p1.is_positive == p2.is_positive;
      }

      friend inline bool ::operator!= (const marley::Parity& p1,
        const marley::Parity& p2)
      {
        return p1.is_positive != p2.is_positive;
      }

      friend inline int ::operator* (const marley::Parity& p,
        const int& i)
      {
        if (p.is_positive) return i;
        else return -i;
      }

      friend inline int ::operator* (const int& i,
        const marley::Parity& p)
      {
        if (p.is_positive) return i;
        else return -i;
      }

      friend inline std::ostream& ::operator<< (std::ostream& out,
        const marley::Parity& p)
      {
        if (p.is_positive) out << "+";
        else out << "-";
        return out;
      }

      friend inline std::istream& ::operator>> (std::istream& in,
        marley::Parity& p)
      {
        char c;
        in >> c;

        if (c == '+') p.is_positive = true;
        else if (c == '-') p.is_positive = false;
        else throw marley::Error(std::string("Invalid parity ")
          + c + " assigned via the >> operator to"
          + " a variable of type marley::Parity");

        return in;
      }

      // Allows implicit or explicit casts of marley::Parity to int
      inline operator int() const {
        if (is_positive) return 1;
        else return -1;
      }

      // Converts the parity to a string
      inline std::string str() const {
        if (is_positive) return std::string("+");
        else return std::string("-");
      }

    private:
      // Boolean value that is true when the parity is +1 and false when it is -1
      bool is_positive;
  };

}

inline marley::Parity operator* (const marley::Parity& p1, const marley::Parity& p2)
{
  return marley::Parity(p1 == p2);
}
