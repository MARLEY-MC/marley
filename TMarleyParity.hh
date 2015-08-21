#pragma once
#include <iostream>
#include <string>

// Class that represents a parity value (either +1 or -1) in a type-safe way
class TMarleyParity {
  public:
    // Default constructor chooses positive parity
    inline TMarleyParity() {
      is_positive = true;
    }

    inline TMarleyParity(bool is_pos) {
      is_positive = is_pos;
    }

    inline TMarleyParity(int i) {
      if (i == 1) is_positive = true;
      else if (i == -1) is_positive = false;
      else throw std::runtime_error(std::string("Invalid parity ")
        + std::to_string(i) + " passed to constructor of TMarleyParity");
    }

    // Copy constructor
    inline TMarleyParity(const TMarleyParity& parity) {
      is_positive = parity.is_positive;
    }

    inline TMarleyParity& operator= (const TMarleyParity& parity) {
      // Do the copy
      is_positive = parity.is_positive;

      // Return the existing object
      return *this;
    }

    inline TMarleyParity& operator= (const int& i) {
      // Do the assignment
      if (i == 1) is_positive = true;
      else if (i == -1) is_positive = false;
      else throw std::runtime_error(std::string("Invalid parity ")
        + std::to_string(i) + " assigned to variable of type TMarleyParity");

      // Return the existing object
      return *this;
    }

    friend inline bool operator== (const TMarleyParity& p1,
      const TMarleyParity& p2)
    {
      return p1.is_positive == p2.is_positive;
    }

    friend inline bool operator!= (const TMarleyParity& p1,
      const TMarleyParity& p2)
    {
      return p1.is_positive != p2.is_positive;
    }

    friend inline TMarleyParity operator* (const TMarleyParity& p1,
      const TMarleyParity& p2)
    {
      return TMarleyParity(p1.is_positive == p2.is_positive);
    }

    friend inline int operator* (const TMarleyParity& p,
      const int& i)
    {
      if (p.is_positive) return i;
      else return -i;
    }

    friend inline int operator* (const int& i,
      const TMarleyParity& p)
    {
      if (p.is_positive) return i;
      else return -i;
    }

    friend inline std::ostream& operator<< (std::ostream& out,
      const TMarleyParity& p)
    {
      if (p.is_positive) out << "+";
      else out << "-";
      return out;
    }

    friend inline std::istream& operator>> (std::istream& in,
      TMarleyParity& p)
    {
      char c;
      in >> c;

      if (c == '+') p.is_positive = true;
      else if (c == '-') p.is_positive = false;
      else throw std::runtime_error(std::string("Invalid parity ")
        + c + " assigned via the >> operator to"
        + " a variable of type TMarleyParity");

      return in;
    }

    // Allows implicit or explicit casts of TMarleyParity to int
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
