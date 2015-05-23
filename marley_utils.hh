#pragma once // Nonstandard but widely-supported (see http://en.wikipedia.org/wiki/Pragma_once)
             // preprocessor directive that prevents this file from being included multiple times.
             // Consider changing to an include guard (http://en.wikipedia.org/wiki/Include_guard)
             // if desired by project collaborators.

#include <algorithm>
#include <complex>
#include <functional>
#include <map>
#include <random>
#include <string>

namespace marley_utils {

  // Conversion factor to use when expressing ENSDF energies (keV) in
  // standard MARLEY energy units (MeV)
  const double MeV = 1e-3;

  // Consistent value of pi to use throughout all of MARLEY
  const double pi = std::acos(-1);

  // The physical constants given here were all taken from
  // the 2014 edition of the Review of Particle Physics
  // published by the Particle Data Group.
  
  // Fine structure constant 
  const double alpha = 7.2973525698e-3;
  // Conversion factor used to switch to natural units (hbar = c = 1)
  const double hbar_c = 197.3269718; // MeV*fm
  // Electron mass
  const double m_e = 0.510998928; // MeV

  // Strings to use for latex table output of ENSDF data
  extern std::string latex_table_1, latex_table_2, latex_table_3, latex_table_4;

  // Random number generator that will be used when selecting gammas for
  // cascade simulations.
  extern unsigned seed;
  extern std::knuth_b rand_gen;

  // Sample a random number uniformly on either [min, max) or [min, max]
  double uniform_random_double(double min, double max, bool inclusive);

  // Create an ENSDF nucid string given a nuclide's atomic number Z
  // and mass number A
  std::string nuc_id(int Z, int A);

  // Take the square root of a number. Assume that a negative argument is
  // due to roundoff error and return zero in such cases rather than NaN.
  double real_sqrt(double num);

  // Compute the complex gamma function using the Lanczos approximation
  std::complex<double> gamma(std::complex<double> z);

  // Numerically integrate a 1D function using the composite trapezoid rule
  double num_integrate(const std::function<double(double)> &f,
    double a, double b, int n);

  // Numerically minimize or maximize a function of one variable using
  // Brent's method (see http://en.wikipedia.org/wiki/Brent%27s_method)
  double minimize(const std::function<double(double)> f, double leftEnd,
    double rightEnd, double epsilon, double& minLoc);

  double maximize(const std::function<double(double)> f, double leftEnd,
    double rightEnd, double epsilon, double& maxLoc);

  // Find both solutions of a quadratic equation while attempting
  // to avoid floating-point arithmetic issues
  void solve_quadratic_equation(double A, double B,
    double C, double &solPlus, double &solMinus);

  // String containing all of the characters that will be
  // considered whitespace by default in the string
  // manipulation functions below 
  const std::string whitespace = " \f\n\r\t\v"; 

  // This version of std::stod will return 0 if it encounters
  // an empty string or an all-whitespace string.
  inline double str_to_double(const std::string& s) {
    size_t endpos = s.find_last_not_of(whitespace);
    if (endpos == std::string::npos) {
      return 0.0; // string was all whitespace
    }
    else {
      return std::stod(s);
    }
  }

  // Function that creates a copy of a std::string object
  // that has been converted to all lowercase
  inline std::string to_lowercase(const std::string& s) {
    std::string new_s = s;
    std::transform(new_s.begin(), new_s.end(), new_s.begin(), ::tolower);
    return new_s;
  }

  // These std::string trimming functions were taken from code
  // presented here: http://www.cplusplus.com/faq/sequences/strings/trim/
  // and here: http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
  // The first three (with the suffix copy) return a trimmed copy of
  // the string without modifying the original.
  inline std::string trim_right_copy(const std::string& s,
    const std::string& delimiters = whitespace)
  {
    size_t endpos = s.find_last_not_of(delimiters);
    return (endpos == std::string::npos) ? "" : s.substr(0, endpos + 1);
  }
  
  inline std::string trim_left_copy(const std::string& s,
    const std::string& delimiters = whitespace)
  {
    size_t startpos = s.find_first_not_of(delimiters);
    return (startpos == std::string::npos) ? "" : s.substr(startpos);
  }
  
  inline std::string trim_copy(const std::string& s,
    const std::string& delimiters = whitespace)
  {
    return trim_left_copy(trim_right_copy(s, delimiters), delimiters);
  }
  
  // The second three alter the original string, returning it
  // after it has been trimmed. 
  inline std::string& trim_right_inplace(std::string& s,
    const std::string& delimiters = whitespace)
  {
    size_t endpos = s.find_last_not_of(delimiters);
    if (endpos == std::string::npos) {
      s.clear();
    }
    else {
      s.erase(endpos + 1);
    }
    return s; 
  }
  
  inline std::string& trim_left_inplace(std::string& s,
    const std::string& delimiters = whitespace)
  {
    size_t startpos = s.find_first_not_of(delimiters);
    if (startpos == std::string::npos) {
      s.clear();
    }
    else {
      s.erase(0, startpos);
    }
    return s; 
  }
  
  inline std::string& trim_inplace(std::string& s,
    const std::string& delimiters = whitespace)
  {
    return trim_left_inplace(trim_right_inplace(s,delimiters), delimiters);
  }

  const std::map<int, std::string> element_symbols = {
    { 1, "H" },
    { 2, "He" },
    { 3, "Li" },
    { 4, "Be" },
    { 5, "B" },
    { 6, "C" },
    { 7, "N" },
    { 8, "O" },
    { 9, "F" },
    { 10, "Ne" },
    { 11, "Na" },
    { 12, "Mg" },
    { 13, "Al" },
    { 14, "Si" },
    { 15, "P" },
    { 16, "S" },
    { 17, "Cl" },
    { 18, "Ar" },
    { 19, "K" },
    { 20, "Ca" },
    { 21, "Sc" },
    { 22, "Ti" },
    { 23, "V" },
    { 24, "Cr" },
    { 25, "Mn" },
    { 26, "Fe" },
    { 27, "Co" },
    { 28, "Ni" },
    { 29, "Cu" },
    { 30, "Zn" },
    { 31, "Ga" },
    { 32, "Ge" },
    { 33, "As" },
    { 34, "Se" },
    { 35, "Br" },
    { 36, "Kr" },
    { 37, "Rb" },
    { 38, "Sr" },
    { 39, "Y" },
    { 40, "Zr" },
    { 41, "Nb" },
    { 42, "Mo" },
    { 43, "Tc" },
    { 44, "Ru" },
    { 45, "Rh" },
    { 46, "Pd" },
    { 47, "Ag" },
    { 48, "Cd" },
    { 49, "In" },
    { 50, "Sn" },
    { 51, "Sb" },
    { 52, "Te" },
    { 53, "I" },
    { 54, "Xe" },
    { 55, "Cs" },
    { 56, "Ba" },
    { 57, "La" },
    { 58, "Ce" },
    { 59, "Pr" },
    { 60, "Nd" },
    { 61, "Pm" },
    { 62, "Sm" },
    { 63, "Eu" },
    { 64, "Gd" },
    { 65, "Tb" },
    { 66, "Dy" },
    { 67, "Ho" },
    { 68, "Er" },
    { 69, "Tm" },
    { 70, "Yb" },
    { 71, "Lu" },
    { 72, "Hf" },
    { 73, "Ta" },
    { 74, "W" },
    { 75, "Re" },
    { 76, "Os" },
    { 77, "Ir" },
    { 78, "Pt" },
    { 79, "Au" },
    { 80, "Hg" },
    { 81, "Tl" },
    { 82, "Pb" },
    { 83, "Bi" },
    { 84, "Po" },
    { 85, "At" },
    { 86, "Rn" },
    { 87, "Fr" },
    { 88, "Ra" },
    { 89, "Ac" },
    { 90, "Th" },
    { 91, "Pa" },
    { 92, "U" },
    { 93, "Np" },
    { 94, "Pu" },
    { 95, "Am" },
    { 96, "Cm" },
    { 97, "Bk" },
    { 98, "Cf" },
    { 99, "Es" },
    { 100, "Fm" },
    { 101, "Md" },
    { 102, "No" },
    { 103, "Lr" },
    { 104, "Rf" },
    { 105, "Db" },
    { 106, "Sg" },
    { 107, "Bh" },
    { 108, "Hs" },
    { 109, "Mt" },
    { 110, "Ds" },
    { 111, "Rg" },
    { 112, "Cn" },
    //{ 113, "Uut" },
    { 114, "Fl" },
    //{ 115, "Uup" },
    { 116, "Lv" },
    //{ 117, "Uus" },
    //{ 118, "Uuo" },
  };
}
