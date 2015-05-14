#pragma once // Nonstandard but widely-supported (see http://en.wikipedia.org/wiki/Pragma_once)
             // preprocessor directive that prevents this file from being included multiple times.
             // Consider changing to an include guard (http://en.wikipedia.org/wiki/Include_guard)
             // if desired by project collaborators.

#include <algorithm>
#include <complex>
#include <functional>
#include <random>
#include <string>

namespace marley_utils {

  // Strings to use for latex table output of ENSDF data
  extern std::string latex_table_1, latex_table_2, latex_table_3, latex_table_4;

  // Random number generator that will be used when selecting gammas for
  // cascade simulations.
  extern unsigned seed;
  extern std::knuth_b rand_gen;

  // Compute the complex gamma function using the Lanczos approximation
  std::complex<double> gamma(std::complex<double> z);

  // Numerically minimize or maximize a function of one variable using
  // Brent's method (see http://en.wikipedia.org/wiki/Brent%27s_method)
  double minimize(std::function<double(double)> f, double leftEnd,
    double rightEnd, double epsilon, double& minLoc);

  double maximize(std::function<double(double)> f, double leftEnd,
    double rightEnd, double epsilon, double& maxLoc);

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

}
