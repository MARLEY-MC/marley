#pragma once // Nonstandard but widely-supported (see http://en.wikipedia.org/wiki/Pragma_once)
             // preprocessor directive that prevents this file from being included multiple times.
             // Consider changing to an include guard (http://en.wikipedia.org/wiki/Include_guard)
             // if desired by project collaborators.

#include <algorithm>
#include <chrono>
#include <complex>
#include <functional>
#include <unordered_map>
#include <random>
#include <regex>
#include <string>

namespace marley_utils {

  // Conversion factor to use when expressing ENSDF energies (keV) in
  // standard MARLEY energy units (MeV)
  const double MeV = 1e-3;

  // Conversion factor to use when expressing atomic masses (micro-amu)
  // in standard MARLEY energy units (MeV)
  const double micro_amu = 0.000931494061;

  // Consistent value of pi to use throughout all of MARLEY
  const double pi = std::acos(-1);

  // The physical constants given here were all taken from
  // the 2014 edition of the Review of Particle Physics
  // published by the Particle Data Group.
  
  // Fine structure constant 
  constexpr double alpha = 7.2973525698e-3;
  // Conversion factor used to switch to natural units (hbar = c = 1)
  constexpr double hbar_c = 197.3269718; // MeV*fm
  constexpr double hbar_c2 = std::pow(hbar_c, 2); // MeV^2 * fm^2
  // Electron mass
  constexpr double m_e = 0.510998928; // MeV
  // Constant to use when converting from mb to MeV^(-2)
  constexpr double mb = 1/3.89379338e5; // MeV^(-2) mb^(-1)
  // Square of the elementary charge
  constexpr double e2 = hbar_c * alpha; // MeV*fm

  // Strings to use for latex table output of ENSDF data
  extern std::string latex_table_1, latex_table_2, latex_table_3, latex_table_4;

  // Random number generator that will be used when selecting gammas for
  // cascade simulations.
  extern uint_fast64_t seed;
  extern std::mt19937_64 rand_gen; // Use the 64-bit Mersenne Twister RNG

  // Sample a random number uniformly on either [min, max) or [min, max]
  double uniform_random_double(double min, double max, bool inclusive);

  // Default tolerance for rejection sampling
  const double rejection_sampling_tolerance = 1e-8;

  // Sample a random number x from the pdf f(x) on the interval [xmin, xmax]
  double rejection_sample(std::function<double(double)> f,
    double xmin, double xmax, double max_search_tolerance
    = rejection_sampling_tolerance);

  // Create an ENSDF nucid string given a nuclide's atomic number Z
  // and mass number A
  std::string nuc_id(int Z, int A);

  // Return the PDG particle ID that corresponds to a ground-state
  // nucleus with atomic number Z and mass number A
  inline int get_nucleus_pid(int Z, int A) {
    return 10000*Z + 10*A + 1000000000;
  }

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

  // Efficiently read in an entire file as a std::string
  std::string get_file_contents(std::string filename);

  // Advance to the next line of an ifstream that either matches (match == true)
  // or does not match (match == false) a given regular expression
  std::string get_next_line(std::ifstream &file_in, const std::regex &rx,
    bool match);

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

  // Function that converts a std::string object to
  // all lowercase in place and returns a reference to
  // it afterwards
  inline std::string& to_lowercase_inplace(std::string& s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
  }

  // Function that converts a std::string object to
  // all uppercase in place and returns a reference to
  // it afterwards
  inline std::string& to_uppercase_inplace(std::string& s) {
    std::transform(s.begin(), s.end(), s.begin(), ::toupper);
    return s;
  }

  // Functions for padding std::string objects in place.  They all return
  // references to the string afterwards.  These functions are based on
  // http://stackoverflow.com/a/667219/4081973
  inline std::string& pad_left_inplace(std::string &str,
    const size_t len, const char pad_char = ' ')
  {
    if(len > str.size())
    str.insert(0, len - str.size(), pad_char);
    return str;
  }

  inline std::string& pad_right_inplace(std::string &str,
    const size_t len, const char pad_char = ' ')
  {
    if(len > str.size())
    str.append(len - str.size(), pad_char);
    return str;
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
  
  // The second three alter the original string, returning a
  // reference to it after it has been trimmed.
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

  // Function that takes a number of bytes and returns a string
  // representing the amount of memory in more readable units
  std::string num_bytes_to_string(double bytes, unsigned precision = 3);

  // Trim an ENSDF nucid string and make two-letter element symbols have a
  // lowercase last letter. Currently, no checking is done to see if the
  // string is a valid nucid.
  std::string nucid_to_symbol(std::string nucid);

  // Similar to nucid_to_symbol, but returns the atomic number as an integer
  // instead
  int nucid_to_Z(std::string nucid);

  // Generalized std::chrono::duration helper types
  template <typename repType> using
    seconds = std::chrono::duration<repType>;
  template <typename repType> using
    minutes = std::chrono::duration<repType, std::ratio<60>>;
  template <typename repType> using
    hours = std::chrono::duration<repType, std::ratio<3600>>;
  template <typename repType> using
    days = std::chrono::duration<repType, std::ratio<86400>>;

  // This function is a generalized version of code taken from the accepted answer at
  // http://stackoverflow.com/questions/15957805/extract-year-month-day-etc-from-stdchronotime-point-in-c
  //
  // Function that returns a string representation (in the
  // format days hours:minutes:seconds) of a std::chrono::duration object
  template <typename repType, typename periodType = std::ratio<1>> std::string duration_to_string(
    std::chrono::duration<repType, periodType> duration)
  {
    int day_count = static_cast<int>(std::chrono::duration_cast
      <marley_utils::days<repType>>(duration) / (marley_utils::days<repType>(1)));
    duration -= marley_utils::days<repType>(day_count);

    int hour_count = static_cast<int>(std::chrono::duration_cast
      <marley_utils::hours<repType>>(duration) / (marley_utils::hours<repType>(1)));
    duration -= marley_utils::hours<repType>(hour_count);

    int minute_count = static_cast<int>(std::chrono::duration_cast
      <marley_utils::minutes<repType>>(duration) / (marley_utils::minutes<repType>(1)));
    duration -= marley_utils::minutes<repType>(minute_count);

    int second_count = static_cast<int>(std::chrono::duration_cast
      <marley_utils::seconds<repType>>(duration) / (marley_utils::seconds<repType>(1)));
    duration -= marley_utils::seconds<repType>(second_count);

    std::ostringstream out;
    if (day_count > 1) {
      out << day_count << " days ";
    }
    if (day_count == 1) {
      out << day_count << " day ";
    }
    if (hour_count < 10) out << "0";
    out << hour_count << ":";
    if (minute_count < 10) out << "0";
    out << minute_count << ":";
    if (second_count < 10) out << "0";
    out << second_count;

    return out.str();
  }

  template <typename durationType> std::string duration_to_string(
    durationType duration)
  {
    return duration_to_string<typename durationType::rep,
      typename durationType::period>(duration);
  }

  // Function that takes two std::system_clock::time_point objects and returns
  // a string (in the same format as marley_utils::duration_to_string)
  // representing the time between them
  std::string elapsed_time_string(
    std::chrono::system_clock::time_point &start_time,
    std::chrono::system_clock::time_point &end_time);

  // Lookup table for element symbols (keys are atomic numbers Z,
  // values are symbols on the periodic table). The symbol "Nn" is
  // used for a neutron to match the ENSDF convention.
  const std::unordered_map<int, std::string> element_symbols = {
    { 0, "Nn"},
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

  // TODO: consider other ways of adding support for fast reverse lookups
  // in the element_symbols map rather than reproducing its contents here
  // Lookup table for atomic numbers (keys are symbols on the periodic table,
  // values are atomic numbers Z). The symbol "Nn" is used for a neutron to
  // match the ENSDF convention.
  const std::unordered_map<std::string, int> atomic_numbers = {
    {"Nn", 0    },
    {"H" , 1    },
    {"He", 2    },
    {"Li", 3    },
    {"Be", 4    },
    {"B" , 5    },
    {"C" , 6    },
    {"N" , 7    },
    {"O" , 8    },
    {"F" , 9    },
    {"Ne", 10   },
    {"Na", 11   },
    {"Mg", 12   },
    {"Al", 13   },
    {"Si", 14   },
    {"P" , 15   },
    {"S" , 16   },
    {"Cl", 17   },
    {"Ar", 18   },
    {"K" , 19   },
    {"Ca", 20   },
    {"Sc", 21   },
    {"Ti", 22   },
    {"V" , 23   },
    {"Cr", 24   },
    {"Mn", 25   },
    {"Fe", 26   },
    {"Co", 27   },
    {"Ni", 28   },
    {"Cu", 29   },
    {"Zn", 30   },
    {"Ga", 31   },
    {"Ge", 32   },
    {"As", 33   },
    {"Se", 34   },
    {"Br", 35   },
    {"Kr", 36   },
    {"Rb", 37   },
    {"Sr", 38   },
    {"Y" , 39   },
    {"Zr", 40   },
    {"Nb", 41   },
    {"Mo", 42   },
    {"Tc", 43   },
    {"Ru", 44   },
    {"Rh", 45   },
    {"Pd", 46   },
    {"Ag", 47   },
    {"Cd", 48   },
    {"In", 49   },
    {"Sn", 50   },
    {"Sb", 51   },
    {"Te", 52   },
    {"I" , 53   },
    {"Xe", 54   },
    {"Cs", 55   },
    {"Ba", 56   },
    {"La", 57   },
    {"Ce", 58   },
    {"Pr", 59   },
    {"Nd", 60   },
    {"Pm", 61   },
    {"Sm", 62   },
    {"Eu", 63   },
    {"Gd", 64   },
    {"Tb", 65   },
    {"Dy", 66   },
    {"Ho", 67   },
    {"Er", 68   },
    {"Tm", 69   },
    {"Yb", 70   },
    {"Lu", 71   },
    {"Hf", 72   },
    {"Ta", 73   },
    {"W" , 74   },
    {"Re", 75   },
    {"Os", 76   },
    {"Ir", 77   },
    {"Pt", 78   },
    {"Au", 79   },
    {"Hg", 80   },
    {"Tl", 81   },
    {"Pb", 82   },
    {"Bi", 83   },
    {"Po", 84   },
    {"At", 85   },
    {"Rn", 86   },
    {"Fr", 87   },
    {"Ra", 88   },
    {"Ac", 89   },
    {"Th", 90   },
    {"Pa", 91   },
    {"U" , 92   },
    {"Np", 93   },
    {"Pu", 94   },
    {"Am", 95   },
    {"Cm", 96   },
    {"Bk", 97   },
    {"Cf", 98   },
    {"Es", 99   },
    {"Fm", 100  },
    {"Md", 101  },
    {"No", 102  },
    {"Lr", 103  },
    {"Rf", 104  },
    {"Db", 105  },
    {"Sg", 106  },
    {"Bh", 107  },
    {"Hs", 108  },
    {"Mt", 109  },
    {"Ds", 110  },
    {"Rg", 111  },
    {"Cn", 112  },
    //{ "Uut", 113 },
    {"Fl", 114  },
    //{ "Uup", 115 },
    {"Lv", 116  },
    //{ "Uus", 117 },
    //{ "Uuo", 118 },
  };

  const std::string marley_logo = "╔╦╗╔═╗╦═╗╦  ╔═╗╦ ╦\n"
                                  "║║║╠═╣╠╦╝║  ║╣ ╚╦╝\n"
                                  "╩ ╩╩ ╩╩╚═╩═╝╚═╝ ╩\n";

  const std::string marley_pic =
    "ZZ77I\?I7777\?\?+I\?\?+=====~=~~~:+=+7I$I$Z$$OZ77Z$O\?\?=\?\?$Z$Z=$=~=~,:::~~:=~===+=+=II\n"
    "Z$7IIIIII+I=++=+==+~=~~~~~===\?ZIZ$OOZZOOI8$ZZO8Z8O$ZOO7Z+:,:::,::~:::::~+~:~=~=+\n"
    "$77I$II\?+I\?+++=+=+====+\?=~=I$7O7$OOZ$8OO8DOO8888O8ZOZZ8Z7I~:::::::::::~~:~==~==~\n"
    "I7I$I7++\?+++=+=++===~~++I\?OO$I$ZOZDO8O888D888OO8888888ZOO$\?I::~,::,::::~:~~=~~==\n"
    "IIII+\?\?I\?=\?+=+=====~~=$$Z$$Z8IOZOZOD888D$=:::,,,,,:Z88888ZZ$=\?:,~::::~~:~,:~~~:~\n"
    "77I$\?\?\?+\?+\?=+=~=~~~=7Z87ZZ$OIZZOO8888DD\?::::,,::::,,::O88O8ZZ=:=:::,::::::::=~:=\n"
    "II+\?+\?\?++++======++7Z8I7OO$$DO$8Z888DZ\?~~:,,,,,,,,,:,:~:$O8ZOZ\?=~,:,,:,:::::~~~:\n"
    "$+I\?+\?+=+===+=+\?+IZ8$DDDIZO8OOD8D88D87+:::=7ZI7~,,,,:::::~8O8ZZ7+~::,:,::,:,:::~\n"
    "7+\?+\?+\?+======\?$OO$IOD8OD88O8OD7D8DD8I:~77$7ZZ$$O+::,:::~~=8O78OI\?~:,:,::,,,::~:\n"
    "I\?+++=\?+=~++$$ZZODZ88ZDDOZZDOZD8DDDD7~~Z\?=~==~Z88O+~,::::~~=88OO8$+~=:,,::,,,:,:\n"
    "+++++==++I$$ZO8$$O8DNDNDOZ8O88O8DDD$+\?$$8O$888O8O$7:~~:~,,~~O8888OO$~,:,,:::,:,:\n"
    "+\?++====\?\?$O88Z88ODD$ZOZ$$OOOOODDD\?~~::~=$Z~Z7887I~::$Z$O\?::I888ZO8$\?~:,,::~::,,\n"
    "\?+++====\?$88DDZ$88OZDDD8$O88O88D8\?:::,:,,:+8OO+7$~,,,O8OOOZ==8888888O\?=~::::::,,\n"
    "++=====+\?OODZ$$DD8DDD8ZZ$$ODO8DD=::~,,:::,,,,:,,,,:,,:OI~:ZOZ88D8$8OO7+:,::::~::\n"
    "+++=+=\?\?8D87DZ8DDDDZ888O88OO88ZD\?\?\?:::~:::~=I::,,,:~,,88D:,7Z8D88O88O+=:,,,:~::,\n"
    "+=+=++78Z7OZZ8DDDD8DDD8$+78$8D8D7I~:~I\?II7$I\?ZI$DD8I\?+==~88II8DD8O88OZ~:,,:::,:,\n"
    "++===+ID$ODDDDNDDDDD7=$O\?$88NDD8I\?++I\?+77\?~~~:=\?IOZ$OO7O8787$8D88$88ZO=~::,:::::\n"
    "=+===+I$8D8=8NDDND=OZZD$\?7ODO8D87=\?++II+Z++\?I\?Z7$Z$I$D8~::~Z7DDD87O$OZ++,,,::::,\n"
    "===~+\?8DZ8N8DNDNNDD+I=77I7D8DDDOZ++==\?I8\?\?$8ZO$7=\?$\?Z$==::,:8DDD8888O$$~,,:::::,\n"
    "====$8OOODN\?DDNDNNDI=+ZIIZD8DDD$7+I+==IZ8Z~~\?~O8$7\?=ZI+$:::~DDDD8O88$87~,::::::,\n"
    "==+8$8O88N8ODNNNDNDND7Z=ZDD8NDNO+7+=\?+\?~Z\?,:,,+,+DDZ7Z+7~~,~DDDD888O8$$\?:,,,:::,\n"
    "+~=+ZO88DD8DNNDNDDND78DI\?O8DDDDZ\?7I=\?==~~ZDDDDO,:,=D87+I==~=DDDD8DDOO$I=::,::::,\n"
    "+=+\?8I8DDNO8NNNNDNNDDD\?IDON8DDD$7$I\?7=~=~~~:~8DDD.~888O7I\?\?8DDD88OD8O$~:::,,:,:,\n"
    "=~\?8$O8DDDZDNDNNDNDNDD78DDNNDDDZI+\?\?7~=+:==77~~:ZO,DOOZZ7$ZDDDD8D8D888$::::,::,,\n"
    "+=\?8IDDNDDODNNNNNNNN$D$DDDNDDDDO$=7=+=~\?+==\?+7Z7I7$I\?Z$$ODDDDD88DDDD88I\?,,::,:,:\n"
    "=\?$OI8NDND8DNNDNNNDD7788DNDD8DDNZ$I+\?~~~==+=$$$II$I~ZIOO8DDDDD8D8D8DDO7=::,:,,:,\n"
    "I$8888DDNDDDNN8NNNNN7Z88NNDDNNDDD8ZO7I=~~:~\?=77II\?\?7O$ZODDDDD8DDD8O88Z\?:,:::,::,\n"
    "$DDZ8DDNND8DNDNNNNND7$8ODN8DNND888D8$$ZI$$\?I~=+I\?+7Z7$ODDDNDD88888OODO:::,,,,:::\n"
    "8DDOD8DDNN$ND7NNNNND8I88NN8DDNDDDZDD$ZZIZ88ZO7\?==78$IZDDDDDDDOD888DD8Z,::,,,,::,\n"
    "DD8DDDDDDNNNO\?NNNNN\?878NNNZOZDDD8ZDNDD888O8DOD87ZOI\?O8DDDDDDDODO$D7OO~:::,:,::::\n"
    "8NZOO8D8NDNNN7D88NNIOID8ND88$NDD8O$DDDDDN88Z8D8DD8Z7D8DDDDDDDOD8\?8OO=::,:::,:,,,\n"
    "DOZ8ZD8DDDNDDDN$8ONDNI88NNZD7DND$$8DNDNND8DNOO$NODODDDNDDDNDDD887IOO=:::::::,,,:\n"
    "NID8ND8DDNN+DNN$DID8D=ODDNZD\?8DN8ODZDDDN8NDOD8DDD8DDDDDDDDNDDD8~IDOO~:::,:,:,:,:\n"
    "IDODDDNDNDNZNNNI77D8+IZDN888=O8DD8DDNDNNNN8ODDDDN8DDDDDDDDDZDZI77+OO$I~:::::,:,:\n"
    "7NODDDD8NNDDN\?I7$INO++ODDI8D+7O8$88NNNDDNNDDDD8D8O$DDDDDND88O$$887$\?OZZ:::,,::,,\n"
    "NNDODDNDDNDNO$7I\?+NZ\?7DD$88=+\?ZZ$OIDDNNDND8DODND8$=DDDNDD8OOI$O87$O8Z8\?+:~::,,::\n"
    "NNNDZD$DN+DNOZZ\?==D7=IDN$DNI=I$ZZZ7$Z$DDNN8D8O$8$::DDNN8D8D8I8O8Z88DD8I\?+:::,,::\n"
    "NNNZ$DDNNZ7NOII+==D8+ZO=8D\?++I$$Z$88O8OZ8OZOIZ8$::NNDNDOD88$$ZIIIZOD88I8O+:::,,,\n"
    "NNNZ8DDND8ON7ZI==~DOZ88OO7$I\?\?7+I$\?$Z8$$Z$ZZZD+::ODNZOD8=7I$D8Z$~I8DDD888+~~,,,,\n"
    "NNN$DD8NDDNN+7I==~7\?+7ZZ8I7I\?7\?77I7+ZZZ$$7ZNI~::~88O8DD+O$$DODZZZZDDODOO+:$\?:~,,\n"
    "NNDD8DNN8ZDN8I\?~~=~Z\?$I8D$7I=\?ZZ77\?\?$II$ZDZ=~::8D8ODZINZO$ODDI78=8DD8D+\?Z$$Z\?::,\n"
    "NN\?\?NNNDD8DD$\?7=NN=OZD8ODZ\?I=7ZI$I7O$$OO$=DDI:~\?7OOIDD8\?$DDDO\?Z\?OZ87I8~OZ888Z+~:\n"
    "NNIO8NND88NN8IIND87ZDDO8$I\?7$ZIII7\?\?Z8:=\?+OZD:IIDIODND7\?O$NDOOI7::O$\?77DZO$D8ZI~\n"
    "NNI7DDNNDONNO\?O$ONDNI\?DN$+\?ZIZ$7I778=++==++ID\?Z$O$ZDD\?\?7$DDO\?$$O~:8+~OO8$\?888OZ=\n"
    "NNO78NNN78NNDDZ7DN$8DZ78$I7I+$I7\?\?Z~+++=++\?II$OOIIZDDD87ODD\?ZO8Z:8Z$7ODDZ78D8O$\?\n"
    "NN=ZONNDDODD77$ZN7D8D=NOO$$I$7I77$8+\?++\?=+\?$7I\?DO7DDDD87ND8Z8$8~:=OO7Z$OI+8DDZZ\?\n"
    "NDOO8NND$ON87NDD\?7\?$D7N8Z7II8ZI+Z8~\?==+==I\?DZID$=ZDD8ZIZNODOZ8D~:~ZD$\?\?\?D7O8DOI\?\n"
    "NNIODNNNDNNI$NO\?+DDNDD8ZZ$78OO$ID~\?+\?=++=IZ$OO8OI$N$$Z8NDZZO$D::~=$OZI+=$$$+O8Z+\n"
    "NI\?O8NNNDD\?ODD8=$$887DOZZ$I7\?\?\?$~+\?=+~=\?\?~O\?ZZN\?I7\?OZ8NDOOIDZI:~==+$OI~~~ZI~78ZI\n"
    "NZ7DNNNNN8DNNOD\?78N8N8$\?7I\?\?++7==+\?\?~~=++Z~OZD7Z=DO7NND7$78DO=~====O$+:~=OZ=\?$O+\n"
    "N8OZDNNNNZNMN7$\?8OO\?OZ\?7$7III7$+++=====++=8+$ZI7\?$7$8D$D\?ODD8O==~=IO8\?\?+=D$I=7I=\n"
    "ND88ODN8NNNNN7\?Z7ZN=D8I$\?\?\?I+7~\?\?=++===+Z8\?+I$\?O\?DD8+8\?DOO8DN7I+==+=Z8=Z~$Z7+\?7\?\n";
}
