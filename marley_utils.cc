#include <algorithm>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
#include "marley_utils.hh"

// Strings to use for latex table output of ENSDF data
std::string marley_utils::latex_table_1 = "\\documentclass[12pt]{article}\n"
  "\n"
  "\\usepackage{amsmath}\n"
  "\\usepackage{booktabs}\n"
  "\\usepackage[justification=justified,\n"
  "width=\\textwidth,\n"
  "labelformat=empty]{caption}\n"
  "\\usepackage[top=1in, bottom=1in, left=0.25in, right=0.25in]{geometry}\n"
  "\\usepackage{isotope}\n"
  "\\usepackage{longtable}\n"
  "\\usepackage{multirow}\n"
  "\\usepackage{siunitx}\n"
  "\n"
  "\\newcommand{\\ExtraRowSpace}{1cm}\n"
  "\n"
  "\\begin{document}\n"
  "\\begin{center}\n"
  "\\begin{longtable}{\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "c\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "%  table-column-width = 2cm]\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "S[table-number-alignment = center,\n"
  "  table-text-alignment = center]\n"
  "}\n"
  "\\caption";

std::string marley_utils::latex_table_2 = "\\toprule\n"
  "%{\\centering\\textbf{Level Energy (keV)}}\n"
  "%& {\\centering\\textbf{Spin-Parity}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}\n"
  "%& {\\centering\\textbf{Final Energy (keV)}} \\\\\n"
  "{\\multirow{3}{2cm}{\\centering\\textbf{Level Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{Spin-Parity}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}}\n"
  "& {\\multirow{3}{2.5cm}{\\centering\\textbf{Final Level Energy (keV)}}} \\\\\n"
  "& & & & \\\\\n"
  "& & & & \\\\\n"
  "\\midrule\n"
  "\\endfirsthead\n"
  "\\caption[]";

std::string marley_utils::latex_table_3 = "\\toprule\n"
  "%{\\centering\\textbf{Level Energy (keV)}}\n"
  "%& {\\centering\\textbf{Spin-Parity}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}\n"
  "%& {\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}\n"
  "%& {\\centering\\textbf{Final Energy (keV)}} \\\\\n"
  "{\\multirow{3}{2cm}{\\centering\\textbf{Level Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{Spin-Parity}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ Energy (keV)}}}\n"
  "& {\\multirow{3}{2cm}{\\centering\\textbf{$\\boldsymbol{\\gamma}$ RI}}}\n"
  "& {\\multirow{3}{2.5cm}{\\centering\\textbf{Final Level Energy (keV)}}} \\\\\n"
  "& & & & \\\\\n"
  "& & & & \\\\\n"
  "\\midrule\n"
  "\\noalign{\\vspace{-\\ExtraRowSpace}}\n"
  "\\endhead\n"
  "\\bottomrule \\multicolumn{4}{r}{\\textit{Continued on next page}} \\\\\n"
  "\\endfoot\n"
  "\\bottomrule\n"
  "\\endlastfoot\n"
  "% Table data\n"
  "\n";

std::string marley_utils::latex_table_4 = "\\end{longtable}\n"
  "\\end{center}\n"
  "\\end{document}";

// Random number generator that will be used when selecting gammas for
// cascade simulations. Seed it using the system time.
// This is an attempt to do a decent job of seeding the random number
// generator, but optimally accomplishing this can be tricky (see, for
// example, http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
uint_fast64_t marley_utils::seed = std::chrono::system_clock::now().time_since_epoch().count();
std::seed_seq seed_sequence{marley_utils::seed};
std::mt19937_64 marley_utils::rand_gen(seed_sequence);

// This implementation of the complex gamma function is based on the
// Lanczos approximation and its Python implementation given
// on Wikipedia (https://en.wikipedia.org/wiki/Lanczos_approximation)
// The C++ version given here is taken almost verbatim from
// http://bytes.com/topic/c/answers/576697-c-routine-complex-gamma-function
std::complex<double> marley_utils::gamma(std::complex<double> z)
{
  // Initialize some constants used in the algorithm. The "static
  // const" keywords ensure that these constants are initialized only
  // once (not reinitialized each time this function is called)
  static const int g=7;
  static const double pi =
  3.1415926535897932384626433832795028841972;
  static const double p[g+2] = {0.99999999999980993, 676.5203681218851,
  -1259.1392167224028, 771.32342877765313, -176.61502916214059,
  12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
  1.5056327351493116e-7};

  if (std::real(z) < 0.5) {
    return pi / (std::sin(pi*z)*gamma(1.0-z));
  }

  z -= 1.0;

  std::complex<double> x = p[0];

  for (int i = 1; i < g + 2; i++) {
    x += p[i]/(z+std::complex<double>(i,0));
  }

  std::complex<double> t = z + (g + 0.5);

  return std::sqrt(2*pi) * std::pow(t, z + 0.5) * std::exp(-t) * x;
}

// This function is a modified version of a public-domain implementation of
// Brent's algorithm for minimizing a function. You can download the original
// source code from http://www.codeproject.com/Articles/30201/Optimizing-a-Function-of-One-Variable

// The return value of minimize is the minimum of the function f.
// The location where f takes its minimum is returned in the variable minLoc.
// Notation and implementation based on Chapter 5 of Richard Brent's book
// "Algorithms for Minimization Without Derivatives".
double marley_utils::minimize(const std::function<double(double)> f, // [in] objective function to minimize
  double leftEnd,     // [in] smaller value of bracketing interval
  double rightEnd,    // [in] larger value of bracketing interval
  double epsilon,     // [in] stopping tolerance
  double& minLoc)     // [out] location of minimum
{
    double d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
    static const double c = 0.5*(3.0 - std::sqrt(5.0));
    static const double SQRT_DBL_EPSILON = std::sqrt(DBL_EPSILON);
    
    double& a = leftEnd;
    double& b = rightEnd;
    double& x = minLoc;

    v = w = x = a + c*(b - a);
    d = e = 0.0;
    fv = fw = fx = f(x);

    // Check stopping criteria
    while (m = 0.5*(a + b),
      tol = SQRT_DBL_EPSILON*fabs(x) + epsilon,
      t2 = 2.0*tol,
      fabs(x - m) > t2 - 0.5*(b - a))
    {
        p = q = r = 0.0;
        if (fabs(e) > tol)
        {
            // fit parabola
            r = (x - w)*(fx - fv);
            q = (x - v)*(fx - fw);
            p = (x - v)*q - (x - w)*r;
            q = 2.0*(q - r);
            (q > 0.0) ? p = -p : q = -q;
            r = e; e = d;
        }
        if (fabs(p) < fabs(0.5*q*r) && p < q*(a - x) && p < q*(b - x))
        {
            // A parabolic interpolation step
            d = p/q;
            u = x + d;
            // f must not be evaluated too close to a or b
            if (u - a < t2 || b - u < t2)
                d = (x < m) ? tol : -tol;
        }
        else
        {
            // A golden section step
            e = (x < m) ? b : a;
            e -= x;
            d = c*e;
        }
        // f must not be evaluated too close to x
        if (fabs(d) >= tol)
            u = x + d;
        else if (d > 0.0)
            u = x + tol;
        else
            u = x - tol;
        fu = f(u);
        // Update a, b, v, w, and x
        if (fu <= fx)
        {
            (u < x) ? b = x : a = x;
            v = w; fv = fw; 
            w = x; fw = fx; 
            x = u; fx = fu;
        }
        else
        {
            (u < x) ? a = u : b = u;
            if (fu <= fw || w == x)
            {
                v = w; fv = fw; 
                w = u; fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u; fv = fu;
            }
        }
    }
    return  fx;
}

// We can maximize a function using the same technique by minimizing its opposite
double marley_utils::maximize(const std::function<double(double)> f, // [in] objective function to maximize
  double leftEnd,     // [in] smaller value of bracketing interval
  double rightEnd,    // [in] larger value of bracketing interval
  double epsilon,     // [in] stopping tolerance
  double& maxLoc)     // [out] location of maximum
{
  double result = minimize([&f](double x) -> double { return -1.0*f(x); },
    leftEnd, rightEnd, epsilon, maxLoc);
  return -1.0*result;
}

// Numerically integrate a given function f (that takes a
// double argument to integrate over and returns a double)
// over the interval [a,b] using the composite trapezoidal
// rule over n subintervals.
// (see http://en.wikipedia.org/wiki/Numerical_integration)
double marley_utils::num_integrate(const std::function<double(double)> &f,
  double a, double b, int n)
{
  double integral = 0;
  for(int k = 1; k < n-1; k++) {
    integral += ((b - a)/n)*f(a + k*(b - a)/n);
  }
  integral += ((b - a)/n)*(f(a)/2 + f(b)/2);
  return integral;
}

// Solves a quadratic equation of the form A*x^2 + B*x + C = 0
// while attempting to minimize errors due to the inherent limitations
// of floating-point arithmetic. The variables solPlus and solMinus
// are loaded with the two solutions x = (-B ± sqrt(B^2 - 4*A*C))/(2*A).
// Unsurprisingly, solPlus corresponds to the choice of the plus sign,
// while solMinus corresponds to the choice of the minus sign.
void marley_utils::solve_quadratic_equation(double A, double B,
  double C, double &solPlus, double &solMinus) {

  // Restructure the calculation to avoid some potentially bad cancellations
  // See, e.g., http://www.petebecker.com/js/js200010.html for details
  double c = C/A;
  double b = B/(2*A);

  double discr = b*b - c;

  // Find both solutions of the quadratic equation while avoiding
  // an extra subtraction (which can potentially lead to catastrophic
  // loss of precision) between -b and the square root of the
  // discriminant.
  if (b > 0) {

    solMinus = -b - marley_utils::real_sqrt(discr);
    solPlus = c/solMinus;
  }
  else {

    solPlus = -b + marley_utils::real_sqrt(discr);
    solMinus = c/solPlus;
    //std::cout << "DEBUG: b = " << b << std::endl;
    //std::cout << "DEBUG: c = " << c << std::endl;
    //std::cout << "DEBUG: b*b - c = " << b*b - c << std::endl;
    //std::cout << "DEBUG: solPlus = " << solPlus << std::endl;
    //std::cout << "DEBUG: solMinus = " << solMinus << std::endl;
  }

}

// Takes the square root of nonnegative arguments. Returns zero otherwise.
// This sqrt implementation is a quick fix for cases where roundoff errors
// give an argument to std::sqrt() that is slightly negative, causing it
// to return NaN.

// TODO: consider other ways of resolving this problem
double marley_utils::real_sqrt(double num) {
  if (num < 0) {
    return 0;
  }
  else {
    return std::sqrt(num);
  }
}

// For a given atomic number Z and mass number A, return a matching ENSDF nucid
std::string marley_utils::nuc_id(int Z, int A) {
  // Check to make sure Z and A have acceptable values
  if (Z < 0 || A < 1 || A > 999) throw std::runtime_error(
    std::string("The atomic number Z = ") + std::to_string(Z)
    + " and the mass number A = " + std::to_string(A)
    + " do not correspond to a valid ENSDF nucid.");

  // Create a three-character string representing the mass number
  std::string atomic_mass_number;
  if (A < 10) {
    atomic_mass_number = "  " + std::to_string(A);
  }
  else if (A < 100) {
    atomic_mass_number = " " + std::to_string(A);
  }
  else {
    atomic_mass_number = std::to_string(A);
  }

  // Get the element symbol as a string
  std::string symbol = element_symbols.at(Z);

  // Make the symbol completely uppercase
  to_uppercase_inplace(symbol);

  // If the symbol is only one character, pad the string so that
  // it is two characters long
  if (symbol.length() == 1) symbol += " ";

  return atomic_mass_number + symbol;
}

// Sample a random double uniformly between min and max using the default
// random number generator marley_utils::rand_gen. The inclusive flag
// determines whether or not max is included in the range. That is,
// when inclusive == false, the sampling is done on the interval [min, max),
// while inclusive == true uses [min, max].
double marley_utils::uniform_random_double(double min, double max, bool inclusive) {
  // Defaults to sampling from [0,1). We will always
  // explicitly supply the upper and lower bounds to
  // this distribution, so we won't worry about the
  // default setting.
  static std::uniform_real_distribution<double> udist;  

  double max_to_use;

  if (inclusive) { // sample from [min, max]

    // Find the double value that comes immediately after max. This allows
    // us to sample uniformly on [min, max] rather than [min,max). This trick comes from
    // http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution/uniform_real_distribution
    max_to_use = std::nextafter(max, std::numeric_limits<double>::max());
  }
  else { // sample from [min, max)
    max_to_use = max;
  }

  std::uniform_real_distribution<double>::param_type params(min, max_to_use);

  // Sample a random double from this distribution
  return udist(rand_gen, params);

}

// Sample from a given 1D probability density function f(x)
// on the interval [xmin, xmax] using a simple rejection method:
// (1) find the maximum of the function on [xmin, xmax]
// (2) sample an x value uniformly over the function f(x)'s domain
// (3) sample a y value uniformly over [0, max(f(x))]
// (4) if y <= f(x), accept the sampled x value
// (5) if y > f(x), reject the sampled x value, and return
// to step 2 to try again
// Note that f(x) does not need to be normalized, but its range
// must be nonnegative
double marley_utils::rejection_sample(std::function<double(double)> f,
  double xmin, double xmax, double max_search_tolerance)
{
  // This variable will be loaded with the value of x
  // that corresponds to the maximum of f(x).
  // We don't actually use this, but currently it's
  // a required parameter of marley_utils::maximize
  double x_at_max;

  // Get the maximum value of f(x). This is needed to
  // correctly apply rejection sampling.
  double fmax = marley_utils::maximize(f, xmin, xmax, max_search_tolerance, x_at_max);

  double x, y;

  do {
    // Sample x value uniformly from [xmin, xmax]
    x = marley_utils::uniform_random_double(xmin, xmax, true);
    // Sample y uniformly from [0, fmax]
    y = marley_utils::uniform_random_double(0, fmax, true);
  }
  // Keep sampling until you get a y value less than f(x) 
  // (the probability density function evaluated at the sampled value of x)
  while (y > f(x));

  return x;
}

// Function that takes a number of bytes and returns a string
// representing the amount of memory in more readable units
std::string marley_utils::num_bytes_to_string(double bytes,
  unsigned precision)
{
  double divisor;
  std::string SI_prefix;

  if (bytes < 1e3) return std::to_string(static_cast<int>(bytes)) + " B";
  else if (bytes < 1e6) {
    divisor = 1e3;
    SI_prefix = "k";
  }
  else if (bytes < 1e9) {
    divisor = 1e6;
    SI_prefix = "M";
  }
  else if (bytes < 1e12) {
    divisor = 1e9;
    SI_prefix = "G";
  }
  else {
    divisor = 1e12;
    SI_prefix = "T";
  }

  std::ostringstream out;
  out << std::fixed;
  out.precision(precision);
  out << bytes/divisor << " " << SI_prefix << "B";
  return out.str();
}

// This function exploits the observation (given in the first answer at
// http://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code)
// that the difference of two std::chrono::system_clock::time_point objects
// can be assigned to a std::chrono::system_clock::duration.
//
// Function that takes two std::system_clock::time_point objects and returns
// a string (in the format days hours:minutes:seconds) representing the time
// between them
std::string marley_utils::elapsed_time_string(
  std::chrono::system_clock::time_point &start_time,
  std::chrono::system_clock::time_point &end_time)
{
  std::chrono::system_clock::duration time_elapsed
    = end_time - start_time;

  return marley_utils::duration_to_string
    <std::chrono::system_clock::duration>(time_elapsed);
}

// Efficiently read in an entire file as a std::string
// This function was taken from
// http://insanecoding.blogspot.in/2011/11/how-to-read-in-file-in-c.html
std::string marley_utils::get_file_contents(std::string filename) {

  // Read from the file using a stream. Ignore windows line ending changes
  std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);

  if (in) {
    std::string contents;

    // Determine the size of the file. Allocate sufficient memory
    // so that the string can hold the entire file without
    // resizing itself
    in.seekg(0, std::ios::end);
    contents.resize(in.tellg());
    in.seekg(0, std::ios::beg);

    // Slurp in the entire file
    in.read(&contents[0], contents.size());
    in.close();
    return contents;
  }
  throw std::runtime_error("Could not read from file " + filename);
}

// Advance to the next line of an ifstream that either matches (match == true)
// or does not match (match == false) a given regular expression
std::string marley_utils::get_next_line(std::ifstream &file_in,
  const std::regex &rx, bool match) {

  std::string line;
  while (!file_in.eof() && file_in.good()) {
    // Get the next line of the file
    std::getline(file_in, line);

    // Check to see if the new line fulfills the search criteria
    if (std::regex_match(line, rx) == match) {
      // If it does, return it
      return line;
    }
    // If not, keep looking
  }

  // If the end of the file is encountered before a suitable
  // line is found, return an empty string
  return std::string("");
}

// Trim an ENSDF nucid string and make two-letter element symbols have a
// lowercase last letter. Currently, no checking is done to see if the
// string is a valid nucid.
std::string marley_utils::nucid_to_symbol(std::string nucid) {
  if (nucid.length() != 5) {
    // Split the string into "A" and "element name" pieces
    std::smatch m;
    std::regex_search(nucid, m, std::regex("[0-9]+"));
    std::string z_str = m.str();
    std::string e_str = m.suffix().str();

    // If the element name has more than one letter,
    // then make the last one lowercase.
    if (e_str.length() > 1)
      e_str.back() = tolower(e_str.back());
    return z_str + e_str;
  }
  // The nucid is 5 characters long, so getting our
  // desired format is a lot easier
  nucid.back() = tolower(nucid.back());
  return marley_utils::trim_copy(nucid);
}

// Converts an ENSDF nucid to an atomic number. Currently, no checking is done
// to see if the string is a valid nucid.
int marley_utils::nucid_to_Z(std::string nucid) {
  // String that will be loaded with the element symbol
  std::string e_str;

  if (nucid.length() != 5) {
    // Extract the element name from the string
    std::smatch m;
    std::regex_search(nucid, m, std::regex("[0-9]+"));
    e_str = m.suffix().str();

    // If the element name has more than one letter,
    // then make the last one lowercase.
    if (e_str.length() > 1)
      e_str.back() = tolower(e_str.back());
  }
  else {
    // The nucid is 5 characters long, so getting our
    // desired format is a lot easier
    nucid.back() = tolower(nucid.back());
    e_str = nucid.substr(nucid.size() - 2); // Get the last two characters
    marley_utils::trim_right_inplace(e_str); // Trims the string if needed
  }
  return atomic_numbers.at(e_str);
}
