#include <algorithm>
#include <cctype>
#include <chrono>
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
// cascade simulations. Seed it using the system time
unsigned marley_utils::seed = std::chrono::system_clock::now().time_since_epoch().count();
std::knuth_b marley_utils::rand_gen(seed);


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
