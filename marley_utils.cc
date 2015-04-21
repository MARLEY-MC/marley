#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
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
