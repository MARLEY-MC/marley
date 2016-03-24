#include <algorithm>
#include <climits>
#include <functional>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "marley_utils.hh"
#include "Generator.hh"
#include "NuclearPhysics.hh"

#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"

int main() {

  int Zi = 19;
  int Ai = 40;
  int twoJi = 2;
  marley::Parity Pi = 1;

  double strength_threshold = 0.2;

  //std::cout << std::setprecision(16) << std::scientific;

  marley::Generator gen("config.txt");

  marley::Reaction r = gen.get_reactions().back();

  std::unordered_map<const marley::Fragment*,
    std::function<double(double)> > funcs;

  std::unordered_map<const marley::Fragment*, double> total_c_widths;
  std::unordered_map<const marley::Fragment*, double> E_c_mins;
  std::unordered_map<const marley::Fragment*, double> Exf_maxes;

  TCanvas canvas;

  for (size_t j = 0; j < r.get_num_levels(); ++j) {

    // Get distributions for HF decay to the continuum for each unbound level
    // accessible by this reaction above the strength (B(F) + B(GT)) threshold.
    bool bound;
    double strength;
    double Exi = r.get_level_energy(j, bound, strength);
    if (bound || strength < strength_threshold) continue;

    std::cout << "Level " << j << " at Ex = " << Exi << " MeV has strength "
      << strength << std::endl;

    marley::Particle initial(marley_utils::get_nucleus_pid(Zi, Ai),
      marley::MassTable::get_atomic_mass(Zi, Ai) + Exi);

    marley::NuclearPhysics::hf_test3(Zi, Ai, initial, Exi, twoJi, Pi,
        gen.get_structure_db(), funcs, total_c_widths, E_c_mins, Exf_maxes);

    int s = marley::NuclearPhysics::get_fragments().size();
    // Prepare plots of the distributions
    for (int i = -1; i < s; ++i) {

      std::string pname_string;
      const marley::Fragment* f;

      if (i == -1) {
        f = nullptr;
        pname_string = marley_utils::particle_symbols.at(marley_utils::PHOTON);
      }

      else {
        f = &marley::NuclearPhysics::get_fragments().at(i);
        pname_string = marley_utils::particle_symbols.at(f->get_pid());
      }

      std::vector<double> xs, ys;
      std::cout << "DEBUG: fragment " << pname_string << " has E_c_min = " << E_c_mins[f]
        << " and Exf_max = " << Exf_maxes[f] << std::endl;
      for (double x = E_c_mins[f]; x < Exf_maxes[f]; x += (Exf_maxes[f] - E_c_mins[f]) / 1000) {
        xs.push_back(x);
        ys.push_back(funcs[f](x));
      }
      if (xs.size() == 0) continue;
      TGraph graph(xs.size(), &xs.front(), &ys.front());
      graph.SetTitle(("Continuum PDF for "
        + pname_string + " emission from "
        + std::to_string(Ai) + marley_utils::element_symbols.at(Zi) + " at Ex = "
        + std::to_string(Exi) + " MeV, B(GT) = "
        + std::to_string(strength)).c_str());
      graph.GetXaxis()->SetTitle("Final Excitation Energy");
      graph.GetXaxis()->CenterTitle();
      graph.GetXaxis()->SetTitleOffset(1.3);
      graph.GetYaxis()->SetTitle("PDF (unnormalized)");
      graph.GetYaxis()->CenterTitle();
      graph.SetLineColor(4);
      graph.SetLineWidth(3);
      graph.Draw("AL");
      canvas.SaveAs(("F_lev" + std::to_string(j) + "_" + pname_string
        + ".pdf").c_str());
      canvas.Clear();
    }
  }
}
