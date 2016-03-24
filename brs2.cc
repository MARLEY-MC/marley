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

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

static const std::unordered_map<int, int> f_br_plot_ids = {
  { marley_utils::PHOTON, 10},
  { marley_utils::NEUTRON, 20},
  { marley_utils::PROTON, 30},
  { marley_utils::DEUTERON, 40},
  { marley_utils::TRITON, 50},
  { marley_utils::HELION, 60},
  { marley_utils::ALPHA, 70},
};

int main() {

  int Zi = 19;
  int Ai = 40;
  int twoJi = 2;
  marley::Parity Pi = 1;

  double strength_threshold = 0.2;

  //std::cout << std::setprecision(16) << std::scientific;

  marley::Generator gen("config.txt");

  marley::Reaction r = gen.get_reactions().back();

  for (size_t j = 0; j < r.get_num_levels(); ++j) {

    // Get HF decay results for each unbound level accessible
    // by this reaction above the strength (B(F) + B(GT)) threshold.
    bool bound;
    double strength;
    double Exi = r.get_level_energy(j, bound, strength);
    if (bound || strength < strength_threshold) continue;

    std::cout << "Level " << j << " at Ex = " << Exi << " MeV has strength "
      << strength << std::endl;

    marley::Particle initial(marley_utils::get_nucleus_pid(Zi, Ai),
      marley::MassTable::get_atomic_mass(Zi, Ai) + Exi);

    std::vector<double> KEs;
    std::vector<int> f_pids;

    size_t num_trials = 1e5;

    //marley::NuclearPhysics::hf_test(Zi, Ai, initial, Exi, twoJi, Pi,
    //  gen.get_structure_db(), gen);
    marley::NuclearPhysics::hf_test2(Zi, Ai, initial, Exi, twoJi, Pi,
        gen.get_structure_db(), gen, num_trials, KEs, f_pids);

    // Initialize fragment type counts
    std::unordered_map<int, int> fcounts;
    std::unordered_map<int, double> f_max_KEs;
    fcounts[marley_utils::PHOTON] = 0;
    f_max_KEs[marley_utils::PHOTON] = 0.;
    for (const auto& f : marley::NuclearPhysics::get_fragments()) {
      fcounts[f.get_pid()] = 0;
      f_max_KEs[f.get_pid()] = 0.;
    }

    // Count fragments and get maximum kinetic energies
    for (size_t i = 0; i < f_pids.size(); ++i) {
      int f_pid = f_pids.at(i);
      ++fcounts.at(f_pid);
      if (f_max_KEs.at(f_pid) <  KEs.at(i))
        f_max_KEs[f_pid] = KEs.at(i);
    }

    // Print branching ratios and prepare spectrum plots
    TCanvas canvas;

    // Set up bin boundaries for the branching ratio plot
    std::vector<double> dpids(1, 0.);
    for (const auto& pair : fcounts) {
      int id = f_br_plot_ids.at(pair.first);
      dpids.push_back(id);
      dpids.push_back(id + 5.);
    }

    // ROOT expects the particle IDs to be in ascending order, so
    // sort them before using them as bin boundaries
    std::sort(dpids.begin(), dpids.end());

    // Add two extra bins on the end
    dpids.push_back(dpids.back() + 5.);
    dpids.push_back(dpids.back() + 10.);

    TH1D frag_br_histogram("fragment_BR_histogram",
     ("De-excitation branching ratios for " + std::to_string(Ai)
     + marley_utils::element_symbols.at(Zi)
     + " at Ex = " + std::to_string(Exi) + " MeV, B(GT) = "
     + std::to_string(strength)
     + "; Particle ID; Branching Ratio").c_str(),
     dpids.size() - 1, &dpids.front());

    for (const auto& pair : fcounts) {
      double br = pair.second / static_cast<double>(num_trials);
      //std::cout << pair.first << " branching ratio: " << br << std::endl;
      // If there were no instances of this fragment being produced, skip
      // ahead to the next fragment.
      if (pair.second == 0) continue;

      double max_KE = f_max_KEs.at(pair.first);

      // Create a histogram to store the emitted fragment kinetic energies.
      // Use 100 bins, and adjust the upper limit of the last bin so that
      // the maximum observed value just barely falls within it.
      std::string pname_string = marley_utils::particle_symbols.at(pair.first);
      TH1D frag_E_histogram((pname_string + "_E_histogram").c_str(),
        ("De-excitation " + pname_string
        + " spectrum (BR = " + std::to_string(br)
        + ") for " + std::to_string(Ai) + marley_utils::element_symbols.at(Zi)
        + " at Ex = " + std::to_string(Exi) + " MeV, B(GT) = "
        + std::to_string(strength)
        + "; Kinetic Energy [MeV]; Partial Branching Ratio").c_str(),
        100, 0, 1.3*std::nextafter(max_KE, DBL_MAX));

      for (size_t i = 0; i < f_pids.size(); ++i)
        if (f_pids.at(i) == pair.first) {
          frag_E_histogram.Fill(KEs.at(i), 1.0 / fcounts.at(pair.first));
          int id = f_br_plot_ids.at(pair.first);
          frag_br_histogram.Fill(id, 1.0 / num_trials);
        }

      frag_E_histogram.Draw();
      frag_E_histogram.GetYaxis()->SetTitleOffset(1.3);
      frag_E_histogram.GetXaxis()->CenterTitle();
      frag_E_histogram.GetYaxis()->CenterTitle();
      frag_E_histogram.GetXaxis()->SetTitleOffset(1.3);
      frag_E_histogram.SetLineWidth(1.8);
      canvas.SaveAs(("E_lev" + std::to_string(j) + "_" + pname_string
        + ".pdf").c_str());
      canvas.Clear();

    }

    //canvas.SetLogx();
    frag_br_histogram.Draw("bar");
    frag_br_histogram.GetYaxis()->SetTitleOffset(1.3);
    frag_br_histogram.GetXaxis()->CenterTitle();
    frag_br_histogram.GetYaxis()->CenterTitle();
    frag_br_histogram.GetXaxis()->SetTitleOffset(1.3);
    frag_br_histogram.SetLineWidth(1.8);

    // Label particles by symbol instead of ID number
    TAxis& Xax = *frag_br_histogram.GetXaxis();
    Xax.SetLabelSize(0.055);
    for (const auto& pair : fcounts) {
      int id = f_br_plot_ids.at(pair.first);
      int bin_index = Xax.FindBin(id);
      Xax.SetBinLabel(bin_index, marley_utils::particle_symbols.at(
        pair.first).c_str());
    }
    canvas.Update();

    canvas.SaveAs(("E_lev" + std::to_string(j) + "_brs.pdf").c_str());
    canvas.Clear();

    TFile test_file("test_plot.root", "UPDATE");
    frag_br_histogram.Write();
    test_file.Close();
  }
}
