#include <algorithm>
#include <climits>
#include <functional>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyGenerator.hh"
#include "TMarleyKinematics.hh"
#include "TMarleyNuclearPhysics.hh"

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
  TMarleyParity Pi = 1;

  double strength_threshold = 0.2;

  //std::cout << std::setprecision(16) << std::scientific;

  TMarleyGenerator gen("config.txt");

  TMarleyReaction r = gen.get_reactions().back();

  for (size_t j = 0; j < r.get_num_levels(); ++j) {

    // Get HF decay results for each unbound level accessible
    // by this reaction above the strength (B(F) + B(GT)) threshold.
    bool bound;
    double strength;
    double Exi = r.get_level_energy(j, bound, strength);
    if (bound || strength < strength_threshold) continue;

    std::cout << "Level " << j << " at Ex = " << Exi << " MeV has strength "
      << strength << std::endl;

    TMarleyParticle initial(marley_utils::get_nucleus_pid(Zi, Ai),
      TMarleyMassTable::get_atomic_mass(Zi, Ai) + Exi);

    //std::vector<double> KEs;
    //std::vector<int> f_pids;

    TMarleyHFTable hftable = TMarleyNuclearPhysics::create_hf_table(Zi, Ai,
      initial, Exi, twoJi, Pi, gen.get_structure_db(), gen);

    // Prepare spectrum plots
    TCanvas canvas;

    std::vector<int> pids(1, marley_utils::PHOTON);
    for (const auto& f : TMarleyNuclearPhysics::get_fragments())
      pids.push_back(f.get_pid());

    // Set up bin boundaries for the branching ratio plot
    std::vector<double> dpids(1, 0.);
    for (const auto& pid : pids) {
      int id = f_br_plot_ids.at(pid);
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

    // Get the decay channel objects from the HFTable
    const auto& dc_vec = hftable.get_channels();

    std::vector<double> KEs;
    std::vector<double> widths;
    std::vector<bool> c_flags;

    double total_width = hftable.get_total_width();

    for (const auto& pid : pids) {

      KEs.clear();
      widths.clear();
      c_flags.clear();

      int id = f_br_plot_ids.at(pid);

      double width = 0.;
      double max_KE = -std::numeric_limits<double>::max();

      std::string pname_string = marley_utils::particle_symbols.at(pid);

      // Number of decay channels for this product particle
      size_t num_channels = 0;

      // Loop over each of the decay channels
      for (size_t j = 0; j < hftable.get_widths().size(); ++j) {
	// If the current decay channel is for the current final particle, then
	// add the partial decay width to the total. Also get the outgoing
	// kinetic energy for the fragment corresponding to the center of the
	// final energy bin for this channel
        if (pid == dc_vec.at(j)->get_fragment_pid()) {

          double dc_width = hftable.get_widths().at(j);
          widths.push_back(dc_width);

          width += dc_width;
 
          TMarleyParticle first_product = TMarleyParticle(pid,
            dc_vec.at(j)->get_fragment_mass());

          double Exf = dc_vec.at(j)->get_bin_center_Exf();

          int Za = dc_vec.at(j)->get_fragment_Z();
          int Zf = Zi - Za;
          int Af = Ai - dc_vec.at(j)->get_fragment_A();
          double Mfgs = TMarleyMassTable::get_atomic_mass(Zf, Af)
            + Za*TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);
          
          TMarleyParticle second_product = TMarleyParticle(
            marley_utils::get_nucleus_pid(Zf, Af), Mfgs + Exf);

          TMarleyKinematics::two_body_decay(initial, first_product,
            second_product, 0., 0.);

          double KE = first_product.get_total_energy()
            - first_product.get_mass();

          KEs.push_back(KE);

          if (max_KE < KE) max_KE = KE;

          c_flags.push_back(dc_vec.at(j)->is_continuum());

          ++num_channels;
        }
      }

      // Load the branching ratio for this final particle into
      // the histogram of branching ratios
      double br = width / total_width;
      frag_br_histogram.Fill(id, br);

      // Skip making a spectrum plot for this final particle if it doesn't
      // have any decay channels available.
      if (num_channels == 0) continue; 

      // Create a histogram to store the emitted fragment kinetic energies.
      // Use 100 bins, and adjust the upper limit of the last bin so that
      // the maximum observed value just barely falls within it.
      TH1D frag_E_histogram((pname_string + "_E_histogram").c_str(),
        ("De-excitation " + pname_string
        + " spectrum (BR = " + std::to_string(br)
        + ") for " + std::to_string(Ai) + marley_utils::element_symbols.at(Zi)
        + " at Ex = " + std::to_string(Exi) + " MeV, B(GT) = "
        + std::to_string(strength)
        + "; Kinetic Energy [MeV]; Partial Branching Ratio").c_str(),
        100, 0, 1.3*std::nextafter(max_KE, DBL_MAX));

      for (size_t i = 0; i < widths.size(); ++i) {
        double partial_br = widths.at(i) / width;
        frag_E_histogram.Fill(KEs.at(i), partial_br);
        std::cout << "Decay via emission of a " << KEs.at(i) << " MeV "
          << pname_string << " has partial BR = " << partial_br;
        if (c_flags.at(i)) std::cout << " (continuum bin)";
        std::cout << std::endl;
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
    for (const auto& pair : f_br_plot_ids) {
      int id = pair.second;
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
