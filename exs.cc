#ifndef USE_ROOT
#error "This MARLEY analysis code currently must be compiled with ROOT support."
#endif

#include <vector>

#include "marley_utils.hh"
//#include "DecayScheme.hh"
//#include "Event.hh"
//#include "Generator.hh"
#include "MassTable.hh"
#include "NeutrinoSource.hh"
//#include "Reaction.hh"

#include "TAxis.h"
//#include "TFile.h"
#include "TGraph.h"
//#include "TH1D.h"
//#include "TTree.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"

std::string nu_type_string(marley::NeutrinoSource::NeutrinoType type) {
  if (type == marley::NeutrinoSource::NeutrinoType::ElectronNeutrino)
    return std::string("#nu_{e}");
  else if (type == marley::NeutrinoSource::NeutrinoType::ElectronAntineutrino)
    return std::string("#bar{#nu}_{e}");
  else if (type == marley::NeutrinoSource::NeutrinoType::MuonNeutrino)
    return std::string("#nu_{#mu}");
  else if (type ==   marley::NeutrinoSource::NeutrinoType::MuonAntineutrino)
    return std::string("#bar{#nu}_{#mu}");
  else if (type ==  marley::NeutrinoSource::NeutrinoType::TauonNeutrino)
    return std::string("#nu_{#tau}");
  else if (type ==  marley::NeutrinoSource::NeutrinoType::TauonAntineutrino)
    return std::string("#bar{#nu}_{#tau}");
  // TODO: add exception here
  else return std::string("");
}

double nu_e_xs(double E_nu_lab, marley::NeutrinoSource::NeutrinoType type) {
  // Fermi coupling constant (MeV^(-2)) 
  static constexpr double GF = 1.16637e-11;

  // sin^2(theta_W) (sine squared of the weak mixing angle)
  // Effective value taken from 2014 PDG Review of Particle Physics,
  // Table 1.1, "Physical Constants"
  static constexpr double sin2thetaw = 0.23155;

  // Mandelstam s (square of the total center of mass energy)
  double me = marley::MassTable::get_particle_mass(marley_utils::ELECTRON);
  double me2 = std::pow(me, 2);
  double s = me2 + 2 * me * E_nu_lab;
  // CM frame neutrino energy
  double E_nu_cm = (s - me2) / (2 * std::sqrt(s));
  // Coupling constants
  double g1, g2;
  if (type == marley::NeutrinoSource::NeutrinoType::ElectronNeutrino) {
    g1 = 0.5 + sin2thetaw;
    g2 = sin2thetaw;
  }
  else if (type == marley::NeutrinoSource::NeutrinoType::ElectronAntineutrino) {
    g1 = sin2thetaw;
    g2 = 0.5 + sin2thetaw;
  }
  else if (type == marley::NeutrinoSource::NeutrinoType::MuonNeutrino
    || type == marley::NeutrinoSource::NeutrinoType::TauonNeutrino)
  {
    g1 = -0.5 + sin2thetaw;
    g2 = sin2thetaw;
  }
  else {
    // type == marley::AntineutrinoSource::AntineutrinoType::MuonAntineutrino
    // || type == marley::AntineutrinoSource::AntineutrinoType::TauonAntineutrino
    g1 = sin2thetaw;
    g2 = -0.5 + sin2thetaw;
  }

  // Helper variables
  double g22_over_three = std::pow(g2, 2) / 3.0;
  double me2_over_s = me2 / s;
  // Total cross section in natural units (MeV^(-2))
  double xs = (4 / marley_utils::pi) * std::pow(GF * E_nu_cm, 2)
    * (std::pow(g1, 2) + (g22_over_three - g1*g2)*me2_over_s
    + g22_over_three*(1 + std::pow(me2_over_s, 2)));
  // Total cross section in 10^-46 cm^2
  double xs_new_units = xs / marley_utils::mb * 1e19;
  return xs_new_units;
}

int main(){

  TCanvas canvas;

  std::vector<TGraph*> graphs;

  std::vector<marley::NeutrinoSource::NeutrinoType>
    nu_types = { marley::NeutrinoSource::NeutrinoType::ElectronNeutrino,
    marley::NeutrinoSource::NeutrinoType::ElectronAntineutrino,
    marley::NeutrinoSource::NeutrinoType::MuonNeutrino,
    marley::NeutrinoSource::NeutrinoType::MuonAntineutrino,
    marley::NeutrinoSource::NeutrinoType::TauonNeutrino,
    marley::NeutrinoSource::NeutrinoType::TauonAntineutrino,
  };

  TLegend xs_legend(0.8, 0.15, 0.9, 0.6);
  xs_legend.SetMargin(0.2);
  xs_legend.SetTextSize(0.03);

  for (size_t j = 0; j < nu_types.size(); ++j) {
    std::vector<double> Es;
    std::vector<double> xs;
    marley::NeutrinoSource::NeutrinoType type = nu_types.at(j);
    for (double Enu = 0.1; Enu <= 100; Enu += 0.01) {
      Es.push_back(Enu);
      xs.push_back(nu_e_xs(Enu, type));
    }
    TGraph* xs_graph = new TGraph(Es.size(), &Es.front(), &xs.front());
    xs_graph->SetLineColor(j + 1);
    xs_graph->SetLineWidth(3);
    if (j == 0) {
      xs_graph->Draw("AL");
      xs_graph->SetTitle("Total Cross Section for  #nu ES on e^{#minus}");
      xs_graph->GetXaxis()->SetTitle("Neutrino Energy (MeV)");
      xs_graph->GetXaxis()->CenterTitle();
      xs_graph->GetXaxis()->SetTitleOffset(1.3);
      xs_graph->GetYaxis()->SetTitle("Cross Section (10^{ -46} cm^{2})");
      xs_graph->GetYaxis()->CenterTitle();
      xs_legend.AddEntry(xs_graph, nu_type_string(type).c_str(), "l");
    }
    else {
      xs_graph->Draw("L");
      xs_legend.AddEntry(xs_graph, nu_type_string(type).c_str(), "l");
    }
    graphs.push_back(xs_graph);
  }

  xs_legend.Draw();

  canvas.SetLogx(1);
  canvas.SetLogy(1);

  canvas.SaveAs("nu_e_xs.pdf");

  canvas.Clear();

  for (const auto g : graphs) delete g;

  return 0;
}
