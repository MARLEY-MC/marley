#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleyMassTable.hh"
#include "TMarleySphericalOpticalModel.hh"
#include "TMarleyNuclearPhysics.hh"
#include "TMarleyIntegrator.hh"
#include "TMarleyGenerator.hh"

int main() {

  int Zi = 19;
  int Ai = 40;
  int twoJi = 2;
  TMarleyParity Pi = 1;

  //std::cout << std::setprecision(16) << std::scientific;

  std::vector<double> widths;
  std::vector<int> pids;

  std::cout << "  Exi";
  for (const auto& f : TMarleyNuclearPhysics::get_fragments())
    std::cout << "  " << f.get_pid();
  std::cout << "  22" << std::endl;

  for (double Exi = 0.1; Exi <= 100; Exi += 0.1) {
    for (const auto& f : TMarleyNuclearPhysics::get_fragments()) {
      int Zf = Zi - f.get_Z();
      int Af = Ai - f.get_A();
      TMarleySphericalOpticalModel om(Zf, Af);
      TMarleyDecayScheme ds(Zf, Af,
      //  //std::string("/home/sjg/Desktop/ensdf_22_Oct_2014/ensdf.0")
      //  //+ std::to_string(Af));
        std::string("/home/sjg/talys/structure/levels/exp/z0")
        + std::to_string(Zf), TMarleyDecayScheme::FileFormat::talys);
      double width = TMarleyNuclearPhysics::hf_fragment_partial_width(Zi, Ai, Exi,
        twoJi, Pi, f, om, ds);
      //std::cout << "Particle ID = " << f.get_pid() << ", width = " << width
      //  << std::endl;
      widths.push_back(width);
    }

    TMarleyDecayScheme ds(Zi, Ai,
    //  //std::string("/home/sjg/Desktop/ensdf_22_Oct_2014/ensdf.0")
    //  //+ std::to_string(Af));
      std::string("/home/sjg/talys/structure/levels/exp/z0")
      + std::to_string(Zi), TMarleyDecayScheme::FileFormat::talys);
    double width = TMarleyNuclearPhysics::hf_gamma_partial_width(Exi, twoJi, Pi, ds);
    //std::cout << "Particle ID = 22, width = " << width << std::endl;
    widths.push_back(width);

    double total_width = 0.;
    std::cout << "  " << Exi;
    for (auto w : widths) total_width += w;
    for (auto w : widths) std::cout << "  " << w / total_width;
    std::cout << std::endl;
    widths.clear();
  }
}
