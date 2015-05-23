#include <iostream>
#include <string>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleyEvent.hh"
#include "TMarleyReaction.hh"

int main(){

  // Select the isotope and ENSDF file to use for the simulation
  std::string nuc_id = marley_utils::nuc_id(19, 40); // 40K
  std::string filename = "ensdf.040";

  // Create a decay scheme object to store data
  // imported from the ENSDF file
  TMarleyDecayScheme ds(nuc_id, filename);

  TMarleyReaction r("ve40ArCC.react", &ds);

  // TODO: debug numerical errors that arise when
  // Ea = E_threshold
  // Incident neutrino energy
  double Ea = 8; // MeV

  // Simulate a reaction at the given
  // incident neutrino energy
  for (int i = 1; i <= 100; i++) {
    TMarleyEvent e = r.create_event(Ea);
    std::cout << std::endl << std::endl;
    e.print_event();
    std::cout << std::endl << std::endl;
  }

  //double E_threshold = r.get_threshold_energy();
  //std::cout << "threshold Ea = "
  //  << E_threshold
  //  << std::endl;
  //std::cout << "max E_level at threshold = "
  //  << r.max_level_energy(E_threshold)
  //  << std::endl;

  //for(int i = 0; i < 200; i++) {
  //  double cos = -1.0 + i*0.01;
  //  std::cout << cos << " " << r.differential_xs(0, Ea, cos) << std::endl;
  //}

  return 0;
}
