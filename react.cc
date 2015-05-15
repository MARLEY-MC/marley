#include <cmath>
#include <iostream>

#include "TMarleyReaction.hh"

int main(){

  // All hard-coded energies are in MeV
  double Ea = 8;
  double E_level = 2.289871;
  double cos_theta_c = -1;

  TMarleyReaction reaction;

  reaction.create_event(8);
  std::cout << "dsigma/dcos(theta) = "
    << reaction.differential_xs(E_level, Ea, cos_theta_c)
    << std::endl;
  std::cout << "sigma = "
    << reaction.total_xs(E_level, Ea)
    << std::endl;
  std::cout << "threshold Ea = "
    << reaction.E_threshold
    << std::endl;
  std::cout << "max E_level = "
    << reaction.max_level_energy(reaction.E_threshold)
    << std::endl;

  return 0;
}
