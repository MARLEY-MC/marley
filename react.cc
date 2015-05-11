#include <cmath>
#include <iostream>

#include "TMarleyReaction.hh"

int main(){

  // All hard-coded energies are in MeV
  double E_level = 2.289871;
  double Ea = 8;
  //double Etot = Ea + mb;

  double theta_c = std::acos(0);
  double cos_theta_c = std::cos(theta_c);

  TMarleyReaction reaction;

  double Ec = reaction.ejectile_energy(E_level, Ea, cos_theta_c);

  // Compute the energy and scattering angle of the residue
  //double Ed = Etot - Ec;
  //double theta_d = std::asin(std::sqrt(Ec*Ec - mc*mc)
  //  * std::sin(theta_c) / std::sqrt(Ed*Ed - md*md));
  //double cos_theta_d = std::cos(theta_d);

  std::cout.precision(15);
  std::cout << std::scientific;
  std::cout << "Ec = " << Ec << std::endl;

  std::cout << "\"Exact:\" " << reaction.fermi_function(18., 20, true) << std::endl;
  std::cout << "Approximate: " << reaction.fermi_approx(18., 20, true) << std::endl; 

  return 0;
}
