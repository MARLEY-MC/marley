#include <algorithm>
#include <climits>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#include "marley_utils.hh"
#include "Generator.hh"
#include "NuclearPhysics.hh"

int main() {

  int Zi = 19;
  int Ai = 40;
  int twoJi = 2;
  marley::Parity Pi = 1;
  double Exi = 12;
  size_t num_trials = 1e5;

  //std::cout << std::setprecision(16) << std::scientific;

  marley::Generator gen("config.txt");

  marley::Particle initial(marley_utils::get_nucleus_pid(Zi, Ai),
    marley::MassTable::get_atomic_mass(Zi, Ai) + Exi);

  marley::HFTable hftable = marley::NuclearPhysics::create_hf_table(Zi, Ai,
    initial, Exi, twoJi, Pi, gen.get_structure_db(), gen);

  std::cout << "Num channels = " << hftable.get_num_channels() << std::endl;

  double Ex;
  int twoJ;
  marley::Parity P;
  for (size_t j = 0; j < num_trials; ++j) {
    Ex = Exi;
    twoJ = twoJi;
    P = Pi;
    hftable.sample_channel(gen).get_post_decay_parameters(Ex, twoJ, P);
    std::cout << "Event " << j << " sampled Ex = " << Ex << ", twoJ = " << twoJ
      << ", P = " << P << std::endl;
  }
}
