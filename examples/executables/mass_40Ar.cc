#include <iostream>

#include "marley/MassTable.hh"

constexpr int Z = 18; // Ar
constexpr int A = 40;

int main() {

  const marley::MassTable& mt = marley::MassTable::Instance();
  double mass_Ar40 = mt.get_atomic_mass( Z, A );

  std::cout << "MARLEY was";
  #ifndef USE_ROOT
    std::cout << " not";
  #endif
  std::cout << " built with ROOT support.\n";

  std::cout << "The atomic mass of 40Ar is " << mass_Ar40 << " MeV/c^2\n";

  return 0;
}
