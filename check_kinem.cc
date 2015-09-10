#include <iostream>

#include "TMarleyParticle.hh"
#include "TMarleyKinematics.hh"
//#include "TMarleyMassTable.hh"

int main() {

  TMarleyParticle p1(30, 2, 8, 7, 50);
  TMarleyParticle p2(20, 15); 
  TMarleyParticle p3(20, 10); 

  TMarleyKinematics::two_body_decay(p1, p2, p3, 0.9, 0.8);
  std::cout << "p2: E = " << p2.get_total_energy()
    << ", px = " << p2.get_px()
    << ", py = " << p2.get_py()
    << ", pz = " << p2.get_pz() << std::endl;
  std::cout << "p3: E = " << p3.get_total_energy()
    << ", px = " << p3.get_px()
    << ", py = " << p3.get_py()
    << ", pz = " << p3.get_pz() << std::endl;
  std::cout << "Energy conservation: E2 + E3 - E1 = "
    << p2.get_total_energy() + p3.get_total_energy() - p1.get_total_energy()
    << std::endl;
  std::cout << "px conservation: px2 + px3 - px1 = "
    << p2.get_px() + p3.get_px() - p1.get_px()
    << std::endl;
  std::cout << "py conservation: py2 + py3 - py1 = "
    << p2.get_py() + p3.get_py() - p1.get_py()
    << std::endl;
  std::cout << "pz conservation: pz2 + pz3 - pz1 = "
    << p2.get_pz() + p3.get_pz() - p1.get_pz()
    << std::endl;
}
