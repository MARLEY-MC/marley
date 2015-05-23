#include <iostream>

#include "TMarleyEvent.hh"

TMarleyEvent::TMarleyEvent(double E_level) {
  E_residue_level = E_level;
  reaction = nullptr;
  //projectile = nullptr;
  //target = nullptr;
  //ejectile = nullptr;
  //residue = nullptr;
}

void TMarleyEvent::add_initial_particle(TMarleyParticle p) {
  initial_particles.push_back(p);
}

void TMarleyEvent::add_final_particle(TMarleyParticle p) {
  final_particles.push_back(p);
}

void TMarleyEvent::set_reaction(TMarleyReaction* r) {
  reaction = r;
}

std::vector<TMarleyParticle>* TMarleyEvent::get_initial_particles() {
  return &initial_particles;
}

std::vector<TMarleyParticle>* TMarleyEvent::get_final_particles() {
  return &final_particles;
}

// Prints information about this event to stdout
void TMarleyEvent::print_event() {
  std::cout << "*** Initial particles ***" << std::endl;
  for (std::vector<TMarleyParticle>::iterator i = initial_particles.begin();
    i != initial_particles.end(); ++i)
  {
    std::cout << "id: " << i->get_id() << "   energy: " << i->get_total_energy()
      << " MeV" << std::endl;
  }
  std::cout << std::endl << std::endl << "*** Final particles ***" << std::endl;
  for (std::vector<TMarleyParticle>::iterator i = final_particles.begin();
    i != final_particles.end(); ++i)
  {
    std::cout << "id: " << i->get_id() << "   energy: " << i->get_total_energy()
      << " MeV" << std::endl;
  }
}
