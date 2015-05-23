#include <iostream>
#include <stdexcept>

#include "TMarleyEvent.hh"

TMarleyEvent::TMarleyEvent(double E_level) {
  E_residue_level = E_level;
  reaction = nullptr;
  projectile = nullptr;
  target = nullptr;
  ejectile = nullptr;
  residue = nullptr;
}

// Sets one of the special particle pointers
// (e.g., TMarleyParticle::projectile) according
// to the supplied ParticleRole
void TMarleyEvent::assign_particle_pointer(TMarleyParticle* p,
  TMarleyEvent::ParticleRole r)
{
  switch(r) {
    case TMarleyEvent::ParticleRole::projectile:
      projectile = p;
      break;
    case TMarleyEvent::ParticleRole::target:
      target = p;
      break;
    case TMarleyEvent::ParticleRole::ejectile:
      ejectile = p;
      break;
    case TMarleyEvent::ParticleRole::residue:
      residue = p;
      break;
    default: // Do nothing if the particle does not have a
      break; // special role in this event
  }
}

TMarleyParticle* TMarleyEvent::get_residue() {
  return residue;
}

void TMarleyEvent::add_initial_particle(TMarleyParticle p,
  TMarleyEvent::ParticleRole r)
{
  initial_particles.push_back(p);

  if (r == TMarleyEvent::ParticleRole::ejectile)
    throw std::runtime_error(std::string("The ejectile")
    + " is not an initial state particle role.");
  if (r == TMarleyEvent::ParticleRole::residue)
    throw std::runtime_error(std::string("The residue")
    + " is not an initial state particle role.");

  assign_particle_pointer(&(initial_particles.back()), r);
}

void TMarleyEvent::add_final_particle(TMarleyParticle p,
  TMarleyEvent::ParticleRole r)
{
  final_particles.push_back(p);

  if (r == TMarleyEvent::ParticleRole::projectile)
    throw std::runtime_error(std::string("The projectile")
    + " is not an final state particle role.");
  if (r == TMarleyEvent::ParticleRole::target)
    throw std::runtime_error(std::string("The target")
    + " is not an final state particle role.");

  assign_particle_pointer(&(final_particles.back()), r);
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
