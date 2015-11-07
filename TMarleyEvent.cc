#include <iostream>
#include <stdexcept>

#include "TMarleyEvent.hh"

//#ifdef USE_ROOT
//#ifndef __CINT__
//ClassImp(TMarleyEvent);
//#endif
//#endif

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
    case TMarleyEvent::ParticleRole::pr_projectile:
      projectile = p;
      break;
    case TMarleyEvent::ParticleRole::pr_target:
      target = p;
      break;
    case TMarleyEvent::ParticleRole::pr_ejectile:
      ejectile = p;
      break;
    case TMarleyEvent::ParticleRole::pr_residue:
      residue = p;
      break;
    default: // Do nothing if the particle does not have a
      break; // special role in this event
  }
}

TMarleyParticle* TMarleyEvent::get_residue() {
  return residue;
}

TMarleyParticle* TMarleyEvent::get_ejectile() {
  return ejectile;
}

TMarleyParticle* TMarleyEvent::get_projectile() {
  return projectile;
}

TMarleyParticle* TMarleyEvent::get_target() {
  return target;
}

void TMarleyEvent::add_initial_particle(const TMarleyParticle& p,
  TMarleyEvent::ParticleRole r)
{
  initial_particles.push_back(p);

  if (r == TMarleyEvent::ParticleRole::pr_ejectile)
    throw std::runtime_error(std::string("The ejectile")
    + " is not an initial state particle role.");
  if (r == TMarleyEvent::ParticleRole::pr_residue)
    throw std::runtime_error(std::string("The residue")
    + " is not an initial state particle role.");

  assign_particle_pointer(&(initial_particles.back()), r);
}

void TMarleyEvent::add_final_particle(const TMarleyParticle& p,
  TMarleyEvent::ParticleRole r)
{
  final_particles.push_back(p);

  if (r == TMarleyEvent::ParticleRole::pr_projectile)
    throw std::runtime_error(std::string("The projectile")
    + " is not an final state particle role.");
  if (r == TMarleyEvent::ParticleRole::pr_target)
    throw std::runtime_error(std::string("The target")
    + " is not an final state particle role.");

  assign_particle_pointer(&(final_particles.back()), r);
}

void TMarleyEvent::set_reaction(TMarleyReaction* r) {
  reaction = r;
}

// Prints information about this event to stdout
void TMarleyEvent::print_event() {
  std::cout << "*** Initial particles ***" << std::endl;
  for (const auto& i : initial_particles)
  {
    std::cout << "id: " << i.get_id() << "   energy: " << i.get_total_energy()
      << " MeV" << std::endl;
  }
  std::cout << std::endl << std::endl << "*** Final particles ***" << std::endl;
  for (const auto& f : final_particles)
  {
    std::cout << "id: " << f.get_id() << "   energy: " << f.get_total_energy()
      << " MeV" << std::endl;
  }
}

double TMarleyEvent::get_E_level() {
  return E_residue_level;
}

std::ostream& operator<< (std::ostream& out, const TMarleyEvent& e) {
  out << *(e.projectile) << std::endl;
  //for (const TMarleyParticle& p : e.initial_particles) {
  //  out << p << std::endl;
  //}
  //out << std::endl << std::endl;
  for (auto& p : e.final_particles) {
    out << p << std::endl;
  }
  return out;
}
