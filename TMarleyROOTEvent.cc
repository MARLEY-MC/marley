#include <iostream>
#include <stdexcept>

#include "TMarleyROOTEvent.hh"

//#ifdef USE_ROOT
//#ifndef __CINT__
//ClassImp(TMarleyROOTEvent);
//#endif
//#endif

TMarleyROOTEvent::TMarleyROOTEvent(double E_level) {
  E_residue_level = E_level;
  projectile = nullptr;
  target = nullptr;
  ejectile = nullptr;
  residue = nullptr;
}

// Sets one of the special particle pointers
// (e.g., TMarleyParticle::projectile) according
// to the supplied ParticleRole
void TMarleyROOTEvent::assign_particle_pointer(TMarleyParticle* p,
  TMarleyROOTEvent::ParticleRole r)
{
  switch(r) {
    case TMarleyROOTEvent::ParticleRole::pr_projectile:
      projectile = p;
      break;
    case TMarleyROOTEvent::ParticleRole::pr_target:
      target = p;
      break;
    case TMarleyROOTEvent::ParticleRole::pr_ejectile:
      ejectile = p;
      break;
    case TMarleyROOTEvent::ParticleRole::pr_residue:
      residue = p;
      break;
    default: // Do nothing if the particle does not have a
      break; // special role in this event
  }
}

void TMarleyROOTEvent::add_initial_particle(const TMarleyParticle& p,
  TMarleyROOTEvent::ParticleRole r)
{
  initial_particles.push_back(p);

  if (r == TMarleyROOTEvent::ParticleRole::pr_ejectile)
    throw std::runtime_error(std::string("The ejectile")
    + " is not an initial state particle role.");
  if (r == TMarleyROOTEvent::ParticleRole::pr_residue)
    throw std::runtime_error(std::string("The residue")
    + " is not an initial state particle role.");

  assign_particle_pointer(&initial_particles.back(), r);
}

void TMarleyROOTEvent::add_final_particle(const TMarleyParticle& p,
  TMarleyROOTEvent::ParticleRole r)
{
  final_particles.push_back(p);

  if (r == TMarleyROOTEvent::ParticleRole::pr_projectile)
    throw std::runtime_error(std::string("The projectile")
    + " is not an final state particle role.");
  if (r == TMarleyROOTEvent::ParticleRole::pr_target)
    throw std::runtime_error(std::string("The target")
    + " is not an final state particle role.");

  assign_particle_pointer(&final_particles.back(), r);
}
