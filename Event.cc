#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "Event.hh"

//#ifdef USE_ROOT
//#ifndef __CINT__
//ClassImp(marley::Event);
//#endif
//#endif

marley::Event::Event(double E_level) {
  E_residue_level = E_level;
  reaction = nullptr;
  projectile = nullptr;
  target = nullptr;
  ejectile = nullptr;
  residue = nullptr;
}

// Sets one of the special particle pointers
// (e.g., marley::Particle::projectile) according
// to the supplied ParticleRole
void marley::Event::assign_particle_pointer(marley::Particle* p,
  marley::Event::ParticleRole r)
{
  switch(r) {
    case marley::Event::ParticleRole::pr_projectile:
      projectile = p;
      break;
    case marley::Event::ParticleRole::pr_target:
      target = p;
      break;
    case marley::Event::ParticleRole::pr_ejectile:
      ejectile = p;
      break;
    case marley::Event::ParticleRole::pr_residue:
      residue = p;
      break;
    default: // Do nothing if the particle does not have a
      break; // special role in this event
  }
}

marley::Particle* marley::Event::get_residue() {
  return residue;
}

marley::Particle* marley::Event::get_ejectile() {
  return ejectile;
}

marley::Particle* marley::Event::get_projectile() {
  return projectile;
}

marley::Particle* marley::Event::get_target() {
  return target;
}

void marley::Event::add_initial_particle(const marley::Particle& p,
  marley::Event::ParticleRole r)
{
  initial_particles.push_back(p);

  if (r == marley::Event::ParticleRole::pr_ejectile)
    throw std::runtime_error(std::string("The ejectile")
    + " is not an initial state particle role.");
  if (r == marley::Event::ParticleRole::pr_residue)
    throw std::runtime_error(std::string("The residue")
    + " is not an initial state particle role.");

  assign_particle_pointer(&(initial_particles.back()), r);
}

void marley::Event::add_final_particle(const marley::Particle& p,
  marley::Event::ParticleRole r)
{
  final_particles.push_back(p);

  if (r == marley::Event::ParticleRole::pr_projectile)
    throw std::runtime_error(std::string("The projectile")
    + " is not an final state particle role.");
  if (r == marley::Event::ParticleRole::pr_target)
    throw std::runtime_error(std::string("The target")
    + " is not an final state particle role.");

  assign_particle_pointer(&(final_particles.back()), r);
}

void marley::Event::set_reaction(marley::Reaction* r) {
  reaction = r;
}

// Prints information about this event to stdout
void marley::Event::print_event() {
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

double marley::Event::get_E_level() {
  return E_residue_level;
}

std::ostream& operator<< (std::ostream& out, const marley::Event& e) {
  out << *(e.projectile) << std::endl;
  //for (const marley::Particle& p : e.initial_particles) {
  //  out << p << std::endl;
  //}
  //out << std::endl << std::endl;
  for (auto& p : e.final_particles) {
    out << p << std::endl;
  }
  return out;
}

// Function that dumps a marley::Particle to an output stream in HEPEVT format. This
// is a private helper function for the publically-accessible write_hepevt.
void marley::Event::dump_hepevt_particle(const marley::Particle& p, std::ostream& os,
  bool track)
{
  if (track) os << "1 ";
  else os << "0 ";

  // TODO: improve this entry to give the user more control over the vertex
  // location and to reflect the parent-daughter relationships between
  // particles.
  // Factors of 1000. are used to convert MeV to GeV for the HEPEvt format
  os << p.get_id() << " 0 0 0 0 " << p.get_px() / 1000. << " " << p.get_py() / 1000. << " "
    << p.get_pz() / 1000. << " " << p.get_total_energy() / 1000. << " " << p.get_mass() / 1000.
    // Spacetime origin is currently used as the initial position 4-vector for
    // all particles
    << " 0. 0. 0. 0." << std::endl;
}

void marley::Event::write_hepevt(size_t event_num, std::ostream& out) {
  out << std::setprecision(16) << std::scientific;
  out << event_num  << " " << final_particles.size() + 1 << std::endl;
  dump_hepevt_particle(*projectile, out, false);
  for (const auto& fp : final_particles) dump_hepevt_particle(fp, out, true);
}
