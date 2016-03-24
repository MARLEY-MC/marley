#include "marley_utils.hh"
#include "Particle.hh"

//#ifdef USE_ROOT
//#ifndef __CINT__
//ClassImp(marley::Particle);
//#endif
//#endif

void marley::Particle::init(int id, double E, double p_x, double p_y, double p_z,
  double m, int q)
{
  particle_id = id;
  charge = q;
  total_energy = E;
  px = p_x;
  py = p_y;
  pz = p_z;
  mass = m;
}

marley::Particle::Particle() {
  init(0, 0., 0., 0., 0., 0., 0);
}

marley::Particle::Particle(int id, double m) {
  int q = marley_utils::get_particle_charge(id);
  init(id, m, 0., 0., 0., m, q);
}

marley::Particle::Particle(int id, double m, int q) {
  init(id, m, 0., 0., 0., m, q);
}

marley::Particle::Particle(int id, double p_x,
  double p_y, double p_z, double m)
{
  int q = marley_utils::get_particle_charge(id);
  total_energy = std::sqrt(std::pow(p_x, 2)
    + std::pow(p_y, 2) + std::pow(p_z, 2)
    + std::pow(m, 2));
  init(id, total_energy, p_x, p_y, p_z, m, q);
}

marley::Particle::Particle(int id, double p_x,
  double p_y, double p_z, double m, int q)
{
  total_energy = std::sqrt(std::pow(p_x, 2)
    + std::pow(p_y, 2) + std::pow(p_z, 2)
    + std::pow(m, 2));
  init(id, total_energy, p_x, p_y, p_z, m, q);
}

marley::Particle::Particle(int id, double E, double p_x,
  double p_y, double p_z, double m)
{
  int q = marley_utils::get_particle_charge(id);
  init(id, E, p_x, p_y, p_z, m, q);
}

marley::Particle::Particle(int id, double E, double p_x,
  double p_y, double p_z, double m, int q)
{
  init(id, E, p_x, p_y, p_z, m, q);
}

double marley::Particle::get_total_energy() const {
  return total_energy;
}

double marley::Particle::get_mass() const {
  return mass;
}

double marley::Particle::get_px() const {
  return px;
}

double marley::Particle::get_py() const {
  return py;
}

double marley::Particle::get_pz() const {
  return pz;
}

int marley::Particle::get_id() const {
  return particle_id;
}

void marley::Particle::add_child(marley::Particle* child) {
  children.push_back(child);
}

std::ostream& operator<< (std::ostream& out, const marley::Particle& p) {
  out << p.particle_id << " " << p.total_energy << " " << p.px << " " << p.py
    << " " << p.pz;
  return out;
}
