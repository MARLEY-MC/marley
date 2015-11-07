#include "marley_utils.hh"
#include "TMarleyParticle.hh"

//#ifdef USE_ROOT
//#ifndef __CINT__
//ClassImp(TMarleyParticle);
//#endif
//#endif

void TMarleyParticle::init(int id, double E, double p_x, double p_y, double p_z,
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

TMarleyParticle::TMarleyParticle() {
  init(0, 0., 0., 0., 0., 0., 0);
}

TMarleyParticle::TMarleyParticle(int id, double m) {
  int q = marley_utils::get_particle_charge(id);
  init(id, m, 0., 0., 0., m, q);
}

TMarleyParticle::TMarleyParticle(int id, double m, int q) {
  init(id, m, 0., 0., 0., m, q);
}

TMarleyParticle::TMarleyParticle(int id, double p_x,
  double p_y, double p_z, double m)
{
  int q = marley_utils::get_particle_charge(id);
  total_energy = std::sqrt(std::pow(p_x, 2)
    + std::pow(p_y, 2) + std::pow(p_z, 2)
    + std::pow(m, 2));
  init(id, total_energy, p_x, p_y, p_z, m, q);
}

TMarleyParticle::TMarleyParticle(int id, double p_x,
  double p_y, double p_z, double m, int q)
{
  total_energy = std::sqrt(std::pow(p_x, 2)
    + std::pow(p_y, 2) + std::pow(p_z, 2)
    + std::pow(m, 2));
  init(id, total_energy, p_x, p_y, p_z, m, q);
}

TMarleyParticle::TMarleyParticle(int id, double E, double p_x,
  double p_y, double p_z, double m)
{
  int q = marley_utils::get_particle_charge(id);
  init(id, E, p_x, p_y, p_z, m, q);
}

TMarleyParticle::TMarleyParticle(int id, double E, double p_x,
  double p_y, double p_z, double m, int q)
{
  init(id, E, p_x, p_y, p_z, m, q);
}

double TMarleyParticle::get_total_energy() const {
  return total_energy;
}

double TMarleyParticle::get_mass() const {
  return mass;
}

double TMarleyParticle::get_px() const {
  return px;
}

double TMarleyParticle::get_py() const {
  return py;
}

double TMarleyParticle::get_pz() const {
  return pz;
}

int TMarleyParticle::get_id() const {
  return particle_id;
}

void TMarleyParticle::add_child(TMarleyParticle* child) {
  children.push_back(child);
}

std::ostream& operator<< (std::ostream& out, const TMarleyParticle& p) {
  out << p.particle_id << " " << p.total_energy << " " << p.px << " " << p.py
    << " " << p.pz;
  return out;
}
