#include "TMarleyParticle.hh"

TMarleyParticle::TMarleyParticle(int id, double E, double p_x,
  double p_y, double p_z, double m)
{
  particle_id = id;
  total_energy = E;
  px = p_x;
  py = p_y;
  pz = p_z;
  mass = m;
  parent = nullptr;
}

double TMarleyParticle::get_total_energy() {
  return total_energy;
}

double TMarleyParticle::get_mass() {
  return mass;
}

double TMarleyParticle::get_px() {
  return px;
}

double TMarleyParticle::get_py() {
  return py;
}

double TMarleyParticle::get_pz() {
  return pz;
}

int TMarleyParticle::get_id() {
  return particle_id;
}

TMarleyParticle* TMarleyParticle::get_parent() {
  return parent;
}

std::vector<TMarleyParticle*>* TMarleyParticle::get_children() {
  return &children;
}

void TMarleyParticle::set_parent(TMarleyParticle* p) {
  parent = p;
}

void TMarleyParticle::add_child(TMarleyParticle* child) {
  children.push_back(child);
}
