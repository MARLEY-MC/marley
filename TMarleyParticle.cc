#include "TMarleyParticle.hh"

//#ifdef USE_ROOT
//#ifndef __CINT__
//ClassImp(TMarleyParticle);
//#endif
//#endif

TMarleyParticle::TMarleyParticle() {
  particle_id = 0;
  total_energy = 0.0;
  px = 0.0;
  py = 0.0;
  pz = 0.0;
  mass = 0.0;
}

TMarleyParticle::TMarleyParticle(int id, double E, double p_x,
  double p_y, double p_z, double m)
{
  particle_id = id;
  total_energy = E;
  px = p_x;
  py = p_y;
  pz = p_z;
  mass = m;
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

std::vector<TMarleyParticle*>* TMarleyParticle::get_children() {
  return &children;
}

void TMarleyParticle::add_child(TMarleyParticle* child) {
  children.push_back(child);
}
