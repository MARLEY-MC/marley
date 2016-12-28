#include <cmath>

#include "marley/marley_utils.hh"
#include "marley/JSON.hh"
#include "marley/Particle.hh"

marley::Particle::Particle() : four_momentum_{0., 0., 0., 0.} {}

marley::Particle::Particle(int pdg_code, double m)
  : four_momentum_{0., 0., 0., 0.}, pdg_code_(pdg_code), mass_(m),
  charge_(marley_utils::get_particle_charge(pdg_code)) {}

marley::Particle::Particle(int pdg_code, double m, int q)
  : four_momentum_{ m, 0., 0., 0. }, pdg_code_(pdg_code), mass_(m),
  charge_(q) {}

marley::Particle::Particle(int pdg_code, double px, double py, double pz,
  double m) : four_momentum_{ 0., px, py, pz }, pdg_code_(pdg_code), mass_(m),
  charge_(marley_utils::get_particle_charge(pdg_code))
{
  double E = std::sqrt(std::pow(m, 2) + std::pow(px, 2) + std::pow(py, 2)
    + std::pow(pz, 2));
  four_momentum_[0] = E;
}

marley::Particle::Particle(int pdg_code, double px, double py, double pz,
  double m, int q) : four_momentum_{ 0., px, py, pz }, pdg_code_(pdg_code),
  mass_(m), charge_(q)
{
  double E = std::sqrt(std::pow(m, 2) + std::pow(px, 2) + std::pow(py, 2)
    + std::pow(pz, 2));
  four_momentum_[0] = E;
}

marley::Particle::Particle(int pdg_code, double E, double px,
  double py, double pz, double m) : four_momentum_{ E, px, py, pz },
  pdg_code_(pdg_code), mass_(m),
  charge_(marley_utils::get_particle_charge(pdg_code)) {}

marley::Particle::Particle(int pdg_code, double E, double px,
  double py, double pz, double m, int q) : four_momentum_{ E, px, py, pz },
  pdg_code_(pdg_code), mass_(m), charge_(q) {}

double marley::Particle::momentum_magnitude() const {
  return std::sqrt(std::pow(four_momentum_[1], 2)
    + std::pow(four_momentum_[2], 2) + std::pow(four_momentum_[3], 2));
}

double marley::Particle::kinetic_energy() const {
  return std::max(four_momentum_[0] - mass_, 0.);
}

/// @note Although marley::Event::print() guarantees full precision for
/// stream output, this function does not if it is called directly. If you
/// need to serialize individual particle objects to text and read them
/// back in again, either use an output stream with the proper precision
/// set in advance (at least std::numeric_limits<double>::max_digits10),
/// or consider using JSON input/output.
void marley::Particle::print(std::ostream& out) const {
  out << pdg_code_ << ' ' << four_momentum_[0] << ' ' << four_momentum_[1]
    << ' ' << four_momentum_[2] << ' ' << four_momentum_[3]
    << ' ' << mass_ << ' ' << charge_;
}

void marley::Particle::read(std::istream& in) {
  in >> pdg_code_ >> four_momentum_[0] >> four_momentum_[1]
    >> four_momentum_[2] >> four_momentum_[3] >> mass_ >> charge_;
}

marley::JSON marley::Particle::to_json() const {
  marley::JSON particle = marley::JSON::object();
  particle["pdg"] = pdg_code_;
  particle["E"] = four_momentum_[0];
  particle["px"] = four_momentum_[1];
  particle["py"] = four_momentum_[2];
  particle["pz"] = four_momentum_[3];
  particle["mass"] = mass_;
  particle["charge"] = charge_;
  return particle;
}
