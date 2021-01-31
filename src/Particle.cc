/// @file
/// @copyright Copyright (C) 2016-2021 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#include <cmath>

#include "marley/marley_utils.hh"
#include "marley/JSON.hh"
#include "marley/Particle.hh"

namespace {
  // Helper functions for converting a JSON object into a marley::Particle
  // TODO: reduce code duplication between read_JSON_double() and
  // read_JSON_integer()
  double read_JSON_double( const std::string& key, const marley::JSON& json ) {

    if ( !json.has_key(key) ) throw marley::Error("Missing key " + key
      + " encountered in an input JSON particle object");

    bool ok = false;
    const auto& temp = json.at( key );
    double result = temp.to_double( ok );
    if ( !ok ) throw marley::Error("Invalid value " + temp.to_string()
      + " associated with the " + key + " key encountered in an input"
      " JSON particle object");
    return result;
  }

  long read_JSON_integer( const std::string& key, const marley::JSON& json ) {

    if ( !json.has_key(key) ) throw marley::Error("Missing key " + key
      + " encountered in an input JSON particle object");

    bool ok = false;
    const auto& temp = json.at( key );
    long result = temp.to_long( ok );
    if ( !ok ) throw marley::Error("Invalid value " + temp.to_string()
      + " associated with the " + key + " key encountered in an input"
      " JSON particle object");
    return result;
  }

}

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

void marley::Particle::clear() {
  for (size_t j = 0u; j < 4u; ++j ) four_momentum_[j] = 0.;
  pdg_code_ = 0;
  mass_ = 0.;
  charge_ = 0;
}

void marley::Particle::from_json(const marley::JSON& json) {
  this->clear();
  pdg_code_ = read_JSON_integer( "pdg", json );
  charge_ = read_JSON_integer( "charge", json );
  four_momentum_[0] = read_JSON_double( "E", json );
  four_momentum_[1] = read_JSON_double( "px", json );
  four_momentum_[2] = read_JSON_double( "py", json );
  four_momentum_[3] = read_JSON_double( "pz", json );
  mass_ = read_JSON_double( "mass", json );
}
