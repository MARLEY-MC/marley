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

#include <algorithm> // std::iter_swap
#include <fstream>
#include <iostream>
#include <iomanip>

#include "marley/Error.hh"
#include "marley/Event.hh"
#include "marley/JSON.hh"
#include "marley/MassTable.hh"
#include "marley/marley_utils.hh"

// Local constants used only within this file
namespace {

  // Indices of the projectile and target in the vector of initial particles
  constexpr size_t PROJECTILE_INDEX = 0u;
  constexpr size_t TARGET_INDEX = 1u;

  // Indices of the ejectile and residue in the vector of final particles
  constexpr size_t EJECTILE_INDEX = 0u;
  constexpr size_t RESIDUE_INDEX = 1u;

  // Conversion factor for converting GeV to MeV (the latter of which
  // is used in MARLEY natural units)
  constexpr double GEV_TO_MEV = 1000.;

  // Helper function for Event::print_human_readable()
  void print_particle_info( std::ostream& out, const marley::Particle& p ) {
    out << "  particle with PDG code = " << p.pdg_code()
      << " has total energy " << p.total_energy() << " MeV,"
      << '\n' << "    3-momentum = (" << p.px() << " MeV, " << p.py()
      << " MeV, " << p.pz() << " MeV)," << '\n'
      << "    mass = " << p.mass() << " MeV, and charge = "
      << p.charge() << " times the proton charge." << '\n';
  }

}

// Creates an empty 2-->2 scattering event with dummy initial and final
// particles. The residue (particle d) has excitation energy Ex and
// spin-parity 0+.
marley::Event::Event(double Ex)
  : initial_particles_{new marley::Particle(), new marley::Particle()},
  final_particles_{new marley::Particle(), new marley::Particle()},
  Ex_(Ex), twoJ_(0), parity_(true) {}

// Creates an 2-->2 scattering event with given initial (a & b) and final
// (c & d) particles. The residue (particle d) has excitation energy Ex,
// spin twoJ divided by two, and parity P immediately following the 2-->2
// reaction.
marley::Event::Event(const marley::Particle& a, const marley::Particle& b,
  const marley::Particle& c, const marley::Particle& d, double Ex, int twoJ,
  const marley::Parity& P)
  : initial_particles_{new marley::Particle(a), new marley::Particle(b)},
  final_particles_{new marley::Particle(c), new marley::Particle(d)},
  Ex_(Ex), twoJ_(twoJ), parity_(P) {}

// Destructor
marley::Event::~Event() {
  this->delete_particles();
}

// Copy constructor
marley::Event::Event(const marley::Event& other_event)
  : initial_particles_(other_event.initial_particles_.size()),
  final_particles_(other_event.final_particles_.size()),
  Ex_(other_event.Ex_), twoJ_(other_event.twoJ_),
  parity_(other_event.parity_)
{
  for (size_t i = 0; i < other_event.initial_particles_.size(); ++i) {
    initial_particles_[i] = new marley::Particle(
      *other_event.initial_particles_[i]);
  }
  for (size_t i = 0; i < other_event.final_particles_.size(); ++i) {
    final_particles_[i] = new marley::Particle(
      *other_event.final_particles_[i]);
  }
}


// Move constructor
marley::Event::Event(marley::Event&& other_event)
  : initial_particles_(other_event.initial_particles_.size()),
  final_particles_(other_event.final_particles_.size()),
  Ex_(other_event.Ex_), twoJ_(other_event.twoJ_),
  parity_(other_event.parity_)
{
  other_event.Ex_ = 0.;
  for (size_t i = 0; i < other_event.initial_particles_.size(); ++i) {
    initial_particles_[i] = other_event.initial_particles_[i];
  }
  for (size_t i = 0; i < other_event.final_particles_.size(); ++i) {
    final_particles_[i] = other_event.final_particles_[i];
  }
  other_event.initial_particles_.clear();
  other_event.final_particles_.clear();
}

// Copy assignment operator
marley::Event& marley::Event::operator=(const marley::Event& other_event) {
  Ex_ = other_event.Ex_;
  twoJ_ = other_event.twoJ_;
  parity_ = other_event.parity_;

  // Delete the old particle objects owned by this event
  this->delete_particles();

  initial_particles_.resize(other_event.initial_particles_.size());
  final_particles_.resize(other_event.final_particles_.size());

  for (size_t i = 0; i < other_event.initial_particles_.size(); ++i) {
    initial_particles_[i] = new marley::Particle(
      *other_event.initial_particles_[i]);
  }
  for (size_t i = 0; i < other_event.final_particles_.size(); ++i) {
    final_particles_[i] = new marley::Particle(
      *other_event.final_particles_[i]);
  }

  return *this;
}

// Move assignment operator
marley::Event& marley::Event::operator=(marley::Event&& other_event) {

  Ex_ = other_event.Ex_;
  other_event.Ex_ = 0.;

  twoJ_ = other_event.twoJ_;
  other_event.twoJ_ = 0;

  parity_ = other_event.parity_;
  other_event.parity_ = marley::Parity( true );

  // Delete the old particle objects owned by this event
  this->delete_particles();

  initial_particles_.resize(other_event.initial_particles_.size());
  final_particles_.resize(other_event.final_particles_.size());

  for (size_t i = 0; i < other_event.initial_particles_.size(); ++i) {
    initial_particles_[i] = other_event.initial_particles_[i];
  }
  for (size_t i = 0; i < other_event.final_particles_.size(); ++i) {
    final_particles_[i] = other_event.final_particles_[i];
  }

  other_event.initial_particles_.clear();
  other_event.final_particles_.clear();

  return *this;
}

marley::Particle& marley::Event::projectile() {
  return *initial_particles_.at(PROJECTILE_INDEX);
}

marley::Particle& marley::Event::target() {
  return *initial_particles_.at(TARGET_INDEX);
}

marley::Particle& marley::Event::ejectile() {
  return *final_particles_.at(EJECTILE_INDEX);
}

marley::Particle& marley::Event::residue() {
  return *final_particles_.at(RESIDUE_INDEX);
}

const marley::Particle& marley::Event::projectile() const {
  return *initial_particles_.at(PROJECTILE_INDEX);
}

const marley::Particle& marley::Event::target() const {
  return *initial_particles_.at(TARGET_INDEX);
}

const marley::Particle& marley::Event::ejectile() const {
  return *final_particles_.at(EJECTILE_INDEX);
}

const marley::Particle& marley::Event::residue() const {
  return *final_particles_.at(RESIDUE_INDEX);
}

void marley::Event::add_initial_particle(const marley::Particle& p)
{
  initial_particles_.push_back(new marley::Particle(p));
}

void marley::Event::add_final_particle(const marley::Particle& p)
{
  final_particles_.push_back(new marley::Particle(p));
}

void marley::Event::clear() {
  this->delete_particles();
  Ex_ = 0.;
  twoJ_ = 0;
  parity_ = marley::Parity( true );
}

void marley::Event::delete_particles() {
  for (auto& p : initial_particles_) if (p) delete p;
  for (auto& p : final_particles_) if (p) delete p;
  initial_particles_.clear();
  final_particles_.clear();
}

void marley::Event::print(std::ostream& out) const {
  // Use an temporary ostringstream object so that we can ensure all
  // floating-point values are output with full precision without disturbing
  // the user's settings in the "out" stream.
  std::ostringstream temp;
  temp << std::scientific;
  temp.precision(std::numeric_limits<double>::max_digits10);

  temp << initial_particles_.size() << ' ' << final_particles_.size()
    << ' ' << Ex_ << ' ' << twoJ_ << ' ' << parity_ << '\n';

  for (const auto i : initial_particles_) temp << *i << '\n';
  for (const auto f : final_particles_) temp << *f << '\n';

  out << temp.str();
}

void marley::Event::read(std::istream& in) {
  this->clear();

  int num_initial;
  int num_final;

  in >> num_initial >> num_final >> Ex_ >> twoJ_ >> parity_;

  // If reading the event header line failed for some
  // reason, just return without trying to do anything else.
  if ( !in ) return;

  // If we have invalid numbers of particles in either the initial or final
  // state (there need to be at least two in each and we can't fill memory when
  // allocating space for them) then complain by throwing an error
  if ( num_initial < 2 || static_cast<size_t>(num_initial)
    > initial_particles_.max_size() )
  {
    throw marley::Error("Invalid number of initial state particles ("
      + std::to_string(num_initial) + ") encountered in marley::Event::read()");
  }
  else if ( num_final < 2 || static_cast<size_t>(num_final)
    > final_particles_.max_size() )
  {
    throw marley::Error("Invalid number of final state particles ("
      + std::to_string(num_final) + ") encountered in marley::Event::read()");
  }

  initial_particles_.resize( num_initial );
  final_particles_.resize( num_final );

  for (int i = 0; i < num_initial; ++i) {
    marley::Particle* p = new marley::Particle;
    in >> *p;
    if ( !in ) throw marley::Error("Parse error while reading initial"
      " particle #" + std::to_string(i) + " from an ASCII-format event"
      " record");
    initial_particles_.at(i) = p;
  }
  for (int f = 0; f < num_final; ++f) {
    marley::Particle* p = new marley::Particle;
    in >> *p;
    if ( !in ) throw marley::Error("Parse error while reading final"
      " particle #" + std::to_string(f) + " from an ASCII-format event"
      " record");
    final_particles_.at(f) = p;
  }
}

// Function that dumps a marley::Particle to an output stream in HEPEVT format.
// This is a private helper function for the publicly-accessible write_hepevt.
void marley::Event::dump_hepevt_particle(const marley::Particle& p,
  std::ostream& os, int status, int jmohep1, int jmohep2) const
{
  // Print the status code to begin the particle entry in the HEPEVT record
  os << status << ' ';

  // TODO: improve this entry to give the user more control over the vertex
  // location and to reflect the parent-daughter relationships between
  // particles.
  // Convert from MARLEY natural units (MeV) to GeV for the HEPEVT format
  os << p.pdg_code() << ' ' << jmohep1 << ' ' << jmohep2 << " 0 0 "
    << p.px() / GEV_TO_MEV << ' ' << p.py() / GEV_TO_MEV
    << ' ' << p.pz() / GEV_TO_MEV << ' ' << p.total_energy() / GEV_TO_MEV
    << ' ' << p.mass() / GEV_TO_MEV
    // The spacetime origin is currently used as the initial position 4-vector
    // for all particles
    << " 0. 0. 0. 0." << '\n';
}

void marley::Event::write_hepevt(size_t event_num, double flux_avg_tot_xsec,
  std::ostream& out) const
{
  // Use an temporary ostringstream object so that we can ensure all
  // floating-point values are output with full precision without disturbing
  // the user's settings in the "out" stream.
  std::ostringstream temp;
  temp.precision( std::numeric_limits<double>::max_digits10 );
  temp << std::scientific;

  // Create a dummy particle that encodes extra MARLEY-specific information in
  // the HEPEVT format. Preserve MARLEY natural units for these quantities
  // (MeV) by pre-multiplying by the conversion factor used in
  // dump_hepevt_particle()
  marley::Particle dummy_particle;
  dummy_particle.set_total_energy( Ex_ * GEV_TO_MEV );
  dummy_particle.set_mass( flux_avg_tot_xsec * GEV_TO_MEV );

  // Add one to the total particle count so that our dummy particle will be
  // included correctly
  size_t num_particles = 1 + initial_particles_.size()
    + final_particles_.size();

  // Write the HEPEVT header line to the event record
  temp << event_num  << ' ' << num_particles << '\n';

  // Write the initial particles to the event record
  for (const auto i : initial_particles_) dump_hepevt_particle(*i, temp,
    HEPEVT_INITIAL_STATE_STATUS_CODE);

  // Write our dummy particle to the event record
  // We use the jmohep1 slot to record the twoJ_ data member and the jmohep2
  // slot to record the parity_ data member (both as integers)
  dump_hepevt_particle( dummy_particle, temp, HEPEVT_MARLEY_INFO_STATUS_CODE,
    twoJ_, static_cast<int>(parity_) );

  // Write the final particles to the event record
  for (const auto f : final_particles_) dump_hepevt_particle(*f, temp,
    HEPEVT_FINAL_STATE_STATUS_CODE);

  // Output the finished HEPEVT format event to the "out" stream
  out << temp.str();
}

marley::JSON marley::Event::to_json() const {
  marley::JSON event = marley::JSON::object();

  event["Ex"] = Ex_;
  event["twoJ"] = twoJ_;
  event["parity"] = static_cast<int>( parity_ );
  event["initial_particles"] = marley::JSON::array();
  event["final_particles"] = marley::JSON::array();

  for (const auto ip : initial_particles_)
    event.at("initial_particles").append(ip->to_json());

  for (const auto fp : final_particles_)
    event.at("final_particles").append(fp->to_json());

  return event;
}

// Reconstructs a marley::Event object from an input HEPEVT-format event record
// in an input stream. Returns a boolean that reflects whether the input stream
// state was still good at the end of the attempt to read the HEPEVT record.
// This can be used as a while loop condition when reading in multiple HEPEVT
// records.
bool marley::Event::read_hepevt(std::istream& in, double* flux_avg_tot_xsec)
{
  // If flux_avg_tot_xsec is not null, then clear any previous value there
  // before continuing. It will be loaded with a new value below if we
  // find a dummy particle encoding extra MARLEY-specific information in
  // the input HEPEVT record.
  if ( flux_avg_tot_xsec ) *flux_avg_tot_xsec = 0.;

  // Remove any pre-existing contents in this event
  this->clear();

  // Check that the event contains exactly two initial-state particles
  size_t initial_state_particles = 0u;
  // Check that one of the initial-state particles is an atom (ion)
  size_t initial_state_ions = 0u;
  // Check that the event contains exactly one final-state lepton
  size_t final_state_leptons = 0u;
  // Check that the event contains at least one final-state ion
  size_t final_state_ions = 0u;

  // Initial indices of the final lepton and largest final ion (as judged
  // by the mass number A) in the final_particles_
  // vector. These will be used below to move them into the correct order.
  size_t final_lepton_idx = 0u;
  size_t largest_final_ion_idx = 0u;
  int largest_A = 0;

  int event_num; // HEPEVT event number (ignored)
  int num_particles; // Total number of particles stored in the event
  in >> event_num  >> num_particles;

  // If reading the first line failed for some
  // reason, just return without trying to do anything else.
  if ( !in ) return false;

  // If the number of particles in the event is negative or extremely
  // huge (so that it would overload memory) then throw a marley::Error
  if ( num_particles < 0 || static_cast<size_t>(num_particles)
    > initial_particles_.max_size() )
  {
    throw marley::Error("Invalid number of particles ("
      + std::to_string(num_particles) + ") encountered in marley::"
      "Event::read_hepevt()");
  }

  // Dummy variables used to skip HEPEVT fields that MARLEY doesn't
  // care about (e.g., the particle production vertex 4-position)
  int dummy_int;
  double dummy_double;

  // Fields that we do care about
  double Etot, px, py, pz, M;
  int status_code, pdg, jmohep1, jmohep2;
  for ( int p = 0; p < num_particles; ++p ) {

    in >> status_code >> pdg >> jmohep1 >> jmohep2;

    // Skip the JDAHEP1 and JDAHEP2 fields (daughter indices)
    for ( int j = 0; j < 2; ++j ) in >> dummy_int;

    // Read in the particle 4-momentum and mass.
    in >> px >> py >> pz >> Etot >> M;

    // Skip the VHEP1 through VHEP4 fields (production vertex 4-position)
    for ( int j = 0; j < 4; ++j ) in >> dummy_double;

    if ( !in ) throw marley::Error("Parse error while reading  particle #"
      + std::to_string(p) + " from a HEPEVT-format event record");

    // If the particle has this status code, it is a dummy particle that
    // contains MARLEY-specific information
    if ( status_code == HEPEVT_MARLEY_INFO_STATUS_CODE ) {
      // The dummy particle total energy stores the nuclear excitation energy.
      // This value can be reconstructed (with some effort) from the
      // particle four-vectors, but just storing it is a lot easier.
      Ex_ = Etot;
      // The dummy particle mass stores the flux-averaged total cross section.
      // If the user passed something other than a nullptr as the
      // flux_avg_tot_xsec argument to this function, then load it
      // with this value, which is not normally stored in the marley::Event
      // itself.
      if ( flux_avg_tot_xsec ) *flux_avg_tot_xsec = M;
      // The JMOHEP1 field contains two times the residue spin immediately
      // following the initial two-two scattering reaction
      twoJ_ = jmohep1;
      // The JMOHEP2 field contains an integer representation of the residue
      // parity immediately following the initial two-two scattering reaction
      parity_ = jmohep2;
    }

    // If the particle has a status code other than the two used by
    // MARLEY, then don't bother to record the particle in the event object
    if ( status_code != HEPEVT_INITIAL_STATE_STATUS_CODE
      && status_code != HEPEVT_FINAL_STATE_STATUS_CODE ) continue;

    // Convert the particle 4-momentum components and mass from GeV (used in
    // the HEPEVT format) to MARLEY natural units (MeV)
    px *= GEV_TO_MEV;
    py *= GEV_TO_MEV;
    pz *= GEV_TO_MEV;
    Etot *= GEV_TO_MEV;
    M *= GEV_TO_MEV;

    // We already checked that the particle has a status code recognized
    // by MARLEY above, so store it in the event object in the appropriate
    // place.
    if ( status_code == HEPEVT_INITIAL_STATE_STATUS_CODE ) {
      initial_particles_.push_back(
        new marley::Particle(pdg, Etot, px, py, pz, M) );
      ++initial_state_particles;
      if ( marley_utils::is_ion(pdg) ) ++initial_state_ions;
    }
    else {
      // status_code == HEPEVT_FINAL_STATE_STATUS_CODE
      final_particles_.push_back(
        new marley::Particle(pdg, Etot, px, py, pz, M) );
      if ( marley_utils::is_lepton(pdg) ) {
        ++final_state_leptons;
        final_lepton_idx = final_particles_.size() - 1;
      }
      if ( marley_utils::is_ion(pdg) ) {
        ++final_state_ions;
        // Update the index of the largest final ion
        int A = marley_utils::get_particle_A( pdg ); // mass number
        if ( A > largest_A ) {
          largest_A = A;
          largest_final_ion_idx = final_particles_.size() - 1;
        }
      }
    }
  }

  // Check that the event satisfies the criteria needed for a reasonable
  // marley::Event
  if ( initial_state_particles != 2u ) throw marley::Error(
    std::to_string(initial_state_particles) + " initial-state particles"
    " encountered while reading HEPEVT input");
  else if ( initial_state_ions != 1u ) throw marley::Error("Multiple"
    " initial-state ions encountered while reading HEPEVT input");
  else if ( final_state_leptons != 1u ) throw marley::Error("Multiple"
    " final-state leptons encountered while reading HEPEVT input");
  else if ( final_state_ions < 1u ) throw marley::Error("Zero"
    " final-state ions encountered while reading HEPEVT input");

  // Reorder (if needed) the initial particle list to ensure that the target
  // (the initial atom) appears in the correct position. This is required
  // because marley::Event objects rely on the ordering of the
  // initial_particles_ array to allow easy retrieval of the projectile and
  // target.
  auto target_iter = initial_particles_.begin() + TARGET_INDEX;
  if ( !marley_utils::is_ion( (*target_iter)->pdg_code() ) ) {
    auto begin_ip = initial_particles_.begin();
    std::iter_swap( begin_ip, begin_ip + 1 );
  }

  // Play the same sort of game with the final particles. This will ensure that
  // the ejectile and residue are stored in the correct places.
  auto begin_fp = final_particles_.begin();
  if ( final_lepton_idx != EJECTILE_INDEX ) {
    std::iter_swap( begin_fp + EJECTILE_INDEX, begin_fp + final_lepton_idx );
  }
  if ( largest_final_ion_idx != RESIDUE_INDEX ) {
    std::iter_swap( begin_fp + RESIDUE_INDEX, begin_fp + largest_final_ion_idx );
  }

  // The marley::Particle constructor that we used above assumes that all ion
  // PDG codes passed to it represent fully ionized particles. Here we will
  // override this behavior and assume instead that the target atom has
  // zero net charge in the initial state.
  this->target().set_charge( 0 );

  // Determine the correct ionization state immediately following the prompt
  // 2-->2 scattering reaction. Do this by assuming that the target atom is
  // neutral (see above). This implies that the final ion charge must be equal
  // and opposite to that of the final lepton.
  int Qf_ion = -this->ejectile().charge();

  // Get the nuclear residue net electric charge by summing all final-state
  // ions besides the residue itself. The difference between this sum and the
  // pre-de-excitation ion charge is the charge that should be used.
  int sum_Q_final_ions = 0;
  for ( size_t f = 0u; f < final_particles_.size(); ++f ) {
    if ( f == RESIDUE_INDEX ) continue; // skip the residue

    const auto& fp = final_particles_.at( f );
    int pdg_f = fp->pdg_code();
    if ( marley_utils::is_ion(pdg_f) ) sum_Q_final_ions += fp->charge();
  }
  this->residue().set_charge( Qf_ion - sum_Q_final_ions );

  // Everything was read in successfully. We're done!
  return true;
}

void marley::Event::from_json(const marley::JSON& json) {

  // TODO: reduce code duplication in this function

  // Remove any existing contents from this event
  this->clear();

  if ( !json.has_key("Ex") ) throw marley::Error("Missing"
    " nuclear excitation energy key in input JSON-format event");

  bool ok = false;
  const auto& temp_Ex = json.at("Ex");
  Ex_ = temp_Ex.to_double( ok );
  if ( !ok ) throw marley::Error("Invalid nuclear excitation energy"
    + temp_Ex.to_string() + " encountered in input JSON-format event");

  if ( !json.has_key("twoJ") ) throw marley::Error("Missing"
    " twoJ key in input JSON-format event");

  ok = false;
  const auto& temp_twoJ = json.at("twoJ");
  twoJ_ = temp_twoJ.to_long( ok );
  if ( !ok ) throw marley::Error("Invalid \"twoJ\" value"
    + temp_twoJ.to_string() + " encountered in input JSON-format event");

  ok = false;
  const auto& temp_parity = json.at("parity");
  parity_ = temp_parity.to_long( ok );
  if ( !ok ) throw marley::Error("Invalid parity value"
    + temp_parity.to_string() + " encountered in input JSON-format event");

  // Retrieve and load the array of initial particles
  if ( !json.has_key("initial_particles") ) throw marley::Error("Missing"
    " initial particle array in input JSON-format event");

  const auto& temp_ip = json.at("initial_particles");
  if ( !temp_ip.is_array() ) throw marley::Error("Invalid initial"
    " particle array " + temp_ip.to_string() + " encountered in input"
    " JSON-format event");

  for ( const auto& p_object : temp_ip.array_range() ) {
    if ( !p_object.is_object() ) throw marley::Error("Invalid particle"
      " object " + p_object.to_string() + " encountered while parsing a"
      " JSON-format particle array");
    initial_particles_.push_back( new marley::Particle );
    initial_particles_.back()->from_json( p_object );
  }

  // Retrieve and load the array of final particles
  if ( !json.has_key("final_particles") ) throw marley::Error("Missing"
    " final particle array in input JSON-format event");

  const auto& temp_fp = json.at("final_particles");
  if ( !temp_fp.is_array() ) throw marley::Error("Invalid final"
    " particle array " + temp_fp.to_string() + " encountered in input"
    " JSON-format event");

  for ( const auto& p_object : temp_fp.array_range() ) {
    if ( !p_object.is_object() ) throw marley::Error("Invalid particle"
      " object " + p_object.to_string() + " encountered while parsing a"
      " JSON-format particle array");
    final_particles_.push_back( new marley::Particle );
    final_particles_.back()->from_json( p_object );
  }

}

void marley::Event::print_human_readable(std::ostream& out, int num) const {
  out << "\n*** MARLEY Event ";
  if ( num >= 0 ) out << num << ' ';
  out << "has " << this->get_initial_particles().size()
    << " initial particles and " << this->get_final_particles().size()
    << " final particles. ***" << '\n';

  int twoJ = this->twoJ();
  bool twoJ_is_odd = ( twoJ % 2 == 1 );
  out << "The residual nucleus initially had excitation energy "
    << this->Ex() << " MeV and spin-parity ";
  if ( twoJ_is_odd ) out << twoJ << "/2";
  else out << twoJ / 2;
  out << this->parity() << '\n';

  out << "Initial particles" << '\n';
  marley::Particle* p = NULL;
  for (size_t i = 0u; i < this->get_initial_particles().size(); ++i) {
    print_particle_info( out, *this->get_initial_particles().at(i) );
  }
  out << "Final particles" << '\n';
  for (size_t i = 0u; i < this->get_final_particles().size(); ++i) {
    print_particle_info( out, *this->get_final_particles().at(i) );
  }
}
