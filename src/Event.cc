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
}

// Creates an empty 2->2 scattering event with dummy initial and final
// particles and residue (particle d) excitation energy Ex.
marley::Event::Event(double Ex)
  : initial_particles_{new marley::Particle(), new marley::Particle()},
  final_particles_{new marley::Particle(), new marley::Particle()},
  Ex_(Ex) {}

// Creates an 2->2 scattering event with given initial (a & b) and final
// (c & d) particles. The residue (particle d) has excitation energy Ex
// immediately following the 2->2 reaction.
marley::Event::Event(const marley::Particle& a, const marley::Particle& b,
  const marley::Particle& c, const marley::Particle& d, double Ex)
  : initial_particles_{new marley::Particle(a), new marley::Particle(b)},
  final_particles_{new marley::Particle(c), new marley::Particle(d)},
  Ex_(Ex) {}

// Destructor
marley::Event::~Event() {
  this->delete_particles();
}

// Copy constructor
marley::Event::Event(const marley::Event& other_event)
  : initial_particles_(other_event.initial_particles_.size()),
  final_particles_(other_event.final_particles_.size()),
  Ex_(other_event.Ex_)
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
  Ex_(other_event.Ex_)
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
    << ' ' << Ex_ << '\n';

  for (const auto i : initial_particles_) temp << *i << '\n';
  for (const auto f : final_particles_) temp << *f << '\n';

  out << temp.str();
}

void marley::Event::read(std::istream& in) {
  this->clear();

  size_t num_initial;
  size_t num_final;

  in >> num_initial >> num_final >> Ex_;

  initial_particles_.resize(num_initial);
  final_particles_.resize(num_final);

  for (size_t i = 0; i < num_initial; ++i) {
    marley::Particle* p = new marley::Particle;
    in >> *p;
    initial_particles_.at(i) = p;
  }
  for (size_t f = 0; f < num_final; ++f) {
    marley::Particle* p = new marley::Particle;
    in >> *p;
    final_particles_.at(f) = p;
  }
}

// Function that dumps a marley::Particle to an output stream in HEPEVT format.
// This is a private helper function for the publicly-accessible write_hepevt.
void marley::Event::dump_hepevt_particle(const marley::Particle& p,
  std::ostream& os, bool track) const
{
  // Print the flag that indicates whether the particle should be tracked
  // (i.e., it is a final particle) or not
  if (track) os << HEPEVT_FINAL_STATE_STATUS_CODE << ' ';
  else os << HEPEVT_INITIAL_STATE_STATUS_CODE << ' ';

  // TODO: improve this entry to give the user more control over the vertex
  // location and to reflect the parent-daughter relationships between
  // particles.
  // Convert from MARLEY natural units (MeV) to GeV for the HEPEVT format
  os << p.pdg_code() << " 0 0 0 0 " << p.px() / GEV_TO_MEV
    << ' ' << p.py() / GEV_TO_MEV << ' ' << p.pz() / GEV_TO_MEV
    << ' ' << p.total_energy() / GEV_TO_MEV << ' ' << p.mass() / GEV_TO_MEV
    // Spacetime origin is currently used as the initial position 4-vector for
    // all particles
    << " 0. 0. 0. 0." << '\n';
}

void marley::Event::write_hepevt(size_t event_num, std::ostream& out) const {

  // Use an temporary ostringstream object so that we can ensure all
  // floating-point values are output with full precision without disturbing
  // the user's settings in the "out" stream.
  std::ostringstream temp;
  temp.precision(std::numeric_limits<double>::max_digits10);
  temp << std::scientific;

  size_t num_particles = initial_particles_.size() + final_particles_.size();
  temp << event_num  << ' ' << num_particles << '\n';

  for (const auto i : initial_particles_) dump_hepevt_particle(*i, temp, false);
  for (const auto f : final_particles_) dump_hepevt_particle(*f, temp, true);

  // Output the finished HEPEVT format event to the "out" stream
  out << temp.str();
}

marley::JSON marley::Event::to_json() const {
  marley::JSON event = marley::JSON::object();

  event["Ex"] = Ex_;
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
bool marley::Event::read_hepevt(std::istream& in)
{
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

  // Dummy variables used to skip HEPEVT fields that MARLEY doesn't
  // care about (e.g., the particle production vertex 4-position)
  int dummy_int;
  double dummy_double;

  // Fields that we do care about
  double Etot, px, py, pz, M;
  int status_code, pdg;
  for ( int p = 0; p < num_particles; ++p ) {

    in >> status_code >> pdg;

    // Skip JMOHEP1, JMOHEP2, JDAHEP1, JDAHEP2 fields
    // (mother and daughter indices)
    for ( int j = 0; j < 4; ++j ) in >> dummy_int;

    // Read in the particle 4-momentum and mass.
    in >> px >> py >> pz >> Etot >> M;

    // Skip the VHEP1 through VHEP4 fields (production vertex 4-position)
    for ( int j = 0; j < 4; ++j ) in >> dummy_double;

    // If the particle has a status code other than the two recognized by
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

  if ( !in ) return false;

  // Check that the event satisfies the criteria needed to compute a
  // nuclear excitation energy value Ex_
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

  // Everything is now in place except for the nuclear excitation energy.
  // We'll invoke 4-momentum conservation to compute it. First, we'll
  // determine the 4-momentum components of the hadronic system immediately
  // following the initial 2->2 scattering reaction:
  double Etot_had = this->projectile().total_energy()
    + this->target().total_energy() - this->ejectile().total_energy();
  double px_had = this->projectile().px() + this->target().px()
    - this->ejectile().px();
  double py_had = this->projectile().py() + this->target().py()
    - this->ejectile().py();
  double pz_had = this->projectile().pz() + this->target().pz()
    - this->ejectile().pz();

  // Assuming that the hadronic system is on-shell allows us to compute
  // its invariant mass:
  double M_had = marley_utils::real_sqrt( Etot_had*Etot_had - px_had*px_had
    - py_had*py_had - pz_had*pz_had );

  // The only thing left to do is to determine the ground-state mass of the
  // final ion before de-excitation. We can determine the correct ionization
  // state by invoking conservation of electric charge and looking at the final
  // lepton charge:
  int Qf_ion = -this->ejectile().charge();

  // The marley::Particle constructor that we used above assumes that all ion
  // PDG codes passed to it represent fully ionized particles. Here we will
  // override this behavior and assume instead that the target atom has
  // zero net charge in the initial state.
  this->target().set_charge( 0 );

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

  // For lepton-nucleus scattering, the initial 2->2 reaction doesn't change the
  // nuclear mass number A. The appropriate proton number Z can be deduced based
  // on the ionization state:
  int Zi = marley_utils::get_particle_Z( this->target().pdg_code() );
  int A = marley_utils::get_particle_A( this->target().pdg_code() );
  int Zf = Zi + Qf_ion;

  // We're now ready to compute the ground-state mass of the final ion
  // (neglecting electron binding energies just as we do elsewhere in MARLEY)
  const auto& mt = marley::MassTable::Instance();
  double Mhad_gs = mt.get_atomic_mass( Zf, A ) - Qf_ion*mt.get_particle_mass(
    marley_utils::ELECTRON );

  // The difference between the invariant mass and the ground-state mass is
  // the nuclear excitation energy. Store it in the event object.
  Ex_ = M_had - Mhad_gs;

  // Everything was read in successfully. We're done!
  return true;
}

void marley::Event::from_json(const marley::JSON& json) {

  // Remove any existing contents from this event
  this->clear();

  if ( !json.has_key("Ex") ) throw marley::Error("Missing"
    " nuclear excitation energy key in input JSON-format event");

  bool ok = false;
  const auto& temp_Ex = json.at("Ex");
  Ex_ = temp_Ex.to_double( ok );
  if ( !ok ) throw marley::Error("Invalid nuclear excitation energy"
    + temp_Ex.to_string() + " encountered in input JSON-format event");

  // Retrieve and load the array of initial particles
  // TODO: reduce code duplication in this function
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
