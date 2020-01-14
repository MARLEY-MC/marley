#include <iostream>

// These functions largely duplicate code that also exists
// in marley::Event::print_human_readable(). They are included
// here as an example of code that accesses all data members
// of the Particle and Event objects.
void print_particle_info(const marley::Particle& p) {
  std::cout << "  particle with PDG code = " << p.pdg_code()
    << " has total energy " << p.total_energy() << " MeV,"
    << '\n' << "    3-momentum = (" << p.px() << " MeV, " << p.py()
    << " MeV, " << p.pz() << " MeV)," << '\n'
    << "    mass = " << p.mass() << " MeV, and charge = "
    << p.charge() << " times the proton charge." << '\n';
}

void print_event_info(const marley::Event& e, const size_t num) {
  std::cout << "\n*** Event " << num << " has "
    << e.get_initial_particles().size()
    << " initial particles and " << e.get_final_particles().size()
    << " final particles. ***" << '\n';

  int twoJ = e.twoJ();
  bool twoJ_is_odd = ( twoJ % 2 == 1 );
  std::cout << "The residual nucleus initially had excitation energy "
    << e.Ex() << " MeV and spin-parity ";
  if ( twoJ_is_odd ) std::cout << twoJ << "/2";
  else std::cout << twoJ / 2;
  std::cout << e.parity() << '\n';

  std::cout << "Initial particles" << '\n';
  marley::Particle* p = NULL;
  for (size_t i = 0; i < e.get_initial_particles().size(); ++i) {
    print_particle_info(*e.get_initial_particles().at(i));
  }
  std::cout << "Final particles" << '\n';
  for (size_t i = 0; i < e.get_final_particles().size(); ++i) {
    print_particle_info(*e.get_final_particles().at(i));
  }
}

void print_events(const std::string& file_name) {

  marley::MacroEventFileReader reader( file_name );
  marley::Event ev;

  int e = 0;

  while ( reader >> ev ) {
    print_event_info( ev, e );
    ++e;
  }
}
