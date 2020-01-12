#include <iostream>

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
  std::cout << "The residual nucleus initially had excitation energy "
    << e.Ex() << " MeV." << '\n';
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
