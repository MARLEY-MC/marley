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

void print_events(const char* filename) {

  TFile* file = new TFile(filename, "read");
  TTree* tree = NULL;
  file->GetObject("MARLEY Event Tree", tree);
  if (!tree) {
    std::cout << "MARLEY event tree not found" << '\n';
    return;
  }

  marley::Event* ev = new marley::Event;
  tree->SetBranchAddress("events", &ev);

  marley::Particle* fp = NULL;

  size_t num_entries = tree->GetEntries();
  for (size_t i = 0; i < num_entries; ++i) {

    tree->GetEntry(i);

    print_event_info(*ev, i);
  }
}
