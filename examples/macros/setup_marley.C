void setup_marley() {

  int return_code = gSystem->Load("libMARLEY_ROOT");
  if (return_code == 0) std::cout << "\nSuccessfully loaded MARLEY_ROOT"
    " dynamic library.\n";
  else std::cout << "\nError loading MARLEY_ROOT dynamic library.\n";

  // Include the appropriate headers if we're using ROOT 6+.
  if (gROOT->GetVersionInt() >= 60000) {
    std::cout << "ROOT 6 or greater detected. Including MARLEY headers...\n";
    gInterpreter->ProcessLine("#include \"marley/Particle.hh\"");
    gInterpreter->ProcessLine("#include \"marley/Event.hh\"");
  }
}
