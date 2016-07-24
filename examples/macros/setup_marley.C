void setup_marley() {

  int return_code = gSystem->Load("libMARLEY_ROOT");
  if (return_code == 0) std::cout << "\nSuccessfully loaded MARLEY_ROOT"
    " dynamic library.\n";
  else std::cout << "\nError loading MARLEY_ROOT dynamic library.\n"
    << "Please add the directory that contains this library to your\n"
    << "LD_LIBRARY_PATH environment variable (or equivalent on\n"
    << "non-Linux systems) and try again.\n";

  // Include the appropriate headers if we're using ROOT 6+.
  if (gROOT->GetVersionInt() >= 60000) {
    std::cout << "ROOT 6 or greater detected. Including MARLEY headers...\n";
    TInterpreter::EErrorCode* ec = new TInterpreter::EErrorCode();
    gInterpreter->ProcessLine("#include \"marley/Particle.hh\"", ec);
    if (*ec != 0) std::cout << "Error loading"
      << " MARLEY header Particle.hh. For MARLEY headers stored in"
      << " /path/to/include/marley/, please add /path/to/include"
      << " to your ROOT_INCLUDE_PATH environment variable and"
      << " try again.\n";
    gInterpreter->ProcessLine("#include \"marley/Event.hh\"", ec);
    if (*ec != 0) std::cout << "Error loading"
      << " MARLEY header Event.hh. For MARLEY headers stored in"
      << " /path/to/include/marley/, please add /path/to/include"
      << " to your ROOT_INCLUDE_PATH environment variable and"
      << " try again.\n";
  }
}
