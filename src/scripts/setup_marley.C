void load_marley_header( const std::string& header_file_name ) {

  TInterpreter::EErrorCode* ec = new TInterpreter::EErrorCode();

  std::string include_statement = "#include \"marley/"
    + header_file_name + '\"';

  gInterpreter->ProcessLine( include_statement.c_str(), ec );
  if ( *ec != 0 ) {
    std::cout << "Error loading MARLEY header " << header_file_name
      << ". For MARLEY headers stored in /path/to/include/marley/,"
      << " please add /path/to/include to your ROOT_INCLUDE_PATH"
      << " environment variable and try again.\n";
  }

  delete ec;
}

void setup_marley() {

  int return_code = 0;

  // Pre-load the 2D graphics library if we're using ROOT 6+.
  // Some versions of ROOT 6 seem to require this in order to load
  // the MARLEY dictionaries correctly.
  if (gROOT->GetVersionInt() >= 60000) {
    return_code = gSystem->Load("libGraf");
    if (return_code != 0) std::cout << "\nError loading ROOT Graf"
      " library.\n";
  }

  // Load the ROOT interface for MARLEY
  return_code = gSystem->Load("libMARLEY_ROOT");
  if (return_code == 0) std::cout << "\nSuccessfully loaded MARLEY_ROOT"
    " dynamic library.\n";
  else std::cout << "\nError loading MARLEY_ROOT dynamic library.\n"
    << "Please add the directory that contains this library to your\n"
    << "LD_LIBRARY_PATH environment variable (or equivalent on\n"
    << "non-Linux systems) and try again.\n";

  // Include the appropriate headers if we're using ROOT 6+.
  if (gROOT->GetVersionInt() >= 60000) {
    std::cout << "ROOT 6 or greater detected. Including MARLEY headers...\n";
    load_marley_header( "Particle.hh" );
    load_marley_header( "Event.hh" );
    load_marley_header( "MacroEventFileReader.hh" );
    load_marley_header( "Parity.hh" );
  }
}
