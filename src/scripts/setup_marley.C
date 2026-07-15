void setup_marley() {

  int return_code = 0;

  // Load the MARLEY shared library
  return_code = gSystem->Load( "libMARLEY" );
  if ( return_code == 0 ) std::cout << "\nSuccessfully loaded MARLEY"
    " shared library.\n";
  else std::cout << "\nError loading MARLEY shared library.\n"
    << "Please add the directory that contains this library to your\n"
    << "LD_LIBRARY_PATH environment variable (or equivalent on\n"
    << "non-Linux systems) and try again.\n";

}
