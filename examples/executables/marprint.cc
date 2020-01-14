// Standard library includes
#include <iostream>
#include <string>

// MARLEY includes
#include "marley/Event.hh"
#ifdef USE_ROOT
  #include "marley/RootEventFileReader.hh"
#else
  #include "marley/EventFileReader.hh"
#endif

// Prints events from an input file in a human-readable format. This function
// also updates an overall event number as it works through the file.
void print_events(const std::string& file_name, int& ev_num) {

  #ifdef USE_ROOT
    marley::RootEventFileReader reader( file_name );
  #else
    marley::EventFileReader reader( file_name );
  #endif

  marley::Event ev;

  while ( reader >> ev ) {
    ev.print_human_readable( std::cout, ev_num );
    ++ev_num;
  }
}

int main(int argc, char* argv[]) {

  // If the user has not supplied any command-line arguments, display the
  // standard help message and exit
  if (argc <= 1) {
    std::cout << "Usage: " << argv[0] << " INPUT_FILE...\n";
    return 0;
  }

  std::vector<std::string> file_names;
  for ( int s = 1; s < argc; ++s ) file_names.push_back( argv[s] );

  int event_number = 0;
  for ( const auto& file_name : file_names ) {
    print_events( file_name, event_number );
  }

  return 0;
}
