#include <iostream>

#include "marley/Event.hh"
#include "marley/EventFileReader.hh"

int main(int argc, char** argv) {

  if ( argc < 2 ) return 1;

  std::string input_file_name( argv[1] );

  marley::EventFileReader efr( input_file_name );
  marley::Event event;

  double avg_xsec = efr.flux_averaged_xsec();
  std::cout << "flux-averaged total cross section = "
    << avg_xsec << " * 10^{-42} cm^2 / atom\n";

  while ( efr >> event ) {
    std::cout << event << '\n';
  }

  return 0;
}
