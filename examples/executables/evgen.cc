// Standard library includes
#include <iostream>
#include <string>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/PrintStreams.h"

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"

constexpr int NUM_EVENTS = 10;

int main( int argc, char* argv[] ) {

  if ( argc != 2 ) {
    std::cout << "Usage: " << argv[0] << " CONFIG_FILE\n";
    return 1;
  }

  std::string config_file( argv[1] );

  marley::JSONConfig cfg( config_file );
  marley::Generator gen = cfg.create_generator();
  gen.set_up_run_info();

  for ( int j = 0; j < NUM_EVENTS; ++j ) {
    auto ev = gen.create_event();
    std::cout << *ev << '\n';
  }

}
