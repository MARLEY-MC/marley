// Standard library includes
#include <iostream>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/PrintStreams.h"

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"

constexpr int NUM_EVENTS = 10;

int main() {

  marley::JSONConfig cfg( "/home/config.js" );
  marley::Generator gen = cfg.create_generator();
  gen.set_up_run_info();

  for ( int j = 0; j < NUM_EVENTS; ++j ) {
    auto ev = gen.create_event();
    std::cout << *ev << '\n';
  }

}
