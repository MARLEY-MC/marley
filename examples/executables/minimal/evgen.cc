// Standard library includes
#include <iostream>

// MARLEY includes
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/JSONConfig.hh"

constexpr int NUM_EVENTS = 10;

int main() {

  marley::JSONConfig cfg( "/home/config.js" );
  marley::Generator gen = cfg.create_generator();

  for ( int j = 0; j < NUM_EVENTS; ++j ) {
    marley::Event ev = gen.create_event();
    std::cout << ev << '\n';
  }

}
