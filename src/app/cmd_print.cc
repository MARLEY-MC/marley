// Standard library includes
#include <iostream>
#include <string>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/PrintStreams.h"

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/EventFileReader.hh"

bool marley::CommandHandler::cmd_print( std::deque< std::string >& args ) {

  std::string first_arg;
  if ( !args.empty() ) first_arg = args.front();

  if ( first_arg.empty() || first_arg == "-h" || first_arg == "--help" ) {
    args.clear();
    args.push_front( "print" );
    return marley::CommandHandler::cmd_help( args );
  }

  size_t num_files = args.size();
  for ( size_t i = 1u; i < num_files; ++i ) {
    marley::EventFileReader reader( args.at(i) );
    HepMC3::GenEvent ev;

    while ( reader >> ev ) {
      std::cout << ev;
    }
  }

  return true;
}
