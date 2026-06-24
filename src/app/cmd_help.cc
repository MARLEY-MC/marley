// Standard library includes
#include <iostream>
#include <string>

// MARLEY includes
#include "marley/CommandHandler.hh"

bool marley::CommandHandler::cmd_help( std::deque< std::string >& args ) {
  if ( args.empty() ) {
    marley::CommandHandler::print_top_level_help();
    return true;
  }

  if ( args.size() > 1u ) {
    std::cerr << "marley help: too many arguments\n";
    return false;
  }

  // If we've made it here, then the help command has exactly one argument
  // which corresponds to a command name or the -h/--help options
  std::string command_name = args.front();

  // These options are considered synonyms of the help command itself
  if ( command_name == "-h" || command_name == "--help" ) {
    command_name = "help";
  }

  auto cmd_iter = command_map_.find( command_name );
  if ( cmd_iter != command_map_.end() ) {
    cmd_iter->second.print_help_();
    return true;
  }

  std::cerr << "marley: unknown command '" << command_name << "'\n";
  std::cerr << "Run 'marley help' for a list of available commands.\n";
  return false;
}
