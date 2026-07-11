/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

// Standard library includes
#include <iomanip>
#include <iostream>

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/marley_utils.hh"

namespace {

  bool print_version( const std::deque< std::string >& /*dummy*/ ) {
    std::cout << "MARLEY (Model of Argon Reaction Low Energy Yields) "
      << MARLEY_VERSION << '\n';
    std::cout << "Copyright (C) 2016-2024 Steven Gardiner\n";
    std::cout << "License: GNU GPL version 3 "
      << "<http://opensource.org/licenses/GPL-3.0>\n";
    std::cout << "This is free software: you are free to change and"
      << " redistribute it.\n";
    return true;
  }

  void print_generate_help() {
    std::cout << "Usage: marley generate CONFIG_FILE\n\n"
      << "Generate Monte Carlo neutrino interaction events according to the\n"
      << "settings in CONFIG_FILE (a MARLEY JSON-based job configuration"
      << " file).\n\n"
      << "This is the default command: 'marley CONFIG_FILE' is equivalent to\n"
      << "'marley generate CONFIG_FILE'.\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n";
  }

  void print_print_help() {
    std::cout << "Usage: marley print [FORMAT] INPUT_FILES...\n\n"
      << "Print one or more MARLEY event files in a human-readable format.\n\n"
      << "  FORMAT         Optional output format selector:\n"
      << "                   hepmc3 (default): stream HepMC3 text to stdout\n"
      << "                   legacy: MARLEY legacy human-readable event"
      << " summary\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n";
  }

  void print_reweight_help() {
    std::cout << "Usage: marley reweight CONFIG_FILE INPUT_FILE\n\n"
      << "Reweight previously generated MARLEY events in INPUT_FILE using the\n"
      << "weight-calculation settings in CONFIG_FILE.\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n";
  }

  void print_dumpxs_help() {
    std::cout << "Usage: marley dump-xs OUTPUT_FILE CONFIG_FILE\n\n"
      << "Tabulate the total cross section versus projectile energy using the\n"
      << "settings in CONFIG_FILE and write the results to OUTPUT_FILE.\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n";
  }

  void print_decay_help() {
    std::cout << "Usage: marley decay CONFIG_FILE\n\n"
      << "Simulate stand-alone nuclear de-excitations according to the decays\n"
      << "settings in CONFIG_FILE.\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n";
  }

  void print_summarize_help() {
    std::cout << "Usage: marley summarize OUTPUT_FILE INPUT_FILE...\n\n"
      << "Convert one or more MARLEY event files into a ROOT TTree"
      << " summary file\n"
      << "suitable for analysis with ROOT C++ macros or Python/PyROOT.\n\n"
      << "  OUTPUT_FILE    Output ROOT file (will be created or overwritten)\n"
      << "  INPUT_FILE...  One or more input MARLEY event files\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n\n"
      << "Note: this command requires a ROOT-enabled build of MARLEY.\n";
  }

  void print_help_command_help() {
    std::cout << "Usage: marley help [COMMAND]\n\n"
      << "Show the top-level help message or detailed help for COMMAND.\n\n"
      << "Options:\n"
      << "  -h, --help     Print this help message\n";
  }

  void print_version_command_help() {
    std::cout << "Usage: marley version\n\n"
      << "Print version information.\n";
  }

  void print_convert_help() {
    std::cout << "Usage: marley convert [--output-format FORMAT]"
      << " -o OUTPUT_FILE INPUT_FILES...\n\n"
      << "Convert MARLEY event files between supported formats.\n\n"
      << "  -o OUTPUT_FILE        Required: path to the output file\n"
      << "  --output-format FORMAT Output format. Supported values:\n"
      << "                          \"ascii\" (default), \"root\",\n"
      << "                          \"legacy\", \"hepevt\"\n"
      << "  -f, --force           Overwrite the output file without"
      << " prompting\n\n"
      << "Options:\n"
      << "  -h, --help            Print this help message\n\n"
      << "Note: ROOT format conversions require a ROOT-enabled build"
      << " of MARLEY.\n"
      << "The \"legacy\" and \"hepevt\" formats are one-way conversions"
      << " from the\n"
      << "current HepMC3-based format to the deprecated output formats"
      << " used in\n"
      << "MARLEY v1.2.1 and earlier.\n";
  }

}

const marley::CommandHandler::CommandMap
  marley::CommandHandler::command_map_ =
{
  { "generate", { "Generate Monte Carlo events",
    print_generate_help, cmd_generate } },

  { "print", { "Print existing events in a human-readable format",
    print_print_help, cmd_print } },

  { "reweight", { "Reweight previously generated events",
    print_reweight_help, cmd_reweight } },

  { "dump-xs", { "Tabulate total cross section vs. projectile kinetic energy",
    print_dumpxs_help, cmd_dumpxs } },

  { "decay",  { "Simulate nuclear de-excitations",
    print_decay_help, cmd_decay } },

  { "summarize", { "Create a ROOT TTree summary of event files",
    print_summarize_help, cmd_summarize, true } },

  // ASCII-to-ASCII conversion does not need ROOT; ROOT I/O is
  // guarded by runtime checks in the implementation.
  { "convert", { "Convert event files between supported formats",
    print_convert_help, cmd_convert, false } },

  { "help", { "Show this help message or help for a specific command",
    print_help_command_help, cmd_help } },

  { "version", { "Print version information",
    print_version_command_help, print_version } },
};

bool marley::CommandHandler::print_top_level_help() {
  std::cout << "Usage: marley <command> [options]\n\n"
    << "Commands:\n";

  for ( const auto& [name, info] : command_map_ ) {
    std::cout << "  " << std::left << std::setw(12) << name
      << info.summary_;
    if ( info.requires_root_ ) std::cout << "  [requires ROOT]";
    std::cout << '\n';
  }

  std::cout << "\nOptions:\n"
    << "  -h, --help     Show top-level help\n"
    << "  -v, --version  Print version information\n";
  std::cout << "\nRun 'marley help <command>' or 'marley <command> --help' for"
    << " details.\n"
    << "MARLEY home page: <https://www.marleygen.org>\n";
  return true;
}

bool marley::CommandHandler::execute() {

  // If we don't have any command-line arguments, then just print the
  // standard top-level help message
  if ( cmds_.empty() ) {
    return this->print_top_level_help();
  }

  // Otherwise, get the first argument which should correspond to a command
  // or one of the top-level options
  std::string subcmd = cmds_.front();

  // Reassign some synonyms to their corresponding commands that appear in the
  // map of accepted values
  if ( subcmd == "--version" || subcmd == "-v" ) {
    subcmd = "version";
  }
  else if ( subcmd == "--help" || subcmd == "-h" ) {
    subcmd = "help";
  }

  auto cmd_iter = command_map_.find( subcmd );
  if ( cmd_iter != command_map_.end() ) {
    // The user provided an explicit command, so drop it from the
    // deque before delegating to the appropriate function. We no longer need to
    // resolve the command name.
    cmds_.pop_front();
    // Call the function corresponding to the chosen command
    return cmd_iter->second.cmd_( cmds_ );
  }

  // The 'marley' command is an easter egg that doesn't appear in the official
  // list
  if ( subcmd == "marley" ) {
    std::cout << marley_utils::marley_pic;
    return true;
  }

  // The default command is 'generate', so if the user didn't explicitly
  // provide it and didn't use an option prefix, then assume that the user
  // intended 'generate'
  if ( !subcmd.empty() && subcmd[0] != '-' ) {
    return this->cmd_generate( cmds_ );
  }

  std::cerr << "marley: unknown option '" << subcmd << "'\n";
  std::cerr << "Run 'marley help' for a list of available commands.\n";
  return false;
}

marley::CommandHandler::CommandHandler( int argc, char* argv[] )
  // Build a deque of the command-line arguments for passing to commands
  : cmds_( argv, argv + argc )
{
  // Strip the executable name off of the front of the deque
  cmds_.pop_front();
}
