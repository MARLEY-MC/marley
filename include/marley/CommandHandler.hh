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

#pragma once

// Standard library includes
#include <deque>
#include <functional>
#include <map>
#include <string>

namespace marley {

  /// @brief  Wrapper for handling commands in the main marley executable
  class CommandHandler {

  public:

    CommandHandler( int argc, char* argv[] );

    struct CommandInfo {
      /// Simple description of the command
      std::string summary_;

      /// Functiont that prints a detailed help message to stdout
      std::function< void() > print_help_;

      /// Function to call to execute the command
      std::function< bool( std::deque< std::string >& ) > cmd_;

      /// Indicates whether MARLEY needs to be linked against ROOT to execute
      /// this command
      bool requires_root_ = false;
    };

    using CommandMap = std::map< std::string, CommandInfo >;

    /// Print the top-level help message by aggregating the command
    /// summaries from the map
    static bool print_top_level_help();

    /// Execute the requested command based on the command-line arguments
    /// @return true if the command was successful, false otherwise
    bool execute();

  protected:

    /// Map containing information about recognized commands handled
    /// by the marley executable
    static const CommandMap command_map_;

    /// Stores the command-line arguments with the executable name removed
    /// from the front
    std::deque< std::string > cmds_;

    // -------------------------------------
    // Implementation of commands
    // -------------------------------------

    /// Generate Monte Carlo events
    static bool cmd_generate( std::deque< std::string >& args );

    /// Simulate nuclear de-excitations
    static bool cmd_decay( std::deque< std::string >& args );

    /// Tabulate energy-dependent total cross section values
    static bool cmd_xsec( std::deque< std::string >& args );

    /// Display top-level or command-specific help messages
    static bool cmd_help( std::deque< std::string >& args );

    /// Print existing MARLEY events in a human-readable format
    static bool cmd_print( std::deque< std::string >& args );

    /// Reweight existing MARLEY events
    static bool cmd_reweight( std::deque< std::string >& args );

    /// Summarize an existing sample of MARLEY events as a ROOT TTree
    static bool cmd_summarize( std::deque< std::string >& args );

    /// Convert MARLEY event files between supported formats
    static bool cmd_convert( std::deque< std::string >& args );

  };

}
