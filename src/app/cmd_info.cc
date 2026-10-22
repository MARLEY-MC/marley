/// @file
/// @copyright Copyright (C) 2016-2026 Steven Gardiner
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
#include <filesystem>
#include <iostream>
#include <string>

// MARLEY includes
#include "marley/CommandHandler.hh"
#include "marley/FileManager.hh"

// "marley info" answers "what will happen", as opposed to marley-config's
// "what was built" (design spec 6.1). It exists because the only diagnostic
// previously available for data discovery was at DEBUG level, filtered out
// by the shipped logging configuration and compiled out of release builds
// entirely -- so a production binary could not be asked where it was
// looking.
//
// Its output is for humans and is NOT contractual. There is deliberately no
// machine-readable mode: nothing in MARLEY's own packaging parses it, and a
// recipe that parsed it would reintroduce exactly the build-time coupling
// that removing marley-config --datadir was meant to prevent. A stable
// key-value mode is a reasonable future request; until it exists, nothing
// should depend on this format.
//
// It answers for the CLI binary, not for an embedder. "marley info" reports
// the roots of the libMARLEY that this executable loaded. An art module
// linked against a different libMARLEY has different roots -- which is
// precisely the situation in which someone reaches for this tool. The
// honest reading is "what would the command-line MARLEY do here"; the
// answer for an embedded one is in that job's log, under io.resolve.

namespace {

  namespace fs = std::filesystem;

  void print_info_help() {
    std::cout << "Usage: marley info [--resolve [--subdir DIR] NAME]\n\n"
      << "Report where this MARLEY executable looks for its data files.\n\n"
      << "With no arguments, print the located library, the data roots in\n"
      << "search order, and the working directory.\n\n"
      << "  --resolve NAME   Show every path MARLEY would try for NAME, and\n"
      << "                   which one it would use\n"
      << "  --subdir DIR     The category directory to resolve NAME under,\n"
      << "                   as the calling code would pass it (e.g. react,\n"
      << "                   structure, config, optical_model). Defaults to\n"
      << "                   none, which is rarely what a real caller does\n\n"
      << "Options:\n"
      << "  -h, --help       Print this message\n\n"
      << "This output is intended for humans and its format is not stable.\n"
      << "Do not parse it from a script or a packaging recipe. The exit\n"
      << "status is stable: --resolve exits 0 if the name resolves and 1 if\n"
      << "it does not.\n";
  }

  /// @brief Print the roots, or say plainly that there are none
  void print_roots( const marley::FileManager& fm ) {
    const auto& roots = fm.roots();

    std::cout << "Library:           ";
    if ( fm.library_path().empty() ) {
      // Self-location failing is not fatal -- the roots list simply loses
      // its shipped entry -- but it is the first thing worth knowing.
      std::cout << "(could not be located)\n";
    }
    else std::cout << fm.library_path().string() << '\n';

    std::cout << "Working directory: " << fm.working_directory().string()
      << "\n\n";

    if ( roots.empty() ) {
      std::cout << "Data roots:        no valid data root was located.\n"
        << "                   Every lookup will fall back to the working\n"
        << "                   directory, and most will fail.\n";
      return;
    }

    std::cout << "Data roots, in search order:\n";
    for ( const auto& r : roots ) {
      std::cout << "  " << r.path.string();
      std::cout << ( r.from_search_path ? "   (MARLEY_SEARCH_PATH)"
                                        : "   (located from libMARLEY)" );
      std::error_code ec;
      if ( !fs::is_directory(r.path, ec) ) std::cout << "   [DOES NOT EXIST]";
      std::cout << '\n';
    }
    std::cout << "  " << fm.working_directory().string()
      << "   (working directory, tried last)\n";
  }

  /// @brief Show the whole trace for one lookup, found or not
  bool resolve_one( const marley::FileManager& fm, const std::string& name,
    const std::string& subdir )
  {
    // The pure resolve() rather than find_file(): this needs the trace for a
    // SUCCESSFUL lookup too, and find_file() returns only the path. Using
    // the same function the real callers use is the point -- a diagnostic
    // that disagrees with the thing it diagnoses is worse than none.
    auto res = marley::FileManager::resolve( name, subdir, {}, fm.roots(),
      fm.working_directory() );

    std::cout << name;
    if ( !subdir.empty() ) std::cout << "   (subdir \"" << subdir << "\")";
    std::cout << "\n\n";

    if ( res.probed.empty() ) {
      std::cout << "  No candidate paths were tried. There are no data roots"
        " and the\n  working directory leg did not apply.\n";
    }
    else {
      std::cout << "  Paths tried, in order:\n";
      for ( const auto& p : res.probed ) {
        const bool hit = ( p == res.path );
        std::cout << "    " << ( hit ? "-> " : "   " ) << p.string()
          << ( hit ? "   USED" : "" ) << '\n';
      }
    }

    std::cout << '\n';

    if ( res.found() ) {
      std::cout << "FOUND  " << res.path.string() << '\n';
      return true;
    }

    std::cout << "NOT FOUND\n";

    // The suggestion is the whole reason this command is worth running on a
    // name that fails. --subdir defaults to none, which is rarely what a
    // real caller passes, so the common mistake is asking about "ES.react"
    // when the generator asks for it under "react" -- and the scan answers
    // that without the user having to know.
    fs::path suggestion = marley::FileManager::suggest_for( name, fm.roots() );
    if ( !suggestion.empty() ) {
      std::cout << "\nDid you mean \"" << suggestion.generic_string()
        << "\"?";
      fs::path parent = suggestion.parent_path();
      if ( !parent.empty() && subdir.empty() ) {
        std::cout << "  Or pass --subdir " << parent.generic_string();
      }
      std::cout << '\n';
    }

    return false;
  }

}

bool marley::CommandHandler::cmd_info( std::deque< std::string >& args ) {

  std::string subdir;
  std::string name;
  bool want_resolve = false;

  while ( !args.empty() ) {
    std::string arg = args.front();
    args.pop_front();

    if ( arg == "-h" || arg == "--help" ) {
      print_info_help();
      return true;
    }
    else if ( arg == "--resolve" ) {
      want_resolve = true;
    }
    else if ( arg == "--subdir" ) {
      if ( args.empty() ) {
        std::cerr << "marley info: --subdir requires a directory name\n";
        return false;
      }
      subdir = args.front();
      args.pop_front();
    }
    else if ( !arg.empty() && arg.front() == '-' ) {
      std::cerr << "marley info: unknown option \"" << arg << "\"\n";
      return false;
    }
    else if ( name.empty() ) {
      name = arg;
    }
    else {
      std::cerr << "marley info: unexpected extra argument \"" << arg
        << "\"\n";
      return false;
    }
  }

  const auto& fm = marley::FileManager::Instance();

  if ( !want_resolve ) {
    if ( !name.empty() ) {
      std::cerr << "marley info: \"" << name << "\" given without --resolve\n";
      return false;
    }
    if ( !subdir.empty() ) {
      std::cerr << "marley info: --subdir only applies to --resolve\n";
      return false;
    }
    print_roots( fm );
    return true;
  }

  if ( name.empty() ) {
    std::cerr << "marley info: --resolve requires a file name\n";
    return false;
  }

  // Always print the roots first: a resolution trace is unreadable without
  // knowing what the roots were.
  print_roots( fm );
  std::cout << '\n';

  // The exit status is the one machine-readable thing this command offers,
  // and the only one that does not involve parsing output whose format is
  // explicitly not contractual. "marley info --resolve X > /dev/null" is a
  // reasonable way to ask whether a name resolves, so a name that is not
  // found exits nonzero even though the command itself did its job.
  return resolve_one( fm, name, subdir );
}
