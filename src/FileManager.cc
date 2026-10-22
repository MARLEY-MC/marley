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

#include <dlfcn.h>

#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>

#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/Logger.hh"
#include "marley/marley_utils.hh"

namespace fs = std::filesystem;

namespace {

  // Delimiter that separates directory names in MARLEY_SEARCH_PATH
  constexpr char SEARCH_PATH_DELIMITER = ':';

  // Probed inside a candidate root to decide whether it really is one
  constexpr char SENTINEL_FILE[] = "mass_table.js";

  // Hop taken from the directory containing libMARLEY to reach the data root
  constexpr char DATA_SUBPATH[] = "share/marley/data";

  /// @brief Anchor for self-location.
  /// @details This must be a function with internal linkage *defined in this
  /// translation unit*. dladdr() reports the object containing whatever
  /// address it is given, so taking the address of anything inline or
  /// header-defined risks naming the consumer's object instead of libMARLEY
  /// — silently, and with a plausible-looking answer.
  void self_location_anchor() {}

  /// @brief Ask the dynamic loader where libMARLEY is
  /// @return The canonical path to the library, or empty on failure
  fs::path locate_self() {

    Dl_info info;
    void (*anchor)() = &self_location_anchor;

    if ( !dladdr(reinterpret_cast<void*>(anchor), &info) || !info.dli_fname ) {
      return {};
    }

    // dli_fname is loader-dependent: glibc reports the path as it was opened
    // (possibly relative, possibly through a symlink), while macOS reports a
    // fully resolved one. Canonicalising makes the two agree, and absolutises
    // a relative result so that a later chdir cannot move the data root.
    std::error_code ec;
    fs::path p = fs::weakly_canonical( fs::path(info.dli_fname), ec );
    if ( ec ) return {};

    return p;
  }

  /// @brief Test whether a directory looks like a MARLEY data root
  bool has_sentinel( const fs::path& dir ) {
    std::error_code ec;
    return fs::is_regular_file( dir / SENTINEL_FILE, ec );
  }

  /// @brief Absolutise and canonicalise a user-supplied directory
  fs::path normalise( const fs::path& p ) {
    std::error_code ec;
    fs::path out = fs::weakly_canonical( p, ec );
    if ( ec ) return p;
    return out;
  }
}

marley::FileManager::FileManager() {

  // Snapshot the working directory. Lookups happen lazily -- structure files
  // are resolved during the de-excitation cascade -- so reading it per call
  // would make resolution depend on when a given nuclide was first needed.
  std::error_code ec;
  cwd_ = fs::current_path( ec );
  if ( ec ) cwd_.clear();

  // Roots supplied by the user come first. They are trusted as given: no
  // sentinel check, because a plain directory legitimately has none.
  const char* msp = std::getenv( "MARLEY_SEARCH_PATH" );
  if ( msp ) {
    for ( const auto& dir : search_path_to_dir_vector( msp ) ) {
      if ( dir.empty() ) continue;
      fs::path p = normalise( dir );
      if ( !fs::is_directory(p, ec) ) {
        MARLEY_LOG( WARN, "io.resolve" ) << "MARLEY_SEARCH_PATH entry \""
          << dir << "\" is not an existing directory. It will be kept, but"
          " nothing will be found there.";
      }
      roots_.push_back( DataRoot{ p, true } );
    }
  }

  // Then the root found by asking the loader where this library is.
  lib_path_ = locate_self();
  if ( lib_path_.empty() ) {
    MARLEY_LOG( WARN, "io.resolve" ) << "Could not determine the location of"
      " the MARLEY shared library. Data files will be found only via"
      " MARLEY_SEARCH_PATH or the working directory.";
  }
  else {
    fs::path root = ( lib_path_.parent_path() / ".." / DATA_SUBPATH )
      .lexically_normal();

    if ( has_sentinel(root) ) {
      roots_.push_back( DataRoot{ root, false } );
    }
    else {
      // Dropped rather than fatal: construction never throws, so the eventual
      // lookup can report the whole trace instead of this one symptom.
      MARLEY_LOG( WARN, "io.resolve" ) << "The expected MARLEY data directory"
        " \"" << root.string() << "\" does not contain " << SENTINEL_FILE
        << ", so it will not be used.";
    }
  }

  if ( roots_.empty() ) {
    MARLEY_LOG( WARN, "io.resolve" ) << "No valid MARLEY data root was"
      " located. Set MARLEY_SEARCH_PATH to the directory containing "
      << SENTINEL_FILE << " if MARLEY's data files are installed elsewhere.";
  }
}

const marley::FileManager& marley::FileManager::Instance() {

  // Function-local static: the singleton is created on first use and,
  // unlike the old implementation, its constructor cannot throw.
  static std::unique_ptr< marley::FileManager >
    the_instance( new marley::FileManager() );

  return *the_instance;
}

std::vector< std::string > marley::FileManager::search_path_to_dir_vector(
  const std::string& search_path )
{
  return marley_utils::split_string( search_path, SEARCH_PATH_DELIMITER );
}

marley::FileManager::Resolution marley::FileManager::resolve(
  const std::string& name, const std::string& subdir,
  const fs::path& base_dir, const std::vector< DataRoot >& roots,
  const fs::path& cwd )
{
  Resolution res;

  auto try_path = [&res]( const fs::path& p ) -> bool {
    res.probed.push_back( p );
    std::error_code ec;
    if ( fs::is_regular_file(p, ec) ) {
      res.path = p;
      return true;
    }
    return false;
  };

  const fs::path name_path( name );

  // Leg 1: a path the caller supplied outright is used as given, and no
  // other candidate is tried. Searching elsewhere for an explicitly named
  // file would be wrong, not helpful.
  if ( name_path.is_absolute() || name.rfind("./", 0) == 0
    || name.rfind("../", 0) == 0 )
  {
    try_path( name_path );
    return res;
  }

  // Leg 2: resolve against a directory the caller supplied. Tested for
  // non-emptiness, not for "was passed": base_dir defaults to an empty path
  // and path("") / "name" == "name", so branching on the latter would probe
  // the working directory ahead of the roots.
  if ( !base_dir.empty() ) {
    if ( try_path(base_dir / name_path) ) return res;
  }

  // Leg 3: the data roots, in order.
  for ( const auto& r : roots ) {
    if ( subdir.empty() ) {
      if ( try_path(r.path / name_path) ) return res;
    }
    else {
      if ( try_path(r.path / subdir / name_path) ) return res;
      // A search-path entry may be a plain directory standing in for a
      // category directory. The self-located root is skipped here: its
      // layout is known, so a hit could only be a file added to the tree.
      if ( r.from_search_path && try_path(r.path / name_path) ) return res;
    }
  }

  // Leg 4: the working directory, last, and without a subdir. Last is what
  // makes it safe -- it can supply a name nothing else provides, but can
  // never replace shipped data.
  if ( !cwd.empty() ) {
    if ( try_path(cwd / name_path) ) {
      res.from_working_directory = true;
      return res;
    }
  }

  return res;
}

fs::path marley::FileManager::find_file( const std::string& name,
  const std::string& subdir, const fs::path& base_dir ) const
{
  Resolution res = resolve( name, subdir, base_dir, roots_, cwd_ );

  if ( !res ) throw marley::Error( failure_message(name, subdir, res) );

  report( name, res );
  return res.path;
}

void marley::FileManager::report( const std::string& name,
  const Resolution& res ) const
{
  // Log each distinct file once. Logging per call would be unusable: a
  // de-excitation cascade resolves a decay-scheme file per nuclide.
  if ( !reported_.insert(res.path).second ) return;

  // A hit on the working-directory leg is worth seeing by default; an
  // ordinary resolution is not. data/config/logger.js ships
  // io.resolve: "notice", which draws exactly that line.
  //
  // The flag comes from resolve(), which knows which leg fired. It used to
  // be inferred here by comparing the result against cwd/name, which is not
  // the same question: an "#include:" of a file beside its configuration,
  // run from that directory, is supplied by the base_dir leg at a path the
  // working-directory leg would also have produced -- and was reported as a
  // working-directory hit.
  if ( res.from_working_directory ) {
    MARLEY_LOG( NOTICE, "io.resolve" ) << '"' << name << "\" resolved from the"
      " working directory: " << res.path.string();
  }
  else {
    MARLEY_LOG( INFO, "io.resolve" ) << '"' << name << "\" -> "
      << res.path.string();
  }
}

std::string marley::FileManager::failure_message( const std::string& name,
  const std::string& subdir, const Resolution& res ) const
{
  std::ostringstream msg;

  msg << "Could not find the MARLEY data file \"" << name << '\"';
  if ( !subdir.empty() ) msg << " (expected under a \"" << subdir
    << "\" directory)";
  msg << ".\n";

  if ( roots_.empty() ) {
    msg << "  No valid data root was located.";
  }
  else {
    msg << "  Data roots, in search order:\n";
    for ( const auto& r : roots_ ) {
      msg << "    " << r.path.string()
        << ( r.from_search_path ? "   (MARLEY_SEARCH_PATH)"
                                : "   (located from libMARLEY)" ) << '\n';
    }
  }

  msg << "\n  Paths tried:\n";
  for ( const auto& p : res.probed ) msg << "    " << p.string() << '\n';

  fs::path suggestion = suggest_for( name, roots_ );
  if ( !suggestion.empty() ) {
    msg << "\n  Did you mean \"" << suggestion.generic_string() << "\"?";
  }

  return msg.str();
}

std::filesystem::path marley::FileManager::suggest_for(
  const std::string& name, const std::vector< DataRoot >& roots )
{
  // Walk the roots looking for a file with this base name. Only reached
  // when a lookup has already failed, so the cost is paid once, in a run
  // that is about to stop -- or in "marley info --resolve", whose entire
  // job is to answer this.
  const std::string base = fs::path( name ).filename().string();

  for ( const auto& r : roots ) {
    std::error_code ec;
    fs::recursive_directory_iterator it( r.path,
      fs::directory_options::skip_permission_denied, ec );
    if ( ec ) continue;

    fs::recursive_directory_iterator end;
    for ( ; it != end; it.increment(ec) ) {
      if ( ec ) { ec.clear(); continue; }
      if ( !it->is_regular_file(ec) ) continue;
      if ( it->path().filename().string() != base ) continue;

      fs::path rel = fs::relative( it->path(), r.path, ec );
      if ( ec ) continue;
      return rel;
    }
  }

  return {};
}
