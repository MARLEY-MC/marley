// Minimal external consumer of an installed MARLEY. See CMakeLists.txt in
// this directory for why it exists.
//
// EventFileReader.hh pulls in HepMC3 types and FileManager.hh does not.
// Including both distinguishes a missing include path from a broken
// imported target.

#include <iostream>

#include "marley/EventFileReader.hh"
#include "marley/FileManager.hh"

int main() {
  const auto& fm = marley::FileManager::Instance();

  // Self-location, from an external program rather than MARLEY's own
  // executable.
  std::cout << "library:  " << fm.library_path().string() << '\n';
  std::cout << "roots:    " << fm.roots().size() << '\n';
  for ( const auto& r : fm.roots() ) {
    std::cout << "  " << r.path.string()
      << ( r.from_search_path ? "  (MARLEY_SEARCH_PATH)" : "  (self-located)" )
      << '\n';
  }

  if ( fm.roots().empty() ) {
    std::cerr << "ERROR: an installed MARLEY located no data root\n";
    return 1;
  }

  // Linking is not enough: check that a shipped data file resolves too.
  auto res = marley::FileManager::resolve( "mass_table.js", "", {},
    fm.roots(), fm.working_directory() );
  if ( !res ) {
    std::cerr << "ERROR: could not resolve mass_table.js\n";
    return 1;
  }
  std::cout << "resolved: " << res.path.string() << '\n';

  return 0;
}
