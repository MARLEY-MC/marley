// Minimal external consumer of an installed MARLEY. See CMakeLists.txt in
// this directory for why it exists.
//
// It deliberately includes a header that pulls in HepMC3 types
// (EventFileReader.hh) as well as one that does not, because the two fail
// in different ways: a missing include path breaks the first, while a
// broken imported target breaks both.

#include <iostream>

#include "marley/EventFileReader.hh"
#include "marley/FileManager.hh"

int main() {
  const auto& fm = marley::FileManager::Instance();

  // Self-location, exercised from an external program rather than from
  // MARLEY's own executable. This is the property the whole data-discovery
  // design exists to provide, and nothing else in CI tests it from outside.
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

  // Resolve a shipped file the way MARLEY's own code does. A consumer that
  // links successfully but cannot find data is still broken.
  auto res = marley::FileManager::resolve( "mass_table.js", "", {},
    fm.roots(), fm.working_directory() );
  if ( !res ) {
    std::cerr << "ERROR: could not resolve mass_table.js\n";
    return 1;
  }
  std::cout << "resolved: " << res.path.string() << '\n';

  return 0;
}
