/// @todo Replace with more portable commands from the C++17 filesystem
/// library (avoided for now to maintain support for old compilers)
#include "dirent.h"

#include "marley/FileManager.hh"

const marley::FileManager& marley::FileManager::Instance() {

  // Create the mass table using a static variable. This ensures
  // that the singleton instance is only created once.
  static std::unique_ptr<marley::FileManager>
    the_instance(new marley::FileManager());

  // Return a reference to the singleton instance
  return *the_instance;
}
