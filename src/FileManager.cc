#include <cstdlib>

/// @todo Replace with more portable commands from the C++17 filesystem
/// library (avoided for now to maintain support for old compilers)
#include "dirent.h"
#include "sys/stat.h"

#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/Logger.hh"
#include "marley/marley_utils.hh"

namespace {
  // Delimiter that separates directory names in search paths used by
  // marley::FileManager
  constexpr char SEARCH_PATH_DELIMITER = ':';
}

// Define the static default search path member variable. Its value will be set
// when the singleton instance of marley::FileManager is constructed.
std::string marley::FileManager::default_search_path_;

marley::FileManager::FileManager() {

  // Check whether the MARLEY environment variable is set
  char* mar = std::getenv("MARLEY");
  if ( !mar ) throw marley::Error("The MARLEY enviornment variable is not set."
    " Please set it (e.g., by sourcing the setup_marley.sh script) and"
    "try again.");

  // If the MARLEY_SEARCH_PATH enviornment variable is set, use that
  // instead of the default search path
  char* msp = std::getenv("MARLEY_SEARCH_PATH");
  if ( msp ) default_search_path_ = std::string( msp );
  else {
    std::string marley_dir( mar );
    default_search_path_ = marley_dir + "/data";
    default_search_path_ += ':' + marley_dir + "/data/react";
    default_search_path_ += ':' + marley_dir + "/data/structure";
  }

  MARLEY_LOG_DEBUG() << "MARLEY search path set to \""
    << default_search_path_ << '\"';
}

const marley::FileManager& marley::FileManager::Instance() {

  // Create the mass table using a static variable. This ensures
  // that the singleton instance is only created once.
  static std::unique_ptr<marley::FileManager>
    the_instance(new marley::FileManager());

  // Return a reference to the singleton instance
  return *the_instance;
}

std::vector<std::string> marley::FileManager::search_path_to_dir_vector(
  const std::string& search_path)
{
  std::vector<std::string> dir_names
    = marley_utils::split_string(search_path, SEARCH_PATH_DELIMITER);
  return dir_names;
}

std::string marley::FileManager::find_file(const std::string& base_name,
  const std::string& search_path) const
{
  auto search_dirs = search_path_to_dir_vector( search_path );
  return find_file(base_name, search_dirs);
}

std::string marley::FileManager::find_file(const std::string& base_name,
  const std::vector<std::string>& search_dirs) const
{
  std::string full_path_to_file;

  // Search the directories in the order listed. In each directory,
  // compare all regular file names to the given base name. If a match
  // is found, stop the search and store the result in full_path_to_file.
  // If no match is found, full_path_to_file will remain empty.
  for (const auto& dir : search_dirs) {
    dir_iterate(dir,
      [&full_path_to_file, &base_name](const std::string& other_full_path,
        const std::string& other_base_name) -> bool
      {
        if ( base_name == other_base_name ) {
          // We've found the desired file, so store its full path
          full_path_to_file = other_full_path;
          // Signal that we don't need to continue the search any longer
          return true;
        }
        return false;
      }
    );
  }

  return full_path_to_file;
}

std::vector<std::string> marley::FileManager::list_all_files(
  const std::string& search_path) const
{
  auto search_dirs = search_path_to_dir_vector( search_path );
  return list_all_files( search_dirs );
}

// Returns a vector of strings loaded with the full paths to all regular files
// in the desired search path. The directories are scanned non-recursively.
std::vector<std::string> marley::FileManager::list_all_files(
  const std::vector<std::string>& search_dirs) const
{
  std::vector<std::string> result;

  for (const auto& dir : search_dirs) {
    dir_iterate(dir,
      [&result](const std::string& full_path,
        const std::string& /*base_name*/) -> bool
      {
        // Add the current file to the vector of full paths
        result.push_back( full_path );
        // Signal that the directory scan should continue
        return false;
      }
    );
  }

  return result;
}


// Iterates non-recursively through all the regular files in a given
// directory. Calls a std::function for each one that it finds.
// If the function returns true, the iteration continues.
// If it returns false, then the iteration will stop before the
// next file is used.
void marley::FileManager::dir_iterate(const std::string& dir_name,
  const std::function<bool(const std::string&, const std::string&)>& func)
{
  // Open the directory for reading
  DIR* directory = opendir( dir_name.c_str() );

  if ( !directory ) {
    MARLEY_LOG_WARNING() << "Could not read from the directory \""
      << dir_name << '\"';
    return;
  }

  // Loop through the files in the directory one by one
  bool stop_iterations = false;
  dirent* file = nullptr;
  while ( file = readdir(directory), !stop_iterations && file ) {

    // Get information about the current file using the stat() function
    struct stat file_stat;
    std::string base_name = file->d_name;
    std::string full_file_name = dir_name + '/' + base_name;

    MARLEY_LOG_DEBUG() << "marley::FileManager found file \""
      << full_file_name << '\"';

    // If we had a problem, complain and try the next file
    if ( stat(base_name.c_str(), &file_stat) ) {
      MARLEY_LOG_DEBUG() << "Couldn't stat the file \""
        << full_file_name << '\"';
      continue;
    }

    // If the file is a regular file (as opposed to, e.g., a subdirectory),
    // then call the function with it as an argument
    if (file_stat.st_mode & S_IFREG) {
      stop_iterations = func( full_file_name, base_name );
    }
  }

  closedir( directory );
}