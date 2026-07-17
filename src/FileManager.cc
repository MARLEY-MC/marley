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

#include <cstdlib>
#include <filesystem>

#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/Logger.hh"
#include "marley/marley_utils.hh"

namespace {
  // Delimiter that separates directory names in search paths used by
  // marley::FileManager
  constexpr char SEARCH_PATH_DELIMITER = ':';

  bool file_exists(const std::string& file_name) {
    return std::filesystem::is_regular_file(file_name);
  }

  //bool file_is_readable(const std::string& file_name) {
  //  std::ifstream test_stream( file_name );
  //  return test_stream.good();
  //}
}

// Define the static default search path member variable. Its value will be set
// when the singleton instance of marley::FileManager is constructed.
std::string marley::FileManager::default_search_path_;

marley::FileManager::FileManager() {

  // Check whether the MARLEY environment variable is set
  char* mar = std::getenv("MARLEY");
  if ( !mar ) throw marley::Error("The MARLEY environment variable is not set."
    " Please set it (e.g., by sourcing the setup_marley.sh script) and"
    " try again.");

  marley_dir_ = std::string( mar );

  // If the MARLEY_SEARCH_PATH environment variable is set, use that
  // instead of the default search path
  char* msp = std::getenv("MARLEY_SEARCH_PATH");
  if ( msp ) default_search_path_ = std::string( msp );
  else {
    default_search_path_ = marley_dir_ + "/data";
    default_search_path_ += ':' + marley_dir_ + "/data/react";
    default_search_path_ += ':' + marley_dir_ + "/data/optical_model";
    default_search_path_ += ':' + marley_dir_ + "/data/structure";
  }

  MARLEY_LOG( DEBUG, "io" ) << "MARLEY search path set to \""
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

  // First check if we can find the file just using base_name. We can
  // avoid the search if it's already accessible (e.g., because the
  // full path was supplied instead of the base name)
  if ( file_exists(base_name) ) return base_name;

  // If the file name contains directory separators, try resolving
  // the relative path against each search directory
  if ( base_name.find('/') != std::string::npos ) {
    for (const auto& dir : search_dirs) {
      std::string candidate = dir + '/' + base_name;
      if ( file_exists(candidate) ) return candidate;
    }
  }

  // Search the directories in the order listed. In each directory,
  // compare all regular file names to the given base name. If a match
  // is found, stop the search and store the result in full_path_to_file.
  // If no match is found, full_path_to_file will remain empty.
  for (const auto& dir : search_dirs) {
    bool file_found = dir_iterate(dir,
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

    if ( file_found ) break;
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
bool marley::FileManager::dir_iterate(const std::string& dir_name,
  const std::function<bool(const std::string&, const std::string&)>& func)
{
  std::error_code ec;
  std::filesystem::directory_iterator it(dir_name,
    std::filesystem::directory_options::skip_permission_denied, ec);

  if (ec) {
    MARLEY_LOG(WARN, "io") << "Could not read from the directory \""
      << dir_name << '\"';
    return false;
  }

  std::filesystem::directory_iterator end;
  bool stop_iterations = false;
  while (!stop_iterations && it != end) {
    const auto& entry = *it;
    std::string full_file_name = entry.path().string();
    std::string base_name = entry.path().filename().string();

    MARLEY_LOG(TRACE, "io") << "marley::FileManager found file \""
      << full_file_name << '\"';

    if (entry.is_regular_file()) {
      stop_iterations = func(full_file_name, base_name);
    }

    it.increment(ec);
    if (ec) {
      MARLEY_LOG(DEBUG, "io") << "Couldn't access a file in directory \""
        << dir_name << '\"';
      ec.clear();
    }
  }

  return stop_iterations;
}
