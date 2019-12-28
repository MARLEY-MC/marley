#pragma once
#include <functional>
#include <vector>

namespace marley {

  /// @brief Singleton class that handles file searches
  class FileManager {

    public:

      /// @brief Deleted copy constructor
      FileManager(const FileManager&) = delete;

      /// @brief Deleted move constructor
      FileManager(FileManager&&) = delete;

      /// @brief Deleted copy assignment operator
      FileManager& operator=(const FileManager&) = delete;

      /// @brief Deleted move assignment operator
      FileManager& operator=(FileManager&&) = delete;

      /// @brief Get a const reference to the singleton instance of the
      /// FileManager
      static const FileManager& Instance();

      /// @brief Converts a ':'-delimited search path (given as a single
      /// string) into a vector of directory names
      /// @param search_path The search path to convert
      /// that should be searched for the file
      /// @return A vector of directory names corresponding to
      /// those given in the search path
      static std::vector<std::string> search_path_to_dir_vector(
        const std::string& search_path = FileManager::default_search_path_);

      /// @brief Searches for a file in the given directories
      /// @details The search is non-recursive, i.e., no subdirectories
      /// are included in the search
      /// @param base_name The base name (file name without any path)
      /// of the desired file
      /// @param search_dirs A vector of strings specifying the directories
      /// that should be searched for the file
      /// @return The full path to the file, or an empty string if
      /// no match could be found
      std::string find_file(const std::string& base_name,
        const std::vector<std::string>& search_dirs) const;

      /// @brief Searches for a file in the given search path
      /// @details The search is non-recursive, i.e., no subdirectories
      /// are included in the search
      /// @param base_name The base name (file name without any path)
      /// of the desired file
      /// @param search_path A string containing a ':'-delimited list
      /// of directories that should be searched for the file
      /// @return The full path to the file, or an empty string if
      /// no match could be found
      std::string find_file(const std::string& base_name,
        const std::string& search_path
        = marley::FileManager::default_search_path_) const;

      /// @brief Get a vector of full paths for all regular files in
      /// the requested directories
      /// @details Files in subdirectories are not included
      /// (the search is non-recursive)
      /// @param search_dirs A vector of strings specifying the directories
      /// that should be scanned
      /// @return A vector of strings giving the full path to each regular file
      /// that was found
      std::vector<std::string> list_all_files(
        const std::vector<std::string>& search_dirs) const;

      /// @brief Get a vector of full paths for all regular files in
      /// the requested search path
      /// @details Files in subdirectories are not included
      /// (the search is non-recursive)
      /// @param search_path A string containing a ':'-delimited list
      /// of directories that should be scanned for files
      /// @return A vector of strings giving the full path to each regular file
      /// that was found
      std::vector<std::string> list_all_files(
        const std::string& search_path
        = FileManager::default_search_path_) const;

      /// @brief Returns the path to the root MARLEY folder
      inline std::string marley_dir() const { return marley_dir_; }

    protected:

      /// @brief Create the singleton FileManager object
      FileManager();

      /// @brief Helper function that executes a std::function on
      /// the name of each regular file found while scanning through
      /// a single directory (non-recursively)
      /// @param dir_name The name of the directory to be scanned
      /// @param func The function to execute. The first argument
      /// is the full file name (including the path). The second is
      /// the base name of the file (no path). The return value should
      /// be true if no further iterations are needed (e.g., because
      /// a file matching a desired name was found) and false otherwise.
      /// @return True if the iterations were stopped early (due to function
      /// returning true) or false if they were completed
      static bool dir_iterate(const std::string& dir_name,
        const std::function<bool(const std::string&, const std::string&)>& func);

      /// @brief By default, use this search path when looking for files
      static std::string default_search_path_;

      /// @brief Stores the value of the MARLEY environment variable
      /// (which points to the root folder of the source code distribution)
      std::string marley_dir_;
  };

}
