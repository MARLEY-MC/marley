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

#pragma once
#include <filesystem>
#include <set>
#include <string>
#include <vector>

namespace marley {

  /// @brief Singleton class that locates MARLEY's data files
  /// @details The data root is found by asking the dynamic loader where
  /// libMARLEY itself is and taking a single fixed hop from there. No
  /// environment variable is required. See the design specification, §2.
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

      /// @brief A directory that may contain MARLEY data files
      struct DataRoot {
        /// @brief Absolute path to the root
        std::filesystem::path path;

        /// @brief True if this root came from MARLEY_SEARCH_PATH
        /// @details Such roots are also probed for <root>/<name> in addition
        /// to <root>/<subdir>/<name>, so that a plain directory may stand in
        /// for a category directory. The self-located root is not probed that
        /// way: its layout is known, so a bare hit there could only match a
        /// file added to the installed tree.
        bool from_search_path = false;
      };

      /// @brief The outcome of a name resolution, including everything that
      /// was tried
      struct Resolution {
        /// @brief The resolved path, or empty if the name was not found
        std::filesystem::path path;

        /// @brief Every candidate that was tested, in the order tested
        std::vector< std::filesystem::path > probed;

        bool found() const { return !path.empty(); }
        explicit operator bool() const { return found(); }
      };

      /// @brief Get a const reference to the singleton instance
      static const FileManager& Instance();

      /// @brief Resolve a name against a set of roots
      /// @details Pure and non-throwing: it touches the filesystem only to
      /// test candidates, and takes everything else as a parameter so that it
      /// can be unit-tested without a singleton. The legs, in order, are
      ///   1. @p name is absolute or begins "./" or "../"  -> used as given
      ///   2. @p base_dir is non-empty                     -> base_dir/name
      ///   3. each root in order                           -> root/subdir/name
      ///                                                   -> root/name
      ///                                        (second probe: search-path
      ///                                         roots only, and only when
      ///                                         subdir is non-empty)
      ///   4. the working directory                        -> cwd/name
      /// @param name The name to resolve, as written by the caller or the user
      /// @param subdir The category directory expected beneath a root, or ""
      /// @param base_dir A directory to resolve against first, or empty
      /// @param roots The data roots, in search order
      /// @param cwd The working directory to use for leg 4
      /// @return A Resolution; check found() before using path()
      static Resolution resolve( const std::string& name,
        const std::string& subdir,
        const std::filesystem::path& base_dir,
        const std::vector< DataRoot >& roots,
        const std::filesystem::path& cwd );

      /// @brief Resolve a name, or throw
      /// @details Calls resolve() against this FileManager's roots and
      /// working directory. On success the resolved path is reported once at
      /// INFO under the "io.resolve" logging category (or NOTICE if it came
      /// from the working directory). On failure a marley::Error is thrown
      /// carrying the whole trace.
      /// @param name The name to resolve
      /// @param subdir The category directory expected beneath a root, or ""
      /// @param base_dir A directory to resolve against first, or empty
      /// @return The resolved path
      std::filesystem::path find_file( const std::string& name,
        const std::string& subdir = "",
        const std::filesystem::path& base_dir = {} ) const;

      /// @brief Get the data roots, in search order
      /// @details May be empty if no valid root could be located. Callers
      /// that report roots to the user must handle that case.
      const std::vector< DataRoot >& roots() const { return roots_; }

      /// @brief Get the working directory used for the final resolution leg
      /// @details Snapshotted when the singleton is constructed, so that a
      /// later chdir cannot change where a lookup looks.
      const std::filesystem::path& working_directory() const { return cwd_; }

      /// @brief Get the path to libMARLEY, as located at construction
      /// @details Empty if self-location failed. Reported by "marley info".
      const std::filesystem::path& library_path() const { return lib_path_; }

    private:

      /// @brief Create the singleton FileManager object
      FileManager();

      /// @brief Converts a ':'-delimited search path into a vector of
      /// directory names
      static std::vector< std::string > search_path_to_dir_vector(
        const std::string& search_path );

      /// @brief Build the message for a failed lookup, including a
      /// "did you mean" suggestion when one can be found
      std::string failure_message( const std::string& name,
        const std::string& subdir, const Resolution& res ) const;

      /// @brief Report a successful resolution once, at INFO (or NOTICE if
      /// the working-directory leg supplied it)
      void report( const std::string& name, const Resolution& res ) const;

      /// @brief The data roots, in search order
      std::vector< DataRoot > roots_;

      /// @brief The working directory, snapshotted at construction
      std::filesystem::path cwd_;

      /// @brief Path to libMARLEY, or empty if self-location failed
      std::filesystem::path lib_path_;

      /// @brief Paths already reported under "io.resolve", so that each
      /// distinct file is logged once rather than once per lookup
      mutable std::set< std::filesystem::path > reported_;
  };

}
