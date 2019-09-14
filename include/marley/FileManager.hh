#pragma once
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

      /// @brief
      /// @param Z Atomic number
      /// @param A Mass number
      /// @return Mass excess (MeV)
      double liquid_drop_model_mass_excess(int Z, int A) const;

    protected:

      /// @brief Create the singleton FileManager object
      FileManager();
  };

}
