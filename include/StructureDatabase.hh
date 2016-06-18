#pragma once
#include <memory>
#include <string>
#include <unordered_map>

#include "BackshiftedFermiGasModel.hh"
#include "ConfigFile.hh"
#include "DecayScheme.hh"
#include "SphericalOpticalModel.hh"

namespace marley {

  /// @brief Container for nuclear structure information organized by nuclide
  /// @details Currently, the StructureDatabase object can hold nuclear
  /// discrete level data (DecayScheme objects), optical models
  /// (SphericalOpticalModel objects), and level density models
  /// (LevelDensityModel objects)
  class StructureDatabase {

    public:

      /// Creates an empty database.
      StructureDatabase();

      /// @brief Initializes the database using structure records from a
      /// ConfigFile object
      StructureDatabase(const marley::ConfigFile& cf);

      /// @brief Add discrete level data for a specific nuclide to the database
      /// @param nucid ENSDF-style nuclide indentifier (e.g., "40K") for the
      /// desired nuclide
      /// @param ds DecayScheme object containing the discrete level data to be
      /// added
      void add_decay_scheme(const std::string nucid,
        marley::DecayScheme ds);

      /// Removes all previously stored data from the database.
      void clear();

      /// @brief Deletes the discrete level data in the database associated
      /// with a given nuclide
      /// @details If a decay scheme does not exist in the database for the
      /// requested nuclide, then this function will return without making any
      /// changes.
      /// @param nucid ENSDF-style nuclide indentifier (e.g., "40Ar")
      void remove_decay_scheme(const std::string nucid);

      /// @brief Retrieves discrete level data from the database
      /// @param particle_id PDG particle ID for the desired nuclide
      /// @return pointer to the requested nuclide's DecayScheme object,
      /// or nullptr if one could not be found
      marley::DecayScheme* get_decay_scheme(const int particle_id);

      /// @brief Retrieves discrete level data from the database
      /// @param nucid ENSDF-style nuclide identifier (e.g., "39Ar")
      /// @return pointer to the requested nuclide's DecayScheme object,
      /// or nullptr if one could not be found
      marley::DecayScheme* get_decay_scheme(const std::string nucid);

      /// @brief Retrieves discrete level data from the database
      /// @param Z atomic number
      /// @param A mass number
      /// @return pointer to the requested nuclide's DecayScheme object,
      /// or nullptr if one could not be found
      marley::DecayScheme* get_decay_scheme(const int Z, const int A);

      /// @brief Retrieves an optical model object from the database, creating
      /// it if one did not already exist
      /// @param particle_id PDG particle ID for the desired nuclide
      marley::SphericalOpticalModel& get_optical_model(int nucleus_pid);

      /// @brief Retrieves an optical model object from the database, creating
      /// it if one did not already exist
      /// @param Z atomic number
      /// @param A mass number
      marley::SphericalOpticalModel& get_optical_model(const int Z,
        const int A);

      /// @brief Retrieves a level density model object from the database,
      /// creating it if one did not already exist
      /// @param Z atomic number
      /// @param A mass number
      marley::LevelDensityModel& get_level_density_model(const int Z,
        const int A);

    private:

      /// @brief Lookup table for marley::DecayScheme objects.
      /// @details Keys are ENSDF-style nucIDs, values are decay schemes.
      std::unordered_map<std::string, marley::DecayScheme> decay_scheme_table;

      /// @brief Table for looking up decay schemes by PDG particle ID
      std::unordered_map<int, marley::DecayScheme*> pid_decay_scheme_table;

      /// @brief Lookup table for marley::SphericalOpticalModel objects.
      /// @details Keys are PDG particle IDs, values are optical models.
      std::unordered_map<int, std::unique_ptr<marley::SphericalOpticalModel> >
        optical_model_table;

      /// @brief Lookup table for marley::LevelDensityModel objects.
      /// @details Keys are PDG particle IDs, values are unique_ptrs to level
      /// density models.
      std::unordered_map<int, std::unique_ptr<marley::LevelDensityModel> >
        level_density_table;

      /// @brief Helper function that loads a DecayScheme object from a single
      /// structure record in a ConfigFile
      void add_from_record(const marley::ConfigFile::StructureRecord& sr);

      /// @brief Helper function that loads DecayScheme objects from every
      /// structure record in a ConfigFile
      void add_all_from_config_file(const marley::ConfigFile& cf);
  };
}
