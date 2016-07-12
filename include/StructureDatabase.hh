#pragma once
#include <memory>
#include <string>
#include <unordered_map>

#include "BackshiftedFermiGasModel.hh"
#include "ConfigurationFile.hh"
#include "DecayScheme.hh"
#include "GammaStrengthFunctionModel.hh"
#include "SphericalOpticalModel.hh"

namespace marley {

  /// @brief Container for nuclear structure information organized by nuclide
  /// @details Currently, the StructureDatabase object can hold nuclear
  /// discrete level data (DecayScheme objects), optical models
  /// (SphericalOpticalModel objects), and level density models
  /// (LevelDensityModel objects)
  class StructureDatabase {

    public:

      /// @brief Creates an empty database.
      StructureDatabase();

      /// @brief Construct and add a DecayScheme object to the database that
      /// contains discrete level data for a specific nuclide
      /// @param pdg PDG code for the desired nuclide
      /// @param filename Name of the file that contains the discrete level data
      /// @param format DecayScheme::FileFormat specifier that indicates
      /// which nuclear data format is used in the file
      void emplace_decay_scheme(int pdg, const std::string& filename,
        DecayScheme::FileFormat format = DecayScheme::FileFormat::talys);

      /// @brief Removes all previously stored data from the database.
      void clear();

      /// @brief Deletes the discrete level data in the database associated
      /// with a given nuclide
      /// @details If a decay scheme does not exist in the database for the
      /// requested nuclide, then this function will return without making any
      /// changes.
      /// @param pdg PDG code for the nuclide to remove
      void remove_decay_scheme(int pdg);

      /// @brief Retrieves discrete level data from the database
      /// @param particle_id PDG code for the desired nuclide
      /// @return pointer to the requested nuclide's DecayScheme object,
      /// or nullptr if one could not be found
      marley::DecayScheme* get_decay_scheme(const int particle_id);

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

      /// @brief Retrieves a gamma-ray strength function model object from the
      /// database, creating it if one did not already exist
      /// @param Z atomic number
      /// @param A mass number
      marley::GammaStrengthFunctionModel& get_gamma_strength_function_model(
        const int Z, const int A);

    private:

      /// @brief Lookup table for marley::DecayScheme objects.
      /// @details Keys are PDG codes, values are unique_ptrs to decay schemes.
      std::unordered_map<int, std::unique_ptr<marley::DecayScheme> >
        decay_scheme_table_;

      /// @brief Lookup table for marley::SphericalOpticalModel objects.
      /// @details Keys are PDG codes, values are unique_ptrs to optical models.
      std::unordered_map<int, std::unique_ptr<marley::SphericalOpticalModel> >
        optical_model_table_;

      /// @brief Lookup table for marley::LevelDensityModel objects.
      /// @details Keys are PDG codes, values are unique_ptrs to level
      /// density models.
      std::unordered_map<int, std::unique_ptr<marley::LevelDensityModel> >
        level_density_table_;

      /// @brief Lookup table for marley::GammaStrengthFunctionModel objects.
      /// @details Keys are PDG codes, values are unique_ptrs to gamma-ray
      /// strength function models.
      std::unordered_map<int, std::unique_ptr<
        marley::GammaStrengthFunctionModel> > gamma_strength_function_table_;
  };
}
