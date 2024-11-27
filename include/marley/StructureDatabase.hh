/// @file
/// @copyright Copyright (C) 2016-2024 Steven Gardiner
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
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>

#include "marley/DecayScheme.hh"
#include "marley/JSON.hh"

namespace marley {

  class GammaStrengthFunctionModel;
  class Fragment;
  class LevelDensityModel;
  class OpticalModel;

  /// @brief Container for nuclear structure information organized by nuclide
  /// @details Currently, the StructureDatabase object can hold nuclear
  /// discrete level data (DecayScheme objects), optical models (OpticalModel
  /// objects), @f$\gamma@f$-ray strength function models
  /// (GammaStrengthFunctionModel objects), and level density models
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

      /// @brief Add a DecayScheme object to the database that
      /// contains discrete level data for a specific nuclide
      /// @param pdg PDG code for the desired nuclide
      /// @param ds A unique_ptr to DecayScheme object to move into the
      /// database.
      void add_decay_scheme(int pdg, std::unique_ptr<marley::DecayScheme>& ds);

      /// @brief Create a set of Particle Data Group codes for every
      /// nuclide in a discrete level data file
      /// @param filename Name of the file that contains the discrete level
      /// data
      /// @param format DecayScheme::FileFormat specifier that indicates
      /// which nuclear data format is used in the file
      std::set<int> find_all_nuclides(const std::string& filename,
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
      marley::OpticalModel& get_optical_model(int nucleus_pid);

      /// @brief Retrieves an optical model object from the database, creating
      /// it if one did not already exist
      /// @param Z atomic number
      /// @param A mass number
      marley::OpticalModel& get_optical_model(const int Z,
        const int A);

      /// @brief Retrieves a level density model object from the database,
      /// creating it if one did not already exist
      /// @param nucleus_pid PDG particle ID for the desired nucleus
      marley::LevelDensityModel& get_level_density_model(const int nucleus_pid);

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

      /// @brief Retrieves a gamma-ray strength function model object from the
      /// database, creating it if one did not already exist
      /// @param nuc_pdg PDG code for the nuclide of interest
      marley::GammaStrengthFunctionModel& get_gamma_strength_function_model(
        const int nuc_pdg);

      /// Retrieves a const reference to the table of DecayScheme objects
      inline const std::unordered_map<int,
        std::unique_ptr<marley::DecayScheme> >& decay_schemes() const
      {
        return decay_scheme_table_;
      }

      /// Retrieves a const reference to the table of Fragment objects
      static inline const std::map<int, marley::Fragment>& fragments() {
        if ( !initialized_gs_spin_parity_table_ ) initialize_jpi_table();
        return fragment_table_;
      }

      /// @brief Retrieves nuclear fragment data from the database
      /// @param fragment_pdg PDG code for the desired fragment
      /// @return Pointer to the requested Fragment object, or nullptr
      /// if it could not be found
      static const marley::Fragment* get_fragment(const int fragment_pdg);

      /// @brief Retrieves nuclear fragment data from the database
      /// @param Z Fragment proton number
      /// @param A Fragment mass number
      /// @return Pointer to the requested Fragment object, or nullptr if
      /// it could not be found
      static const marley::Fragment* get_fragment(const int Z, const int A);

      /// @brief Returns the maximum orbital angular momentum to consider
      /// when simulating fragment emission to the continuum
      inline int get_fragment_l_max() const { return fragment_l_max_; }

      /// @brief Returns the maximum multipolarity to consider
      /// when simulating gamma-ray emission to the continuum
      inline int get_gamma_l_max() const { return gamma_l_max_; }

      /// @brief Sets the maximum orbital angular momentum to consider
      /// when simulating fragment emission to the continuum
      inline void set_fragment_l_max( int ell ) { fragment_l_max_ = ell; }

      /// @brief Sets the maximum multipolarity to consider when simulating
      /// gamma-ray emission to the continuum
      inline void set_gamma_l_max( int ell ) { gamma_l_max_ = ell; }

      /// @brief Looks up the ground-state spin-parity for a particular nuclide
      /// @param[in] nuc_pdg PDG code for the nuclide of interest
      /// @param[out] twoJ Two times the ground-state nuclear spin
      /// @param[out] Pi Ground-state nuclear parity
      static void get_gs_spin_parity(int nuc_pdg, int& twoJ,
        marley::Parity& Pi);

      /// @brief Looks up the ground-state spin-parity for a particular nuclide
      /// @param[in] Z Proton number for the nuclide of interest
      /// @param[in] A Mass number for the nuclide of interest
      /// @param[out] twoJ Two times the ground-state nuclear spin
      /// @param[out] Pi Ground-state nuclear parity
      static void get_gs_spin_parity(const int Z, const int A, int& twoJ,
        marley::Parity& Pi);

      /// @brief Helper function that initializes the map of JSON settings
      /// for the optical model parameters
      void load_optical_model_params( const marley::JSON* om_config = nullptr );

    private:

      /// @brief Lookup table for marley::DecayScheme objects.
      /// @details Keys are PDG codes, values are unique_ptrs to decay schemes.
      static std::unordered_map<int, std::unique_ptr<marley::DecayScheme> >
        decay_scheme_table_;

      /// @brief Lookup table for marley::OpticalModel objects.
      /// @details Keys are PDG codes, values are unique_ptrs to optical models.
      std::unordered_map<int, std::unique_ptr<marley::OpticalModel> >
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

      /// @brief Lookup table for nuclear fragments that will be considered
      /// when modeling de-excitations in the unbound continuum
      static std::map<int, marley::Fragment> fragment_table_;

      /// @brief Default value of fragment_l_max_
      static constexpr int DEFAULT_FRAGMENT_L_MAX = 5;

      /// @brief Default value of gamma_l_max_
      static constexpr int DEFAULT_GAMMA_L_MAX = DEFAULT_FRAGMENT_L_MAX;

      /// @brief Maximum orbital angular momentum @f$ \ell @f$ to consider
      /// when computing differential widths (via a FragmentContinuumExitChannel
      /// object) for decays to the unbound continuum via fragment emission
      int fragment_l_max_ = DEFAULT_FRAGMENT_L_MAX;

      /// @brief Maximum multipolarity @f$ \ell @f$ to consider
      /// when computing differential widths (via a GammaContinuumExitChannel
      /// object) for decays to the unbound continuum via gamma-ray emission
      int gamma_l_max_ = DEFAULT_GAMMA_L_MAX;

      /// @brief Flag that indicates whether the ground-state spin-parities
      /// have already been loaded from the relevant data file
      static bool initialized_gs_spin_parity_table_;

      /// @brief Name of the file used as an index for nuclear structure data
      const std::string structure_index_filename_ = "nuclide_index.txt";

      /// @brief Flag that indicates whether the index to the structure
      /// data files has been loaded or not
      bool loaded_structure_index_ = false;

      /// @brief Lookup table for decay scheme data files. Keys are nuclide
      /// PDG codes, values are data file names
      std::map< int, std::string > decay_scheme_filenames_;

      /// @brief Name of the data file which will be used to read in the
      /// ground-state nuclear spin-parity values
      static const std::string jpi_data_file_name_;

      /// @brief Lookup table for ground-state nuclear spin-parities
      static std::map< int, std::pair<int, marley::Parity> > jpi_table_;

      /// @brief Helper function that initializes the table of ground-state
      /// nuclear spin-parities
      static void initialize_jpi_table();

      /// @brief Helper function that initializes the file index for
      /// loading nuclear structure data
      void load_structure_index();

      /// @brief Storage for JSON configuration settings for optical model
      /// parameters
      std::map< std::string, marley::JSON > om_config_map_;
  };

}
