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

// standard library includes
#include <string>

// MARLEY includes
#include "marley/CoulombCorrector.hh"
#include "marley/Generator.hh"
#include "marley/InterpolationGrid.hh"
#include "marley/JSON.hh"

namespace marley {

  enum class SubContinuumMode;

  class JSONConfig {

    public:

      using InterpMethod = InterpolationGrid<double>
        ::InterpolationMethod;

      explicit JSONConfig( const JSON& object );
      explicit JSONConfig( const std::string& json_filename );

      Generator create_generator() const;

      void prepare_direction( Generator& gen ) const;
      void prepare_neutrino_source( Generator& gen ) const;
      void prepare_reactions( Generator& gen,
        CoulombCorrector::CoulombMode coulomb_mode,
        const JSON& ff_config, const SubContinuumMode sc_mode ) const;
      void prepare_structure( Generator& gen ) const;
      void prepare_target( Generator& gen ) const;
      void prepare_weights( Generator& gen ) const;

      InterpMethod get_interpolation_method( const std::string& rule ) const;
      int neutrino_pdg( const std::string& nu ) const;

      /// @brief Helper function used to define ROOT-based neutrino source
      /// types
      /// @detail This function is a no-op when MARLEY is built without
      /// ROOT support
      bool process_extra_source_types( const std::string& type,
        const JSON& source_spec, int pdg_code,
        std::unique_ptr< NeutrinoSource >& source ) const;

      /// Helper function that checks whether an input JSON object representing
      /// the form factor configuration corresponds to a valid request for the
      /// allowed approximation to be used
      static bool check_for_allowed_approximation(
        const JSON& ff_config );

      inline const JSON& get_json() const;
      inline void set_json( const JSON& json );

      static void handle_json_error( const std::string& name,
        const JSON& json );

    protected:

      /// @brief Helper function for loading strings from the JSON
      /// configuration
      std::string source_get( const char* name, const JSON& source_spec,
        const char* description, const char* default_str ) const;

      /// @brief JSON object describing this configuration
      JSON json_;
  };

  inline const JSON& JSONConfig::get_json() const { return json_; }

  inline void JSONConfig::set_json( const JSON& json ) { json_ = json; }

}
