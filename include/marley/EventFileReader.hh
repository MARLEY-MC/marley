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

// Standard library includes
#include <fstream>
#include <memory>
#include <string>

// MARLEY includes
#include "marley/OutputFile.hh"

#ifdef USE_ROOT
// HepMC3 includes
#include "HepMC3/Data/GenEventData.h"
#include "HepMC3/Data/GenRunInfoData.h"
#include "HepMC3/GenRunInfo.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#endif

// Forward-declare some HepMC3 classes
namespace HepMC3 {
  class GenEvent;
  class GenRunInfo;
  class ReaderAscii;
}

namespace marley {

  /// @brief Object that parses MARLEY output files
  class EventFileReader {

    public:

      EventFileReader( const std::string& file_name );

      virtual ~EventFileReader() = default;

      /// @brief Read the next MARLEY event record from the file
      /// @param[out] ev Reference to the object that will be filled with
      /// the next event record
      /// @return True if reading the next event was successful, or false
      /// otherwise. This behavior is designed to be used as a while loop
      /// condition for iterating over events in the output file
      virtual bool next_event( HepMC3::GenEvent& ev );

      /// @brief Returns the flux-averaged total cross section
      /// used to produce the events in the file
      /// @details For file formats which do not include this information,
      /// this function will return zero
      /// @param[in] natural_units If true, then this function will
      /// return the flux-averaged total cross section in natural units
      /// (MeV<sup> -2</sup>). If false (default), then 10<sup>-42</sup>
      /// cm<sup>2</sup> will be used.
      double flux_averaged_xsec( bool natural_units = false );

      /// @brief Stream operator for reading in the next event
      inline EventFileReader& operator>>( HepMC3::GenEvent& ev ) {
        next_event( ev );
        return *this;
      }

      /// @brief Implicit boolean conversion allows the state of the
      /// input stream (or ROOT file) to be tested for readiness to
      /// read in another event
      virtual operator bool() const;

    protected:

      /// @brief Name of the file (with any needed path specification) to be
      /// read
      std::string file_name_;

      /// @brief Format of the output file being read
      /// @details This format will be determined automatically by
      /// deduce_file_format() and does not need to be specified by the user
      OutputFile::Format format_;

      /// @brief Input stream used to read from textual output formats
      std::ifstream in_;

      /// @brief Helper object used to interpret ASCII-format HepMC3 files
      std::shared_ptr< HepMC3::ReaderAscii > reader_;

      #ifdef USE_ROOT
      /// @brief Helper object used to interpret ROOT-format HepMC3 files
      std::unique_ptr< TFile > tfile_;

      // @brief Pointer to the TTree containing the MARLEY events to be loaded
      TTree* ttree_ = nullptr;

      /// @brief The index of the last TTree entry that was read
      /// @details If no events have been read from the TTree yet, then
      /// this variable will have the value -1
      long event_num_ = -1;

      /// @brief Temporary storage for reading back HepMC3 events from
      /// a ROOT-format file
      std::unique_ptr< HepMC3::GenEventData > temp_event_data_;

      /// @brief Bare pointer to the event data needed for use when setting
      /// the branch address for ROOT I/O
      HepMC3::GenEventData* temp_event_data_ptr_ = nullptr;

      /// @brief Temporary storage for reading back HepMC3 run information from
      /// a ROOT-format file
      HepMC3::GenRunInfoData* temp_run_info_data_ = nullptr;

      /// @brief Stores the retrieved run information
      std::shared_ptr< HepMC3::GenRunInfo > run_info_;
      #endif

      /// @brief Flux-averaged total cross section (MeV<sup> -2</sup>) used to
      /// produce the events in the file, or zero if that information is not
      /// included in a particular format
      double flux_avg_tot_xs_ = 0.;

      /// @brief Flag that indicates whether initialize() has been called or
      /// not @details To avoid problems with using virtual functions in the
      /// constructor, we defer nearly all of the initialization to the first
      /// call to one of the other public member functions.
      bool initialized_ = false;

      /// @brief Helper function that auto-detects which of the available
      /// output formats is appropriate for the requested file
      virtual bool deduce_file_format();

      /// @brief Prepares the file for reading the events
      virtual void initialize();

      /// @brief This function should be called at the beginning of all public
      /// member functions of EventFileReader that interact with data in the
      /// file
      /// @details It provides necessary initialization as a workaround to
      /// calling virtual functions in the constructor
      void ensure_initialized();

      /// Helper function for loading the flux-averaged total cross section
      /// information from the HepMC3 run information
      void get_flux_averaged_xsec( const HepMC3::GenRunInfo& run_info );
  };

}
