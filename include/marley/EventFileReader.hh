/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/LICENSE or
// visit http://opensource.org/licenses/GPL-3.0

#pragma once
#include <deque>
#include <fstream>
#include <string>

#include "marley/JSON.hh"
#include "marley/OutputFile.hh"

namespace marley {

  // Forward-declare the Event class
  class Event;

  /// @brief Object that parses MARLEY output files written in any of the
  /// available formats, except for ROOT format
  /// @details For a version of this class that can also handle ROOT-format
  /// output files, see the RootEventFileReader class
  class EventFileReader {

    public:

      EventFileReader( const std::string& file_name );

      virtual ~EventFileReader() = default;

      /// @brief Read the next MARLEY event record from the file
      /// @param ev Reference to the object that will be filled with
      /// the next event record
      /// @return True if reading the next event was successful, or false
      /// otherwise. This behavior is designed to be used as a while loop
      /// condition for iterating over events in the output file
      virtual bool next_event( marley::Event& ev );

      /// @brief Returns the flux-averaged total cross section
      /// (MeV<sup> -2</sup>) used to produce the events in the file
      /// @details For file formats which do not include this information,
      /// this function will return zero
      inline double flux_averaged_xsec() {
        this->ensure_initialized();
        return flux_avg_tot_xs_;
      }

      /// @brief Stream operator for reading in the next event
      inline EventFileReader& operator>>( marley::Event& ev ) {
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

      /// @brief Used to parse events from JSON-format files
      marley::JSON json_event_array_;
      /// @brief Used to iterate over events from JSON-format files
      marley::JSON::JSONWrapper< std::deque< marley::JSON > >
        json_event_array_wrapper_;
      /// @brief Iterator to the next JSON event in json_event_array
      std::deque<marley::JSON>::iterator json_event_iter_;

      /// @brief Flux-averaged total cross section
      /// (MeV<sup> -2</sup>) used to produce the events in the file,
      /// or zero if that information is not included in a particular
      /// format
      double flux_avg_tot_xs_ = 0.;

      /// @brief Flag that indicates whether initialize() has been called or not
      /// @details To avoid problems with using virtual functions in the constructor,
      /// we defer nearly all of the initialization to the first call to
      /// one of the other public member functions.
      bool initialized_ = false;

      /// @brief Helper function that auto-detects which of the available output
      /// formats is appropriate for the requested file
      virtual bool deduce_file_format();

      /// @brief Prepares the file for reading the events
      virtual void initialize();

      /// @brief This function should be called at the beginning of all public
      /// member functions of EventFileReader that interact with data in
      /// the file
      /// @details It provides necessary initialization as a workaround to
      /// calling virtual functions in the constructor
      void ensure_initialized();
  };

}
