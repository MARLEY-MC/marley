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

// Standard library includes
#include <iterator>

// HepMC3 includes
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/EventFileReader.hh"

#ifdef USE_ROOT
  #include "TError.h"
  #include "TFile.h"
  #include "TTree.h"

  #include "HepMC3/Data/GenEventData.h"
#endif


marley::EventFileReader::EventFileReader( const std::string& file_name )
  : file_name_( file_name )
{
}

// Try to read an event from the file using each possible format. If we
// succeed, set the appropriate format code and return true. If all fail,
// return false.
bool marley::EventFileReader::deduce_file_format() {

  // Before we bother to check anything else, see if the file exists and is
  // readable. Skip the MARLEY search path in this case (the event file should
  // have been passed to the constructor with any needed path specification).
  // Complain if the file cannot be read.
  const auto& fm = marley::FileManager::Instance();
  if ( fm.find_file(file_name_, "").empty() ) throw marley::Error( "Could"
    " not read from the file \"" + file_name_ + '\"' );

  // Create a temporary event object to use for the following format checks
  HepMC3::GenEvent temp_event;

  #ifdef USE_ROOT
  // Before checking if the file is in ROOT format, completely suppress any
  // error messages from ROOT.
  auto temp_error_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  // Try to read in a HepMC3::GenEvent from the file assuming that the ROOT
  // output format was used
  tfile_ = std::unique_ptr<TFile>( TFile::Open(file_name_.c_str(), "read") );

  // We've completed the ROOT format check, so restore the old error messaging
  // behavior
  gErrorIgnoreLevel = temp_error_level;

  if ( tfile_ ) {
    format_ = marley::OutputFile::Format::ROOT;
    return true;
  }
  #endif

  // Try to read in a HepMC3::GenEvent from the file assuming that the ASCII
  // output format was used
  HepMC3::ReaderAscii temp_reader( file_name_ );
  bool read_ok = temp_reader.read_event( temp_event );
  if ( read_ok ) {

    // Save the flux-averaged total xsec from the run information for easy
    // retrieval
    auto run_info = temp_event.run_info();
    if ( !run_info ) {
      throw marley::Error( "Missing run information while parsing an"
        " ASCII-format HepMC3 file" );
    }

    this->get_flux_averaged_xsec( *run_info );

    format_ = marley::OutputFile::Format::ASCII;
    return true;
  }

  // TODO: add other formats here as needed

  // If everything else failed, then complain that events could not be read
  return false;
}

void marley::EventFileReader::initialize() {
  switch ( format_ ) {

    case marley::OutputFile::Format::ASCII:
    {
      in_.open( file_name_ );
      reader_ = std::make_shared< HepMC3::ReaderAscii >( in_ );
      break;
    }

    #ifdef USE_ROOT
    case marley::OutputFile::Format::ROOT:
    {
      tfile_->GetObject( "MARLEY_event_tree", ttree_ );
      if ( !ttree_ ) throw marley::Error( "Failed to load MARLEY event TTree"
        " from the ROOT file \"" + file_name_ + '\"' );

      temp_event_data_ = std::make_unique< HepMC3::GenEventData >();
      temp_event_data_ptr_ = temp_event_data_.get();
      ttree_->SetBranchAddress( "event", &temp_event_data_ptr_ );

      tfile_->GetObject( "MARLEY_run_info", temp_run_info_data_ );
      if ( !temp_run_info_data_ ) throw marley::Error( "Failed to load"
        " MARLEY run information from the ROOT file \"" + file_name_ + '\"' );

      run_info_ = std::make_shared< HepMC3::GenRunInfo >();
      run_info_->read_data( *temp_run_info_data_ );

      this->get_flux_averaged_xsec( *run_info_ );

      break;
    }
    #endif

    default:
      throw marley::Error( "Unrecognized file format encountered in"
        " marley::EventFileReader::initialize()" );
  }
}

bool marley::EventFileReader::next_event( HepMC3::GenEvent& ev )
{
  this->ensure_initialized();
  switch ( format_ ) {

    case marley::OutputFile::Format::ASCII:
    {
      bool read_ok = reader_->read_event( ev );
      return read_ok && in_;
      break;
    }

    #ifdef USE_ROOT
    case marley::OutputFile::Format::ROOT:
    {
      ++event_num_;
      if ( event_num_ < ttree_->GetEntries() ) {

        temp_event_data_->particles.clear();
        temp_event_data_->vertices.clear();
        temp_event_data_->links1.clear();
        temp_event_data_->links2.clear();
        temp_event_data_->attribute_id.clear();
        temp_event_data_->attribute_name.clear();
        temp_event_data_->attribute_string.clear();

        ttree_->GetEntry( event_num_ );
        ev.read_data( *temp_event_data_ );
        ev.set_run_info( run_info_ );
        return true;
      }

      ev = HepMC3::GenEvent();
      return false;

      break;
    }
    #endif

    default:
      throw marley::Error( "Unrecognized file format encountered in"
        " marley::EventFileReader::next_event()" );
  }

  ev = HepMC3::GenEvent();
  return false;
}

marley::EventFileReader::operator bool() const {

  switch ( format_ ) {

    case marley::OutputFile::Format::ASCII:
    {
      return static_cast< bool >( in_ );
      break;
    }

    #ifdef USE_ROOT
    case marley::OutputFile::Format::ROOT:
    {
      return ( tfile_ && ttree_ && event_num_ < ttree_->GetEntries() );
      break;
    }
    #endif

    default:
      throw marley::Error( "Unrecognized file format encountered in"
        " marley::EventFileReader::operator bool()" );
  }

  return false;
}

void marley::EventFileReader::ensure_initialized() {
  if ( !initialized_ ) {

    if ( !this->deduce_file_format() ) throw marley::Error( "Could not"
      " read MARLEY events from the file \"" + file_name_ + '\"' );

    this->initialize();
    initialized_ = true;
  }
}

double marley::EventFileReader::flux_averaged_xsec( bool natural_units ) {
  this->ensure_initialized();
  if ( natural_units ) return flux_avg_tot_xs_;
  double result = flux_avg_tot_xs_ * marley_utils::hbar_c2
    * marley_utils::fm2_to_minus40_cm2 * 1e2; // 10^{-42} cm^2
  return result;
}

void marley::EventFileReader::get_flux_averaged_xsec(
  const HepMC3::GenRunInfo& run_info )
{
  auto avg_xsec_attr = run_info.attribute< HepMC3::DoubleAttribute >(
    "NuHepMC.FluxAveragedTotalCrossSection" );
  if ( !avg_xsec_attr ) {
    throw marley::Error( "Missing flux-averaged total cross section while"
      " parsing a HepMC3 file" );
  }

  // Retrieve the flux-averaged total cross section and convert it back to
  // natural units (MeV^{-2}) from picobarn
  flux_avg_tot_xs_ = avg_xsec_attr->value()
    / ( marley_utils::hbar_c2 * marley_utils::fm2_to_picobarn );
}
