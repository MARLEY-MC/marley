// ROOT includes
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"

// MARLEY includes
#include "marley/OutputFileReader.hh"

marley::OutputFileReader::OutputFileReader(const std::string& file_name)
  : marley::TextOutputFileReader(file_name)
{
}

bool marley::OutputFileReader::deduce_file_format() {

  // First delegate simpler format checks to TextOutputFileReader
  bool set_format_already =  marley::TextOutputFileReader::deduce_file_format();
  if ( set_format_already ) return true;

  // Check if the file is in ROOT format. Before doing so, completely suppress
  // any error messages from ROOT.
  auto temp_error_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;

  tfile_ = std::unique_ptr<TFile>( TFile::Open(file_name_.c_str(), "read") );

  // We've completed the ROOT format check, so restore the old error messaging
  // behavior
  gErrorIgnoreLevel = temp_error_level;

  // If the file was in ROOT format, then set the format_ data member
  // accordingly
  if ( tfile_ ) {
    format_ = marley::OutputFile::Format::ROOT;
    return true;
  }
  else return false;
}

void marley::OutputFileReader::initialize() {

  if ( format_ == marley::OutputFile::Format::ROOT ) {

    tfile_->GetObject( "MARLEY_event_tree", ttree_ );
    if ( !ttree_ ) throw marley::Error("Failed to load MARLEY event TTree"
      " from the ROOT file \"" + file_name_ + '\"');

    event_ = std::make_unique<marley::Event>();
    event_ptr_ = event_.get();
    ttree_->SetBranchAddress( "event", &event_ptr_ );

    TParameter<double>* temp_param = nullptr;
    tfile_->GetObject( "MARLEY_flux_avg_xsec", temp_param );
    if ( !temp_param ) throw marley::Error("Failed to load the flux-averaged"
      " total cross section from the ROOT file \"" + file_name_ + '\"');

    flux_avg_tot_xs_ = temp_param->GetVal();

    delete temp_param;
  }
  else return marley::TextOutputFileReader::initialize();
}

bool marley::OutputFileReader::next_event( marley::Event& ev )
{
  this->ensure_initialized();

  if ( format_ == marley::OutputFile::Format::ROOT ) {
    if ( event_num_ < ttree_->GetEntries() ) {
      ttree_->GetEntry( event_num_ );
      ev = *event_;
      ++event_num_;
      return true;
    }

    ev = marley::Event();
    return false;
  }
  else return marley::TextOutputFileReader::next_event( ev );
}

marley::OutputFileReader::operator bool() const {
  if ( format_ == marley::OutputFile::Format::ROOT ) {
    return ( tfile_ && ttree_ && event_num_ < ttree_->GetEntries() );
  }
  else return marley::TextOutputFileReader::operator bool();
}
