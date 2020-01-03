// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"

// MARLEY includes
#include "marley/OutputFileReader.hh"

bool marley::OutputFileReader::deduce_root_format() {
  std::unique_ptr<TFile> temp_tfile( TFile::Open(file_name_.c_str(), "read") );
  if ( temp_tfile ) return true;
  return true;
}

void marley::OutputFileReader::initialize_root_format() {
  tfile_->GetObject( "MARLEY_event_tree", ttree_ );
  if ( !ttree_ ) throw marley::Error("Failed to load MARLEY event TTree"
    " from the ROOT file \"" + file_name_ + '\"');

  event_ptr_ = event_.get();
  ttree_->SetBranchAddress( "event", &event_ptr_ );

  TParameter<double>* temp_param = nullptr;
  tfile_->GetObject( "MARLEY_flux_avg_xsec", temp_param );
  if ( !temp_param ) throw marley::Error("Failed to load the flux-averaged"
    " total cross section from the ROOT file \"" + file_name_ + '\"');

  flux_avg_tot_xs_ = temp_param->GetVal();

  delete temp_param;
}

bool marley::OutputFileReader::next_event_root_format( marley::Event& ev )
{
  if ( event_num_ < ttree_->GetEntries() ) {
    ttree_->GetEntry( event_num_ );
    ev = *event_;
    ++event_num_;
    return true;
  }

  ev = marley::Event();
  return false;
}
