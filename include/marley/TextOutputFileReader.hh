#pragma once
#include <string>

#include "marley/Event.hh"
#include "marley/Error.hh"
#include "marley/JSON.hh"
#include "marley/OutputFile.hh"

namespace marley {

  /// @brief Object that parses MARLEY output files written in any of the
  /// available formats
  class TextOutputFileReader {

    public:

      TextOutputFileReader( const std::string& file_name );

      /// @brief Read the next MARLEY event record from the file
      /// @param ev Reference to the object that will be filled with
      /// the next event record
      /// @return True if reading the next event was successful, or false
      /// otherwise. This behavior is designed to be used as a while loop
      /// condition for iterating over events in the output file
      bool next_event( marley::Event& ev );

      inline double flux_averaged_xsec() const { return flux_avg_tot_xs_; }

    protected:

      std::string file_name_;
      OutputFile::Format format_;

      std::ifstream in_;

      marley::JSON json_;

      //json_config_ = json_.at("gen_state").at("config");
      //seed_ = json_.at("gen_state").at("seed").to_long();
      //state_string_ = json_.at("gen_state").at("generator_state_string")
      //  .to_string();

      //std::string* temp_str = nullptr;
      //tfile_->GetObject( "MARLEY_config", temp_str );
      //if ( !temp_str ) throw marley::Error("Failed to load JSON configuration"
      //  " from the ROOT file \"" + file_name_ + '\"');

      //json_config_.load( *temp_str );

      //delete temp_str;
      //tfile_->GetObject( "MARLEY_state", temp_str );
      //if ( !temp_str ) throw marley::Error("Failed to load the generator state"
      //  " string from the ROOT file \"" + file_name_ + '\"');

      //state_string_ = *temp_str;

      //delete temp_str;
      //tfile_->GetObject( "MARLEY_seed", temp_str );
      //if ( !temp_str ) throw marley::Error("Failed to load the random number"
      //  " seed from the ROOT file \"" + file_name_ + '\"');

      //std::stringstream temp_ss( *temp_str );
      //temp_ss >> seed_;

      //delete temp_str;

      double flux_avg_tot_xs_ = 0.;

      void deduce_file_format();
      void initialize();

      inline virtual bool deduce_root_format() { return false; }
      inline virtual void initialize_root_format() {}
      virtual bool next_event_root_format( marley::Event& ev );
  };


}
