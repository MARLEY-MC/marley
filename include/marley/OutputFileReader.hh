#pragma once
#include <memory>
#include <string>

#include "marley/Event.hh"
#include "marley/Error.hh"
#include "marley/TextOutputFileReader.hh"

#include "TFile.h"
#include "TTree.h"

namespace marley {

  /// @brief Object that parses MARLEY output files written in any of the
  /// available formats, including ROOT format
  class OutputFileReader : public TextOutputFileReader {

    public:

      inline OutputFileReader( const std::string& file_name )
        : TextOutputFileReader( file_name ) {}

    protected:

      // TFile used to read ROOT format output files
      std::unique_ptr<TFile> tfile_;
      // Pointer to the TTree containing the MARLEY events to be loaded
      TTree* ttree_;

      /// @brief Temporary storage for reading events in from a TFile
      std::unique_ptr<marley::Event> event_;
      marley::Event* event_ptr_;

      long event_num_ = 0; ///< The index of the next TTree entry to be read

      /// @return True if file format deduction was successful, or false
      /// otherwise
      virtual bool deduce_root_format() override;
      virtual void initialize_root_format() override;
      virtual bool next_event_root_format(marley::Event& ev) override;
  };


}
