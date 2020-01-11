#pragma once
#include <memory>
#include <string>

#include "marley/Event.hh"
#include "marley/Error.hh"
#include "marley/EventFileReader.hh"

#include "TFile.h"
#include "TTree.h"

namespace marley {

  /// @brief Object that parses MARLEY output files written in any of the
  /// available formats, including ROOT format
  class RootEventFileReader : public EventFileReader {

    public:

      RootEventFileReader( const std::string& file_name );

      virtual ~RootEventFileReader() = default;

      virtual bool next_event(marley::Event& ev) override;

      virtual operator bool() const override;

    protected:

      /// @brief TFile used to read ROOT format output files
      std::unique_ptr<TFile> tfile_;

      // @brief Pointer to the TTree containing the MARLEY events to be loaded
      TTree* ttree_;

      /// @brief Temporary storage for reading events in from a TFile
      std::unique_ptr<marley::Event> event_;

      /// @brief Bare pointer used to interface the event_ data member
      /// with a branch of ttree_ without requiring a manual delete statement
      /// anywhere
      marley::Event* event_ptr_;

      /// @brief The index of the last TTree entry that was read
      /// @details If no events have been read from the TTree yet, then
      /// this variable will have the value -1
      long event_num_ = -1;

      virtual bool deduce_file_format() override;
      virtual void initialize() override;
  };

}
