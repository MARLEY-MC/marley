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
#include <memory>

// MARLEY includes
#include "marley/OutputFile.hh"

namespace HepMC3 {
  class GenEvent;
  class ReaderAscii;
  class WriterAscii;
}

namespace marley {

  /// @brief Class that represents an output file in the standard ASCII format
  /// for HepMC3 events
  class OutputFileAscii : public OutputFile {

    public:

      OutputFileAscii( const JSON& output_config );

      virtual ~OutputFileAscii() = default;

      virtual bool resume( std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events ) override;

      virtual int_fast64_t bytes_written() override;

      virtual void write_event( HepMC3::GenEvent* event ) override;

    protected:

      /// Opens the owned std::fstream with the correct settings
      virtual void open();

      /// Stream used to read and write from the output file as needed
      std::fstream stream_;

      /// Helper object used to read back HepMC3 events for restoring
      /// the generator state, etc.
      std::shared_ptr< HepMC3::ReaderAscii > reader_;

      /// Helper object used to produce standard ASCII HepMC3 output
      std::shared_ptr< HepMC3::WriterAscii > writer_;

      /// Storage for the number of bytes written to disk
      int_fast64_t byte_count_ = 0;

      /// Stream position at which to truncate on the first write after
      /// resume, removing the stale GeneratorState and HepMC3 footer
      std::streampos pending_truncate_pos_ = std::streampos( -1 );

      /// Last event without its GeneratorState, to be flushed on the first
      /// write() call after resume
      std::shared_ptr< HepMC3::GenEvent > pending_flush_event_;
  };

}
