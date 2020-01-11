#pragma once

// Standard library includes
#include <string>

namespace marley {

  // Forward-declare the Event class
  class Event;

  class R5EFR {

    public:

      R5EFR( const std::string& file_name );

      ~R5EFR();

      bool next_event(marley::Event& ev);

      operator bool() const;

      /// @brief Stream operator for reading in the next event
      inline R5EFR& operator>>( marley::Event& ev ) {
        next_event( ev );
        return *this;
      }

    protected:

      void* event_file_reader_;
  };

}
