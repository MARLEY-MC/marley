// MARLEY includes
#include "marley/FileManager.hh"
#include "marley/EventFileReader.hh"

#include <iterator>

marley::EventFileReader::EventFileReader(
  const std::string& file_name) : file_name_(file_name),
  json_event_array_(),
  json_event_array_wrapper_( json_event_array_.array_range() ),
  json_event_iter_( json_event_array_wrapper_.end() )
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
  if ( fm.find_file(file_name_, "").empty() ) throw marley::Error("Could"
    " not read from the file \"" + file_name_ + '\"');

  // Temporarily turn off logging of marley::Error messages for the
  // try/catch blocks below. Otherwise, we'll end up with lots of noise
  // in the Logger from this function.
  bool log_marley_errors = marley::Error::logging_status();
  marley::Error::set_logging_status( false );

  // Create a temporary event object to use for the following format checks
  marley::Event temp_event;


  // If the first character in the file is '{', then
  // assume that the file is a JSON output file
  in_.open( file_name_ );
  char temp_char;
  if ( in_ >> temp_char && temp_char == '{' ) {
    format_ = marley::OutputFile::Format::JSON;

    // Reset the position of the input stream to before the initial '{'
    in_.seekg(0);
    return true;
  }

  // If we can read in a marley::Event from the file via the input stream
  // operator, then assume that the file is in ASCII format
  // Reopen the file to reset the state flags and the stream cursor position
  in_ = std::ifstream( file_name_ );
  try {
    if ( in_ >> flux_avg_tot_xs_ >> temp_event ) {
      format_ = marley::OutputFile::Format::ASCII;
      in_.seekg(0);
      return true;
    }
  }
  // Ignore any marley::Error exceptions that were thrown. These can
  // be expected if the file is not in ASCII format
  catch (const marley::Error&) { }

  // Now try reading in a marley::Event assuming that the file is in HEPEVT
  // format.
  in_ = std::ifstream( file_name_ );
  try {
    if ( temp_event.read_hepevt(in_, &flux_avg_tot_xs_) ) {
      format_ = marley::OutputFile::Format::HEPEVT;
      // Return to the start of the file so that we don't skip the first event
      in_.seekg(0);
      return true;
    }
  }
  // Ignore any marley::Error exceptions that were thrown. These can
  // be expected if the file is not in ASCII format
  catch (const marley::Error&) { }

  // Now that we're past the blocks, restore the old marley::Error logging
  // behavior
  marley::Error::set_logging_status( log_marley_errors );

  // If everything else failed, then complain that events could not be read
  return false;
}

void marley::EventFileReader::initialize() {
  switch ( format_ ) {

    case marley::OutputFile::Format::ASCII:
      in_ >> flux_avg_tot_xs_;
      break;

    case marley::OutputFile::Format::HEPEVT:
      // No further preparation is needed to read the HEPEVT format
      break;

    case marley::OutputFile::Format::JSON: {

      // Turn off auto-logging of marley::Error objects so that
      // we can add extra information via the try/catch block below
      bool log_marley_errors = marley::Error::logging_status();
      marley::Error::set_logging_status( false );

      try {
        auto json = marley::JSON::load( in_ );

        flux_avg_tot_xs_ = json.at("gen_state").at("flux_avg_xsec")
          .to_double();

        json_event_array_ = json.at("events");
        json_event_array_wrapper_ = json_event_array_.array_range();
        json_event_iter_ = json_event_array_wrapper_.begin();
        // Move the iterator to "one before the beginning" so that
        // we get the correct behavior in next_event() and
        // operator bool()
        --json_event_iter_;
      }
      catch (const std::exception& err) {
        // Rethrow the error after adding commentary
        MARLEY_LOG_ERROR() << "Parsing of a JSON-format output file"
          << " failed with error \"" << err.what() << '\"';
        throw err;
      }

      // Restore the previous marley::Error logging behavior now that
      // we're done with our mischief
      marley::Error::set_logging_status( log_marley_errors );
      break;
    }

    default:
      throw marley::Error("Unrecognized file format encountered in"
        " marley::EventFileReader::initialize()");
  }
}

bool marley::EventFileReader::next_event( marley::Event& ev )
{
  this->ensure_initialized();
  switch ( format_ ) {

    case marley::OutputFile::Format::ROOT:
      MARLEY_LOG_WARNING() << "Cannot read from ROOT file."
        << "Please use the marley::EventFileReader class and build"
        << " MARLEY with ROOT support.";
      return false;
      break;

    case marley::OutputFile::Format::ASCII:
      in_ >> ev;
      if ( in_ ) return true;
      break;

    case marley::OutputFile::Format::HEPEVT:
      if ( ev.read_hepevt(in_, &flux_avg_tot_xs_) ) return true;
      break;

    case marley::OutputFile::Format::JSON:
      ++json_event_iter_;
      if ( json_event_iter_ != json_event_array_wrapper_.end() ) {
        ev.from_json( *json_event_iter_ );
        return true;
      }
      break;

    default:
      throw marley::Error("Unrecognized file format encountered in"
        " marley::EventFileReader::next_event()");
  }

  ev = marley::Event();
  return false;
}

marley::EventFileReader::operator bool() const {

  switch ( format_ ) {

    case marley::OutputFile::Format::ROOT:
      return false;
      break;

    case marley::OutputFile::Format::ASCII:
    case marley::OutputFile::Format::HEPEVT:
      return static_cast<bool>( in_ );
      break;

    case marley::OutputFile::Format::JSON:
      return ( json_event_iter_ != json_event_array_wrapper_.end() );
      break;

    default:
      throw marley::Error("Unrecognized file format encountered in"
        " marley::EventFileReader::operator bool()");
  }

  return false;
}

void marley::EventFileReader::ensure_initialized() {
  if ( !initialized_ ) {
    if ( !this->deduce_file_format() ) throw marley::Error("Could not"
      " read MARLEY events from the file \"" + file_name_ + '\"');

    this->initialize();
    initialized_ = true;
  }
}
