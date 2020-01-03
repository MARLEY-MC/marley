// MARLEY includes
#include "marley/TextOutputFileReader.hh"

marley::TextOutputFileReader::TextOutputFileReader(
  const std::string& file_name) : file_name_(file_name)
{
  this->deduce_file_format();
  this->initialize();
}

// Try to read an event from the file using each possible format. If we
// succeed, set the appropriate format code. If all fail, complain by throwing
// a marley::Error.
void marley::TextOutputFileReader::deduce_file_format() {

  // If the first character in the file is '{', then
  // assume that the file is a JSON output file

  std::ifstream in_stream( file_name_ );
  char temp_char;
  if ( in_stream >> temp_char && temp_char == '{' ) {
    format_ = marley::OutputFile::Format::JSON;
    return;
  }

  // If we can read in a marley::Event from the file via the input stream
  // operator, then assume that the file is in ASCII format
  in_stream.close();
  double flux_avg_xsec;
  marley::Event temp_event;
  in_stream.open( file_name_ );
  if ( in_stream >> flux_avg_xsec >> temp_event ) {
    format_ = marley::OutputFile::Format::ASCII;
    return;
  }

  // Now try reading in a marley::Event assuming that the file is in HEPEVT
  // format.
  in_stream.close();
  in_stream.open( file_name_ );
  if ( temp_event.read_hepevt(in_stream) ) {
    format_ = marley::OutputFile::Format::HEPEVT;
    return;
  }

  // None of the other formats worked, so try ROOT format if it is available
  if ( this->deduce_root_format() ) {
    format_ = marley::OutputFile::Format::ROOT;
    return;
  }

  throw marley::Error("Could not read MARLEY events from the file "
    + file_name_ );
}

void marley::TextOutputFileReader::initialize() {
  switch ( format_ ) {

    case marley::OutputFile::Format::ROOT:
      this->initialize_root_format();
      break;

    case marley::OutputFile::Format::ASCII:
      in_ >> flux_avg_tot_xs_;
      break;

    case marley::OutputFile::Format::HEPEVT:
      // TODO: add xsec extraction here
      // in_ >> flux_avg_tot_xs_;
      break;

    case marley::OutputFile::Format::JSON:
      // TODO: add error handling here
      json_ = marley::JSON::load_file( file_name_ );
      flux_avg_tot_xs_ = json_.at("gen_state").at("flux_avg_xsec")
        .to_double();
      break;

    default:
      throw marley::Error("Unrecognized file format encountered in"
        " marley::TextOutputFileReader::initialize()");
  }
}

bool marley::TextOutputFileReader::next_event( marley::Event& ev )
{
  switch ( format_ ) {

    case marley::OutputFile::Format::ROOT:
      this->next_event_root_format( ev );
      break;

    case marley::OutputFile::Format::ASCII:
      in_ >> ev;
      if ( in_ ) return true;
      break;

    case marley::OutputFile::Format::HEPEVT:
      if ( ev.read_hepevt(in_) ) return true;
      break;

    case marley::OutputFile::Format::JSON:
      // TODO: make this
      if ( true ) return true;
      break;

    default:
      throw marley::Error("Unrecognized file format encountered in"
        " marley::TextOutputFileReader::next_event()");
  }

  ev = marley::Event();
  return false;
}

bool marley::TextOutputFileReader::next_event_root_format(
  marley::Event& /*ev*/)
{
  MARLEY_LOG_WARNING() << "Cannot read from ROOT file."
    << "Please use the marley::OutputFileReader class and build"
    << " MARLEY with ROOT support.";
  return false;
}
