/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
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

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/Error.hh"
#include "marley/JSONConfig.hh"
#include "marley/OutputFile.hh"

#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#endif

marley::OutputFile::OutputFile(const std::string& name,
  const std::string& format, const std::string& mode,
  bool force) : name_(name), force_(force)
{
  if (format == "root") format_ = Format::ROOT;
  else if (format == "hepevt") format_ = Format::HEPEVT;
  else if (format == "json") format_ = Format::JSON;
  else if (format == "ascii") format_ = Format::ASCII;
  else throw marley::Error("Invalid output file format \"" + format
    + "\" given in an output file specification");

  if (mode == "overwrite") mode_ = Mode::OVERWRITE;
  else if (mode == "append") {
    if (format_ == Format::HEPEVT || format_ == Format::ASCII)
      mode_ = Mode::APPEND;
    else throw marley::Error("The output mode \"" + mode + "\" is not"
      " allowed for the file format \"" + format + '\"');
  }
  else if (mode == "resume") {
    if (format_ == Format::ROOT || format_ == Format::JSON)
      mode_ = Mode::RESUME;
    else throw marley::Error("The output mode \"" + mode + "\" is not"
      " allowed for the file format \"" + format + '\"');
  }
  else throw marley::Error("Invalid output mode \"" + mode
    + "\" given in an output file specification");
}

std::unique_ptr<marley::Generator> marley::OutputFile::restore_generator(
  const marley::JSON& config)
{
  #ifdef USE_ROOT
    marley::RootJSONConfig jc( config );
  #else
    marley::JSONConfig jc( config );
  #endif
  auto gen = std::make_unique<marley::Generator>( jc.create_generator() );
  return gen;
}

marley::TextOutputFile::TextOutputFile(const std::string& name,
  const std::string& format, const std::string& mode, bool force, int indent)
  : marley::OutputFile(name, format, mode, force), indent_(indent)
{
  this->open();
}

void marley::TextOutputFile::start_json_output(bool start_array) {
  if (format_ != Format::JSON) throw marley::Error("TextOutputFile"
    "::start_json_output() called for a non-JSON file format");

  stream_ << '{';
  if (indent_ >= 0) {
    stream_ << '\n';
    for (int k = 0; k < indent_; ++k) stream_ << ' ';
  }
  stream_ << "\"events\"";
  if (indent_ >= 0) stream_ << ' ';
  stream_ << ':';
  if (indent_ >= 0) stream_ << ' ';

  if (!start_array) return;

  stream_ << '[';
  if (indent_ >= 0) {
    stream_ << '\n';
    for (int k = 0; k < 2*indent_; ++k) stream_ << ' ';
  }
}

void marley::TextOutputFile::open() {
  bool file_exists = check_if_file_exists(name_);

  auto open_mode_flag = std::ios::out | std::ios::trunc;

  if (mode_ == Mode::OVERWRITE && file_exists && !force_) {
    bool overwrite = marley_utils::prompt_yes_no("Overwrite file "
      + name_);
    if (!overwrite) {
      MARLEY_LOG_INFO() << "Cancelling overwrite of output file \""
        << name_ << '\"';
      open_mode_flag = std::ios::app | std::ios::out;
      if (format_ == Format::JSON) mode_ = Mode::RESUME;
      else mode_ = Mode::APPEND;
    }
  }

  if (mode_ == Mode::RESUME) {
    if (format_ != Format::JSON) throw marley::Error("The"
      " \"resume\" mode may only be used with the ROOT and JSON output"
      " formats");
    if (!file_exists) throw marley::Error("Cannot resume run. Could"
      " not open the JSON file \"" + name_ + '\"');
    else open_mode_flag = std::ios::app | std::ios::out;
  }
  else if (mode_ == Mode::APPEND)
    open_mode_flag = std::ios::app | std::ios::out;
  else if (mode_ != Mode::OVERWRITE)
    throw marley::Error("Unrecognized file mode encountered in"
      " TextOutputFile::open()");

  stream_.open(name_, open_mode_flag);

  // Get the event array started if we're writing a fresh JSON file
  if (format_ == Format::JSON && mode_ == Mode::OVERWRITE) {
    start_json_output(true);
  }

}

bool marley::TextOutputFile::resume(std::unique_ptr<marley::Generator>& gen,
  long& num_previous_events)
{
  if (mode_ != Mode::RESUME) {
    throw marley::Error("Cannot call TextOutput"
      "File::resume() for an output mode other than \"resume\"");
    return false;
  }

  if (format_ != Format::JSON) {
    throw marley::Error("Cannot call"
      " TextOutputFile::resume() for an output format other"
      " than \"json\"");
    return false;
  }

  MARLEY_LOG_INFO() << "Continuing previous run from JSON file "
    << name_;

  // Get the JSON objects from the file
  stream_.close();
  stream_.open(name_, std::ios::in);
  marley::JSON temp_json = marley::JSON::load(stream_);
  stream_.close();

  if (!temp_json.has_key("gen_state")) {
    throw marley::Error("Missing generator configuration in JSON"
      " file \"" + name_ + "\": could not restore previous state");
    return false;
  }
  const marley::JSON& gen_state = temp_json.at("gen_state");

  if (!gen_state.has_key("config")) {
    throw marley::Error("Failed to load previous configuration from"
      " the JSON file \"" + name_ + '\"');
    return false;
  }
  const marley::JSON& config = gen_state.at("config");

  if (!gen_state.has_key("generator_state_string")) {
    throw marley::Error("Failed to load previous generator state from"
      " the JSON file \"" + name_ + '\"');
    return false;
  }
  std::string state_string
    = gen_state.at("generator_state_string").to_string();

  if (!gen_state.has_key("seed")) {
    throw marley::Error("Failed to load previous random number"
      " generator seed from the JSON file \"" + name_ + '\"');
    return false;
  }
  std::string seed = gen_state.at("seed").to_string();

  bool count_ok = true;
  if (!gen_state.has_key("event_count")) count_ok = false;
  else num_previous_events = gen_state.at(
    "event_count").to_long(count_ok);

  if (!count_ok) {
    throw marley::Error("Failed to load previous event count"
      " from the JSON file \"" + name_ + '\"');
    return false;
  }

  gen = this->restore_generator( config );
  gen->seed_using_state_string( state_string );

  MARLEY_LOG_INFO() << "The previous run was initialized using"
    << " the random number generator seed " << seed;

  // We've loaded all the metadata we need, so erase the file,
  // and write out all the previous events to it again.
  stream_.open(name_, std::ios::out | std::ios::trunc);

  start_json_output(true);

  const auto& evt_array = temp_json.at("events").array_range();

  // Pretty-print the events one-by-one as needed
  auto begin = evt_array.begin();
  auto end = evt_array.end();

  for (auto iter = begin; iter != end; ++iter) {
    if (iter != begin) stream_ << ',';
    if (indent_ < 0) stream_ << iter->dump_string();
    else {
      stream_ << '\n';
      for (int i = 0; i < 2*indent_; ++i) stream_ << ' ';
      iter->print(stream_, indent_, true, 2*indent_);
    }
  }

  // Unless the event array is empty, add a comma before continuing
  // to write events to the file
  if (begin != end) needs_comma_ = true;
  else needs_comma_ = false;

  // Close the file and re-open it, this time appending to the end
  stream_.close();
  stream_.open(name_, std::ios::out | std::ios::app);

  return true;
}

int_fast64_t marley::TextOutputFile::bytes_written() {
  // If the stream is open, then update the byte count. Otherwise, just
  // use the saved value.
  if (stream_.is_open()) {
    stream_.flush();
    byte_count_ = static_cast<int_fast64_t>( stream_.tellp() );
  }
  return byte_count_;
}

void marley::TextOutputFile::write_event(const marley::Event* event) {
  if (!event) throw marley::Error("Null pointer passed to"
    " TextOutputFile::write_event()");

  switch (format_) {
    case Format::ASCII:
      stream_ << *event;
      break;
    case Format::JSON:
      if (needs_comma_) {
        stream_ << ',';
        if (indent_ > 0) {
          stream_ << '\n';
          for (int k = 0; k < 2*indent_; ++k) stream_ << ' ';
        }
      }
      else needs_comma_ = true;
      if (indent_ == -1) stream_ << event->to_json();
      else {
        event->to_json().print(stream_, indent_, true, 2*indent_);
      }
      break;
    case Format::HEPEVT:
      // TODO: consider incrementing event numbers each time instead of
      // just writing a zero
      event->write_hepevt(0, flux_avg_tot_xsec_, stream_);
      break;
    case Format::ROOT:
      throw marley::Error("ROOT format encountered in TextOutputFile::"
        "write_event()");
      break;
    default:
      throw marley::Error("Invalid format value encountered in"
        " TextOutputFile::write_event()");
  }
}

void marley::TextOutputFile::write_generator_state(
  const marley::JSON& json_config, const marley::Generator& gen,
  const long num_events)
{
  if (format_ != Format::JSON) throw marley::Error("TextOutputFile::"
    "write_generator_state() should only be used with the JSON format");

  stream_ << ',';
  if (indent_ > 0) {
    stream_ << '\n';
    for (int k = 0; k < indent_; ++k) stream_ << ' ';
  }
  stream_ << "\"gen_state\"";
  if (indent_ > 0) stream_ << ' ';
  stream_ << ':';
  if (indent_ > 0) stream_ << ' ';

  marley::JSON temp = marley::JSON::object();

  temp["config"] = json_config;
  temp["generator_state_string"] = gen.get_state_string();
  temp["seed"] = std::to_string(gen.get_seed());
  temp["event_count"] = num_events;
  temp["flux_avg_xsec"] = gen.flux_averaged_total_xs();

  if (indent_ < 0) stream_ << temp.dump_string();
  else temp.print(stream_, indent_, true, indent_);
}

void marley::TextOutputFile::close(const marley::JSON& json_config,
  const marley::Generator& gen, const long num_events)
{
  if (format_ == Format::JSON) {
    // End the JSON array of event objects
    if (indent_ > 0) {
      stream_ << '\n';
      for (int k = 0; k < indent_; ++k) stream_ << ' ';
    }
    stream_ << ']';

    // Save the current state of the generator to the JSON file in case
    // we want to resume a run later
    write_generator_state(json_config, gen, num_events);

    // Terminate the JSON file with a closing curly brace
    if (indent_ > 0) stream_ << '\n';
    stream_ << '}';
  }

  stream_.close();
}

void marley::TextOutputFile::write_flux_avg_tot_xsec(double avg_tot_xsec)
{
  // If we're at the beginning of an ASCII-format file, write
  // the flux-averaged total cross section before the events. This will
  // be written upon closing the file in the case of the JSON output
  // format. If we're not at the start of the file, don't bother.
  // We probably are appending to an existing file. We'll have to trust
  // the user not to mix events with different flux-averaged cross sections
  // in these file formats.
  /// @todo Consider other ways of handling this
  if (format_ == Format::ASCII) {
    bool at_start_of_file = stream_.tellp() == 0;
    if ( !at_start_of_file ) return;

    // Use the same trick as in marley::Event::print() to preserve
    // full numerical precision in the output text
    std::ostringstream temp;
    temp << std::scientific;
    temp.precision(std::numeric_limits<double>::max_digits10);

    temp << avg_tot_xsec;
    stream_ << temp.str() << '\n';
  }

  // Store the value for later (it is needed for the HEPEVT format)
  flux_avg_tot_xsec_ = avg_tot_xsec;
}
