#include <chrono>
#include <csignal>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "marley/marley_utils.hh"

#ifdef USE_ROOT
  #include "marley/RootJSONConfig.hh"
#else
  #include "marley/JSONConfig.hh"
#endif

#include "marley/Generator.hh"
#include "marley/Event.hh"
#include "marley/Logger.hh"

#ifdef USE_ROOT
  #include "TFile.h"
  #include "TInterpreter.h"
  #include "TROOT.h"
  #include "TTree.h"
#endif

// Flag that will alert the main event generation loop that
// the user has interrupted the program's execution (probably
// by pressing ctrl+c). This is a global variable (out of
// necessity, since it must be modified by our signal handler),
// but its scope is limited to this source file by the use
// of the static keyword in its declaration.
volatile static std::sig_atomic_t interrupted = false;

// Function that will be used to handle a SIGINT signal in our event generation
// loop. The method used here is taken from the second answer at
// http://tinyurl.com/jm6zmwg
// Note that, since we don't use the signal code, we leave the integer argument
// of the function nameless. This prevents "unused parameter" warnings from the
// compiler.
void signal_handler(int)
{
  interrupted = true;
}

namespace {

  // Show a number using one decimal digit without scientific notation.
  // Used to print certain numbers in this way without affecting the settings
  // currently in use for std::cout.
  std::string format_number(double number) {
    static std::stringstream temp_stream;
    static bool configured = false;
    if (!configured) {
      temp_stream << std::fixed << std::setprecision(1);
    }

    // Clear the std::stringstream of any of its previous contents
    // This can be done more elegantly by move-assigning a default-constructed
    // std::stringstream object like this:
    // temp_stream = std::stringstream();
    // However, g++ versions < 5 don't implement the move assignment
    // operator for std::stringstream
    // (see discussion at http://stackoverflow.com/q/27152263/4081973).
    // The approach below is a reasonable workaround.
    temp_stream.str("");
    temp_stream.clear();

    temp_stream << number;

    return temp_stream.str();
  }

  // Helper class used to parse the output file specifications from the
  // configuration file and deliver output to the files.
  class OutputFile {

    public:

      OutputFile(const std::string& name, const std::string& format,
        const std::string& mode, bool force = false) : name_(name),
        force_(force)
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

      virtual ~OutputFile() = default;

      const std::string& name() const { return name_; }

      // Load a marley::Generator object whose configuration and state was
      // saved to the metadata in a ROOT or JSON format file.
      // Returns true if the generator is successfully restored using the file's
      // saved metadata, or false otherwise. The parameter num_previous_events
      // is loaded with the number of events saved in the file.
      virtual bool resume(std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events) = 0;

      // Closes the output file, performing any necessary cleanup. Also saves
      // information about the generator configuration and state to those
      // output file formats that support it.
      virtual void close(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) = 0;

      virtual int_fast64_t bytes_written() = 0;

      // Write a new marley::Event to this output file
      virtual void write_event(const marley::Event* event) = 0;

      bool mode_is_resume() const { return mode_ == Mode::RESUME; }

    protected:

      bool check_if_file_exists(const std::string& filename) {
        static std::ifstream test_stream;
        test_stream.open(filename);
        bool file_exists = test_stream.good();
        test_stream.close();
        return file_exists;
      }

      // Opens the output file for writing, making whatever preparations are
      // necessary before the first call to write_event()
      virtual void open() = 0;

      // For the file formats that support it, writes metadata containing the
      // generator state and configuration to the output file.
      virtual void write_generator_state(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) = 0;

      // Add more formats to the enum class as needed. This should include
      // every event format that MARLEY knows how to write. The "ASCII" format
      // is MARLEY's native format for textual input and output of
      // marley::Event objects (via the << and >> operators on std::ostream and
      // std::istream objects).
      enum class Format { ROOT, HEPEVT, JSON, ASCII };

      // Modes to use when writing output files.
      // OVERWRITE = removes any previous contents of the file, then writes new
      //   events in the requested format. This mode is allowed for all output
      //   formats.
      // APPEND = adds new events to the end of the file without altering any
      //   previous contents. This mode is only allowed for formats that don't
      //   save metadata (random number generator state, etc.) to the file,
      //   i.e., the HEPEVT and ASCII formats.
      // RESUME = uses saved metadata in a file to restore a previous MARLEY
      //   configuration (including the random number generator state), then
      //   appends new events after those currently saved in the file. This
      //   mode is only allowed for output formats that include such metadata,
      //   i.e., the ROOT and JSON formats.
      enum class Mode { OVERWRITE, APPEND, RESUME };

      std::string name_; // Name of the file to receive output
      Format format_;
      Mode mode_;

      // Whether to force the current output mode without prompting the user.
      // Currently only used when setting the OVERWRITE output mode.
      bool force_;
  };

  #ifdef USE_ROOT
    class RootOutputFile : public OutputFile {
      public:

        RootOutputFile(const std::string& name, const std::string& format,
          const std::string& mode, bool force = false)
          : OutputFile(name, format, mode, force)
        {
          load_marley_headers();
          open();
        }

        virtual ~RootOutputFile() = default;

        int_fast64_t bytes_written() override {
          return file_->GetBytesWritten();
        }

      private:

        virtual void open() override {
          bool file_exists = check_if_file_exists(name_);

          std::string tfile_open_mode("recreate");
          if (mode_ == Mode::OVERWRITE && file_exists && !force_) {
            bool overwrite = marley_utils::prompt_yes_no(
              "Overwrite ROOT file " + name_);
            if (!overwrite) {
              MARLEY_LOG_INFO() << "Cancelling overwrite of ROOT file \""
                << name_ << '\"';
              tfile_open_mode = std::string("update");
              mode_ = Mode::RESUME;
            }
          }
          else if (mode_ == Mode::RESUME) {
            if (!file_exists) throw marley::Error("Cannot resume run. Could"
              " not open the ROOT file \"" + name_ + '\"');
            else tfile_open_mode = std::string("update");
          }
          else if (mode_ == Mode::APPEND) throw marley::Error("Cannot use the"
            " \"append\" mode with the ROOT output format");

          file_ = std::make_unique<TFile>(name_.c_str(),
            tfile_open_mode.c_str());

          // Check if there was a problem opening the file (e.g., pre-existing
          // file that has the wrong format, etc.). If so, complain and quit.
          if (file_->IsZombie()) throw marley::Error("Invalid format or other"
            " error encountered while opening the ROOT file \"" + name_ + '\"');

          if (mode_ == Mode::OVERWRITE) {
            // Create a ROOT tree to store the events
            tree_ = new TTree("MARLEY_event_tree",
              "Neutrino events generated by MARLEY");

            // We create a branch to store the events here, but set the branch
            // address to nullptr. This will be fixed later when write_event()
            // is called.
            tree_->Branch("event", "marley::Event", nullptr);
          }

          else if (mode_ == Mode::RESUME) {
            file_->GetObject("MARLEY_event_tree", tree_);
            if (!tree_) throw marley::Error("Cannot resume run. Could not find"
              " a valid MARLEY event tree in the ROOT file \"" + name_ + '\"');

            MARLEY_LOG_INFO() << "Continuing previous run from ROOT file \""
              << name_ << "\"\nwhich contains " << tree_->GetEntries()
              << " events.";

            // Get previous RNG seed from the ROOT file
            std::string* seed = nullptr;
            file_->GetObject("MARLEY_seed", seed);

            if (!seed) MARLEY_LOG_WARNING() << "Unable to read"
              << " random number generator seed saved to the ROOT file from"
              << " the previous run.";
            else MARLEY_LOG_INFO() << "The previous run was initialized using"
              << " the random number generator seed " << *seed;
          }

          else throw marley::Error("Unrecognized file mode encountered in"
            " RootOutputFile::open()");
        }

        void load_marley_headers() {
          // Current (24 July 2016) versions of ROOT 6 require runtime loading
          // of headers for custom classes in order to use dictionaries
          // correctly. If we're running ROOT 6+, do the loading here, and give
          // the user guidance if there are any problems.
          // TODO: see if you can use the -inlineInputHeader to include the
          // headers in the libMARLEY_ROOT library and then load them at
          // runtime from there. This isn't very well documented, so your early
          // efforts at doing this didn't work out.
          static bool already_loaded = false;
          if (!already_loaded) {
            already_loaded = true;
            if (gROOT->GetVersionInt() >= 60000) {
              MARLEY_LOG_INFO() << "ROOT 6 or greater detected. Loading class"
                << " information\nfrom headers \"marley/Particle.hh\""
                << " and \"marley/Event.hh\"\n";
              TInterpreter::EErrorCode* ec = new TInterpreter::EErrorCode();
              gInterpreter->ProcessLine("#include \"marley/Particle.hh\"", ec);
              if (*ec != 0) throw marley::Error("Error loading MARLEY header"
                " Particle.hh. For MARLEY headers stored in /path/to/include/"
                "marley/, please add /path/to/include to your"
                " ROOT_INCLUDE_PATH environment variable and try again.");
              gInterpreter->ProcessLine("#include \"marley/Event.hh\"");
              if (*ec != 0) throw marley::Error("Error loading MARLEY header"
                " Event.hh. For MARLEY headers stored in"
                " /path/to/include/marley/, please add /path/to/include to"
                " your ROOT_INCLUDE_PATH environment variable and try again.");
            }
          }
        }

        // Write the internal state string of the random number generator to
        // disk, as well as the generator's JSON configuration. This will allow
        // MARLEY to resume event generation from where it left off with no
        // loss of consistency. This trick is based on
        // http://tinyurl.com/hb7rqsj
        void write_generator_state(const marley::JSON& json_config,
          const marley::Generator& gen, const long /*num_events*/) override
        {
          // Use std::string objects to store MARLEY configuration and
          // random number generator state information to the ROOT file.
          // You can retrieve these strings using the TFile::GetObject()
          // function (see elsewhere in this file for an example).
          file_->cd();
          std::string config( json_config.dump_string() );
          std::string state( gen.get_state_string() );
          std::string seed( std::to_string(gen.get_seed()) );

          file_->WriteObject(&config, "MARLEY_config", "WriteDelete");
          file_->WriteObject(&state, "MARLEY_state", "WriteDelete");
          file_->WriteObject(&seed, "MARLEY_seed", "WriteDelete");
        }

        // Clean up once we're done with the ROOT file. Save some information
        // about the current generator configuration as we clean things up.
        virtual void close(const marley::JSON& json_config,
          const marley::Generator& gen, const long dummy) override
        {

          // Write the event tree to the ROOT file, replacing the previous
          // version if one exists. Avoid data loss by not deleting the
          // previous version until the new version is completely written to
          // disk.
          file_->cd();
          tree_->Write(tree_->GetName(), TTree::kWriteDelete);

          // Save the current state of the generator to the ROOT file in case
          // we want to resume a run later
          write_generator_state(json_config, gen, dummy);

          file_->Close();
        }

        virtual bool resume(std::unique_ptr<marley::Generator>& gen,
          long& num_previous_events) override
        {
          if (mode_ != Mode::RESUME) {
            throw marley::Error("Cannot call Output"
              "File::resume() for an output mode other than"
              " \"resume\"");
            return false;
          }

          std::string* conf = nullptr;
          file_->GetObject("MARLEY_config", conf);
          if (!conf) {
            throw marley::Error("Cannot resume run. Failed to load"
              " JSON configuration from the ROOT file \"" + name_ + '\"');
            return false;
          }

          marley::JSON config = marley::JSON::load(*conf);

          std::string* state = nullptr;
          file_->GetObject("MARLEY_state", state);
          if (!state) {
            throw marley::Error("Cannot resume run. Failed to load"
              " random number generator state string from the ROOT file \""
              + name_ + '\"');
            return false;
          }

          marley::RootJSONConfig jc(config);
          gen = std::make_unique<marley::Generator>(jc.create_generator());
          gen->seed_using_state_string(*state);

          if (!tree_) throw marley::Error("Cannot resume run. Error"
            " accessing the MARLEY event tree stored in the ROOT file \""
            + name_ + '\"');

          num_previous_events = tree_->GetEntries();

          return true;
        }

        virtual void write_event(const marley::Event* event) override {
          if (!event) throw marley::Error("Null pointer passed to"
            " RootOutputFile::write_event()");
          tree_->SetBranchAddress("event", &event);
          tree_->Fill();
        }

        // TFile object that will be used to access the ROOT file
        std::unique_ptr<TFile> file_ = nullptr;

        // This is a bare pointer, but ROOT will associate it with file_, so we
        // don't want to delete it ourselves or let a smart pointer do it.
        TTree* tree_ = nullptr;
    };
  #endif

  class TextOutputFile : public OutputFile {
    private:

      // Formats the beginning of a JSON output file properly based on the
      // current indent level. Writes the result to a std::ostream.
      void start_json_output(bool start_array) {
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

      // Stream used to read and write from the output file as needed
      std::fstream stream_;

      // Flag used to see if we need a comma in front of the current
      // JSON event or not. Unused by the other formats.
      bool needs_comma_ = false;

      // Indent level used when writing JSON files. Unused by other
      // formats.
      int indent_ = -1; // -1 gives the most compact JSON file possible

      // Storage for the number of bytes written to disk
      int_fast64_t byte_count_ = 0;

    public:

      TextOutputFile(const std::string& name, const std::string& format,
        const std::string& mode, bool force = false, int indent = -1)
        : OutputFile(name, format, mode, force), indent_(indent)
      {
        open();
      }

      virtual ~TextOutputFile() = default;

      void set_indent(int indent) { indent_ = indent; }

      virtual void open() override {
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
        //MARLEY_LOG_INFO() << "Events for this run will be ";
        //if (hepevt_file_exists && open_mode_flag == std::ios:ate)
        //  std::cout << "appended ";
        //else std::cout << "written ";
        //std::cout << "to the HEPEvt format file " << hepevt_file_name << '\n';
      }

      // TODO: consider a better way of doing this
      // Returns true if the generator was successfully restored from the
      // saved metadata, or false otherwise.
      virtual bool resume(std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events) override
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

        #ifdef USE_ROOT
          marley::RootJSONConfig jc(config);
        #else
          marley::JSONConfig jc(config);
        #endif
        gen = std::make_unique<marley::Generator>(jc.create_generator());
        gen->seed_using_state_string(state_string);

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

      int_fast64_t bytes_written() override {
        // If the stream is open, then update the byte count. Otherwise, just
        // use the saved value.
        if (stream_.is_open()) {
          stream_.flush();
          byte_count_ = static_cast<int_fast64_t>(stream_.tellp());
        }
        return byte_count_;
      }

      // Write a new marley::Event to this output file
      virtual void write_event(const marley::Event* event) override
      {
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
            // TODO: fix event number here
            event->write_hepevt(0, stream_);
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

      void write_generator_state(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) override
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

        if (indent_ < 0) stream_ << temp.dump_string();
        else temp.print(stream_, indent_, true, indent_);
      }

      virtual void close(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) override
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
  };

  void print_help(std::string executable_name) {
    constexpr char help_message1[] = "Usage: ";
    constexpr char help_message2[] = " [OPTION...] FILE\n"
      "\n"
      "  -h, --help     Print this help message\n"
      "  -v, --version  Print version and exit\n";

    std::cout << help_message1 + executable_name + help_message2;
    std::cout << "\nMARLEY home page: <http://www.marleygen.org>\n";
    std::cout << "E-mail bug reports to: <support@marleygen.org>\n";
    exit(0);
  }

  void print_version() {
    std::cout << "MARLEY (Model of Argon Reaction Low Energy Yields) "
      << marley_utils::MARLEY_VERSION << '\n';

    std::cout << "Copyright (C) 2016-2019 Steven Gardiner\n";
    std::cout << "License: GNU GPL version 3 "
      << "<http://opensource.org/licenses/GPL-3.0>\n";
    std::cout << "This is free software: you are free to change and"
      << " redistribute it.\n";
    exit(0);
  }

  // Use this instead of std::put_time to allow this executable to be built
  // with g++ 4.9 (issue fixed in 5.0). See discussion here:
  // http://stackoverflow.com/a/14137287
  std::string put_time(std::tm* time, const char* format)
  {
    constexpr size_t TIME_STR_SIZE = 100;
    std::string time_str(TIME_STR_SIZE, ' ');
    // Pass a pointer to the first character in the time string as the char*
    // argument of std::strftime(). This is well-defined behavior because the
    // storage of std::string is guaranteed by the C++11 standard to be
    // contiguous, and we have pre-allocated storage inside the string of the
    // necessary size using the constructor above.
    std::strftime(&time_str.front(), TIME_STR_SIZE, format, time);
    marley_utils::trim_right_inplace(time_str);
    return time_str;
  }

  // Formats the lines of the status display shown at the bottom of the
  // screen when this executable is running
  std::string makeStatusLines(long ev_count, long num_events, long num_old_events,
    std::chrono::system_clock::time_point start_time_point,
    const std::vector<std::unique_ptr<OutputFile> >& output_files)
  {
    std::ostringstream temp_oss;
    // Print a status message showing the current number of events
    temp_oss << "\nEvent Count = " << ev_count << "/" << num_events
      << " (" << format_number(ev_count*100
        / static_cast<double>(num_events))
      << "% complete)\033[K\n";

    // Print timing information
    std::chrono::system_clock::time_point current_time_point
      = std::chrono::system_clock::now();
    temp_oss << "Elapsed time: "
      << marley_utils::elapsed_time_string(start_time_point,
      current_time_point) << " (Estimated total run time: ";

    marley_utils::seconds<float> estimated_total_time =
      (current_time_point - start_time_point)
      * (static_cast<float>(num_events - num_old_events)
      / (ev_count - num_old_events));

    temp_oss << marley_utils::duration_to_string
      <marley_utils::seconds<float>>(estimated_total_time)
      << ")\033[K\n";

    for (const auto& file : output_files) {
      temp_oss << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string(file->bytes_written(), 2)
        << "\033[K\n";
    }

    std::time_t estimated_end_time = std::chrono::system_clock::to_time_t(
      start_time_point + std::chrono::duration_cast
      <std::chrono::system_clock::duration>(estimated_total_time));

    temp_oss << "MARLEY is estimated to terminate on "
      << put_time(std::localtime(&estimated_end_time), "%c %Z") << '\n';

    // Move up an extra line in the terminal for each output file because
    // we're displaying the amount of data written to disk for each of them
    // on a separate line.
    for (size_t i = 0; i < output_files.size(); ++i) temp_oss << "\033[F";

    // Move up four lines
    temp_oss << "\033[F\033[F\033[F\033[F";

    return temp_oss.str();
  }

  // Inserts an updated copy of the status lines after each
  // newline is written to std::cout (e.g., by the marley::Logger
  // class). Based on https://stackoverflow.com/a/22043916
  class StatusInserter : public std::streambuf {
    public:
      StatusInserter(std::streambuf* dest, const long& ev_count,
        const long& num_events, const long& num_old_events,
        const std::chrono::system_clock::time_point& start_time_point,
        const std::vector<std::unique_ptr<OutputFile> >& output_files)
        : std::streambuf(), myDest_( dest ), myIsAtStartOfLine_(true),
        do_status_(true), ev_count_( ev_count ), num_events_( num_events ),
        num_old_events_( num_old_events ), start_time_point_( start_time_point ),
        output_files_( output_files ) {}

      inline void set_do_status(bool do_it) { do_status_ = do_it; }

    protected:
      std::streambuf* myDest_;
      bool myIsAtStartOfLine_;
      bool do_status_;
      const long& ev_count_;
      const long& num_events_;
      const long& num_old_events_;
      const std::chrono::system_clock::time_point& start_time_point_;
      const std::vector<std::unique_ptr<OutputFile> >& output_files_;

      int overflow( int ch ) override {
        int retval = 0;
        if ( ch != traits_type::eof() ) {
          if ( do_status_ && myIsAtStartOfLine_ ) {
            std::string status = makeStatusLines(ev_count_,
              num_events_, num_old_events_, start_time_point_, output_files_);
            myDest_->sputn( status.data(), status.size() );
          }
          myIsAtStartOfLine_ = ch == '\n';
          if ( myIsAtStartOfLine_ ) {
            std::string erase( "\033[K" );
            myDest_->sputn( erase.data(), erase.size() );
          }
          retval = myDest_->sputc( ch );
        }
        return retval;
      }
  };

}

int main(int argc, char* argv[]) {

  std::streambuf* cout_default_buf = std::cout.rdbuf();
  std::streambuf* cerr_default_buf = std::cerr.rdbuf();

  try {

    // Initialize the logger. Use a default log level of INFO. This will
    // quickly be overwritten while parsing the JSON configuration file.
    marley::Logger::Instance().add_stream(std::cout,
      marley::Logger::LogLevel::INFO);

    // TODO: if you need command line parsing beyond the trivial stuff used
    // here, consider using a header-only option parser library like cxxopts
    // (https://github.com/jarro2783/cxxopts)
    std::string config_file_name;

    // If the user has not supplied any command-line
    // arguments, display the standard help message
    // and exit
    if (argc <= 1) {
      print_help(argv[0]);
    }
    // The first command-line argument does not begin with a
    // hyphen, so assume that it is the configuration file name.
    else if (std::string(argv[1]).substr(0,1) != "-") {
      config_file_name = argv[1];
    }
    // The first command-line argument begins with a hyphen,
    // so treat it as an option and react accordingly.
    // All of the current options cannot be combined, so
    // just parse the first one and ignore the others.
    else {
      std::string option = argv[1];
      if (option == "-h" || option == "--help") {
        print_help(argv[0]);
      }
      else if (option == "-v" || option == "--version") {
        print_version();
      }
      else if (option == "--marley") {
        std::cout << marley_utils::marley_pic;
        exit(0);
      }
      else {
       std::cout << argv[0] << ": unrecognized "
         << "command line option '" <<  option << "'\n";
       print_help(argv[0]);
      }
    }

    // Parse the objects from the JSON-based configuration file
    marley::JSON json = marley::JSON::load_file(config_file_name);

    // Process them to get the generator configuration. Enable ROOT support if
    // it is available.
    #ifdef USE_ROOT
      marley::RootJSONConfig jc(json);
    #else
      marley::JSONConfig jc(json);
    #endif

    // Get the time that the program was started
    std::chrono::system_clock::time_point start_time_point
      = std::chrono::system_clock::now();

    std::time_t start_time = std::chrono::system_clock::to_time_t(
      start_time_point);

    std::cout << "\nMARLEY started on "
      << put_time(std::localtime(&start_time), "%c %Z") << '\n';

    long num_old_events = 0;

    marley::JSON ex_set = json.get_object("executable_settings");

    // Desired number of events to be generated in this run. This
    // will be read back from the ROOT file and overwritten
    // if we're doing a continuation run.
    long num_events = ex_set.get_long("events", 1e3);

    std::vector<std::unique_ptr<OutputFile> > output_files;

    if (ex_set.has_key("output")) {
      marley::JSON output_set = ex_set.at("output");
      if (!output_set.is_array()) throw marley::Error("The"
        " \"output\" key in the executable settings must have a value that"
        " is a JSON array.");
      else for (const auto& el : output_set.array_range()) {
        if (!el.has_key("file")) throw marley::Error("Missing file name"
          " for an output file specification in the configuration file.");
        else if (!el.has_key("format")) throw marley::Error("Missing format"
          " for an output file specification in the configuration file.");

        std::string filename = el.at("file").to_string();
        std::string format = el.at("format").to_string();

        std::string mode("overwrite"); // default mode is "overwrite"
        if (el.has_key("mode")) mode = el.at("mode").to_string();

        bool force = false; // default behavior is to prompt before overwriting
        if (el.has_key("force")) force = el.at("force").to_bool();

        // Set the indent level (if needed) for JSON format files
        int indent = -1;
        if (format == "json" && el.has_key("indent")) {
          bool ok = false;
          const marley::JSON& idt = el.at("indent");
          indent = static_cast<int>(idt.to_long(ok));
          if (indent < 0) ok = false;
          if (!ok) {
            throw marley::Error("Invalid indent value \""
              + idt.dump_string() + "\" for JSON file \""
              + filename + '\"');
            indent = -1;
          }
        }

        #ifdef USE_ROOT
          if (format == "root") output_files.push_back(
            std::make_unique<RootOutputFile>(filename, format, mode, force));
          else output_files.push_back(std::make_unique<TextOutputFile>(
            filename, format, mode, force, indent));
        #else
          output_files.push_back(std::make_unique<TextOutputFile>(filename,
            format, mode, force, indent));
        #endif
      }
    }
    else {
      // If the user didn't specify anything for the output key, then
      // by default write to a single JSON-format file.
      output_files.push_back(std::make_unique<TextOutputFile>("events.json",
        "json", "overwrite", false));
    }

    // This std::unique_ptr to a Generator object will be initialized below
    // with the proper configuration.
    std::unique_ptr<marley::Generator> gen;

    // Check if we need to resume from one of the files. If so, go ahead and
    // do it. If more than one file is requested for the resume, complain.
    bool need_to_resume = false;
    for (auto& file : output_files) {
      if (file->mode_is_resume()) {
        if (need_to_resume) throw marley::Error("Only one file may be used"
          " to resume a previous run.");
        else {
          need_to_resume = true;
          bool resume_ok = file->resume(gen, num_old_events);
          if (!resume_ok) throw marley::Error("Failed to resume previous run"
            " from the file \"" + file->name() + '\"');
        }
      }
    }

    // If we didn't resume a run from any of the files, use the current
    // configuration
    if (!need_to_resume) gen = std::make_unique<marley::Generator>(
      jc.create_generator());

    // Use the signal handler defined above to deal with
    // SIGINT signals (e.g., ctrl+c interruptions initiated
    // by the user). This will allow us to terminate the
    // loop gracefully, leaving a valid event tree file, etc.
    std::signal(SIGINT, signal_handler);

    // Generate all of the requested events. End the loop early
    // if the user interrupts execution (e.g., via ctrl+C)
    std::cout << '\n';
    auto event = std::make_unique<marley::Event>();
    long ev_count = 1 + num_old_events;

    // Make std::cout use our "status inserter" std::streambuf
    // object so that the status lines get automatically updated
    // with every newline
    StatusInserter my_status_inserter(cout_default_buf, ev_count,
      num_events, num_old_events, start_time_point, output_files);
    std::cout.rdbuf( &my_status_inserter );
    std::cerr.rdbuf( &my_status_inserter );

    for (; ev_count <= num_events && !interrupted; ++ev_count) {

      // Create an event using the generator object
      *event = gen->create_event();

      for (const auto& file : output_files) {
        file->write_event(event.get());
      }

      // Print status messages about simulation progress after every 100
      // events have been generated
      if ((ev_count - num_old_events) % 100 == 1 || ev_count == num_events) {

        // Temporarily disable the auto-printing of the status lines by
        // the StatusInserter class. This will avoid duplication of the
        // information.
        my_status_inserter.set_do_status( false );

        // Print a status message showing the current number of events
        std::cout << makeStatusLines(ev_count, num_events, num_old_events,
          start_time_point, output_files);

        // Re-enable the auto-printing of the status lines now that we've
        // printed them manually
        my_status_inserter.set_do_status( true );
      }
    }

    // Restore the default std::streambuf to std::cout
    std::cout.rdbuf( cout_default_buf );
    std::cerr.rdbuf( cerr_default_buf );

    std::cout << "\033[E";

    for (const auto& file : output_files) {
      file->close(json, *gen, ev_count - 1);
      std::cout << "Data written to " << file->name() << ' '
        << marley_utils::num_bytes_to_string(file->bytes_written())
        << "\033[K\n";
    }

    // Display the time that the program terminated
    std::chrono::system_clock::time_point end_time_point
      = std::chrono::system_clock::now();
    std::time_t end_time = std::chrono::system_clock::to_time_t(
      end_time_point);

    if (!interrupted) std::cout << "MARLEY terminated normally on ";
    else std::cout << "MARLEY was interrupted by the user on ";
    std::cout << put_time(std::localtime(&end_time), "%c %Z")
      << "\033[K\033[E\033[K\n\033[K";

    return 0;
  }

  catch (const marley::Error& error) {
    // Restore the old std::streambuf for std::cout
    std::cout.rdbuf( cout_default_buf );
    std::cerr.rdbuf( cerr_default_buf );

    // Flush the Logger, then rethrow the error for the system to handle.
    // This will probably result in std::terminate() being called.
    auto& log = marley::Logger::Instance();
    log.newline();
    log.flush();
    throw error;
  }

  return 1;
}
