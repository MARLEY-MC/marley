#pragma once

// standard library includes
#include <fstream>
#include <string>

namespace marley {

  class Generator;

  /// @brief Abstract base class for objects that deliver output to a file
  /// opened by the marley command-line executable
  class OutputFile {

    public:

      // Arguments to the constructor correspond to fields in
      // the JSON objects that appear under the "output" key
      // in the "executable_settings" section of a MARLEY job
      // configuration file
      OutputFile(const std::string& name, const std::string& format,
        const std::string& mode, bool force = false);

      virtual ~OutputFile() = default;

      const std::string& name() const { return name_; }

      /// @brief Load a marley::Generator object whose configuration and state
      /// were saved to the metadata in a ROOT or JSON format output file.
      /// @param[out] num_previous_events The number of events previously
      /// saved to the output file
      /// @return True if the generator is successfully restored using the file's
      /// saved metadata, or false otherwise.
      virtual bool resume(std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events) = 0;

      /// @brief Close the output file and perform any necessary cleanup
      /// @details This function also saves information about the generator
      /// configuration and state to those output file formats that support it.
      virtual void close(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) = 0;

      /// @brief The number of bytes that have been written during the
      /// current MARLEY session to this file
      virtual int_fast64_t bytes_written() = 0;

      /// @brief Write a new marley::Event to this output file
      /// @details If a nullptr is passed to this function, then
      /// a marley::Error will be thrown
      virtual void write_event(const marley::Event* event) = 0;

      bool mode_is_resume() const { return mode_ == Mode::RESUME; }

      /// @brief If needed (for the HEPEVT and ASCII formats), write the
      /// flux-averaged total cross section to the file.
      /// @details This function is a no-op unless we're
      /// at the beginning of the output stream.
      virtual void write_flux_avg_tot_xsec(double avg_tot_xsec) = 0;

      // Add more formats to the enum class as needed. This should include
      // every event format that MARLEY knows how to write. The "ASCII" format
      // is MARLEY's native format for textual input and output of
      // marley::Event objects (via the << and >> operators on std::ostream and
      // std::istream objects).
      enum class Format { ROOT, HEPEVT, JSON, ASCII };

    protected:

      /// @brief Helper function for resume() that instantiates a
      /// marley::Generator object given the previous JSON configuration
      /// object
      std::unique_ptr<marley::Generator> restore_generator(
        const marley::JSON& config);

      /// @brief Checks that a given file exists (and is readable)
      inline bool check_if_file_exists(const std::string& filename) {
        std::ifstream test_stream( filename );
        return test_stream.good();
      }

      /// @brief Open the output file for writing
      /// @details This function will also make  whatever preparations are
      /// necessary before the first call to write_event()
      virtual void open() = 0;

      /// @brief For the file formats that support it, write metadata containing the
      /// generator state and configuration to the output file.
      virtual void write_generator_state(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) = 0;

      // Modes to use when writing output to files that are not initially empty
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

      std::string name_; ///< Name of the file to receive output
      Format format_; ///< Format to use when writing events to the file
      Mode mode_;

      // Whether to force the current output mode without prompting the user.
      // Currently only used when setting the OVERWRITE output mode.
      bool force_;
  };

  class TextOutputFile : public OutputFile {
    private:

      // Formats the beginning of a JSON output file properly based on the
      // current indent level. Writes the result to a std::ostream.
      void start_json_output(bool start_array);

      // Stream used to read and write from the output file as needed
      std::fstream stream_;

      // Flag used to see if we need a comma in front of the current
      // JSON event or not. Unused by the other formats.
      bool needs_comma_ = false;

      /// @brief Indent level used when writing JSON files. Unused by other
      /// formats.
      int indent_ = -1; // -1 gives the most compact JSON file possible

      /// @brief Storage for the number of bytes written to disk
      int_fast64_t byte_count_ = 0;

      /// @brief Persistent storage for the flux-averaged total cross
      /// section value (needed for the HEPEVT output format, which includes
      /// it in every event)
      double flux_avg_tot_xsec_ = 0.;

    public:

      TextOutputFile(const std::string& name, const std::string& format,
        const std::string& mode, bool force = false, int indent = -1);

      virtual ~TextOutputFile() = default;

      inline void set_indent(int indent) { indent_ = indent; }

      virtual void open() override;

      // TODO: consider a better way of doing this
      // Returns true if the generator was successfully restored from the
      // saved metadata, or false otherwise.
      virtual bool resume(std::unique_ptr<marley::Generator>& gen,
        long& num_previous_events) override;

      int_fast64_t bytes_written() override;

      // Write a new marley::Event to this output file
      virtual void write_event(const marley::Event* event) override;

      void write_generator_state(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) override;

      virtual void close(const marley::JSON& json_config,
        const marley::Generator& gen, const long num_events) override;

      // This function is a no-op for the JSON format (we will write the
      // flux-averaged cross section to the output file when saving the
      // generator state)
      virtual void write_flux_avg_tot_xsec(double avg_tot_xsec) override;
  };

}
