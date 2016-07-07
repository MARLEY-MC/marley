#pragma once
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <unordered_set>

#include "DecayScheme.hh"
#include "Error.hh"
#include "NeutrinoSource.hh"

namespace marley {

  class ConfigFile {
    public:

      // Type used to organize nuclear structure file
      // load requests within this configuration file
      struct StructureRecord {
        std::string filename;
        marley::DecayScheme::FileFormat format;
        std::unordered_set<int> nucleus_pdg_codes;
      };

      ConfigFile();
      ConfigFile(const std::string file_name);
      inline uint_fast64_t get_seed() const { return seed; }
      inline void set_seed(uint_fast64_t s) { seed = s; };
      inline const std::unordered_set<std::string>&
        get_reaction_filenames() const
      {
        return reaction_filenames;
      }
      inline void add_reaction_filename(std::string rfile) {
        reaction_filenames.insert(rfile);
      }
      inline void remove_reaction_filename(std::string rfile) {
        reaction_filenames.erase(rfile);
      }
      inline void clear_reaction_filenames() {
        reaction_filenames.clear();
      }
      inline const std::vector<StructureRecord>& get_structure_records() const {
        return structure_records;
      }
      #ifdef USE_ROOT
      inline std::string get_root_filename() const { return root_filename; }
      inline bool check_write_root() const { return writeroot; }
      inline bool check_overwrite_root() const {
        return check_before_root_file_overwrite;
      }
      #endif
      void print_summary(std::ostream& os = std::cout);

      inline std::vector<std::unique_ptr<marley::NeutrinoSource> >& get_sources() {
        return sources;
      }

      inline size_t get_num_events() { return num_events; }

      inline std::string get_hepevt_filename() const { return hepevt_filename; }
      inline bool check_write_hepevt() const { return writehepevt; }
      inline bool check_overwrite_hepevt() const {
        return check_before_hepevt_file_overwrite;
      }

    protected:

      // This will be overridden by the derived class marley::RootConfigFile to
      // allow sources that require ROOT libraries to load. The function will
      // return true if a source was successfully created from one of the extra
      // types, or false otherwise. The input values (unused in the base class)
      // are the string from the file indicating what type of source is being
      // requested and the neutrino type particle ID
      inline virtual bool process_extra_source_types(const std::string&, int)
        { return false; }

      // Helper function for the constructor. We avoid polymorphism in the
      // constructor by delegating the parsing to this method.
      void parse();

      std::string filename;
      uint_fast64_t seed;
      std::unordered_set<std::string> reaction_filenames;

      std::vector<StructureRecord> structure_records;

      // The number of events to generate in a run
      size_t num_events;
      static constexpr size_t DEFAULT_NUM_EVENTS = 1000;

      // Neutrino sources will be created based on the specifications given in
      // this file. When a marley::Generator object is constructed using this
      // configuration file, ownership will be transferred to the generator.
      std::vector<std::unique_ptr<marley::NeutrinoSource> > sources;

      std::string root_filename;
      bool writeroot;
      bool check_before_root_file_overwrite;

      std::string hepevt_filename;
      bool writehepevt;
      bool check_before_hepevt_file_overwrite;

      // Convert a lowercase, trimmed string file format
      // description for the structure data to a
      // marley::DecayScheme::FileFormat
      marley::DecayScheme::FileFormat string_to_format(
        const std::string& string);

      // Convert a marley::DecayScheme::FileFormat value
      // to a std::string
      std::string format_to_string(const marley::DecayScheme::FileFormat ff);

      // Get the next word from a parsed line. If errors occur, complain.
      // The last argument determines whether the next word should be
      // converted to all lowercase or left as is.
      bool next_word_from_line(std::string& word,
        bool enable_exceptions = true, bool make_lowercase = true);

      // Helper variables used when parsing the file
      std::ifstream file_in;
      std::string line; // String to store the current line
                        // of the configuration file during parsing

      std::string keyword; // String to store the current keyword read in
                           // from the configuration file

      std::istringstream iss; // Stream used to parse each line
      int line_num;  // Keeps track of the current line of the file
  };

}
