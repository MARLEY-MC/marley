#pragma once
#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

class TMarleyConfigFile {
  public:

    // Type used to organize nuclear structure file
    // load requests within this configuration file
    struct StructureRecord {
      std::string filename;
      std::string format;
      std::unordered_set<std::string> nucids;
    };

    TMarleyConfigFile();
    TMarleyConfigFile(const std::string file_name);
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

  private:

    std::string filename;
    uint_fast64_t seed;
    std::unordered_set<std::string> reaction_filenames;

    std::vector<StructureRecord> structure_records;

    #ifdef USE_ROOT
    std::string root_filename;
    bool writeroot;
    bool check_before_root_file_overwrite;
    #endif

    // Matches comment lines and empty lines
    static const std::regex rx_comment_or_empty;
    // Matches positive integers
    static const std::regex rx_num;
    // Matches trimmed ENSDF-style nucids
    static const std::regex rx_nucid;


    // Get the next word from a parsed line. If errors occur, complain.
    bool next_word_from_line(std::istringstream& iss, std::string& word,
      const std::string& keyword, const int line_number,
      bool enable_exceptions = true);
};
