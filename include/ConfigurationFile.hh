#pragma once
#include <array>
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
#include "StructureDatabase.hh"

namespace marley {

  /// @brief Parser for a MARLEY Generator configuration file
  /// @details For details about the format of the configuration file, please
  /// see chapter 1 of the MARLEY user manual.
  class ConfigurationFile {

    public:

      /// @brief Create a ConfigurationFile object with all options set to
      /// their default values
      ConfigurationFile();

      /// @brief Parse a file to create this ConfigurationFile object
      /// @param file_name The name of the file to parse
      ConfigurationFile(const std::string& file_name);

      /// @brief Get the random number seed that will be used to initialize a
      /// Generator
      inline uint_fast64_t get_seed() const;

      /// @brief Set the random number seed that will be used to initialize a
      /// Generator
      inline void set_seed(uint_fast64_t s);

      /// @brief Get the name(s) of the reaction data files that will be used
      /// to create the Reaction object(s)
      inline const std::unordered_set<std::string>&
        get_reaction_filenames() const;

      /// @brief Add a new reaction data file name
      /// @param rfile Name of the file to add
      inline void add_reaction_filename(const std::string& rfile);

      /// @brief Remove an existing reaction data file name
      /// @param rfile Name of the file to remove
      inline void remove_reaction_filename(const std::string& rfile);

      /// @brief Clear the list of reaction data files
      inline void clear_reaction_filenames();

      /// @brief Gets a three-vector pointing in the direction of the incident
      /// neutrinos
      inline const std::array<double, 3>& get_neutrino_direction();

      /// @brief Sets the incident neutrino direction
      /// @param dir_vec Three-vector that points in the direction of the
      /// incident neutrinos. It does not need to be normalized, but it must be
      /// nonzero, or a marley::Error will be thrown by this function
      inline void set_neutrino_direction(const std::array<double, 3>& dir_vec);

      /// @brief Gets a three-vector pointing in the default incident neutrino
      /// direction
      static inline std::array<double, 3> get_default_neutrino_direction();

      /// @brief Get a unique_ptr reference to the StructureDatabase
      /// that will be used to help create a Generator object
      inline std::unique_ptr<marley::StructureDatabase>& get_structure_db();

      /// @brief Print a summary of the configuration file settings to
      /// a std::ostream
      void print_summary(std::ostream& os = std::cout);

      /// @brief Get a unique_ptr reference to the NeutrinoSource
      /// that will be used to help create a Generator object
      inline std::unique_ptr<marley::NeutrinoSource>& get_source();

      /// @brief Transfers ownership of a NeutrinoSource to this
      /// ConfigurationFile, deleting the previous source if one exists
      void set_source(std::unique_ptr<marley::NeutrinoSource>&& new_source);

      /// @name Executable option accessors
      /// @brief Functions that access configuration file options that are
      /// only used if MARLEY is being run as an executable
      //@{

      /// @brief Check whether a pre-existing
      /// <a href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEvt</a> file will be silently overwritten (true) or whether the
      /// executable will prompt the user to confirm an overwrite (false)
      inline bool check_overwrite_hepevt() const;

      /// @brief Check whether a
      /// <a href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEvt</a> format file will be written (true) or not (false)
      inline bool check_write_hepevt() const;

      /// @brief Get the name of the
      /// <a href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEvt</a> format file that will receive Event output
      inline std::string get_hepevt_filename() const;

      /// @brief Get the number of events that should be generated during a
      /// single run of the executable
      inline size_t get_num_events();

      //@}

    protected:

      /// @brief Helper function for the constructor that handles file parsing.
      void parse();

      /// @brief Helper function that allows extra NeutrinoSource types
      /// to be created by derived classes (like RootConfigurationFile)
      /// without needing to override parse().
      /// @details This method is overridden by the derived class
      /// RootConfigurationFile to allow NeutrinoSource objects to be created
      /// that require ROOT libraries to load.
      /// @return If a NeutrinoSource was successfully created using one of the
      /// extra types, then this function returns true. Otherwise, it returns
      /// false.
      /// @param type String from the configuration file indicating
      /// what NeutrinoSource type is being requested
      /// @param int neutrino_pdg PDG code for the neutrino type produced by
      /// the new source
      inline virtual bool process_extra_source_types(const std::string& type,
        int neutrino_pdg);

      /// @brief Helper function that allows extra configuration file keywords
      /// to be processed by derived classes (like RootConfigurationFile)
      /// without needing to override parse().
      /// @details This method is overridden by the derived class
      /// RootConfigurationFile to handle keywords specifying details
      /// of the ROOT output when MARLEY is run as an executable.
      /// @return If an extra keyword was successfully processed,
      /// then this function returns true. Otherwise, it returns false.
      virtual bool process_extra_keywords();

      /// @brief Helper function to add DecayScheme objects to the structure
      /// database
      /// @param file_name Name of the file containing the discrete level data
      /// @param format Format of the file
      /// @param nucleus_pdg_codes PDG codes for each nuclide whose discrete
      /// level data should be loaded from the file
      void add_decay_schemes(const std::string& file_name,
        const marley::DecayScheme::FileFormat format,
        const std::set<int>& nucleus_pdg_codes);

      /// @brief Helper function for parse() that gets the next word from a
      /// parsed line of the configuration file
      /// @param[out] word The next word from the line
      /// @param enable_exceptions Whether to throw exceptions (true) or not
      /// (false) if a parsing problem occurs
      /// @param make_lowercase Whether to convert the word to all-lowercase
      /// (true) or leave it as is (false) before returning it
      /// @return true if another word was found on the line, or false if the
      /// end of the line has been reached
      bool next_word_from_line(std::string& word,
        bool enable_exceptions = true, bool make_lowercase = true);

      /// @brief Name of the configuration file to parse
      std::string filename_;

      /// @brief Random number seed to use when creating a Generator
      uint_fast64_t seed_;

      /// @brief Reaction data files to parse when creating Reaction objects
      std::unordered_set<std::string> reaction_filenames_;

      /// @brief StructureDatabase that holds the discrete level data specified
      /// in the configuration file
      std::unique_ptr<marley::StructureDatabase> structure_db_;

      /// @brief Default incident neutrino direction to use for generating
      /// events
      static const std::array<double, 3> DEFAULT_INCIDENT_NEUTRINO_DIRECTION_;

      /// @brief Three-vector that points in the direction of the incident
      /// neutrinos
      std::array<double, 3> dir_vec_ = DEFAULT_INCIDENT_NEUTRINO_DIRECTION_;

      /// @brief A NeutrinoSource object that is created based on the
      /// specifications given in the configuration file.
      /// @details When a marley::Generator object is constructed using this
      /// ConfigurationFile, ownership of the NeutrinoSource will be
      /// transferred to the Generator. After this occurs, source_.get() will
      /// return nullptr.
      std::unique_ptr<marley::NeutrinoSource> source_;

      /// @name Executable options
      /// @brief Data members that store configuration file options that are
      /// only used if MARLEY is being run as an executable
      //@{

      /// @brief Whether to prompt the user (true) or not (false) before
      /// overwriting a previously-existing <a
      /// href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEvt</a> format file
      bool check_hepevt_overwrite_;

      /// @brief The default number of events to generate
      static constexpr size_t DEFAULT_NUM_EVENTS_ = 1000;

      /// @brief Name of the <a
      /// href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEvt</a> format file that will receive Event output
      std::string hepevt_filename_;

      /// @brief The number of events to generate
      size_t num_events_;

      /// @brief Whether to create a <a
      /// href="http://home.fnal.gov/~mrenna/lutp0613man2/node49.html">
      /// HEPEvt</a> format output file (true) or not (false)
      bool write_hepevt_;

      //@}

      /// @name File parsing helper variables
      //@{

      /// @brief std::ifstream used to read the configuration file
      std::ifstream file_in_;

      std::istringstream iss_; ///< Stream used to parse each line

      /// @brief String that stores the current keyword read in from the
      /// configuration file
      std::string keyword_;

      /// @brief String that stores the current line of the configuration file
      /// during parsing
      std::string line_;

      size_t line_num_;  ///< The current line number being processed
      //@}
  };

  // Inline function definitions
  inline uint_fast64_t ConfigurationFile::get_seed() const { return seed_; }

  inline void ConfigurationFile::set_seed(uint_fast64_t s) { seed_ = s; }

  inline const std::unordered_set<std::string>&
    ConfigurationFile::get_reaction_filenames() const
    { return reaction_filenames_; }

  inline void ConfigurationFile::add_reaction_filename(const std::string& rfile)
    { reaction_filenames_.insert(rfile); }

  inline void ConfigurationFile::remove_reaction_filename(
    const std::string& rfile) { reaction_filenames_.erase(rfile); }

  inline void ConfigurationFile::clear_reaction_filenames()
    { reaction_filenames_.clear(); }

  inline std::unique_ptr<marley::StructureDatabase>&
    ConfigurationFile::get_structure_db() { return structure_db_; }

  inline std::unique_ptr<marley::NeutrinoSource>&
    ConfigurationFile::get_source() { return source_; }

  inline size_t ConfigurationFile::get_num_events() { return num_events_; }

  inline std::string ConfigurationFile::get_hepevt_filename() const
    { return hepevt_filename_; }

  inline bool ConfigurationFile::check_write_hepevt() const
    { return write_hepevt_; }

  inline bool ConfigurationFile::check_overwrite_hepevt() const
    { return check_hepevt_overwrite_; }

  inline bool ConfigurationFile::process_extra_source_types(const std::string&,
    int) { return false; }

  inline const std::array<double, 3>&
    ConfigurationFile::get_neutrino_direction() { return dir_vec_; }

  inline std::array<double, 3>
    ConfigurationFile::get_default_neutrino_direction()
    { return DEFAULT_INCIDENT_NEUTRINO_DIRECTION_; }
}
