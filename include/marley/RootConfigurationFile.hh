#pragma once
#include <string>
#include <vector>

#include "marley/ConfigurationFile.hh"

namespace marley {

  /// @brief ConfigurationFile that adds support for <a
  /// href="http://root.cern.ch">ROOT</a> input and output
  class RootConfigurationFile : public marley::ConfigurationFile {

    public:

      /// @brief Create a RootConfigurationFile object with all options set to
      /// their default values
      RootConfigurationFile();

      /// @brief Parse a file to create this RootConfigurationFile object
      /// @param file_name The name of the file to parse
      RootConfigurationFile(const std::string& file_name);

      /// @name Executable option accessors
      /// @brief Functions that access configuration file options that are
      /// only used if MARLEY is being run as an executable
      //@{

      /// @brief Check whether a pre-existing <a
      /// href="http://root.cern.ch/">ROOT</a> file will be silently
      /// overwritten (true) or whether the executable will prompt the user to
      /// confirm an overwrite (false)
      inline bool check_overwrite_root() const;

      /// @brief Check whether a <a href="http://root.cern.ch">ROOT</a> format
      /// file will be written (true) or not (false)
      inline bool check_write_root() const;

      /// @brief Get the name of the <a href="http://root.cern.ch">ROOT</a>
      /// format file that will receive Event output
      inline std::string get_root_filename() const;

      //@}

    protected:

      /// @name Executable options
      /// @brief Data members that store configuration file options that are
      /// only used if MARLEY is being run as an executable
      //@{

      /// @brief Whether to prompt the user (true) or not (false) before
      /// overwriting a previously-existing <a href="http://root.cern.ch">
      /// ROOT</a> format file
      bool check_root_overwrite_;

      /// @brief Name of the <a href="http://root.cern.ch">ROOT</a> format file
      /// that will receive Event output
      std::string root_filename_;

      /// @brief Whether to create a <a href="http://root.cern.ch"> ROOT</a>
      /// format output file (true) or not (false)
      bool write_root_;

      //@}

      virtual bool process_extra_source_types(const std::string& type,
        int neutrino_pid) override;

      virtual bool process_extra_keywords() override;

    private:

      /// @brief Helper function that checks whether (E, PDF) pairs represent a
      /// valid probability density function
      /// @param Es vector of energy values (MeV)
      /// @param PDFs vector of probability density values (MeV<sup> -1</sup>)
      void check_pdf_pairs(const std::vector<double>& Es,
        const std::vector<double>& PDFs);

      /// @brief Helper function for loading objects from a <a
      /// href="http://root.cern.ch">ROOT</a> <a
      /// href="https://root.cern.ch/doc/master/classTFile.html">TFile</a>
      /// @param tfile_name Name of the ROOT file to read
      /// @param namecycle Namecycle of the object to load from the file
      /// @return A pointer to the object after it has been loaded from
      /// the file, or nullptr if the load failed
      template<typename T> T* get_root_object(const std::string& tfile_name,
        const std::string& namecycle);
  };

  // Inline function definitions
  inline std::string RootConfigurationFile::get_root_filename() const
    { return root_filename_; }

  inline bool RootConfigurationFile::check_write_root() const
    { return write_root_; }

  inline bool RootConfigurationFile::check_overwrite_root() const
    { return check_root_overwrite_; }

}
