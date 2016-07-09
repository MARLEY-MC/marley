#pragma once
#include <string>
#include <vector>

#include "ConfigFile.hh"

namespace marley {

  class RootConfigFile : public marley::ConfigFile {
    public:
      RootConfigFile();
      RootConfigFile(const std::string& file_name);

      inline std::string get_root_filename() const { return root_filename; }
      inline bool check_write_root() const { return writeroot; }
      inline bool check_overwrite_root() const {
        return check_before_root_file_overwrite;
      }

    protected:

      std::string root_filename;
      bool writeroot;
      bool check_before_root_file_overwrite;

      virtual bool process_extra_source_types(const std::string& type,
        int neutrino_pid) override;

      virtual bool process_extra_keywords() override;

      // Helper function for loading objects from a ROOT TFile
      template<typename T> T* get_root_object(const std::string& tfile_name,
        const std::string& namecycle);

      // Helper function to check that (E, PDF) pairs represent a valid
      // probability density function
      void check_pdf_pairs(const std::vector<double>& Es,
        const std::vector<double>& PDFs);
  };

}
