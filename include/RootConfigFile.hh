#pragma once
#include <string>

#include "ConfigFile.hh"

namespace marley {

  class RootConfigFile : public marley::ConfigFile {
    public:
      inline RootConfigFile() : ConfigFile() {};
      inline RootConfigFile(const std::string file_name) : ConfigFile() {
        filename = file_name;
        parse();
      };

    protected:
      virtual bool process_extra_source_types(const std::string& type,
        int neutrino_pid, double weight) override;

      // Helper function for loading objects from a ROOT TFile
      template<typename T> T* get_root_object(const std::string& tfile_name,
        const std::string& namecycle);
  };

}
