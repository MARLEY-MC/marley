#pragma once
#include <memory>
#include <string>

#include "ConfigFile.hh"
#include "Generator.hh"

namespace marley {

  class Factory {
    public:
      inline Factory() : cf(nullptr) {};
      std::unique_ptr<marley::Generator> make_generator(
        const std::string& file_name);
    protected:
      std::unique_ptr<marley::ConfigFile> cf;
  };

}
