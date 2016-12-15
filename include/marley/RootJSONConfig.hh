#pragma once

// MARLEY includes
#include "marley/JSONConfig.hh"

namespace marley {

  class RootJSONConfig : public JSONConfig {

    public:

      explicit RootJSONConfig(const marley::JSON& object)
        : JSONConfig(object) {}

      virtual bool process_extra_source_types(const std::string& type,
        const marley::JSON& source_spec, int pdg_code,
        std::unique_ptr<marley::NeutrinoSource>& source) const override;
  };

}
