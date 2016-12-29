#pragma once

// standard library includes
#include <string>

// MARLEY includes
#include "marley/Generator.hh"
#include "marley/InterpolationGrid.hh"
#include "marley/JSON.hh"

namespace marley {

  class JSONConfig {

    public:

      using InterpMethod = marley::InterpolationGrid<double>
        ::InterpolationMethod;

      explicit JSONConfig(const marley::JSON& object);

      marley::Generator create_generator() const;

      void prepare_direction(marley::Generator& gen) const;
      void prepare_neutrino_source(marley::Generator& gen) const;
      void prepare_reactions(marley::Generator& gen) const;
      void prepare_structure(marley::Generator& gen) const;

      void update_logger_settings() const;

      InterpMethod get_interpolation_method(const std::string& rule) const;
      int neutrino_pdg(const std::string& nu) const;

      inline virtual bool process_extra_source_types(
        const std::string& /*type*/, const marley::JSON& /*source_spec*/,
        int /*pdg_code*/, std::unique_ptr<marley::NeutrinoSource>& /*source*/)
        const;

      inline const marley::JSON& get_json() const;
      inline void set_json(const marley::JSON& json);

    protected:

      /// @brief Helper function for loading strings from the JSON
      /// configuration
      std::string source_get(const char* name, const marley::JSON& source_spec,
        const char* description, const char* default_str) const;

      /// @brief JSON object describing this configuration
      marley::JSON json_;
  };

  inline const marley::JSON& marley::JSONConfig::get_json() const
    { return json_; }

  inline void marley::JSONConfig::set_json(const marley::JSON& json)
    { json_ = json; }

  inline bool marley::JSONConfig::process_extra_source_types(
    const std::string& /*type*/, const marley::JSON& /*source_spec*/,
    int /*pdg_code*/, std::unique_ptr<marley::NeutrinoSource>& /*source*/)
    const
  { return false; }

}
