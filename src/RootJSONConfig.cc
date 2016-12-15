// MARLEY includes
#include "marley/marley_root.hh"
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/NeutrinoSource.hh"
#include "marley/RootJSONConfig.hh"
#include "marley/Logger.hh"

bool marley::RootJSONConfig::process_extra_source_types(
  const std::string& type, const marley::JSON& source_spec, int pdg_code,
  std::unique_ptr<marley::NeutrinoSource>& source) const
{
  if (type != "th1" && type != "tgraph") return false;

  std::string tfile = source_get("tfile", source_spec, type.c_str(),
    nullptr);
  std::string namecycle = source_get("namecycle", source_spec,
    type.c_str(), nullptr);

  if (type == "th1") {
    auto th1 = marley_root::get_root_object<TH1>(tfile, namecycle);
    source = marley_root::make_root_neutrino_source(pdg_code, th1);
    MARLEY_LOG_INFO() << "Created a TH1 "
      << marley_utils::neutrino_pdg_to_string(pdg_code)
      << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << source->get_Emin() << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << source->get_Emax() << " MeV";
    return true;
  }

  else if (type == "tgraph") {
    auto tg = marley_root::get_root_object<TGraph>(tfile, namecycle);
    source = marley_root::make_root_neutrino_source(pdg_code, tg);
    MARLEY_LOG_INFO() << "Created a TGraph "
      << marley_utils::neutrino_pdg_to_string(pdg_code)
      << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << source->get_Emin() << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << source->get_Emax() << " MeV";
    return true;
  }

  // We shouldn't ever get here, but I'll include something just in case.
  return false;
}
