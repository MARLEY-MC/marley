// standard library includes
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

// MARLEY includes
#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/FileManager.hh"
#include "marley/JSONConfig.hh"
#include "marley/NeutrinoSource.hh"
#include "marley/NuclearReaction.hh"
#include "marley/Logger.hh"
#include "marley/StructureDatabase.hh"

using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;
using ProcType = marley::Reaction::ProcessType;

// anonymous namespace for helper functions, etc.
namespace {

  void source_check_positive(double x, const char* description,
    const char* source_type)
  {
    if (x <= 0.) throw marley::Error(std::string("Non-positive ") + description
      + " value defined for a " + source_type + " neutrino source");
  }

  void source_check_nonnegative(double x, const char* description,
    const char* source_type)
  {
    if (x < 0.) throw marley::Error(std::string("Negative ") + description
      + " value defined for a " + source_type + " neutrino source");
  }

  double source_get_double(const char* name, const marley::JSON& source_spec,
    const char* description)
  {
    if (!source_spec.has_key(name)) throw marley::Error(
      std::string("Missing source.") + name + " key for " + description
      + " source");
    bool ok;
    double result = source_spec.at(name).to_double(ok);
    if (!ok) throw marley::Error(std::string("Invalid value given for source.")
      + name + " key for " + description + " source");
    return result;
  }

  std::vector<double> get_vector(const char* name, const marley::JSON& spec,
    const char* description)
  {
    if (!spec.has_key(name)) throw marley::Error(std::string("The")
      + " specification for a " + description + " source should include"
      + " a " + name + " key.");

    const marley::JSON& vec = spec.at(name);

    if (!vec.is_array()) throw marley::Error(std::string("The value given")
      + " for the " + name + " key for a " + description + " source should"
      + " be an array.");

    std::vector<double> result;

    auto elements = vec.array_range();
    if (elements.begin() != elements.end()) {
      bool ok;
      for (const auto& el : elements) {
        double dub = el.to_double(ok);
        if (ok) result.push_back(dub);
        else throw marley::Error(std::string("Invalid array entry '")
          + el.dump_string() + "' given for the " + name + " key for a "
          + description + " source specification.");
      }
    }

    return result;
  }

  void handle_json_error(const char* name, const marley::JSON& json) {
    std::ostringstream message;
    message << "The JSON parameter \"" << name << "\" was set to the"
      << " invalid value " << json;
    throw marley::Error(message.str());
  }

}

marley::JSONConfig::JSONConfig(const marley::JSON& json) : json_(json)
{
  update_logger_settings();
}

marley::JSONConfig::JSONConfig(const std::string& json_filename)
{
  json_ = marley::JSON::load_file( json_filename );
  update_logger_settings();
}

int marley::JSONConfig::neutrino_pdg(const std::string& nu) const {

  int pdg = 0;

  bool bad = false;

  // Matches integers
  static const std::regex rx_int = std::regex("[-+]?[0-9]+");
  if (std::regex_match(nu, rx_int)) {
    pdg = std::stoi(nu);
    if (!marley::NeutrinoSource::pdg_is_allowed(pdg)) bad = true;
  }
  else if (!marley_utils::string_to_neutrino_pdg(nu, pdg)) {
    bad = true;
  }

  if (bad) throw marley::Error(std::string("Invalid neutrino type")
    + " specification '" + nu + "' given for the MARLEY"
    + " neutrino source.");

  return pdg;
}

marley::Generator marley::JSONConfig::create_generator() const
{
  uint_fast64_t seed;
  if (json_.has_key("seed")) {
    bool ok;
    seed = static_cast<uint_fast64_t>(json_.at("seed").to_long(ok));
    if (!ok) handle_json_error("seed", json_.at("seed"));
  }
  else seed = std::chrono::system_clock::now().time_since_epoch().count();

  // Start with a default-constructed Generator seeded with either a
  // user-supplied seed or the current number of seconds since the Unix epoch.
  marley::Generator gen(seed);

  // Turn off calls to Generator::normalize_E_pdf() until we
  // have set up all the needed pieces
  gen.dont_normalize_E_pdf_ = true;

  // Use the JSON settings to update the generator's parameters
  prepare_direction(gen);
  prepare_structure(gen);
  prepare_neutrino_source(gen);
  prepare_reactions(gen);

  // Now that the reactions and source are both prepared, check that a neutrino
  // from the source can interact via at least one of the enabled reactions
  bool found_matching_pdg = false;
  int source_pdg = gen.get_source().get_pid();
  for ( const auto& react : gen.get_reactions() ) {
    if ( source_pdg == react->pdg_a() ) found_matching_pdg = true;
  }
  // If neutrinos from the source can never interact, then complain about it
  if ( !found_matching_pdg ) throw marley::Error("The neutrino source"
    " produces " + marley_utils::get_particle_symbol(source_pdg)
    + ", which cannot participate in any of the configured reactions.");

  // Before returning the newly-created Generator object, print logging
  // messages describing the reactions that are active.
  MARLEY_LOG_INFO() << "Generator configuration complete. Active reactions:";
  for ( const auto& r : gen.get_reactions() ) {
    if ( r->pdg_a() == source_pdg ) {

      std::string proc_type_str;
      if ( r->process_type() == ProcType::NeutrinoCC
        || r->process_type() == ProcType::AntiNeutrinoCC )
      {
        proc_type_str = "CC";
      }
      else if ( r->process_type() == ProcType::NC )
      {
        proc_type_str = "NC";
      }
      else if ( r->process_type() == ProcType::NuElectronElastic )
      {
        proc_type_str = "ES";
      }
      else throw marley::Error("Unrecognized process type encountered in"
        " marley::JSONConfig::prepare_reactions()");

      // Show the threshold in red if it's above the maximum energy
      // produced by the source
      std::ostringstream temp_oss;
      double threshold_KE = r->threshold_kinetic_energy();
      bool no_flux = ( threshold_KE > gen.get_source().get_Emax() );
      if ( no_flux ) temp_oss << "\u001b[31m";
      temp_oss << threshold_KE << " MeV";
      if ( no_flux ) temp_oss << "\u001b[30m";
      temp_oss << ')';

      MARLEY_LOG_INFO() << "  " << proc_type_str << ": "
        << r->get_description() << " (KE @ threshold: "
        << temp_oss.str();
    }
  }

  // We've prepared the Generator and checked that at least one cross section
  // should be non-vanishing. Now actually integrate the cross section to check
  // that there is flux above threshold (and normalize the energy PDF at the
  // same time). An exception will be thrown if no neutrinos can interact.
  gen.dont_normalize_E_pdf_ = false;
  gen.normalize_E_pdf();

  return gen;
}

//------------------------------------------------------------------------------
void marley::JSONConfig::prepare_direction(marley::Generator& gen) const {
  // Get the incident neutrino direction if the user has specified one
  if (json_.has_key("direction")) {

    const marley::JSON& direction = json_.at("direction");
    bool ok;

    std::array<double, 3> dir_vec = gen.neutrino_direction();

    if (direction.has_key("x")) {
      dir_vec.at(0) = direction.at("x").to_double(ok);
      if (!ok) handle_json_error("direction.x", direction.at("x"));
    }

    if (direction.has_key("y")) {
      dir_vec.at(1) = direction.at("y").to_double(ok);
      if (!ok) handle_json_error("direction.y", direction.at("y"));
    }

    if (direction.has_key("z")) {
      dir_vec.at(2) = direction.at("z").to_double(ok);
      if (!ok) handle_json_error("direction.z", direction.at("z"));
    }

    gen.set_neutrino_direction(dir_vec);
  }
}

void marley::JSONConfig::prepare_reactions(marley::Generator& gen) const {

  const auto& fm = marley::FileManager::Instance();

  if ( json_.has_key("reactions") ) {

    const marley::JSON& rs = json_.at("reactions");

    if ( rs.is_array() ) {

      auto reactions = rs.array_range();
      if ( reactions.begin() != reactions.end() ) {

        // Create a vector to cache the ProcessType values for
        // reactions that have been loaded. Complain if there is
        // duplication.
        std::vector<ProcType> loaded_proc_types;

        for (const auto& r : reactions) {

          std::string filename = r.to_string();

          // Find the reaction data file using the MARLEY search path
          std::string full_file_name = fm.find_file( filename );
          if ( full_file_name.empty() ) {
            throw marley::Error("Could not locate the reaction data file "
              + filename + ". Please check that the file name is spelled"
              " correctly and that the file is in a folder"
              " on the MARLEY search path.");
          }

          auto nuc_reacts = marley::NuclearReaction::load_from_file(
            full_file_name, gen.get_structure_db());

          // All of the NuclearReaction objects loaded from a single file
          // will have the same process type, so just save the first one
          auto temp_pt = nuc_reacts.front()->process_type();

          // If we have a duplicate, warn the user that we'll ignore it
          auto begin = loaded_proc_types.cbegin();
          auto end = loaded_proc_types.cend();
          if ( std::find(begin, end, temp_pt) != end ) {
            MARLEY_LOG_WARNING() << "Reaction matrix elements for the "
              << marley::Reaction::proc_type_to_string(temp_pt) << " process"
              " were already loaded. To avoid duplication, those in "
              << full_file_name << " will be ignored.";
            continue;
          }
          // Otherwise, save the process type for later checks of this kind
          else {
            MARLEY_LOG_INFO() << "Loaded "
              << marley::Reaction::proc_type_to_string(temp_pt)
              << " reaction data from " << full_file_name;
            loaded_proc_types.push_back( temp_pt );
          }

          // Transfer ownership of the new reactions to the generator
          for ( auto& nr : nuc_reacts ) gen.add_reaction( std::move(nr) );
        }

        return;
      }
      else {
        throw marley::Error("At least one reaction matrix data file must be"
          " specified using the \"reactions\" parameter");
      }
    }

    handle_json_error("reactions", rs);
  }

  throw marley::Error("Missing \"reactions\" key in the MARLEY configuration"
    " file.");
}

void marley::JSONConfig::prepare_structure(marley::Generator& gen) const
{
  // Copy of the structure portion of the JSON configuration.
  // If the user didn't provide one, we'll build the default
  // one ourselves.
  marley::JSON st;

  // If the "structure" key is not present, default to
  // loading data from all files in ${MARLEY}/data/structure
  // (not including files in subdirectories). We don't have
  // to bother to exclude README.md or any other file
  // with the wrong format. All such files will be ignored.
  if ( !json_.has_key("structure") ) {

    // Construct the directory name ${MARLEY}/data/structure
    char* mar = std::getenv("MARLEY");
    if ( !mar ) throw marley::Error("The MARLEY environment"
      " variable was not set in call to marley::JSONConfig::"
      "prepare_structure(). Please set it and try again.");
    std::string structure_dir( mar );
    structure_dir += "/data/structure";

    // Get a vector containing all the file names in that directory
    // (each one with its full path included)
    const auto& fm = marley::FileManager::Instance();
    auto file_vec = fm.list_all_files( structure_dir );

    MARLEY_LOG_INFO() << "Loading default nuclear structure data from "
      << structure_dir;

    // Load these file names into our temporary JSON object before continuing
    st = marley::JSON::array();
    for (const std::string& file_name : file_vec) st.append( file_name );
  }
  // If the user supplied a nuclear structure data configuration, then just use
  // that
  else st = json_.at("structure");

  if (!st.is_array()) return;
  auto structure = st.array_range();

  for (const auto& s : structure) {

    // If the user specifies the structure file as an object, then it
    // will have a format specifier. If it's not a JSON object, assume that
    // it is a string giving the file name, with the default TALYS input
    // format being implied.
    std::string filename;
    auto format = marley::DecayScheme::FileFormat::talys;
    if ( s.is_object() ) {
      if ( s.has_key("file") && s.has_key("format") ) {
        filename = s.at("file").to_string();
        std::string form = s.at("format").to_string();
        if ( form == "native" || form == "marley" ) {
          format = marley::DecayScheme::FileFormat::native;
        }
        else if ( form == "talys" ) {
          format = marley::DecayScheme::FileFormat::talys;
        }
        else handle_json_error("structure.format", s.at("format"));
      }
      else handle_json_error("structure", s);
    }
    else {
      // No structure file format was specified, so assume that it is in MARLEY's
      // native format
      filename = s.to_string();
      format = marley::DecayScheme::FileFormat::native;
    }

    // Load data for all nuclides in the structure data file.
    auto& sdb = gen.get_structure_db();

    int file_nuclide_count = 0;
    if ( format == marley::DecayScheme::FileFormat::talys ) {
      // TALYS format data
      std::set<int> nucleus_PDGs = sdb.find_all_nuclides(filename);

      for (int pdg : nucleus_PDGs) {
        int Z = marley_utils::get_particle_Z(pdg);
        int A = marley_utils::get_particle_A(pdg);

        MARLEY_LOG_DEBUG() << "Loading nuclear structure data for "
          << A << marley_utils::element_symbols.at(Z) << " from file "
          << filename;

        sdb.emplace_decay_scheme(pdg, filename);
        ++file_nuclide_count;
      }
    }
    else if ( format == marley::DecayScheme::FileFormat::native ) {
      // Data in MARLEY's native format
      std::ifstream in_file( filename );

      // Make a dummy decay scheme
      auto temp_ds = std::make_unique<marley::DecayScheme>(0, 0);
      auto& log = marley::Logger::Instance();
      // Use it to read one or more decay scheme objects from the input file.
      // Temporarily disable the logger in an ugly hack to avoid a spurious
      // error message about an invalid parity value when reaching EOF.
      /// @\todo Rewrite this in a better way
      while ( log.disable(), in_file >> *temp_ds )
      {
        log.enable(true);
        // Load each one into the structure database
        int A = temp_ds->A();
        int Z = temp_ds->Z();
        MARLEY_LOG_DEBUG() << "Loading nuclear structure data for "
          << A << marley_utils::element_symbols.at(Z) << " from file "
          << filename;
        int pdg = marley_utils::get_nucleus_pid(Z, A);
        sdb.add_decay_scheme(pdg, temp_ds);
        ++file_nuclide_count;
        if ( in_file.eof() ) break;
        temp_ds.reset( new marley::DecayScheme(0, 0) );
      }
      log.enable(true);
    }

    if ( file_nuclide_count > 0 ) {
      MARLEY_LOG_INFO() << "Loaded structure data for " << file_nuclide_count
        << " nuclides from " << filename;
    }
  }
}

//------------------------------------------------------------------------------
InterpMethod marley::JSONConfig::get_interpolation_method(
  const std::string& rule) const
{
  // Try using the ENDF-style numerical codes first
  static const std::regex rx_nonneg_int("[0-9]+");

  if (std::regex_match(rule, rx_nonneg_int)) {
    int endf_interp_code = std::stoi(rule);
    if (endf_interp_code == 1) return InterpMethod::Constant;
    else if (endf_interp_code == 2) return InterpMethod::LinearLinear;
    else if (endf_interp_code == 3) return InterpMethod::LinearLog;
    else if (endf_interp_code == 4) return InterpMethod::LogLinear;
    else if (endf_interp_code == 5) return InterpMethod::LogLog;
  }

  // Interpolation rules may also be given as strings
  else if (rule == "const" || rule == "constant")
    return InterpMethod::Constant;
  else if (rule == "lin" || rule == "linlin")
    return InterpMethod::LinearLinear;
  else if (rule == "log" || rule == "loglog")
    return InterpMethod::LogLog;
  // linear in energy, logarithmic in probability density
  else if (rule == "linlog")
    return InterpMethod::LinearLog;
  // logarithmic in energy, linear in probability density
  else if (rule == "loglin")
    return InterpMethod::LogLinear;
  else throw marley::Error(std::string("Invalid interpolation rule '")
    + rule + "' given in the neutrino source specification");

  // We shouldn't ever end up here, but return something just in case
  return InterpMethod::Constant;
}

//------------------------------------------------------------------------------
void marley::JSONConfig::prepare_neutrino_source(marley::Generator& gen) const
{
  // Check whether the user provided their own estimate of the source PDF
  // maximum value. If they did, adopt that before building the source.
  // This is useful when automatic searches for the maximum don't work.
  // The user can put in a value manually to get rejection sampling to work.
  if ( json_.has_key("energy_pdf_max" ) ) {
    bool ok;
    const marley::JSON& max_spec = json_.at("energy_pdf_max");
    double user_max = max_spec.to_double( ok );
    if ( !ok ) handle_json_error("energy_pdf_max", max_spec);
    else gen.set_default_E_pdf_max( user_max );
  }

  // Check whether the JSON configuration includes a neutrino source
  // specification
  if ( !json_.has_key("source") ) return;
  const marley::JSON& source_spec = json_.at("source");

  // Complain if the user didn't specify a source type
  if ( !source_spec.has_key("type") ) {
    throw marley::Error(std::string("Missing \"type\" key in")
      + " neutrino source specification.");
    return;
  }

  // Get the neutrino source type
  bool ok;
  std::string type = source_spec.at("type").to_string(ok);
  if ( !ok ) handle_json_error("source.type", source_spec.at("type"));

  // Complain if the user didn't specify a neutrino type
  if (!source_spec.has_key("neutrino")) {
    throw marley::Error(std::string("Missing \"neutrino\" key in")
      + " neutrino source specification.");
    return;
  }
  // Get the neutrino type
  std::string nu = source_spec.at("neutrino").to_string(ok);
  if (!ok) handle_json_error("source.neutrino", source_spec.at("neutrino"));

  // Particle Data Group code for the neutrino type produced by this source
  int pdg = neutrino_pdg(nu);

  std::unique_ptr<marley::NeutrinoSource> source;

  if (type == "mono" || type == "monoenergetic") {
    double energy = source_get_double("energy", source_spec,
      "monoenergetic");
    source_check_positive(energy, "energy", "monoenergetic");
    source = std::make_unique<marley::MonoNeutrinoSource>(pdg, energy);
    MARLEY_LOG_INFO() << "Created monoenergetic "
      << marley_utils::get_particle_symbol(pdg) << " source with"
      << " neutrino energy = " << energy << " MeV";
  }
  else if (type == "dar" || type == "decay-at-rest") {
    source = std::make_unique<marley::DecayAtRestNeutrinoSource>(pdg);
     MARLEY_LOG_INFO() << "Created muon decay-at-rest "
       << marley_utils::get_particle_symbol(pdg) << " source";
  }
  else if (type == "fd" || type == "fermi-dirac" || type == "fermi_dirac") {
    double Emin = source_get_double("Emin", source_spec, "Fermi-Dirac");
    double Emax = source_get_double("Emax", source_spec, "Fermi-Dirac");
    double temp = source_get_double("temperature", source_spec,
      "Fermi-Dirac");

    double eta = 0.;
    if (source_spec.has_key("eta")) eta = source_get_double("eta", source_spec,
      "Fermi-Dirac");

    source_check_nonnegative(Emin, "Emin", "Fermi-Dirac");
    source_check_positive(temp, "temperature", "Fermi-Dirac");

    if (Emax <= Emin) throw marley::Error(std::string("Emax <= Emin")
      + " for a Fermi-Dirac neutrino source");

    source = std::make_unique<marley::FermiDiracNeutrinoSource>(pdg, Emin,
      Emax, temp, eta);
    MARLEY_LOG_INFO() << "Created Fermi-Dirac "
      << marley_utils::get_particle_symbol(pdg) << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << Emin << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << Emax << " MeV";
    MARLEY_LOG_INFO() << "  temperature = " << temp << " MeV";
    MARLEY_LOG_INFO() << "  eta = " << eta;
  }
  else if (type == "bf" || type == "beta" || type == "beta-fit") {
    double Emin = source_get_double("Emin", source_spec, "beta-fit");
    double Emax = source_get_double("Emax", source_spec, "beta-fit");
    double Emean = source_get_double("Emean", source_spec, "beta-fit");

    double beta = 4.5;
    if (source_spec.has_key("beta")) beta = source_get_double("beta", source_spec,
      "beta-fit");

    source_check_nonnegative(Emin, "Emin", "beta-fit");
    source_check_positive(Emean, "Emean", "beta-fit");

    if (Emax <= Emin) throw marley::Error(std::string("Emax <= Emin")
      + " for a beta-fit neutrino source");

    source = std::make_unique<marley::BetaFitNeutrinoSource>(pdg, Emin,
      Emax, Emean, beta);
    MARLEY_LOG_INFO() << "Created beta-fit "
      << marley_utils::get_particle_symbol(pdg) << " source with parameters";
    MARLEY_LOG_INFO() << "  Emin = " << Emin << " MeV";
    MARLEY_LOG_INFO() << "  Emax = " << Emax << " MeV";
    MARLEY_LOG_INFO() << "  average energy = " << Emean << " MeV";
    MARLEY_LOG_INFO() << "  beta = " << beta;
  }
  else if (type == "hist" || type == "histogram") {

    std::vector<double> Es = get_vector("E_bin_lefts", source_spec,
      "histogram");
    std::vector<double> weights = get_vector("weights", source_spec,
      "histogram");

    if (Es.size() != weights.size()) throw marley::Error(std::string("The")
      + " sizes of the arrays of energy bin left edges and weights given"
      + " for a histogram neutrino source are unequal.");

    double Emax = source_get_double("Emax", source_spec, "histogram");
    source_check_positive(Emax, "Emax", "histogram");

    // Add Emax to the grid
    Es.push_back(Emax);

    // Set the probability density at E = Emax to be zero (this ensures
    // that no energies outside of the histogram will be sampled)
    weights.push_back(0.);

    // Convert from bin weights to probability densities by dividing by the
    // width of each bin
    int jmax = Es.size() - 1;
    for (int j = 0; j < jmax; ++j) {

      double width = Es.at(j + 1) - Es.at(j);
      if (width <= 0) throw marley::Error(std::string("Invalid bin width")
        + std::to_string(width) + " encountered when creating a histogram"
        + " neutrino source");

      weights.at(j) /= width;
    }

    // Create the source
    source = std::make_unique<marley::GridNeutrinoSource>(Es, weights, pdg,
      InterpMethod::Constant);
    MARLEY_LOG_INFO() << "Created histogram "
      << marley_utils::get_particle_symbol(pdg) << " source";
  }
  else if (type == "grid") {
    std::vector<double> energies = get_vector("energies", source_spec, "grid");
    std::vector<double> PDs = get_vector("prob_densities", source_spec, "grid");
    std::string rule = source_get("rule", source_spec, "grid", "linlin");

    InterpMethod method = get_interpolation_method(rule);

    source = std::make_unique<marley::GridNeutrinoSource>(energies, PDs, pdg,
      method);
    MARLEY_LOG_INFO() << "Created grid "
      << marley_utils::get_particle_symbol(pdg) << " source";
  }
  else if (!process_extra_source_types(type, source_spec, pdg, source)) {
    throw marley::Error(std::string("Unrecognized MARLEY neutrino source")
      + " type '" + type + "'");
  }

  // If the user has specified whether to weight the incident neutrino spectrum
  // by the reaction cross section(s), then set the weight_flux_ flag in the
  // new Generator object accordingly
  if ( source_spec.has_key("weight_flux") ) {
    bool ok = false;
    bool should_we_weight = source_spec.at("weight_flux").to_bool(ok);
    if (!ok) handle_json_error("source.weight_flux",
      source_spec.at("weight_flux"));
    gen.set_weight_flux(should_we_weight);
  }

  // Load the generator with the new source object
  gen.set_source(std::move(source));

}

std::string marley::JSONConfig::source_get(const char* name,
  const marley::JSON& source_spec, const char* description,
  const char* default_str) const
{
  if (!source_spec.has_key(name)) {
    if (default_str) return default_str;
    else throw marley::Error(std::string("Missing source.") + name
      + " key for " + description + " source");
  }
  bool ok;
  std::string result = source_spec.at(name).to_string(ok);
  if (!ok) throw marley::Error(std::string("Invalid value given for source.")
    + name + " key for " + description + " source");
  return result;
}

void marley::JSONConfig::update_logger_settings() const {

  using LogLevel = marley::Logger::LogLevel;
  marley::Logger& logger = marley::Logger::Instance();

  logger.clear_streams();

  if (!json_.has_key("log")) {
    // If the user hasn't specified a logger configuration, use the default,
    // which is logging at the INFO level to stdout.
    logger.add_stream(std::cout, LogLevel::INFO);
    logger.enable();
    return;
  }
  else {
    const marley::JSON& log_config = json_.at("log");

    if (!log_config.is_array()) {
      throw marley::Error(std::string("The configuration given")
        + " for the \"log\" key should be an array of JSON objects.");
    }

    auto elements = log_config.array_range();
    bool ok;
    // If the user has specified an empty list of logger files, then
    // disable the Logger and return immediately.
    if (elements.begin() == elements.end()) {
      logger.disable();
      return;
    }
    else logger.enable();

    // Loop over the list of log files and add them one-by-one to the Logger.
    for (const auto& el : elements) {
      // Get the file name for the new log file
      if (!el.has_key("file")) throw marley::Error(std::string("Missing")
        + " file name in a log file specification");
      std::string file_name = el.at("file").to_string(ok);
      if (!ok) throw marley::Error("Invalid log file name \""
        + file_name + '\"');

      // Use "info" as the default logging level
      LogLevel level = LogLevel::INFO;
      // Set the logging level for the current file to a non-default value
      // if the user has specified one. Complain if you get confused.
      if (el.has_key("level")) {
        std::string level_str = el.at("level").to_string();
        if (level_str == "error") level = LogLevel::ERROR;
        else if (level_str == "warning") level = LogLevel::WARNING;
        else if (level_str == "info") level = LogLevel::INFO;
        else if (level_str == "debug") level = LogLevel::DEBUG;
        else if (level_str == "disabled") level = LogLevel::DISABLED;
        else throw marley::Error("Invalid logging level \""
          + el.dump_string() + '\"');
      }

      // If the file name is "stdout", then add std::cout as a logging
      // stream. Do the same sort of thing for std::cerr. Otherwise, open the
      // requested file and add it to the logger streams.
      if (file_name == "stdout")
        logger.add_stream(std::cout, level);
      else if (file_name == "stderr")
        logger.add_stream(std::cerr, level);
      else {
        // If the user specified a value for the "overwrite" key, use it
        // to determine whether we should append to the file (false) or
        // overwrite it (true). Otherwise, assume we want to append to it.
        auto file_mode = std::ios::out;
        if (el.has_key("overwrite")) {

          marley::JSON ow = el.at("overwrite");

          bool overwrite = ow.to_bool(ok);
          if (!ok) throw marley::Error("Invalid log file overwrite setting \""
            + ow.dump_string() + '\"');

          if (overwrite) file_mode |= std::ios::trunc;
          else file_mode |= std::ios::app;
        }
        else file_mode |= std::ios::app;

        auto outfile = std::make_shared<std::ofstream>(file_name, file_mode);
        if (!outfile || (!outfile->good())) throw marley::Error("Unable"
          " to open the log file \"" + file_name + "\"");
        else logger.add_stream(outfile, level);
      }
    }
  }
}
