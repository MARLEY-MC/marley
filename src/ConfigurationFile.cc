#include "marley_utils.hh"
#include "ConfigurationFile.hh"
#include "InterpolationGrid.hh"
#include "Logger.hh"

// Static constant definitions
const std::array<double, 3>
  marley::ConfigurationFile::DEFAULT_INCIDENT_NEUTRINO_DIRECTION_
  = { 0., 0., 1.};

// Define type abbreviations, constant regular expressions, and helper
// functions that are only visible within this file
namespace {

  using InterpMethod = marley::InterpolationGrid<double>::InterpolationMethod;

  // Matches comment lines and empty lines
  const std::regex rx_comment_or_empty = std::regex("#.*|\\s*");

  // Matches non-negative integers
  const std::regex rx_nonneg_int = std::regex("[0-9]+");

  // Matches integers
  const std::regex rx_int = std::regex("[-+]?[0-9]+");

  // Matches all numbers, including floats (see
  // http://www.regular-expressions.info/floatingpoint.html, but note that
  // we've modified the version given there to allow for a trailing decimal
  // point [e.g., '1.' matches])
  const std::regex rx_num
    = std::regex("[-+]?[0-9]*\\.?[0-9]*([eE][-+]?[0-9]+)?");

  // Matches trimmed ENSDF-style nucids
  const std::regex rx_nucid = std::regex("[1-9][0-9]{0,2}[A-Za-z]{1,2}");

  /// @brief Convert a string to a marley::DecayScheme::FileFormat
  marley::DecayScheme::FileFormat string_to_format(const std::string& string)
  {
    if (string == "talys")
      return marley::DecayScheme::FileFormat::talys;
    else throw marley::Error(std::string("Unrecognized file")
      + " format '" + string + "' passed to"
      + " string_to_format()");
  }

}

// The default constructor assigns default values to all
// configuration file parameters
marley::ConfigurationFile::ConfigurationFile()
  : seed_(std::chrono::system_clock::now().time_since_epoch().count()),
  structure_db_(new marley::StructureDatabase),
  check_hepevt_overwrite_(true), hepevt_filename_("events.hepevt"),
  num_events_(DEFAULT_NUM_EVENTS_), write_hepevt_(false) {}

// Call the default constructor first to set all configuration
// file parameters to their default values
marley::ConfigurationFile::ConfigurationFile(const std::string& file_name)
  : marley::ConfigurationFile()
{
  filename_ = file_name;
  parse();
}

void marley::ConfigurationFile::parse() {

  file_in_.open(filename_);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in_.good()) {
    throw marley::Error(std::string("Could not read from the ") +
      "file " + filename_);
  }

  // Stores the number of lines checked by each call to
  // marley_utils::get_next_line
  int lines_checked;

  // Loop through each line in the file. Search for keywords,
  // and react appropriately whenever they are encountered.
  for (line_num_ = 0;
    line_ = marley_utils::get_next_line(file_in_, rx_comment_or_empty, false,
      lines_checked), line_num_ += lines_checked,
    line_ != "";)
  {
    // Load the stream with the new line
    iss_.clear();
    iss_.str(line_);

    // Read the first word (the keyword for that line) from the stream
    iss_ >> keyword_;
    // Convert the keyword to all lowercase (this is an easy way of
    // getting case-insensitivity). We do the case conversion word-by-word
    // rather than as a whole line because case is important for certain
    // items in the input (particularly file names).
    marley_utils::to_lowercase_inplace(keyword_);

    // String to store the keyword argument currently being processed
    std::string arg;

    // Respond accordingly depending on which keyword it is
    if (keyword_ == "seed") {
      next_word_from_line(arg);
      // If "time" is used as the seed argument, use the system time as the
      // seed
      if (arg == "time") ; // default value, so do nothing
      // If "device" is used, seed the random number generator using
      // random_device
      // TODO: implement this
      //else if (arg == "device");
      else if (!std::regex_match(arg, rx_nonneg_int))
        throw marley::Error(std::string("Invalid random number seed")
          + "encountered on line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);
      else {
        // Convert the seed string to a 64-bit unsigned integer
        seed_ = static_cast<uint_fast64_t>(std::stoull(arg));
      }
    }
    else if (keyword_ == "reaction") {
      next_word_from_line(arg, true, false);
      reaction_filenames_.insert(arg);
    }
    else if (keyword_ == "hepevtfile") {
      next_word_from_line(arg, true, false);
      hepevt_filename_ = arg;
    }
    else if (keyword_ == "writehepevt") {
      next_word_from_line(arg, true, true);
      if (arg == "yes") write_hepevt_ = true;
      else if (arg == "no") write_hepevt_ = false;
      else if (arg == "overwrite") {
        write_hepevt_ = true;
        check_hepevt_overwrite_ = false;
      }
      else {
        throw marley::Error(std::string("Invalid")
          + " HEPEvt file write flag '" + arg
          + "' encountered on line" + std::to_string(line_num_)
          + " of the configuration file " + filename_);
      }
    }
    else if (keyword_ == "structure") {

      std::string structure_file_name;
      next_word_from_line(structure_file_name, true, false);

      std::string format_string;
      marley::DecayScheme::FileFormat format;
      next_word_from_line(format_string);
      try {
        format = string_to_format(format_string);
      }
      // If the format is invalid, catch the exception thrown
      // by the converter and throw a different one so that
      // the error message is more easily understood
      catch (const marley::Error& e) {
         throw marley::Error(std::string("Unknown")
          + " nuclear structure format '" + format_string
          + "' specified on line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);
      }

      std::set<int> nucleus_pdg_codes;

      // Require at least one nuclide specifier on this line
      // by using next_word_from_line here with exceptions enabled
      // and then using it with exceptions disabled in the loop
      next_word_from_line(arg);
      do {
        // If a numeric specifier is given for the nuclide,
        // assume that it has the format Z*1000 + A.
        if (std::regex_match(arg, rx_nonneg_int)) {
          int dummy = std::stoi(arg);
          int Z = dummy / 1000;
          int A = dummy % 1000;
          nucleus_pdg_codes.insert(marley_utils::get_nucleus_pid(Z, A));
        }
        // If the nuclide specifier matches the ENSDF-style
        // nucid format (A + element symbol without a
        // space in between them) then parse it, adding
        // spaces if needed.
        else if (std::regex_match(arg, rx_nucid)) {
          // Valid ENSDF nucids must be all uppercase
          marley_utils::to_uppercase_inplace(arg);

          // Enforce ENSDF's requirement that nucids are 5 characters
          // long by adding spaces as needed
          if (arg.length() != 5) {
            // Split the string into "A" and "element name" pieces
            std::smatch m;
            std::regex_search(arg, m, rx_nonneg_int);
            std::string a_str = m.str();
            std::string e_str = m.suffix().str();

            // Pad each piece so that it has the proper length based
            // on ENSDF conventions
            if (a_str.size() < 3) marley_utils::pad_left_inplace(a_str, 3);
            if (e_str.size() < 2) marley_utils::pad_right_inplace(e_str, 2);
            arg = a_str + e_str;
          }
          int Znuc = marley_utils::nucid_to_Z(arg);
          int Anuc = std::stoi(arg.substr(0,3));
          // Add the nucid to the set of requested nucids
          // for this structure file
          nucleus_pdg_codes.insert(marley_utils::get_nucleus_pid(Znuc, Anuc));
        }
        // Other nuclide specifiers are not allowed
        else throw marley::Error(std::string("Invalid")
          + " nuclide specifier '" + arg
          + "' given on line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);

      } while (next_word_from_line(arg, false));

      // Create DecayScheme objects for the requested nuclides and add them to
      // the database
      add_decay_schemes(structure_file_name, format, nucleus_pdg_codes);
    }
    else if (keyword_ == "source") {

      // PDG particle ID for this source
      int neutrino_pid;

      // Get the particle ID for the neutrino type produced by this source
      next_word_from_line(arg, true, false);
      if (std::regex_match(arg, rx_int)) {
        neutrino_pid = stoi(arg);
        if (!marley::NeutrinoSource::pid_is_allowed(neutrino_pid)) {
          throw marley::Error(std::string("Unallowed")
            + " neutrino source particle ID '" + arg
            + "' given on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        }
      }
      // Anything other than an integer entry here is not allowed
      else throw marley::Error(std::string("Invalid")
        + " neutrino source particle ID '" + arg
        + "' given on line " + std::to_string(line_num_)
        + " of the configuration file " + filename_);

      next_word_from_line(arg, true, true);

      // Determine which source style was requested, process the remaining
      // arguments appropriately, and construct the new marley::NeutrinoSource
      if (arg == "mono" || arg == "monoenergetic") {
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          double nu_energy = std::stod(arg);
          if (nu_energy <= 0.) throw marley::Error(
            std::string("Non-positive") + " energy '" + arg
            + "' given for a monoenergetic neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
          // Create the monoenergetic neutrino source
          else set_source(std::make_unique<marley::MonoNeutrinoSource>(
            neutrino_pid, nu_energy));
        }
        else throw marley::Error(std::string("Invalid")
          + " energy '" + arg
          + "' given for a monoenergetic neutrino source specification on line "
          + std::to_string(line_num_) + " of the configuration file "
          + filename_);
      }

      else if (arg == "dar" || arg == "decay-at-rest") {
        set_source(std::make_unique<marley::DecayAtRestNeutrinoSource>(
          neutrino_pid));
      }

      else if (arg == "fd" || arg == "fermi-dirac" || arg == "fermi_dirac") {
        double Emin, Emax, temperature, eta;
        // Process the minimum energy for this source
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emin = std::stod(arg);
          if (Emin < 0.) throw marley::Error(
            std::string("Negative") + " minimum energy " + std::to_string(Emin)
            + " given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        }
        else throw marley::Error(
            std::string("Invalid") + " minimum energy '" + arg
            + "' given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        // Process the maximum energy for this source
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emax = std::stod(arg);
          if (Emax < Emin) throw marley::Error(
            std::string("Maximum") + " energy " + std::to_string(Emax)
            + " given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_
            + " is less than the minimum energy " + std::to_string(Emin));
        }
        else throw marley::Error(
            std::string("Invalid") + " maximum energy '" + arg
            + "' given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        // Process the temperature (in MeV) for this source
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          temperature = std::stod(arg);
          if (temperature <= 0.) throw marley::Error(
            std::string("Non-positive") + " temperature "
            + std::to_string(temperature)
            + " given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        }
        else throw marley::Error(
            std::string("Invalid") + " temperature '" + arg
            + "' given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        // Process the dimensionless pinching parameter eta for this source, or
        // set it to zero if one is not given. Any real number is allowed for
        // this parameter.
        if (next_word_from_line(arg, false, false)) {
          if (std::regex_match(arg, rx_num)) eta = std::stod(arg);
          else throw marley::Error(
              std::string("Invalid") + " pinching parameter '" + arg
              + "' given for a Fermi-Dirac neutrino source"
              + " specification on line " + std::to_string(line_num_)
              + " of the configuration file " + filename_);
        }
        else eta = 0.;

        // Now that we have all of the necessary parameters, create the new
        // Fermi-Dirac neutrino source
        set_source(std::make_unique<marley::FermiDiracNeutrinoSource>(
          neutrino_pid, Emin, Emax, temperature, eta));
      }

      else if (arg == "bf" || arg == "beta" || arg == "beta-fit") {
        double Emin, Emax, Eavg, beta;
        // Process the minimum energy for this source
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emin = std::stod(arg);
          if (Emin < 0.) throw marley::Error(
            std::string("Negative") + " minimum energy " + std::to_string(Emin)
            + " given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        }
        else throw marley::Error(
            std::string("Invalid") + " minimum energy '" + arg
            + "' given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        // Process the maximum energy for this source
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emax = std::stod(arg);
          if (Emax < Emin) throw marley::Error(
            std::string("Maximum") + " energy " + std::to_string(Emax)
            + " given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_
            + " is less than the minimum energy " + std::to_string(Emin));
        }
        else throw marley::Error(
            std::string("Invalid") + " maximum energy '" + arg
            + "' given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        // Process the temperature (in MeV) for this source
        next_word_from_line(arg, true, false);
        if (std::regex_match(arg, rx_num)) {
          Eavg = std::stod(arg);
          if (Eavg <= 0.) throw marley::Error(
            std::string("Non-positive") + " average energy "
            + std::to_string(Eavg)
            + " given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        }
        else throw marley::Error(
            std::string("Invalid") + " average energy '" + arg
            + "' given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        // Process the fit parameter beta for this source, or
        // set it to 4.5 if it is not given. Any real number is allowed for
        // this parameter.
        if (next_word_from_line(arg, false, false)) {
          if (std::regex_match(arg, rx_num)) beta = std::stod(arg);
          else throw marley::Error(
              std::string("Invalid") + " fit parameter '" + arg
              + "' given for a beta-fit neutrino source"
              + " specification on line " + std::to_string(line_num_)
              + " of the configuration file " + filename_);
        }
        else beta = 4.5;

        // Now that we have all of the necessary parameters, create the new
        // Fermi-Dirac neutrino source
        set_source(std::make_unique<marley::BetaFitNeutrinoSource>(
          neutrino_pid, Emin, Emax, Eavg, beta));
      }

      // The histogram and grid source types are both implemented using an
      // interpolation grid object and share many similarities, so use
      // the same parsing code for both
      else if (arg == "hist" || arg == "histogram" || arg == "grid") {

        bool source_is_histogram = (arg == "hist" || arg == "histogram");
        std::string description;
        if (source_is_histogram) description = "histogram";
        else description = "grid";

        InterpMethod method;

        // If this is a histogram, the interpolation method is already known.
        if (source_is_histogram) method = InterpMethod::Constant;

        // If this is a grid, then parse the interpolation method to use;
        else {
          next_word_from_line(arg, true, true);
          // Specifying an interpolation method using the standard ENDF codes
          // is allowed.
          if (std::regex_match(arg, rx_nonneg_int)) {
            int endf_interp_code = std::stoi(arg);
            if (endf_interp_code == 1) method = InterpMethod::Constant;
            else if (endf_interp_code == 2) method = InterpMethod::LinearLinear;
            else if (endf_interp_code == 3) method = InterpMethod::LinearLog;
            else if (endf_interp_code == 4) method = InterpMethod::LogLinear;
            else if (endf_interp_code == 5) method = InterpMethod::LogLog;
            else throw marley::Error(
              std::string("Invalid") + " interpolation method '" + arg
              + "' given for a grid neutrino source "
              + " specification on line " + std::to_string(line_num_)
              + " of the configuration file " + filename_);
          }
          // Specifying it using a string is also allowed.
          // Note that "hist" should probably be used to enter grids with a
          // constant interpolation rule, but we can include this option for
          // flexibility.
          else if (arg == "const" || arg == "constant")
            method = InterpMethod::Constant;
          else if (arg == "lin" || arg == "linlin")
            method = InterpMethod::LinearLinear;
          else if (arg == "log" || arg == "loglog")
            method = InterpMethod::LogLog;
          // linear in energy, logarithmic in probability density
          else if (arg == "linlog")
            method = InterpMethod::LinearLog;
          // logarithmic in energy, linear in probability density
          else if (arg == "loglin")
            method = InterpMethod::LogLinear;
          else throw marley::Error(
            std::string("Invalid") + " interpolation method '" + arg
            + "' given for a grid neutrino source "
            + " specification on line " + std::to_string(line_num_)
            + " of the configuration file " + filename_);
        }

        // Store the grid point energies and probability densities
        // in these vectors for now.
        std::vector<double> energies, prob_densities;

        // Counter to keep track of the number of bins or grid points parsed so
        // far
        int points_parsed = 0;

        // Flag to tell us when we're done reading in bin or grid point entries
        bool not_finished = true;

        // Get the bin or grid point specifications using these strings
        std::string E_str, PD_str;

	// Previous energy value (used to ensure that all energy values given
	// are strictly increasing). It's set to negative infinity to ensure
	// that the first energy value always passes the strictly increasing
	// test.
        double previous_E = marley_utils::minus_infinity;

        // Parse the bin or grid point specifications for this source
        while(not_finished) {

          // Get the first energy
          bool found_E = next_word_from_line(E_str, false, false);
          if (found_E) {
            if (std::regex_match(E_str, rx_num)) {
              double E = std::stod(E_str);
              if (E < 0.) throw marley::Error(
                std::string("Negative") + " energy "
                + std::to_string(E) + " given in a "
                + description + " neutrino source specification on line "
                + std::to_string(line_num_) + " of the configuration file "
                + filename_);
              if (E <= previous_E) throw marley::Error(
                std::string("Energy") + " values "
                + " given in a " + description
                + " neutrino source specification on line "
                + std::to_string(line_num_) + " of the configuration file "
                + filename_ + " are not strictly increasing");
              energies.push_back(E);
              previous_E = E;
            }
            else throw marley::Error(std::string("Invalid") + " energy '"
              + E_str + "' given in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num_) + " of the configuration file "
              + filename_);
          }
          // Get the first probability density
          bool found_PD = next_word_from_line(PD_str, false, false);
          if (found_PD) {
            if (std::regex_match(PD_str, rx_num)) {
              double pd = std::stod(PD_str);
              if (pd < 0.) throw marley::Error(
                std::string("Negative") + " probability density "
                + std::to_string(pd) + " given in a "
                + description + " neutrino source specification on line "
                + std::to_string(line_num_) + " of the configuration file "
                + filename_);
              prob_densities.push_back(pd);
            }
            else throw marley::Error(std::string("Invalid")
              + " probability density '"
              + PD_str + "' given in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num_) + " of the configuration file "
              + filename_);
          }

          // Stop the loop for a histogram if we read an energy without
          // a matching probability density and at least one point has
          // already been parsed.
          if (found_E && !found_PD && points_parsed >= 1
            && source_is_histogram)
          {
            not_finished = false;
          }
          // Stop the loop for a grid if we have run out of entries
          // to read and at least two points have already been parsed.
          else if (!found_E && !found_PD && points_parsed >= 2
            && !source_is_histogram)
          {
            not_finished = false;
          }
          else if (found_E && !found_PD) {
            std::string weight_desc;
            if (source_is_histogram) weight_desc = "bin weight";
            else weight_desc = "probability density";
            throw marley::Error(std::string("Missing ")
              + weight_desc + " entry in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num_) + " of the configuration file "
              + filename_);
          }
          else if (!found_E) {
            throw marley::Error(std::string("Missing ")
              + "energy entry in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num_) + " of the configuration file "
              + filename_);
          }

          // Increment the entry counter now that we've successfully parsed
          // another point
          ++points_parsed;
        }

        // If this is a histogram, we still need to convert the bin weights
        // into probability density values, so we'll divide each weight by its
        // bin width.
        if (source_is_histogram) {
          int jmax = energies.size() - 1;
          for (int j = 0; j < jmax; ++j) {
            // TODO: The earlier check for monotonically increasing energies
            // should ensure that we don't divide by zero here, but this might
            // be worth looking at more carefully in the future.
            prob_densities.at(j) /= (energies.at(j + 1) - energies.at(j));
          }
          // Set the probability density at E = Emax to be zero (this ensures
          // that no energies outside of the histogram will be sampled and
          // makes the energy and probability densities vector lengths equal).
          prob_densities.push_back(0.);
        }

        // Now that we've processed grid points, create the grid neutrino
        // source
        set_source(std::make_unique<marley::GridNeutrinoSource>(
          energies, prob_densities, neutrino_pid, method));
      }

      else if (!process_extra_source_types(arg, neutrino_pid))
        throw marley::Error(std::string("Unrecognized")
        + " neutrino source type '" + arg
        + "' given on line " + std::to_string(line_num_)
        + " of the configuration file " + filename_);
    }

    else if (keyword_ == "events") {
      next_word_from_line(arg, true, false);
      if (!std::regex_match(arg, rx_num))
        throw marley::Error(std::string("Invalid")
        + " number of events '" + arg
        + "' given on line " + std::to_string(line_num_)
        + " of the configuration file " + filename_);
      double d_n_events = std::stod(arg);
      int n_events = static_cast<int>(d_n_events);
      // Allow n_events == 0 for testing purposes
      if (n_events < 0)
        throw marley::Error(std::string("Number")
        + " of events '" + arg + "' given on line "
        + std::to_string(line_num_) + " of the configuration file "
        + filename_ + " must be non-negative.");

      num_events_ = static_cast<size_t>(n_events);
    }
    else if (keyword_ == "direction") {
      for (size_t i = 0; i < 3; ++i) {
        next_word_from_line(arg, true, false);
        if (!std::regex_match(arg, rx_num))
          throw marley::Error(std::string("Invalid component #")
          + std::to_string(i + 1) + " of the"
          + " incident neutrino direction '" + arg
          + "' given on line " + std::to_string(line_num_)
          + " of the configuration file " + filename_);
        dir_vec_[i] = std::stod(arg);
      }
      static constexpr std::array<double, 3> null_vec = { 0., 0., 0. };
      if (dir_vec_ == null_vec) {
        throw marley::Error(std::string("Null")
          + " vector used with the 'direction' keyword on line "
          + std::to_string(line_num_) + " of the configuration file "
          + filename_);
      }
    }

    else if (!process_extra_keywords()) {
      MARLEY_LOG_WARNING() << "Ignoring unrecognized keyword '" << keyword_
        << "' on line " << line_num_ << " of the configuration file "
        << filename_ << '\n';
    }
  }

  // Check that required keywords were found
  if (reaction_filenames_.empty())
    throw marley::Error(std::string("Configuration file")
    + " must contain at least one use of the reaction keyword.");
  if (!source_)
    throw marley::Error(std::string("Configuration file")
    + " must contain at least one use of the source keyword.");
}

// Get the next word from a parsed line. If errors occur, complain.
// The last argument determines whether the next word should be
// converted to all lowercase or left as is.
bool marley::ConfigurationFile::next_word_from_line(std::string& word,
  bool enable_exceptions, bool make_lowercase)
{
  if (iss_.eof()) {
    if (!enable_exceptions) return false;
    throw marley::Error(std::string("Missing argument ")
      + "for keyword '" + keyword_ + "' encountered while"
      + " parsing line " + std::to_string(line_num_)
      + " of the configuration file " + filename_);
  }
  iss_ >> word;
  if (iss_.fail()) {
    if (!enable_exceptions) return false;
    std::string snum = std::to_string(line_num_);
    throw marley::Error(std::string("Error")
      + " occurred while parsing arguments for keyword '"
      + keyword_ + "' on line "
      + std::to_string(line_num_)
      + " of the configuration file " + filename_);
  }
  else if (word == "&") {
    // A single "&" by itself is the line continuation flag. Anything after it
    // on the current line is ignored, and the next line is treated as a
    // continuation of the previous line
    int lines_checked = 0;
    line_ = marley_utils::get_next_line(file_in_, rx_comment_or_empty, false,
      lines_checked);
    line_num_ += lines_checked,
    iss_.clear();
    iss_.str(line_);
    return next_word_from_line(word, enable_exceptions, make_lowercase);
  }
  else if (make_lowercase) marley_utils::to_lowercase_inplace(word);
  return true;
}

bool marley::ConfigurationFile::process_extra_keywords() {
  if (keyword_ == "writeroot" || keyword_ == "rootfile") {
    MARLEY_LOG_WARNING() << "Ignoring keyword '" << keyword_
      << "' found on line " << std::to_string(line_num_)
      << " of the configuration file. To enable ROOT"
      << " TFile output, use a marley::RootConfigurationFile object instead"
      << " of a marley::ConfigurationFile object.\n";
    return true;
  }
  return false;
}

void marley::ConfigurationFile::add_decay_schemes(const std::string& file_name,
  const marley::DecayScheme::FileFormat format,
  const std::set<int>& nucleus_pdg_codes)
{
  /// @todo Consider altering the parsing process used here so
  /// that all nuclides are loaded from the file in one pass.

  // Add a decay scheme for each nuclear pdg code to the database
  for (const auto& pdg : nucleus_pdg_codes) {

    int Z = marley_utils::get_particle_Z(pdg);
    int A = marley_utils::get_particle_A(pdg);

    std::string trimmed_nucid = marley_utils::trim_copy(
      marley_utils::nuc_id(Z, A));
    MARLEY_LOG_INFO() << "Loading nuclear structure data for "
      << trimmed_nucid << " from file " << file_name;

    structure_db_->emplace_decay_scheme(pdg, file_name, format);
  }
}

void marley::ConfigurationFile::set_source(
  std::unique_ptr<marley::NeutrinoSource>&& new_source)
{
  if (source_) MARLEY_LOG_WARNING() << "Replacing existing NeutrinoSource"
    << " in a ConfigurationFile object.";
  source_.reset(new_source.release());
}

void marley::ConfigurationFile::set_neutrino_direction(
  const std::array<double, 3>& dir_vec)
{
  static constexpr std::array<double, 3> null_vector = { 0., 0., 0. };
  if (dir_vec == null_vector) {
    throw marley::Error(std::string("Null")
      + " three-vector passed to"
      + " marley::ConfigurationFile::set_neutrino_direction()");
  }
  else dir_vec_ = dir_vec;
}
