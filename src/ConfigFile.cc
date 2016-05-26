#include <thread>

#include "marley_utils.hh"
#include "ConfigFile.hh"
#include "InterpolationGrid.hh"

// Define type abbreviations and constant regular expressions that are only
// visible within this file
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
  const std::regex rx_num = std::regex("[-+]?[0-9]*\\.?[0-9]*([eE][-+]?[0-9]+)?");
  // Matches trimmed ENSDF-style nucids
  const std::regex rx_nucid = std::regex("[1-9][0-9]{0,2}[A-Za-z]{1,2}");
}

// The default constructor assigns default values to all
// configuration file parameters
marley::ConfigFile::ConfigFile() {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  writeroot = false;
  check_before_root_file_overwrite = true;
  root_filename = "events.root";
  writehepevt = false;
  check_before_hepevt_file_overwrite = true;
  hepevt_filename = "events.hepevt";
/***
  contbin_width = DEFAULT_CONTINUUM_BIN_RESOLUTION;
  contbin_num_subs = DEFAULT_CONTINUUM_BIN_SUBINTERVALS;
  num_threads = 1;
***/
  num_events = DEFAULT_NUM_EVENTS;
}

// Call the default constructor first to set all configuration
// file parameters to their default values
marley::ConfigFile::ConfigFile(std::string file_name)
  : marley::ConfigFile()
{

  filename = file_name;
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw marley::Error(std::string("Could not read from the ") +
      "file " + filename);
  }

  std::string line; // String to store the current line
                    // of the configuration file during parsing

  std::string keyword; // String to store the current keyword read in
                       // from the configuration file

  std::istringstream iss; // Stream used to help parse each line

  // Stores the number of lines checked by each call to
  // marley_utils::get_next_line
  int lines_checked;

  // Loop through each line in the file. Search for keywords,
  // and react appropriately whenever they are encountered. 
  for (int line_num = 0;
    line = marley_utils::get_next_line(file_in, rx_comment_or_empty, false,
      lines_checked), line_num += lines_checked,
    line != ""; /*++line_num*/)
  {
    // Load the stream with the new line
    iss.clear();
    iss.str(line);

    // Read the first word (the keyword for that line) from the stream
    iss >> keyword;
    // Convert the keyword to all lowercase (this is an easy way of
    // getting case-insensitivity). We do the case conversion word-by-word
    // rather than as a whole line because case is important for certain
    // items in the input (particularly file names).
    marley_utils::to_lowercase_inplace(keyword);

    // String to store the keyword argument currently being processed
    std::string arg;

    // Respond accordingly depending on which keyword it is
    if (keyword == "seed") {
      next_word_from_line(iss, arg, keyword, line_num);
      // If "time" is used as the seed argument, use the system time as the seed
      if (arg == "time") ; // default value, so do nothing
      // If "device" is used, seed the random number generator using random_device 
      // TODO: implement this
      //else if (arg == "device");
      else if (!std::regex_match(arg, rx_nonneg_int))
        throw marley::Error(std::string("Invalid random number seed")
          + "encountered on line " + std::to_string(line_num)
          + " of the configuration file " + filename);
      else {
        // Convert the seed string to a 64-bit unsigned integer
        seed = static_cast<uint_fast64_t>(std::stoull(arg));
      }
    }
    else if (keyword == "reaction") {
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      reaction_filenames.insert(arg);
    }
    else if (keyword == "rootfile") {
      #ifdef USE_ROOT
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      root_filename = arg;
      #else
      // TODO: change this to a warning that this will be ignored
      throw marley::Error(std::string("Cannot set")
        + " ROOT filename. MARLEY must be compiled"
        + " with ROOT support.");
      #endif
    }
    else if (keyword == "hepevtfile") {
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      hepevt_filename = arg;
    }
    else if (keyword == "writeroot") {
      #ifdef USE_ROOT
      next_word_from_line(iss, arg, keyword, line_num, true, true);
      if (arg == "yes") writeroot = true;
      else if (arg == "no") writeroot = false;
      else if (arg == "overwrite") {
        writeroot = true;
        check_before_root_file_overwrite = false;
      }
      else {
        throw marley::Error(std::string("Invalid")
          + " ROOT file write flag '" + arg
          + "' encountered on line" + std::to_string(line_num)
          + " of the configuration file " + filename);
      }
      #else
      throw marley::Error(std::string("The writeroot")
        + " keyword may only be used when MARLEY is compiled"
        + " with ROOT support.");
      #endif
    }
    else if (keyword == "writehepevt") {
      next_word_from_line(iss, arg, keyword, line_num, true, true);
      if (arg == "yes") writehepevt = true;
      else if (arg == "no") writehepevt = false;
      else if (arg == "overwrite") {
        writehepevt = true;
        check_before_hepevt_file_overwrite = false;
      }
      else {
        throw marley::Error(std::string("Invalid")
          + " HEPEvt file write flag '" + arg
          + "' encountered on line" + std::to_string(line_num)
          + " of the configuration file " + filename);
      }
    }
    else if (keyword == "structure") {
      StructureRecord sr;
      next_word_from_line(iss, sr.filename, keyword, line_num, true, false);
      std::string format_string;
      next_word_from_line(iss, format_string, keyword, line_num);
      try {
        sr.format = string_to_format(format_string);
      }
      // If the format is invalid, catch the exception thrown
      // by the converter and throw a different one so that
      // the error message is more easily understood
      catch (const marley::Error& e) {
         throw marley::Error(std::string("Unknown")
          + " nuclear structure format '" + format_string
          + "' specified on line " + std::to_string(line_num)
          + " of the configuration file " + filename);
      }
      // Require at least one nuclide specifier on this line
      // by using next_word_from_line here with exceptions enabled
      // and then using it with exceptions disabled in the loop
      next_word_from_line(iss, arg, keyword, line_num);
      do {
        // If a numeric specifier is given for the nuclide,
        // assume that it has the format Z*1000 + A, and
        // create an ENSDF nucid accordingly
        if (std::regex_match(arg, rx_nonneg_int)) {
          int dummy = std::stoi(arg);
          int Z = dummy / 1000;
          int A = dummy % 1000;
          sr.nucids.insert(marley_utils::nuc_id(Z, A));
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
            std::string z_str = m.str();
            std::string e_str = m.suffix().str();

            // Pad each piece so that it has the proper length based
            // on ENSDF conventions
            if (z_str.size() < 3) marley_utils::pad_left_inplace(z_str, 3); 
            if (e_str.size() < 2) marley_utils::pad_right_inplace(e_str, 2);
            arg = z_str + e_str;
          }
          // Add the nucid to the set of requested nucids
          // for this structure file
          sr.nucids.insert(arg);
        }
        // Other nuclide specifiers are not allowed
        else throw marley::Error(std::string("Invalid")
          + " nuclide specifier '" + arg
          + "' given on line " + std::to_string(line_num)
          + " of the configuration file " + filename);

      } while (next_word_from_line(iss, arg,
        keyword, line_num, false));

      // Add this structure record to the master vector
      structure_records.push_back(sr);

    }
    else if (keyword == "source") {

      // Weight and PDG particle ID for this source
      double weight;
      int neutrino_pid;

      // Get the particle ID for the neutrino type produced by this source
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      if (std::regex_match(arg, rx_int)) {
        neutrino_pid = stoi(arg);
        if (!marley::NeutrinoSource::pid_is_allowed(neutrino_pid)) {
          throw marley::Error(std::string("Unallowed")
            + " neutrino source particle ID '" + arg
            + "' given on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        }
      }
      // Anything other than an integer entry here is not allowed
      else throw marley::Error(std::string("Invalid")
        + " neutrino source particle ID '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename);

/***  Do not allow the user to specify a source weight for now, since
 * only one neutrino source is currently allowed. Re-enable this code
 * if and when you add the capability for marley::Generator to handle
 * multiple sources concurrently.
      // Get the weight for this source if it is specified. Otherwise,
      // set the weight to 1.
      next_word_from_line(iss, arg, keyword, line_num, true, true);
      bool weight_specified = false;
      if (std::regex_match(arg, rx_num)) {
        weight = std::stod(arg);
        weight_specified = true;
        // Only positive weights are allowed
        if (weight <= 0.) {
          throw marley::Error(std::string("Non-positive")
          + " neutrino source weight '" + arg
          + "' given on line " + std::to_string(line_num)
          + " of the configuration file " + filename);
        }
      }
      else weight = 1.;

      // Advance to the next argument if the weight was given
      if (weight_specified) next_word_from_line(iss, arg, keyword, line_num,
        true, true);
*/

/* Replacement code (delete when you re-enable the block above) */

      weight = 1.;
      next_word_from_line(iss, arg, keyword, line_num, true, true);
/* End of replacement code */

      // Determine which source style was requested, process the remaining
      // arguments appropriately, and construct the new marley::NeutrinoSource
      if (arg == "mono" || arg == "monoenergetic") {
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          double nu_energy = std::stod(arg);
          if (nu_energy <= 0.) throw marley::Error(
            std::string("Non-positive") + " energy '" + arg
            + "' given for a monoenergetic neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
          // Create the monoenergetic neutrino source
          else sources.push_back(std::make_unique<marley::MonoNeutrinoSource>(
            neutrino_pid, weight, nu_energy));
        }
        else throw marley::Error(std::string("Invalid")
          + " energy '" + arg
          + "' given for a monoenergetic neutrino source specification on line "
          + std::to_string(line_num) + " of the configuration file "
          + filename);
      }

      else if (arg == "dar" || arg == "decay-at-rest") {
        sources.push_back(std::make_unique<marley::DecayAtRestNeutrinoSource>(
          weight, neutrino_pid));
      }

      else if (arg == "fd" || arg == "fermi-dirac" || arg == "fermi_dirac") {
        double Emin, Emax, temperature, eta;
        // Process the minimum energy for this source
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emin = std::stod(arg);
          if (Emin < 0.) throw marley::Error(
            std::string("Negative") + " minimum energy " + std::to_string(Emin)
            + " given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        }
        else throw marley::Error(
            std::string("Invalid") + " minimum energy '" + arg
            + "' given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        // Process the maximum energy for this source
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emax = std::stod(arg);
          if (Emax < Emin) throw marley::Error(
            std::string("Maximum") + " energy " + std::to_string(Emax)
            + " given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename
            + " is less than the minimum energy " + std::to_string(Emin));
        }
        else throw marley::Error(
            std::string("Invalid") + " maximum energy '" + arg
            + "' given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        // Process the temperature (in MeV) for this source
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          temperature = std::stod(arg);
          if (temperature <= 0.) throw marley::Error(
            std::string("Non-positive") + " temperature "
            + std::to_string(temperature)
            + " given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        }
        else throw marley::Error(
            std::string("Invalid") + " temperature '" + arg
            + "' given for a Fermi-Dirac neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        // Process the dimensionless pinching parameter eta for this source, or
        // set it to zero if one is not given. Any real number is allowed for
        // this parameter.
        if (next_word_from_line(iss, arg, keyword, line_num, false, false)) {
          if (std::regex_match(arg, rx_num)) eta = std::stod(arg);
          else throw marley::Error(
              std::string("Invalid") + " pinching parameter '" + arg
              + "' given for a Fermi-Dirac neutrino source"
              + " specification on line " + std::to_string(line_num)
              + " of the configuration file " + filename);
        }
        else eta = 0.;

        // Now that we have all of the necessary parameters, create the new
        // Fermi-Dirac neutrino source
        sources.push_back(std::make_unique<marley::FermiDiracNeutrinoSource>(
          neutrino_pid, weight, Emin, Emax, temperature, eta));
      }

      else if (arg == "bf" || arg == "beta" || arg == "beta-fit") {
        double Emin, Emax, Eavg, beta;
        // Process the minimum energy for this source
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emin = std::stod(arg);
          if (Emin < 0.) throw marley::Error(
            std::string("Negative") + " minimum energy " + std::to_string(Emin)
            + " given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        }
        else throw marley::Error(
            std::string("Invalid") + " minimum energy '" + arg
            + "' given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        // Process the maximum energy for this source
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          Emax = std::stod(arg);
          if (Emax < Emin) throw marley::Error(
            std::string("Maximum") + " energy " + std::to_string(Emax)
            + " given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename
            + " is less than the minimum energy " + std::to_string(Emin));
        }
        else throw marley::Error(
            std::string("Invalid") + " maximum energy '" + arg
            + "' given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        // Process the temperature (in MeV) for this source
        next_word_from_line(iss, arg, keyword, line_num, true, false);
        if (std::regex_match(arg, rx_num)) {
          Eavg = std::stod(arg);
          if (Eavg <= 0.) throw marley::Error(
            std::string("Non-positive") + " average energy "
            + std::to_string(Eavg)
            + " given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        }
        else throw marley::Error(
            std::string("Invalid") + " average energy '" + arg
            + "' given for a beta-fit neutrino source"
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
        // Process the fit parameter beta for this source, or
        // set it to 4.5 if it is not given. Any real number is allowed for
        // this parameter.
        if (next_word_from_line(iss, arg, keyword, line_num, false, false)) {
          if (std::regex_match(arg, rx_num)) beta = std::stod(arg);
          else throw marley::Error(
              std::string("Invalid") + " fit parameter '" + arg
              + "' given for a beta-fit neutrino source"
              + " specification on line " + std::to_string(line_num)
              + " of the configuration file " + filename);
        }
        else beta = 4.5;

        // Now that we have all of the necessary parameters, create the new
        // Fermi-Dirac neutrino source
        sources.push_back(std::make_unique<marley::BetaFitNeutrinoSource>(
          neutrino_pid, weight, Emin, Emax, Eavg, beta));
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
          next_word_from_line(iss, arg, keyword, line_num, true, true);
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
              + " specification on line " + std::to_string(line_num)
              + " of the configuration file " + filename);
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
            + " specification on line " + std::to_string(line_num)
            + " of the configuration file " + filename);
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
          bool found_E = next_word_from_line(iss, E_str, keyword, line_num,
            false, false);
          if (found_E) {
            if (std::regex_match(E_str, rx_num)) {
              double E = std::stod(E_str);
              if (E < 0.) throw marley::Error(
                std::string("Negative") + " energy "
                + std::to_string(E) + " given in a "
                + description + " neutrino source specification on line "
                + std::to_string(line_num) + " of the configuration file "
                + filename);
              if (E <= previous_E) throw marley::Error(
                std::string("Energy") + " values "
                + " given in a " + description
                + " neutrino source specification on line "
                + std::to_string(line_num) + " of the configuration file "
                + filename + " are not strictly increasing");
              energies.push_back(E);
              previous_E = E;
            }
            else throw marley::Error(std::string("Invalid") + " energy '"
              + E_str + "' given in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num) + " of the configuration file "
              + filename);
          }
          // Get the first probability density
          bool found_PD = next_word_from_line(iss, PD_str, keyword, line_num,
            false, false);
          if (found_PD) {
            if (std::regex_match(PD_str, rx_num)) {
              double pd = std::stod(PD_str);
              if (pd < 0.) throw marley::Error(
                std::string("Negative") + " probability density "
                + std::to_string(pd) + " given in a "
                + description + " neutrino source specification on line "
                + std::to_string(line_num) + " of the configuration file "
                + filename);
              prob_densities.push_back(pd);
            }
            else throw marley::Error(std::string("Invalid")
              + " probability density '"
              + PD_str + "' given in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num) + " of the configuration file "
              + filename);
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
              + std::to_string(line_num) + " of the configuration file "
              + filename);
          }
          else if (!found_E) {
            throw marley::Error(std::string("Missing ")
              + "energy entry in a " + description
              + " neutrino source specification on line "
              + std::to_string(line_num) + " of the configuration file "
              + filename);
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
          // that no energies outside of the histogram will be sampled and makes
          // the energy and probability densities vector lengths equal).
          prob_densities.push_back(0.);
        }

        // Now that we've processed grid points, create the grid neutrino
        // source
        sources.push_back(std::make_unique<marley::GridNeutrinoSource>(
          energies, prob_densities, neutrino_pid, weight, method));
      }

      else throw marley::Error(std::string("Unrecognized")
      + " neutrino source type '" + arg
      + "' given on line " + std::to_string(line_num)
      + " of the configuration file " + filename);
    }

    else if (keyword == "events") {
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      if (!std::regex_match(arg, rx_num))
        throw marley::Error(std::string("Invalid")
        + " number of events '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename);
      double d_n_events = std::stod(arg);
      int n_events = static_cast<int>(d_n_events);
      // Allow n_events == 0 for testing purposes
      if (n_events < 0)
        throw marley::Error(std::string("Number")
        + " of events '" + arg + "' given on line "
        + std::to_string(line_num) + " of the configuration file "
        + filename + " must be non-negative.");

      num_events = static_cast<size_t>(n_events);
    }

/*** Continuum binning and multiple thread options removed 02/01/2016.
 * Neither of these features is fully implemented yet, so they're not
 * necessary to have as part of the configuration file format.
    else if (keyword == "contbin") {
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      if (!std::regex_match(arg, rx_num))
        throw marley::Error(std::string("Non-numeric")
        + " continuum bin width '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename);
      double width = std::stod(arg);
      if (width <= 0.)
        throw marley::Error(std::string("Continuum")
        + " bin width '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename
        + " must be positive.");
      contbin_width = width;
    }
    else if (keyword == "contbinsubs") {
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      if (!std::regex_match(arg, rx_nonneg_int))
        throw marley::Error(std::string("Invalid")
        + " number of continuum bin subintervals '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename);
      int subs = std::stoi(arg);
      if (subs <= 0)
        throw marley::Error(std::string("Number")
        + " of continuum bin subintervals '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename
        + " must be positive.");
      contbin_num_subs = subs;
    }
    else if (keyword == "threads") {
      next_word_from_line(iss, arg, keyword, line_num, true, false);
      if (!std::regex_match(arg, rx_nonneg_int))
        throw marley::Error(std::string("Invalid")
        + " number of parallel threads '" + arg
        + "' given on line " + std::to_string(line_num)
        + " of the configuration file " + filename);
      int threads = std::stoi(arg);
      if (threads <= 0)
        throw marley::Error(std::string("Number")
        + " of threads '" + arg + "' given on line "
        + std::to_string(line_num) + " of the configuration file "
        + filename + " must be positive.");

      num_threads = static_cast<size_t>(threads);

      std::string force;
      bool forced = next_word_from_line(iss, force, keyword, line_num, false,
        true);
      if (forced && force != "force") throw marley::Error(
        std::string("Trailing") + " argument that is not 'force'"
        + " given for keyword " + keyword + " on line "
        + std::to_string(line_num) + " of the configuration file "
        + filename);

      size_t num_allowed_threads = std::thread::hardware_concurrency();

      // If the requested number of threads exceeds the recommendation for this
      // machine, then complain. Allow the user to override this precautionary
      // message using the word "force" on this line of the input file.
      if (num_threads > 1 && num_threads > num_allowed_threads && !forced)
        throw marley::Error(std::string("Number")
          + " of threads '" + arg + "' given on line "
          + std::to_string(line_num) + " of the configuration file "
          + filename + " exceeds the number of concurrent threads supported"
          + " by this machine (" + std::to_string(num_allowed_threads)
          + ") according to std::thread::hardware_concurrency()."
          + " To force usage of " + arg + " threads, add the word "
          + "'force' (without quotes) at the end of line "
          + std::to_string(line_num) + " of the configuration file "
          + filename);
    }
*******/
    else {
      std::cerr << "Warning: Ignoring unrecognized keyword '"
        << keyword << "' on line " << line_num
        << " of the configuration file "
        << filename << std::endl;
    }
  }

  // Check that required keywords were found
  if (reaction_filenames.empty())
    throw marley::Error(std::string("Configuration file")
    + " must contain at least one use of the reaction keyword.");
  if (sources.empty())
    throw marley::Error(std::string("Configuration file")
    + " must contain at least one use of the source keyword.");
}

// Get the next word from a parsed line. If errors occur, complain.
// The last argument determines whether the next word should be
// converted to all lowercase or left as is.
bool marley::ConfigFile::next_word_from_line(std::istringstream& iss,
  std::string& word, const std::string& keyword, const int line_number,
  bool enable_exceptions, bool make_lowercase)
{
  if (iss.eof()) {
    if (!enable_exceptions) return false;
    throw marley::Error(std::string("Missing argument ")
      + "for keyword '" + keyword + "' encountered while"
      + " parsing line " + std::to_string(line_number)
      + " of the configuration file " + filename);
  }
  iss >> word;
  if (iss.fail()) {
    if (!enable_exceptions) return false;
    std::string snum = std::to_string(line_number);
    throw marley::Error(std::string("Error")
      + " occurred while parsing arguments for keyword '"
      + keyword + "' on line "
      + std::to_string(line_number)
      + " of the configuration file " + filename);
  }
  if (make_lowercase) marley_utils::to_lowercase_inplace(word);
  return true;
}


void marley::ConfigFile::print_summary(std::ostream& os) {
  os << "seed: " << seed << std::endl;
  #ifdef USE_ROOT
  os << "writeroot: ";
  if (writeroot && check_before_root_file_overwrite)
    os << "yes";
  else if (writeroot) os << "overwrite";
  else os << "no";
  os << std::endl;
  os << "rootfile: " << root_filename << std::endl;
  #endif
  for (const auto& rfile: reaction_filenames)
    os << "reaction: " << rfile << std::endl;
  for (const auto& sr: structure_records) {
    os << "structure:";
    // TODO: consider using an ordered set that keeps the
    // nucids in periodic table order. This will make the
    // output here look nicer
    for (const auto& id: sr.nucids) {
      std::string trimmed_id = marley_utils::trim_copy(id);
      // If the element symbol is not a single character,
      // then make its second letter lowercase. Also change
      // the ENSDF code for a neutron ("NN") to "n"
      std::string element_symbol = trimmed_id.substr(trimmed_id.size() - 2);
      if (!std::regex_match(element_symbol.substr(0,1), rx_nonneg_int)) {
        element_symbol.back() = tolower(element_symbol.back());
        if (element_symbol == "Nn") element_symbol = "n";
        trimmed_id = trimmed_id.erase(trimmed_id.size() - 2) + element_symbol;
      }
      os << " " << trimmed_id;
    }
    os << " from " << sr.filename << " ("
      << marley::ConfigFile::format_to_string(sr.format)
      << " format)" << std::endl;
  }
}

marley::DecayScheme::FileFormat
  marley::ConfigFile::string_to_format(
  const std::string& string)
{
  if (string == "ensdf")
    return marley::DecayScheme::FileFormat::ensdf;
  else if (string == "talys")
    return marley::DecayScheme::FileFormat::talys;
  else throw marley::Error(std::string("Unrecognized file")
    + " format '" + string + "' passed to"
    + " marley::ConfigFile::string_to_format()");
}

std::string marley::ConfigFile::format_to_string(
  const marley::DecayScheme::FileFormat ff)
{
  if (ff == marley::DecayScheme::FileFormat::ensdf)
    return std::string("ensdf");
  else if (ff == marley::DecayScheme::FileFormat::talys)
    return std::string("talys");
  else throw marley::Error(std::string("Unrecognized")
    + " marley::DecayScheme::FileFormat value passed to"
    + " marley::ConfigFile::format_to_string()");
}
