#include "marley_utils.hh"
#include "TMarleyConfigFile.hh"

// Matches comment lines and empty lines
const std::regex TMarleyConfigFile::rx_comment_or_empty = std::regex("#.*|\\s*");
// Matches positive integers
const std::regex TMarleyConfigFile::rx_num = std::regex("[0-9]+");
// Matches trimmed ENSDF-style nucids
const std::regex TMarleyConfigFile::rx_nucid = std::regex("[1-9][0-9]{0,2}[A-Za-z]{1,2}");

// The default constructor assigns default values to all
// configuration file parameters
TMarleyConfigFile::TMarleyConfigFile() {
  seed = std::chrono::system_clock::now().time_since_epoch().count();
  #ifdef USE_ROOT
  writeroot = false;
  check_before_root_file_overwrite = true;
  root_filename = "events.root";
  #endif
}

// Call the default constructor first to set all configuration
// file parameters to their default values
TMarleyConfigFile::TMarleyConfigFile(std::string file_name)
  : TMarleyConfigFile()
{

  filename = file_name;
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw std::runtime_error(std::string("Could not read from the ") +
      "file " + filename);
  }

  std::string line; // String to store the current line
                    // of the configuration file during parsing

  std::string keyword; // String to store the current keyword read in
                       // from the configuration file

  std::istringstream iss; // Stream used to help parse each line

  // Loop through each line in the file. Search for keywords,
  // and react appropriately whenever they are encountered. 
  for (int line_num = 1;
    line = marley_utils::get_next_line(file_in, rx_comment_or_empty, false),
    line != ""; ++line_num)
  {
    // Convert the line to all lowercase (this is an easy way of
    // getting case-insensitivity)
    marley_utils::to_lowercase_inplace(line);

    // Load the stream with the new line
    iss.clear();
    iss.str(line);

    // Read the first word (the keyword for that line) from the stream
    iss >> keyword;

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
      else if (!std::regex_match(arg, rx_num))
        throw std::runtime_error(std::string("Invalid random number seed")
          + "encountered on line " + std::to_string(line_num)
          + " of the configuration file " + filename);
      else {
        // Convert the seed string to a 64-bit unsigned integer
        seed = static_cast<uint_fast64_t>(std::stoull(arg));
      }
    }
    else if (keyword == "reaction") {
      next_word_from_line(iss, arg, keyword, line_num);
      reaction_filenames.insert(arg);
    }
    else if (keyword == "rootfile") {
      #ifdef USE_ROOT
      next_word_from_line(iss, arg, keyword, line_num);
      root_filename = arg;
      #else
      // TODO: change this to a warning that this will be ignored
      throw std::runtime_error(std::string("Cannot set")
        + " ROOT filename. MARLEY must be compiled"
        + "with ROOT support.");
      #endif
    }
    else if (keyword == "writeroot") {
      #ifdef USE_ROOT
      next_word_from_line(iss, arg, keyword, line_num);
      if (arg == "yes") writeroot = true;
      else if (arg == "no") writeroot = false;
      else if (arg == "overwrite") {
        writeroot = true;
        check_before_root_file_overwrite = false;
      }
      else {
        throw std::runtime_error(std::string("Invalid")
          + " ROOT file write flag '" + arg
          + "' encountered on line" + std::to_string(line_num)
          + " of the configuration file " + filename);
      }
      #else
      throw std::runtime_error(std::string("Cannot write")
        + " events to a ROOT file. MARLEY must be compiled"
        + "with ROOT support.");
      #endif
    }
    else if (keyword == "structure") {
      StructureRecord sr;
      next_word_from_line(iss, sr.filename, keyword, line_num);
      std::string format_string;
      next_word_from_line(iss, format_string, keyword, line_num);
      try {
        sr.format = string_to_format(format_string);
      }
      // If the format is invalid, catch the exception thrown
      // by the converter and throw a different one so that
      // the error message is more easily understood
      catch (const std::runtime_error& e) {
         throw std::runtime_error(std::string("Unknown")
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
        if (std::regex_match(arg, rx_num)) {
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
            std::regex_search(arg, m, rx_num);
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
        else throw std::runtime_error(std::string("Invalid")
          + " nuclide specifier '" + arg
          + "' given on line " + std::to_string(line_num)
          + " of the configuration file " + filename);

      } while (next_word_from_line(iss, arg,
        keyword, line_num, false));

      // Add this structure record to the master vector
      structure_records.push_back(sr);

    }
    else {
      std::cerr << "Warning: Ignoring unrecognized keyword '"
        << keyword << "' on line " << line_num
        << " of the configuration file "
        << filename << std::endl;
    }
  }

  // Check that required keywords were found
  if (reaction_filenames.empty())
    throw std::runtime_error(std::string("Configuration file")
    + " must contain at least one use of the reaction keyword.");
}

// Get the next word from a parsed line. If errors occur, complain.
bool TMarleyConfigFile::next_word_from_line(std::istringstream& iss, std::string& word,
  const std::string& keyword, const int line_number,
  bool enable_exceptions)
{
  if (iss.eof()) {
    if (!enable_exceptions) return false;
    throw std::runtime_error(std::string("Missing argument ")
      + "for keyword '" + keyword + "' encountered while"
      + " parsing line " + std::to_string(line_number)
      + " of the configuration file " + filename);
  }
  iss >> word;
  if (iss.fail()) {
    if (!enable_exceptions) return false;
    std::string snum = std::to_string(line_number);
    throw std::runtime_error(std::string("Error")
      + " occurred while parsing arguments for keyword '"
      + keyword + "' on line "
      + std::to_string(line_number)
      + " of the configuration file " + filename);
  }
  return true;
}


void TMarleyConfigFile::print_summary(std::ostream& os) {
  os << "seed: " << seed << std::endl;
  os << "writeroot: ";
  if (writeroot && check_before_root_file_overwrite)
    os << "yes";
  else if (writeroot) os << "overwrite";
  else os << "no";
  os << std::endl;
  os << "rootfile: " << root_filename << std::endl;
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
      if (!std::regex_match(element_symbol.substr(0,1), rx_num)) {
        element_symbol.back() = tolower(element_symbol.back());
        if (element_symbol == "Nn") element_symbol = "n";
        trimmed_id = trimmed_id.erase(trimmed_id.size() - 2) + element_symbol;
      }
      os << " " << trimmed_id;
    }
    os << " from " << sr.filename << " ("
      << TMarleyConfigFile::format_to_string(sr.format)
      << " format)" << std::endl;
  }
}

TMarleyDecayScheme::FileFormat
  TMarleyConfigFile::string_to_format(
  const std::string& string)
{
  if (string == "ensdf")
    return TMarleyDecayScheme::FileFormat::ensdf;
  else if (string == "talys")
    return TMarleyDecayScheme::FileFormat::talys;
  else throw std::runtime_error(std::string("Unrecognized file")
    + " format '" + string + "' passed to"
    + " TMarleyConfigFile::string_to_format()");
}

std::string TMarleyConfigFile::format_to_string(
  const TMarleyDecayScheme::FileFormat ff)
{
  if (ff == TMarleyDecayScheme::FileFormat::ensdf)
    return std::string("ensdf");
  else if (ff == TMarleyDecayScheme::FileFormat::talys)
    return std::string("talys");
  else throw std::runtime_error(std::string("Unrecognized")
    + " TMarleyDecayScheme::FileFormat value passed to"
    + " TMarleyConfigFile::format_to_string()");
}
