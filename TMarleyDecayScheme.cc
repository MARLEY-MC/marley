#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "marley_utils.hh"
#include "TMarleyDecayScheme.hh"

// ***** Constants used for ENSDF file parsing *****
//const std::regex TMarleyDecayScheme::ensdf_generic_nuc_id("^[[:alnum:] ]{5}");
const std::string TMarleyDecayScheme::ensdf_primary_record = "[ 1]";
const std::string TMarleyDecayScheme::ensdf_continuation_record = "[^ 1]";
const std::string TMarleyDecayScheme::ensdf_record_data = ".{71}";

const std::regex TMarleyDecayScheme::rx_ensdf_end_record("\\s*"); // Matches blank lines
//const std::regex TMarleyDecayScheme::rx_generic_primary_identification_record(
//  ensdf_generic_nuc_id + ensdf_primary_record + "   " + ensdf_record_data);

// ***** Constants used for TALYS file parsing *****
// Each line describing a level has fortran
// format (i4, f11.6, f6.1, 3x, i2, i3, 19x, e9.3, 1x, 2a1, a18)
//const std::regex TMarleyDecayScheme::rx_talys_level_line(std::string("[0-9 ]{8}\\.[0-9]{6}")
//  + "(?:[0-9 ]{4}\\.[0-9]| {6}) {3}(?:[0-9 ]{2}| {2})[0-9 ]{3} {19}(?:[0-9]\\."
//  + "[0-9]{3}[Ee][+-][0-9]{2}| {9}) .{20}");

// Each line describing a gamma ray transition has
// fortran format (29x, i3, f10.6, e10.3, 5x, a1)
//const std::regex TMarleyDecayScheme::rx_talys_gamma_line(std::string(" {29}[0-9 ]{3}")
//  + "[0-9 ]{3}\\.[0-9]{6}[0-9 ]{2}\\.[0-9]{3}"
//  + "[Ee][+-][0-9]{2} {5}.");

//void TMarleyDecayScheme::do_cascade(std::string initial_energy) {
//  std::map<std::string, TMarleyLevel>::iterator it = levels.find(initial_energy);
//  if (it == levels.end()) {
//    throw std::range_error("Could not do cascade. Level with energy "
//      + initial_energy + " MeV not found.");
//  }
//  else {
//    TMarleyLevel* plevel = &(it->second);
//    do_cascade(plevel);
//  }
//}

// Returns a pointer to the level owned by this decay scheme object
// that has the closest excitation energy to E_level (E_level
// has units of MeV).
TMarleyLevel* TMarleyDecayScheme::get_pointer_to_closest_level(double E_level) {
  // Search for the level whose energy is closest to the given value of E_level
  std::vector<double>::iterator it = std::lower_bound(
    sorted_level_energies.begin(), sorted_level_energies.end(),
    E_level);

  // Determine the index of the level energy appropriately
  unsigned int e_index = std::distance(sorted_level_energies.begin(), it);
  if (e_index == sorted_level_energies.size()) {
    // The given energy is greater than every level energy in our
    // decay scheme. We will therefore assume that the initial level is the
    // highest level.  Its energy's index is given by one less than the
    // number of elements in the sorted vector, so subtract one from our
    // previous result.
    --e_index;
  }
  else if (e_index > 0) {
    // If the calculated index does not correspond to the
    // first element, we still need to check which of the
    // two levels found (one on each side) is really the
    // closest. Do so and reassign the index if needed.
    if (std::abs(E_level - sorted_level_energies[e_index])
      > std::abs(E_level - sorted_level_energies[e_index - 1]))
    {
      --e_index;
    }
  }

  // Get a pointer to the selected level object
  TMarleyLevel* plevel = pv_sorted_levels[e_index];

  return plevel;

}

//void TMarleyDecayScheme::do_cascade(double initial_energy) {
//
//  // Since we were given a numerical initial energy in this
//  // version of the function, search for the level whose energy is
//  // closest to the given value of initial_energy. Then handle
//  // the cascade in the usual way.
//  TMarleyLevel* plevel = get_pointer_to_closest_level(initial_energy);
//  do_cascade(plevel);
//}

void TMarleyDecayScheme::do_cascade(TMarleyLevel* initial_level,
  TMarleyEvent* p_event)
{
  //std::cout << "Beginning gamma cascade at level with energy "
  //  << initial_level->get_string_energy() << " keV" << std::endl;

  bool cascade_finished = false;

  TMarleyLevel* p_current_level = initial_level;

  while (!cascade_finished) {
    // Randomly select a gamma to produce
    TMarleyGamma* p_gamma = p_current_level->sample_gamma();
    if (p_gamma == nullptr) {
      //std::cout << "  this level does not have any gammas" << std::endl;
      cascade_finished = true;
    }
    else {
      p_current_level = p_gamma->get_end_level();
      if (p_current_level == nullptr) {
        throw std::runtime_error(std::string("This gamma does not have an end level. ")
          + "Cannot continue cascade.");
      }
      //std::cout << "  emitted gamma with energy " << p_gamma->get_energy()
      //  << " MeV." New level has energy " << p_current_level->get_string_energy()
      //  << " keV." << std::endl;
      //std::cout.precision(15);
      //std::cout << std::scientific;
      //std::cout << "gamma energy = " << p_gamma->get_energy() << std::endl;
      
      // Add the gamma to the event object as a final particle. Note that
      // photons have a particle id number of 22
      double gamma_energy = p_gamma->get_energy();

      // Sample a direction assuming that the gammas are emitted
      // isotropically in the nucleus's rest frame. Also don't bother
      // to do a Lorentz boost to the lab frame (the nucleus is moving
      // nonrelativistically, and it will likely collide with other
      // nuclei before emitting some of the gammas anyway)
      // sample from [-1,1]
      double gamma_cos_theta = marley_utils::uniform_random_double(-1, 1, true);
      double gamma_theta = std::acos(gamma_cos_theta);

      // sample from [0,2*pi)
      double gamma_phi = marley_utils::uniform_random_double(0, 2*marley_utils::pi, false);

      // Compute 3-momentum components using the sampled angles
      double gamma_px = std::sin(gamma_theta)*std::cos(gamma_phi)*gamma_energy;
      double gamma_py = std::sin(gamma_theta)*std::sin(gamma_phi)*gamma_energy;
      double gamma_pz = gamma_cos_theta*gamma_energy;
      
      // Add the new gamma to the event
      p_event->add_final_particle(TMarleyParticle(22, gamma_energy, gamma_px,
        gamma_py, gamma_pz, 0));

      // Add this gamma to the residue's children
      TMarleyParticle* p_residue = p_event->get_residue();
      p_residue->add_child(&(p_event->get_final_particles()->back()));
    }
  }

  //std::cout << "Finished gamma cascade at level with energy "
  //  << p_current_level->get_string_energy() << std::endl;
}


TMarleyDecayScheme::TMarleyDecayScheme(std::string nucid, std::string filename,
  TMarleyDecayScheme::FileFormat ff)
{

  this->nuc_id = nucid;
  this->file_name = filename;

  // Get the atomic mass number from the first 3 characters
  // of the ENSDF nucid
  this->atomic_mass_number = std::stoi(nucid.substr(0,3));

  // Parse the data file using the appropriate format
  switch (ff) {

    case FileFormat::ensdf:
      parse_ensdf(); 
      break;

    case FileFormat::talys:
      parse_talys();
      break;

    default:
      throw std::runtime_error(std::string("Invalid file format ")
        + " supplied to TMarleyDecayScheme constructor.");
  }

}

void TMarleyDecayScheme::parse_ensdf() {

  // Regular expressions for identifying ensdf record types
  const std::regex rx_primary_identification_record(nuc_id + ensdf_primary_record
    + "   ADOPTED LEVELS.{57}"); //"   ADOPTED LEVELS(, GAMMAS.{49}|.{57})");
  const std::regex rx_primary_level_record(nuc_id + ensdf_primary_record + " L " + ensdf_record_data);
  const std::regex rx_continuation_level_record(nuc_id + ensdf_continuation_record + " L " + ensdf_record_data);
  const std::regex rx_primary_gamma_record(nuc_id + ensdf_primary_record + " G " + ensdf_record_data);
  const std::regex rx_continuation_gamma_record(nuc_id + ensdf_continuation_record + " G " + ensdf_record_data);
 
  // Open the ENSDF file for parsing
  std::ifstream file_in(file_name);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw std::runtime_error(std::string("Could not read from the ") +
      "ENSDF data file " + file_name);
  }

  std::string line; // String to store the current line
                    // of the ENSDF file during parsing

  std::string record; // String that will store the full text
                      // (all lines) of the current ENSDF record

  bool found_decay_scheme = false; // Flag that indicates whether or not
                                   // gamma decay scheme data were found for the
                                   // current nuc_id.

  //std::cout << "DEBUG: Searching for nucid '" << nuc_id << "'" << std::endl;
  while (!file_in.eof()) {
    std::getline(file_in, line);
    if (std::regex_match(line, rx_primary_identification_record)) {
      found_decay_scheme = true;
      //std::cout << "DEBUG: Found nucid '" << nuc_id << "'" << std::endl;
      break;
    }
  }

  if (!found_decay_scheme) {
    throw std::runtime_error(std::string("Gamma decay scheme data ")
      + "(adopted levels, gammas) for " + marley_utils::nucid_to_symbol(nuc_id)
      + " could not be found in the ENSDF data file " + file_name);
  }
  //else {
  //  std::cout << "Gamma decay scheme data for " + nuc_id << " found. Using ENSDF dataset" << std::endl;
  //  std::cout << line << std::endl;
  //}

  bool no_advance = false; // Flag that prevents advancing through
                           // the ENSDF file when a continuation record
                           // is found

  TMarleyLevel* p_current_level = nullptr; // Pointer to the current level object
                                          // being filled with gamma ray data 

  //std::cout << "DEBUG: Parsing data for nucid '" << nuc_id << "'" << std::endl;
  while (!file_in.eof()) {
    // Get the next line of the file
    // unless we already did
    if (!no_advance) std::getline(file_in, line);
    no_advance = false;

    // Level Record
    if (std::regex_match(line, rx_primary_level_record)) {
      //std::cout << "DEBUG:   Parsing level" << std::endl;
      record = line;
      line = this->process_continuation_records(file_in, record, rx_continuation_level_record);
      no_advance = true;
        
      // Extract the level energy (in keV) as a trimmed string from the ENSDF level record
      std::string level_energy = marley_utils::trim_copy(record.substr(9,10)); 

      // Also extract the spin-parity of the level
      std::string spin_parity = marley_utils::trim_copy(record.substr(21,18));

      // Add a new level object to this decay scheme object using these data
      this->add_level(TMarleyLevel(level_energy, spin_parity));

      // Update the current level pointer
      p_current_level = this->get_level(level_energy);
    }

    // Gamma Record
    else if (std::regex_match(line, rx_primary_gamma_record)) {
      //std::cout << "DEBUG:   Parsing gamma" << std::endl;
      record = line;
      line = this->process_continuation_records(file_in,
        record, rx_continuation_gamma_record);
      no_advance = true;

      // Extract the gamma ray's energy and relative (photon)
      // intensity from the ENSDF gamma record. Convert its
      // energy from keV to MeV.
      double gamma_energy = std::stod(record.substr(9,10)) * marley_utils::MeV;
      double gamma_ri = marley_utils::str_to_double(record.substr(21,8));

      // If this gamma belongs to a level record, then add its
      // data to the corresponding level object. Gammas that
      // have no assigned level appear before any level records,
      // so p_current_level will be a null pointer for them.
      if (p_current_level != nullptr) {
        p_current_level->add_gamma(TMarleyGamma(gamma_energy, gamma_ri, p_current_level)); 
      }

    }

    else if (std::regex_match(line, rx_ensdf_end_record)) {
      //std::cout << "DEBUG:   Found end of data" << std::endl;
      break;
    }

  }

  // Now that we've loaded our decay scheme with all of the data,
  // we can fill in the end level pointers for each of our gamma
  // ray objects. This will enable the generator to descend through
  // the decay chains after only looking up the starting level.

  // We will use the sorted_level_energies vector (whose entries were created
  // when we parsed each level record) to find the closest final level's index.

  // Cycle through each of the level objects. We will assign end level pointers
  // to each gamma owned by each level.
  for(std::vector<TMarleyLevel*>::iterator j = this->pv_sorted_levels.begin();
    j != this->pv_sorted_levels.end(); ++j)
  {

    // Calculate Weisskopf Estimates here if gamma data is not known
    if( (*j)->get_gamma_status() == false ) {
      //std::cout << "DEBUG: Computing WE" << std::endl;
      do_weisskopf(j - this->pv_sorted_levels.begin());
    }

    std::vector<TMarleyGamma>* p_gammas = (*j)->get_gammas();
    double initial_level_energy = (*j)->get_numerical_energy();

    for(std::vector<TMarleyGamma>::iterator k = p_gammas->begin();
      k != p_gammas->end(); ++k) 
    {
      double gamma_energy = k->get_energy();

      // Approximate the final level energy so we can search for the final level
      double final_level_energy = initial_level_energy - gamma_energy; 

      // Search for the corresponding energy using the vector of level energies
      std::vector<double>::iterator p_final_level_energy = std::lower_bound(
        this->sorted_level_energies.begin(), this->sorted_level_energies.end(),
        final_level_energy); 

      // Determine the index of the final level energy appropriately
      unsigned int e_index = std::distance(this->sorted_level_energies.begin(),
        p_final_level_energy);
      if (e_index == this->sorted_level_energies.size()) {
        // The calculated final level energy is greater than
        // every energy in the vector. We will therefore assume
        // that the gamma decay takes us to the highest level in
        // the vector. Its index is given by one less than the
        // number of elements in the vector, so subtract one from
        // our previous result.
        --e_index;
      }
      else if (e_index > 0) {
        // If the calculated index does not correspond to the
        // first element, we still need to check which of the
        // two levels found (one on each side) is really the
        // closest. Do so and reassign the index if needed.
        if (std::abs(final_level_energy - this->sorted_level_energies[e_index])
          > std::abs(final_level_energy - this->sorted_level_energies[e_index - 1]))
        {
          --e_index;
        }
      }

      // Use the index to assign the appropriate end level pointer to this gamma
      k->set_end_level(this->pv_sorted_levels[e_index]);

    }

  }

  file_in.close();

}

void TMarleyDecayScheme::parse_talys() {
  // First line in a TALYS level dataset has fortran
  // format (2i4, 2i5, 56x, i4, a2)
  // General regex for this line:
  // std::regex nuclide_line("[0-9 ]{18} {56}[0-9 ]{4}.{2}");
  const std::regex nuclide_line("[0-9 ]{18} {57}" + nuc_id);

  // Open the TALYs level data file for parsing
  std::ifstream file_in(file_name);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw std::runtime_error(std::string("Could not read from the ") +
      "TALYS data file " + file_name);
  }

  std::string line; // String to store the current line
                    // of the TALYS file during parsing

  bool found_decay_scheme = false;

  while (!file_in.eof()) {
    std::getline(file_in, line);
    if (std::regex_match(line, nuclide_line)) {
      found_decay_scheme = true; 
      break;
    }
  }

  if (!found_decay_scheme) {
    throw std::runtime_error(std::string("Gamma decay scheme data ")
      + "(adopted levels, gammas) for " + marley_utils::nucid_to_symbol(nuc_id)
      + " could not be found in the TALYS data file " + file_name);
  }
  //else {
  //  std::cout << "Gamma decay scheme data for " + nuc_id << " found. Using TALYS dataset" << std::endl;
  //  std::cout << line << std::endl;
  //}

  // Dummy integer and number of excited levels for this nuclide
  int dummy, num_excited_levels; 

  // Read in the number of excited levels from the first line of data
  std::istringstream iss(line);
  iss >> dummy >> dummy >> dummy >> num_excited_levels;

  // Pointer to the current level object being filled with gamma ray data
  TMarleyLevel* p_current_level = nullptr; 

  // Temporary vectors to store level energies and pointers to newly
  // created level objects. These are used to calculate gamma ray energies
  // and assign level pointers
  std::vector<double> level_energies; std::vector<TMarleyLevel*> level_ps;

  for (int i = 0; i <= num_excited_levels; i++) {

    // Get the next line of the file. This will be a discrete level record
    std::getline(file_in, line);

    // Load the new line into our istringstream object for parsing. Reset
    // the stream so that we start parsing from the beginning of the string.
    iss.str(line);
    iss.clear();

    // Read in this level's index, energy, spin, parity, and
    // number of gamma transitions
    int level_num, parity, num_gammas;
    double level_energy, spin;
    iss >> level_num >> level_energy >> spin >> parity >> num_gammas;

    // Convert the level energy and spin parity to suitable strings.
    // These are needed for the current version of the TMarleyLevel constructor.
    // If you encounter a nonsensical parity value, complain.
    // When creating the level energy string, convert to keV, then to a string.
    // This will keep the string consistent with the conventions we use
    // for parsing ENSDF format data
    std::string string_level_energy = std::to_string(level_energy/marley_utils::MeV);
    std::string string_spin_parity;

    // Get the fractional part of the spin. If it is zero, construct
    // the spin string using an integer spin
    double ipart;
    if (std::modf(spin, &ipart) == 0.0) {
      string_spin_parity = std::to_string(static_cast<int>(ipart));
    }

    // If the fractional part of the spin is nonzero, assume that
    // the spin is half-integer, and construct the spin string
    // appropriately
    else {
      string_spin_parity = std::to_string(static_cast<int>(ipart));
      string_spin_parity += "/2";
    }

    // Add the sign of the parity to the spin-parity string
    if (parity == 1) {
      string_spin_parity += "+";
    }
    else if (parity == -1) {
      string_spin_parity += "-";
    }
    else throw std::runtime_error(std::string("Invalid parity value ")
      + std::to_string(parity) + " encountered in TALYS dataset "
      + file_name); 

    // Construct a new level object and add it to the decay scheme
    this->add_level(TMarleyLevel(string_level_energy, string_spin_parity));

    // Get a pointer to the newly-added level
    p_current_level = this->get_level(string_level_energy);

    // Add the numerical energy for this level and a pointer
    // to the newly-created level object to our temporary arrays.
    // These will be used when creating gammas for each level.
    level_energies.push_back(level_energy);
    level_ps.push_back(p_current_level);

    for (int j = 0; j < num_gammas; j++) {

      // Get the next line of the file. This will be a gamma record
      std::getline(file_in, line);

      // Load the new line into our istringstream object for parsing. Reset
      // the stream so that we start parsing from the beginning of the string.
      iss.str(line);
      iss.clear();

      // Read in the index of the final level and branching ratio
      // for this gamma transition
      int gamma_final_level_num;
      double br;
      iss >> gamma_final_level_num >> br;

      // Process this gamma if it has a nonvanishing
      // branching ratio and is associated with a discrete level
      if (br > 0 && p_current_level != nullptr) {

        // Compute the gamma ray's energy in MeV by subtracting the energy
        // of the final level from the energy of the initial level
        double gamma_energy = level_energy - level_energies.at(gamma_final_level_num);

        // Create a gamma ray object and add it to the current level object
        p_current_level->add_gamma(TMarleyGamma(gamma_energy, br, p_current_level));

        // Get a pointer to the newly-created gamma ray object
        TMarleyGamma* p_current_gamma = &(p_current_level->get_gammas()->back());

        // Set the gamma ray's final level pointer to point to
        // the appropriate level object
        p_current_gamma->set_end_level(level_ps.at(gamma_final_level_num));
      }
    }
  } 

  file_in.close();
}

std::string TMarleyDecayScheme::process_continuation_records(std::ifstream &file_in,
  std::string &record, const std::regex &rx_cont_record) const {

  std::string line;
  while (!file_in.eof()) {
    // Get the next line of the file
    std::getline(file_in, line);
 
    // Check to see if the next line is a continuation record 
    if (std::regex_match(line, rx_cont_record)) {
      // If it is, add the next line to the current
      // record text.
      record += "\n" + line;
    }

    // If a non-continuation-record line is found,
    // stop reading from the file.
    else return line; 
  
  }

  // If the end of the file is encountered while processing
  // the continuation records, return an empty string
  return std::string("");
}

void TMarleyDecayScheme::print_report(std::ostream& ostr) {
  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for(std::vector<TMarleyLevel*>::iterator j = this->pv_sorted_levels.begin();
    j != this->pv_sorted_levels.end(); ++j)
  {

    //std::string sp = (*j)->get_spin_parity();
    int spin = (*j)->get_ispin();
    int parity = (*j)->get_iparity();
    
    //if (sp.empty()) sp = "UNKNOWN";
    
    ostr << "Level at " << (*j)->get_string_energy()
	 << " keV has spin-parity " << spin << " " << parity << std::endl; // changed from sp to spin and parity
    std::vector<TMarleyGamma>* p_gammas = (*j)->get_gammas();

    // Cycle through each of the gammas owned by the current level
    // (according to the ENSDF specification, these will already be
    // sorted in order of increasing energy)
    for(std::vector<TMarleyGamma>::iterator k = p_gammas->begin();
      k != p_gammas->end(); ++k) 
    {
      ostr << "  has a gamma with energy " << k->get_energy() << " MeV";
      ostr << " (transition to level at "
        << k->get_end_level()->get_string_energy() << " keV)" << std::endl;
      ostr << "    and relative photon intensity " << k->get_ri() << std::endl; 
    }

  }

}

void TMarleyDecayScheme::print_latex_table(std::ostream& ostr) {

  std::string caption_beginning =
    std::string("{\\textbf{Levels") +
    " and $\\boldsymbol{\\gamma}$ transitions \n for " +
    "\\isotope[\\boldsymbol{" +
    marley_utils::trim_copy(this->nuc_id.substr(0,3)) +
    "}]{\\textbf{" + this->nuc_id.substr(3,1) +
    marley_utils::trim_copy(
      marley_utils::to_lowercase(nuc_id.substr(4,1))) + 
    "}} \n from file " + this->file_name + " ";

  ostr << marley_utils::latex_table_1;

  ostr << caption_beginning + "}}\\\\\n";

  ostr << marley_utils::latex_table_2;

  ostr << caption_beginning + " -- \\textit{continued}}} \\\\\n";  

  ostr << marley_utils::latex_table_3;
  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for(std::vector<TMarleyLevel*>::iterator j = this->pv_sorted_levels.begin();
    j != this->pv_sorted_levels.end(); ++j)
  {

    std::string sp = (*j)->get_spin_parity();
    if (sp.empty()) sp = "?";
    if (sp == "UNNATURAL") sp = "unnat.";

    ostr << (*j)->get_string_energy() << " & "
      << sp  << " & ";

    std::vector<TMarleyGamma>* p_gammas = (*j)->get_gammas();

    // If there aren't any gammas for this level, finish writing
    // the current row of the table. Add extra space between this
    // level and the next one.
    if (p_gammas->empty()) {
      ostr << " &  &";
      // If this is the last row of the table, don't add extra space.
      if (j == this->pv_sorted_levels.end() - 1) {
        ostr << "" << std::endl;
      }
      else {
        ostr << " \\\\ \\addlinespace[\\ExtraRowSpace]" << std::endl;
      }
    }

    // Cycle through each of the gammas owned by the current level
    // (according to the ENSDF specification, these will already be
    // sorted in order of increasing energy)
    for(std::vector<TMarleyGamma>::iterator k = p_gammas->begin();
      k != p_gammas->end(); ++k) 
    {
      // If this is not the first gamma, add empty columns
      // for the level energy and spin-parity
      if (k != p_gammas->begin()) ostr << " & & ";
      // Output information about the current gamma
      ostr << k->get_energy() << " & "
        << k->get_ri() << " & "
        << k->get_end_level()->get_string_energy();
      // Add vertical space after the final gamma row. Also prevent page breaks
      // in the middle of a list of gammas by outputting a star at the end of
      // each row except the final gamma row.
      if (k == p_gammas->end() - 1) {
	// Don't add the extra row space for the very last row in the table
	if (j >= this->pv_sorted_levels.end() - 1) {
          ostr << std::endl;
        }
        else {
          ostr << " \\\\ \\addlinespace[\\ExtraRowSpace]" << std::endl;
        }
      }
      else {
        ostr << " \\\\*" << std::endl;
      }
    }

  }

  ostr << marley_utils::latex_table_4 << std::endl;
}


std::string TMarleyDecayScheme::get_nuc_id() const {
  return nuc_id;
}

void TMarleyDecayScheme::set_nuc_id(std::string id) {
  nuc_id = id;
}

bool TMarleyDecayScheme::compare_level_energies(TMarleyLevel* first,
  TMarleyLevel* second)
{
  return first->get_numerical_energy() < second->get_numerical_energy();
}

void TMarleyDecayScheme::add_level(TMarleyLevel level) {
  // Add the level to the std::map of level objects. Use the
  // string version of its energy in keV as the key.
  std::string energy_string = level.get_string_energy();
  levels[energy_string] = level;

  // Get a pointer to the just-added level object
  TMarleyLevel* p_level = &(levels[energy_string]);

  // Get this level's numerical energy
  double l_energy = p_level->get_numerical_energy();

  // Figure out where this level should go in the
  // vector of sorted level energies 
  std::vector<double>::iterator
    insert_point = std::lower_bound(sorted_level_energies.begin(),
    sorted_level_energies.end(), l_energy);

  // Compute the numerical index for where we will
  // insert the energy of the new level
  int index = std::distance(sorted_level_energies.begin(), insert_point);

  // Add an entry in the vector of sorted level
  // energies for this new level
  sorted_level_energies.insert(insert_point, l_energy);

  // Also add a pointer to the new level object at the
  // appropriate place in our energy-sorted vector of pointers
  pv_sorted_levels.insert(pv_sorted_levels.begin() + index, p_level);
}


std::vector<TMarleyLevel*>* TMarleyDecayScheme::get_sorted_level_pointers() {
  return &pv_sorted_levels;
}


std::map<std::string, TMarleyLevel>* TMarleyDecayScheme::get_levels() {
  return &levels;
}

TMarleyLevel* TMarleyDecayScheme::get_level(std::string energy) {
  return &(levels.at(energy));
}

void TMarleyDecayScheme::do_weisskopf(int i){
  TMarleyLevel* init = this->pv_sorted_levels[i];
  TMarleyLevel* current = nullptr;
  double deltaE;
  int deltaJ; // Difference in spin
  int parProd; // Product of parities
  double trans_rate;
  double min_rate = 1e+16; // Minimum rate for transitions. Those below this are neglected. Come up with a good way to determine this rate.
  int lambda;
  double BE, BM, f;
  double E_i, E_f;

  for(int j = i-1; j > 0; j--){
      current = this->pv_sorted_levels[j];
      E_i = init->get_numerical_energy();
      E_f = current->get_numerical_energy();
	
      deltaE = E_i - E_f;
	    
      //Determine which type of EM transition takes place
      deltaJ = std::abs(init->get_ispin() - current->get_ispin());
      parProd = init->get_iparity() * current->get_iparity();

      if(parProd == 1) //Meaning parity stays the same
	switch (deltaJ){
	  case 0:
	  case 1:
	    lambda = 0;
	    BM = calcBM(lambda + 1);
	    f = calcf(lambda + 1);
	    trans_rate = calcTM(lambda, deltaE, BM, f);
	    if ( trans_rate > min_rate )
	      init->add_weiss(E_f, trans_rate);
	    break;
		  
	  case 2:
	    lambda = 1;
	    BE = calcBE(lambda + 1);
	    f = calcf(lambda + 1);
	    trans_rate = calcTE(lambda, deltaE, BE, f);
	    if ( trans_rate > min_rate )
	      init->add_weiss(E_f, trans_rate);
	    break;
		  
	  case 3:
	    lambda = 2;
	    BM = calcBM(lambda + 1);
	    f = calcf(lambda + 1);
	    trans_rate = calcTM(lambda, deltaE, BM, f);
	    if ( trans_rate > min_rate )
	      init->add_weiss(E_f, trans_rate);
	    break;
		  
	  default:
	    if( deltaJ % 2 == 0 ){
	      BE = calcBE(deltaJ);
	      f = calcf(deltaJ);
	      trans_rate = calcTE(deltaJ, deltaE, BE, f);
	      if ( trans_rate > min_rate )
		init->add_weiss(E_f, trans_rate);
	    }
	    else if( (deltaJ + 1) % 2 == 0 ){
		BM = calcBM(deltaJ);
		f = calcf(deltaJ);
		trans_rate = calcTM(deltaJ, deltaE, BM, f);
		if ( trans_rate > min_rate )
		  init->add_weiss(E_f, trans_rate);
	      }
	    else
	      std::cout << "Warning! Something bad happened! deltaJ = " << deltaJ << std::endl;
	    break;
	  }
	    
      else //Meaning parity changes
	switch ( deltaJ )
	  {
	  case 0:
	  case 1:
	    lambda = 0;
	    BE = calcBE(deltaJ + 1);
	    f = calcf(deltaJ + 1);
	    trans_rate = calcTE(lambda, deltaE, BE, f);
	    if ( trans_rate > min_rate )
	      init->add_weiss(E_f, trans_rate);
	    break;
		  
	  case 2:
	    lambda = 1;
	    BM = calcBM(deltaJ + 1);
	    f = calcf(deltaJ + 1);
	    trans_rate = calcTM(lambda, deltaE, BM, f);
	    if ( trans_rate > min_rate )
	      init->add_weiss(E_f, trans_rate);
	    break;
		  
	  case 3:
	    lambda = 2;
	    BE = calcBE(deltaJ + 1);
	    f = calcf(deltaJ + 1);
	    trans_rate = calcTE(lambda, deltaE, BE, f);
	    if ( trans_rate > min_rate )
	      init->add_weiss(E_f, trans_rate);
	    break;
		  
	  default:
	    if( deltaJ % 2 == 0 ){
	      BM = calcBM(deltaJ);
	      f = calcf(deltaJ);
	      trans_rate = calcTM(deltaJ, deltaE, BM, f);
	      if ( trans_rate > min_rate )
		init->add_weiss(E_f, trans_rate);
	    }
	    else if( (deltaJ + 1) % 2 == 0 ){
	      BE = calcBM(deltaJ);
	      f = calcf(deltaJ);
	      trans_rate = calcTE(deltaJ, deltaE, BE, f);
	      if ( trans_rate > min_rate )
		init->add_weiss(E_f, trans_rate);
	    }
	    else
	      std::cout << "Warning! Something bad happened! deltaJ = " << deltaJ << std::endl;
	    break;
	  }  
    }
  init->calc_ri();
}

double TMarleyDecayScheme::doubleFact(int i){
  if (i == 0 || i == 1)
    return 1;
  else
    return i * doubleFact(i - 2);
}

double TMarleyDecayScheme::calcBE (int i){
  return pow(1.2, 2*i) * 9 / (4*marley_utils::pi) / (i + 3) / (i + 3) * pow(atomic_mass_number, 2*i/3.0);
}

double TMarleyDecayScheme::calcBM (int i){
  return pow(1.2, 2*i - 2) * 90/marley_utils::pi / (i + 3) / (i + 3) * pow(atomic_mass_number, (2*i - 2)/3.0);
}

double TMarleyDecayScheme::calcf (int i){
  double x = doubleFact(2*i + 1);
  return (i + 1)/(i*x*x);
}

double TMarleyDecayScheme::calcTE(int i, double dE, double BE_i, double f_i){
  return 5.498e+22 * pow (dE/197330, 2*i + 1) * BE_i * f_i;
}

double TMarleyDecayScheme::calcTM(int i, double dE, double BM_i, double f_i){
  return 6.080e+20 * pow (dE/197330, 2*i + 1) * BM_i * f_i;
}
