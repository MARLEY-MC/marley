#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "marley_utils.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleyGenerator.hh"
#include "TMarleyNuclearPhysics.hh"

// Default level parity is +
const TMarleyParity TMarleyDecayScheme::DEFAULT_PARITY = 1;

// ***** Constants used for ENSDF file parsing *****
const std::string TMarleyDecayScheme::ensdf_primary_record = "[ 1]";
const std::string TMarleyDecayScheme::ensdf_continuation_record = "[^ 1]";
const std::string TMarleyDecayScheme::ensdf_record_data = ".{71}";
const std::regex TMarleyDecayScheme::rx_ensdf_end_record("\\s*"); // Matches blank lines
const std::regex TMarleyDecayScheme::rx_paren_group("\\(([^)]+)\\)([+-])");
const std::regex TMarleyDecayScheme::rx_spin("([0-9]+(?:/2)?)");
const std::regex TMarleyDecayScheme::rx_jpi("([0-9]+(?:/2)?)([+-])?");

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

void TMarleyDecayScheme::do_cascade(TMarleyLevel* initial_level,
  TMarleyEvent* p_event, TMarleyGenerator& gen) const
{
  //std::cout << "Beginning gamma cascade at level with energy "
  //  << initial_level->get_energy() << " MeV" << std::endl;

  bool cascade_finished = false;

  TMarleyLevel* p_current_level = initial_level;

  while (!cascade_finished) {
    // Randomly select a gamma to produce
    TMarleyGamma* p_gamma = p_current_level->sample_gamma(gen);
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
      //  << " MeV." New level has energy " << p_current_level->get_energy()
      //  << " MeV." << std::endl;
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
      double gamma_cos_theta = gen.uniform_random_double(-1, 1, true);
      double gamma_theta = std::acos(gamma_cos_theta);

      // sample from [0,2*pi)
      double gamma_phi = gen.uniform_random_double(0, 2*marley_utils::pi, false);

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
  //  << p_current_level->get_energy() << std::endl;
}

TMarleyDecayScheme::TMarleyDecayScheme(int z, int a, std::string filename,
  TMarleyDecayScheme::FileFormat ff)
{
  Z = Z;
  A = a;
  nuc_id = marley_utils::nuc_id(z, a);
  parse(filename, ff);
}


TMarleyDecayScheme::TMarleyDecayScheme(std::string nucid, std::string filename,
  TMarleyDecayScheme::FileFormat ff)
{

  nuc_id = nucid;

  // Get the atomic mass number from the first 3 characters
  // of the ENSDF nucid
  A = std::stoi(nucid.substr(0,3));
  // Get the atomic number from the nucid using the element symbol
  Z = marley_utils::nucid_to_Z(nucid);
  // Parse the data file
  parse(filename, ff);
}

// Load a level with TMarleyGamma objects based on theoretical gamma-ray
// transition relative intensities
void TMarleyDecayScheme::assign_theoretical_RIs(TMarleyLevel* level_i) {
  // Remove any existing gamma objects from the level
  level_i->clear_gammas();

  // Loop over the decay scheme's levels in ascending order of energy until you
  // reach the initial level.
  int l;
  bool initial_spin_is_zero = level_i->get_two_J() == 0;
  for (auto level_f : pv_sorted_levels) {

    // If we've reached the initial level, end the loop
    if (level_f == level_i) break;

    // 0->0 EM transitions aren't allowed due to angular momentum conservation
    // (photons are spin 1), so if the initial and final spins are zero, skip
    // ahead to the next final level.
    if (initial_spin_is_zero && level_f->get_two_J() == 0) continue;

    // Determine the type (electric or magnetic) and multipolarity of the gamma
    // transition between these two levels
    auto type = TMarleyNuclearPhysics::determine_gamma_transition_type(level_i,
      level_f, l);

    // Approximate the gamma energy by the energy difference between the two levels
    // TODO: consider adding a nuclear recoil correction here
    double e_gamma = level_i->get_energy() - level_f->get_energy();

    // TODO: allow the user to choose which gamma-ray transition model to use 
    // (Weisskopf single-particle estimates, Brink-Axel strength functions, etc.)
    // Note that we don't need to normalize the relative intensities since the
    // std::discrete_distribution already takes care of that for us
    double ri = TMarleyNuclearPhysics::weisskopf_partial_decay_width(A, type,
      l, e_gamma);

    // Add a new gamma object representing this transition to the initial level
    level_i->add_gamma(TMarleyGamma(e_gamma, ri, level_i));
  }
}

void TMarleyDecayScheme::parse_ensdf(std::string filename) {

  // Regular expressions for identifying ensdf record types
  const std::regex rx_primary_identification_record(nuc_id + ensdf_primary_record
    + "   ADOPTED LEVELS.{57}"); //"   ADOPTED LEVELS(, GAMMAS.{49}|.{57})");
  const std::regex rx_primary_level_record(nuc_id + ensdf_primary_record + " L " + ensdf_record_data);
  const std::regex rx_continuation_level_record(nuc_id + ensdf_continuation_record + " L " + ensdf_record_data);
  const std::regex rx_primary_gamma_record(nuc_id + ensdf_primary_record + " G " + ensdf_record_data);
  const std::regex rx_continuation_gamma_record(nuc_id + ensdf_continuation_record + " G " + ensdf_record_data);
 
  // Open the ENSDF file for parsing
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw std::runtime_error(std::string("Could not read from the ") +
      "ENSDF data file " + filename);
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
      + " could not be found in the ENSDF data file " + filename);
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
      line = process_continuation_records(file_in, record, rx_continuation_level_record);
      no_advance = true;
        
      // Extract the level energy (in keV) as a trimmed string from the ENSDF level record
      std::string level_energy_str = marley_utils::trim_copy(record.substr(9,10)); 

      // Convert the energy to MeV and double-precision
      double level_energy = std::stod(level_energy_str) * marley_utils::MeV;

      // Also extract the spin-parity of the level
      std::string jpi = marley_utils::trim_copy(record.substr(21,18));

      // Process the spin-parity string
      // TODO: Improve the processing procedure used here. Right now, all that
      // is done is the following:
      // 1. Preprocess any parenthesized groups with multiple J values and
      // a single parity so that "(1,2,3)+" becomes "(1+,2+,3+)"
      // 2. Search for the first number plus sign or number without a sign.
      // If one is found, use it to determine the level spin-parity. If not,
      // use the default values.
      std::smatch match;
    
      // If a parity is applied to a parenthesized group in the ENSDF file, remove
      // the parentheses and apply that parity to every spin in the group
      if (std::regex_match(jpi, match, rx_paren_group)) {
        std::string spins = match[1];
        std::string parity = match[2];
        jpi = std::regex_replace(spins, rx_spin, "$1" + parity);
      }
    
      int two_J;
      TMarleyParity parity;
    
      // Search for a spin-parity entry like "3" or "1/2-"
      if (std::regex_search(jpi, match, rx_jpi)) {
        // Determine the level parity
        if (match.size() > 2) {
          if (match[2] == "+") parity = TMarleyParity(1);
          else if (match[2] == "-") parity = TMarleyParity(-1);
          else if (marley_utils::trim_copy(match[2]) == "")
            parity = DEFAULT_PARITY;
          else throw std::runtime_error(std::string("Invalid parity ")
            + match[2].str() + " encountered while setting a TMarleyLevel "
            + "parity value.");
        }
        // Use the default parity value if the parity could not be determined
        else parity = DEFAULT_PARITY;
    
        std::string spin = match[1]; 
    
        // Check for half-integer spin and compute 2J appropriately
        // based on the result of the check
        size_t spin_size = spin.size();
        if (spin_size > 2 && spin.substr(spin_size - 2) == "/2")
          two_J = std::stoi(spin.substr(0, spin_size - 2));
        else two_J = 2 * std::stoi(spin);
      }
      else {
        parity = DEFAULT_PARITY;
        if (A % 2) two_J = DEFAULT_FERMION_TWOJ;
        else two_J = DEFAULT_BOSON_TWOJ;
      }

      // Add a new level object to this decay scheme object using these data
      // and update the current level pointer to point to it
      p_current_level = add_level(TMarleyLevel(level_energy,
        two_J, parity));
    }

    // Gamma Record
    else if (std::regex_match(line, rx_primary_gamma_record)) {
      //std::cout << "DEBUG:   Parsing gamma" << std::endl;
      record = line;
      line = process_continuation_records(file_in,
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
  for(std::vector<TMarleyLevel*>::iterator j = pv_sorted_levels.begin();
    j != pv_sorted_levels.end(); ++j)
  {
    // Calculate theoretical gamma-ray relative intensities here if they are
    // unknown for the current level
    if( (*j)->get_gamma_status() == false ) {
      assign_theoretical_RIs(*j);
    }

    std::vector<TMarleyGamma>* p_gammas = (*j)->get_gammas();
    double initial_level_energy = (*j)->get_energy();

    for(std::vector<TMarleyGamma>::iterator k = p_gammas->begin();
      k != p_gammas->end(); ++k) 
    {
      double gamma_energy = k->get_energy();

      // Approximate the final level energy so we can search for the final level
      double final_level_energy = initial_level_energy - gamma_energy; 

      // Search for the corresponding energy using the vector of level energies
      std::vector<double>::iterator p_final_level_energy = std::lower_bound(
        sorted_level_energies.begin(), sorted_level_energies.end(),
        final_level_energy); 

      // Determine the index of the final level energy appropriately
      unsigned int e_index = std::distance(sorted_level_energies.begin(),
        p_final_level_energy);
      if (e_index == sorted_level_energies.size()) {
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
        if (std::abs(final_level_energy - sorted_level_energies[e_index])
          > std::abs(final_level_energy - sorted_level_energies[e_index - 1]))
        {
          --e_index;
        }
      }

      // Use the index to assign the appropriate end level pointer to this gamma
      k->set_end_level(pv_sorted_levels[e_index]);

    }

  }

  file_in.close();

}

void TMarleyDecayScheme::parse_talys(std::string filename) {
  // First line in a TALYS level dataset has fortran
  // format (2i4, 2i5, 56x, i4, a2)
  // General regex for this line:
  // std::regex nuclide_line("[0-9 ]{18} {56}[0-9 ]{4}.{2}");
  std::string nuc_id_copy = nuc_id;
  // Make the last character of the nuc_id lowercase to follow the TALYS convention
  nuc_id_copy.back() = tolower(nuc_id_copy.back());
  
  const std::regex nuclide_line("[0-9 ]{18} {57}" + nuc_id_copy);

  // Open the TALYs level data file for parsing
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw std::runtime_error(std::string("Could not read from the ") +
      "TALYS data file " + filename);
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
      + " could not be found in the TALYS data file " + filename);
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
    int level_num, pi, num_gammas;
    double level_energy, spin;
    iss >> level_num >> level_energy >> spin >> pi >> num_gammas;

    // Compute two times the spin so that we can represent half-integer
    // nuclear level spins as integers
    int twoJ = std::round(2 * spin);

    // Create a parity object to use when constructing the level
    TMarleyParity parity = pi;

    // Construct a new level object and add it to the decay scheme. Get
    // a pointer to the newly-added level
    p_current_level = add_level(TMarleyLevel(level_energy, twoJ, parity));

    // Add the energy for this level and a pointer to the newly-created level
    // object to our temporary arrays.  These will be used when creating gammas
    // for each level.
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

void TMarleyDecayScheme::print_report(std::ostream& ostr) const {
  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for(auto j : pv_sorted_levels)
  {

    int twoj = j->get_two_J();
    std::string spin = std::to_string(twoj / 2);
    // If 2*J is odd, then the level has half-integer spin
    if (twoj % 2) spin += "/2";
    TMarleyParity parity = j->get_parity();
    
    //if (sp.empty()) sp = "UNKNOWN";
    
    ostr << "Level at " << j->get_energy()
	 << " MeV has spin-parity " << spin << parity << std::endl;
    std::vector<TMarleyGamma>* p_gammas = j->get_gammas();

    // Cycle through each of the gammas owned by the current level
    // (according to the ENSDF specification, these will already be
    // sorted in order of increasing energy)
    for(auto &k : *p_gammas)
    {
      ostr << "  has a gamma with energy " << k.get_energy() << " MeV";
      ostr << " (transition to level at "
        << k.get_end_level()->get_energy() << " MeV)" << std::endl;
      ostr << "    and relative photon intensity " << k.get_ri() << std::endl; 
    }

  }

}

void TMarleyDecayScheme::print_latex_table(std::ostream& ostr) {

  std::string caption_beginning =
    std::string("{\\textbf{Levels") +
    " and $\\boldsymbol{\\gamma}$ transitions \n for " +
    "\\isotope[\\boldsymbol{" +
    marley_utils::trim_copy(nuc_id.substr(0,3)) +
    "}]{\\textbf{" + nuc_id.substr(3,1) +
    marley_utils::trim_copy(
      marley_utils::to_lowercase(nuc_id.substr(4,1))) + 
    "}} \n";// from file " + filename + " ";

  ostr << marley_utils::latex_table_1;

  ostr << caption_beginning + "}}\\\\\n";

  ostr << marley_utils::latex_table_2;

  ostr << caption_beginning + " -- \\textit{continued}}} \\\\\n";  

  ostr << marley_utils::latex_table_3;
  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for(std::vector<TMarleyLevel*>::iterator j = pv_sorted_levels.begin();
    j != pv_sorted_levels.end(); ++j)
  {

    std::string sp = (*j)->get_spin_parity_string();
    //if (sp.empty()) sp = "?";
    //if (sp == "UNNATURAL") sp = "unnat.";

    ostr << (*j)->get_energy() << " & "
      << sp  << " & ";

    std::vector<TMarleyGamma>* p_gammas = (*j)->get_gammas();

    // If there aren't any gammas for this level, finish writing
    // the current row of the table. Add extra space between this
    // level and the next one.
    if (p_gammas->empty()) {
      ostr << " &  &";
      // If this is the last row of the table, don't add extra space.
      if (j == pv_sorted_levels.end() - 1) {
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
        << k->get_end_level()->get_energy();
      // Add vertical space after the final gamma row. Also prevent page breaks
      // in the middle of a list of gammas by outputting a star at the end of
      // each row except the final gamma row.
      if (k == p_gammas->end() - 1) {
	// Don't add the extra row space for the very last row in the table
	if (j >= pv_sorted_levels.end() - 1) {
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

// Adds a new level to the decay scheme and returns a pointer to it
TMarleyLevel* TMarleyDecayScheme::add_level(TMarleyLevel level) {

  // Add the level to the list of level objects.
  levels.push_back(level);

  // Get a pointer to the just-added level object
  TMarleyLevel* p_level = &levels.back();

  // Get this level's energy
  double l_energy = p_level->get_energy();

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

  // Return a pointer to the newly-added level
  return p_level;
}

std::vector<TMarleyLevel*>* TMarleyDecayScheme::get_sorted_level_pointers() {
  return &pv_sorted_levels;
}

std::list<TMarleyLevel>* TMarleyDecayScheme::get_levels() {
  return &levels;
}

std::ostream& operator<< (std::ostream& out,
  const TMarleyDecayScheme& ds)
{
  size_t num_levels = ds.pv_sorted_levels.size();
  out << ds.Z << " " << ds.A << " " << num_levels << std::endl;
  for (size_t i = 0; i < num_levels; ++i) {
    TMarleyLevel* l = ds.pv_sorted_levels.at(i);
    out << "  " << l->get_energy() << " " << l->get_two_J() << " "
      << l->get_parity() << " " << l->get_gamma_status()
      << " " << l->get_gammas()->size() << std::endl;
    for (auto &g : *l->get_gammas()) {
      out << "    " << g.get_energy() << " " << g.get_ri();
      std::vector<TMarleyLevel*>::const_iterator pv_end
        = ds.pv_sorted_levels.cbegin() + i;
      std::vector<TMarleyLevel*>::const_iterator cit
        = std::find(ds.pv_sorted_levels.cbegin(), pv_end, g.get_end_level());
      int level_f_idx = -1;
      if (cit != pv_end) {
        level_f_idx = std::distance(ds.pv_sorted_levels.cbegin(), cit);
      }
      out << " " << level_f_idx << std::endl;
    }
  }
  return out;
}

std::istream& operator>> (std::istream& in,
  TMarleyDecayScheme& ds)
{
  size_t num_levels;
  in >> ds.Z >> ds.A >> num_levels;

  double energy, ri;
  int two_j;
  TMarleyParity pi;
  bool g_status;
  size_t num_gammas;
  int level_f_idx;
  TMarleyLevel* l;
  TMarleyGamma* g;

  for (size_t i = 0; i < num_levels; ++i) {
    in >> energy >> two_j >> pi >> g_status >> num_gammas;
    l = ds.add_level(TMarleyLevel(energy, two_j, pi));
    for (size_t j = 0; j < num_gammas; ++j) {
      in >> energy >> ri >> level_f_idx;
      g = l->add_gamma(TMarleyGamma(energy, ri, l));
      if (level_f_idx >= 0)
        g->set_end_level(ds.pv_sorted_levels.at(level_f_idx));
    }
  }
  return in;
}
