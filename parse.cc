#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
#include "ensdf_utils.hh"

class TEnsdfLevel;
class TEnsdfGamma;
class TEnsdfDecayScheme;

/// A class that represents an ENSDF nuclear energy level

class TEnsdfLevel {
  public:
    /// @param energy a string containing the
    /// level energy in keV
    /// @param jpi a string containing the spin
    /// and parity of the level (e.g., 0+)
    TEnsdfLevel(std::string energy = "0", std::string jpi = "");
    void add_gamma(const TEnsdfGamma& gamma);
    void clear_gammas();
    std::vector<TEnsdfGamma>* get_gammas();
    double get_numerical_energy() const;
    std::string get_string_energy() const;
    std::string get_spin_parity() const;
    void set_energy(std::string energy);
    void set_spin_parity(std::string jpi);

  private:
    std::string sEnergy;
    double fEnergy;
    std::string spin_parity;
    std::vector<TEnsdfGamma> gammas;
};

TEnsdfLevel::TEnsdfLevel(std::string energy, std::string jpi) {
  sEnergy = energy;
  fEnergy = std::stod(energy); // Converts the energy string into a double
  spin_parity = jpi;
}

void TEnsdfLevel::add_gamma(const TEnsdfGamma& gamma) {
  gammas.push_back(gamma);
}

double TEnsdfLevel::get_numerical_energy() const {
  return fEnergy; 
}

std::string TEnsdfLevel::get_string_energy() const {
  return sEnergy; 
}

void TEnsdfLevel::set_energy(std::string energy) {
  sEnergy = energy;
  fEnergy = std::stod(energy);
}

std::string TEnsdfLevel::get_spin_parity() const {
  return spin_parity;
}

void TEnsdfLevel::set_spin_parity(std::string jpi) {
  spin_parity = jpi;
}

void TEnsdfLevel::clear_gammas() {
  gammas.clear();
}

std::vector<TEnsdfGamma>* TEnsdfLevel::get_gammas() {
  return &gammas;
}


class TEnsdfGamma {
  public:
    TEnsdfGamma(double energy = 0, double ri = 0, TEnsdfLevel* start_level = nullptr);
    void set_start_level(TEnsdfLevel* start_level);
    void set_end_level(TEnsdfLevel* end_level);
    TEnsdfLevel* get_end_level() const;
    TEnsdfLevel* get_start_level() const;
    double get_energy() const;
    double get_ri() const;

  private:
    double fEnergy;
    double fRI;
    double fCC;
    double fTI; 
    TEnsdfLevel* pStartLevel;
    TEnsdfLevel* pEndLevel;
};

TEnsdfGamma::TEnsdfGamma(double energy, double ri, TEnsdfLevel* start_level) {
  fEnergy = energy;
  fRI = ri;
  fCC = 0;
  fTI = 0;
  pStartLevel = start_level;
  pEndLevel = nullptr;
}

void TEnsdfGamma::set_start_level(TEnsdfLevel* start_level) {
  pStartLevel = start_level;
}

void TEnsdfGamma::set_end_level(TEnsdfLevel* end_level) {
  pEndLevel = end_level;
}

TEnsdfLevel* TEnsdfGamma::get_start_level() const {
  return pStartLevel;
}

TEnsdfLevel* TEnsdfGamma::get_end_level() const {
  return pEndLevel;
}

double TEnsdfGamma::get_energy() const {
  return fEnergy;
}

double TEnsdfGamma::get_ri() const {
  return fRI;
}


class TEnsdfDecayScheme {
  public:
    TEnsdfDecayScheme(std::string nucid, std::string filename);
    std::string get_nuc_id() const;
    void set_nuc_id(std::string id);
    void add_level(const TEnsdfLevel level);
    std::map<std::string, TEnsdfLevel>* get_levels();
    TEnsdfLevel* get_level(std::string energy);
    std::vector<TEnsdfLevel*>* get_sorted_level_pointers();
    void print_report();
    void do_cascade(std::string initial_energy);
    void do_cascade(double initial_energy);
    void do_cascade(TEnsdfLevel* initial_level);

  private:
    std::string nuc_id;
    std::map<std::string, TEnsdfLevel> levels;
    std::vector<TEnsdfLevel*> pv_sorted_levels;
    std::vector<double> sorted_level_energies;
    static bool compare_level_energies(TEnsdfLevel* first,
      TEnsdfLevel* second);
    std::string process_continuation_records(std::ifstream &file_in,
      std::string &record, std::regex &rx_cont_record) const;

};

void TEnsdfDecayScheme::do_cascade(std::string initial_energy) {
  std::map<std::string, TEnsdfLevel>::iterator it = levels.find(initial_energy);
  if (it == levels.end()) {
    throw std::range_error("Could not do cascade. Level with energy "
      + initial_energy + " keV not found.");
  }
  else {
    TEnsdfLevel* plevel = &(it->second);
    do_cascade(plevel);
  }
}

void TEnsdfDecayScheme::do_cascade(double initial_energy) {

  // Since we were given a numerical initial energy in this
  // version of the function, search for the level whose energy is
  // closest to the given value of initial_energy 
  std::vector<double>::iterator it = std::lower_bound(
    sorted_level_energies.begin(), sorted_level_energies.end(),
    initial_energy); 

  // Determine the index of the initial level energy appropriately
  int e_index = std::distance(sorted_level_energies.begin(), it);
  if (e_index == sorted_level_energies.size()) {
    // The given initial energy is greater than every level energy in our
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
    if (std::abs(initial_energy - sorted_level_energies[e_index])
      > std::abs(initial_energy - sorted_level_energies[e_index - 1]))
    {
      --e_index;
    }
  }

  TEnsdfLevel* plevel = pv_sorted_levels[e_index];
    
  do_cascade(plevel);
}


void TEnsdfDecayScheme::do_cascade(TEnsdfLevel* initial_level) {
  std::cout << "Beginning gamma cascade at level with energy "
    << initial_level->get_string_energy() << std::endl;

  bool cascade_finished = false;

  TEnsdfLevel* current_level = initial_level;

  //while (!cascade_finished) {

  //}
}


TEnsdfDecayScheme::TEnsdfDecayScheme(std::string nucid, std::string filename) {

  this->nuc_id = nucid;

  // Regular expressions for identifying ensdf record types
  //std::string generic_nuc_id = "^[[:alnum:] ]{5}";
  std::string primary_record = "[ 1]";
  std::string continuation_record = "[^ 1]";
  std::string record_data = ".{71}";
  
  std::regex rx_end_record("\\s*"); // Matches blank lines
  //std::regex rx_generic_primary_identification_record(nuc_id + primary_record + "   " + record_data);
  std::regex rx_primary_identification_record(nuc_id + primary_record
    + "   ADOPTED LEVELS, GAMMAS.{49}");
  std::regex rx_primary_level_record(nuc_id + primary_record + " L " + record_data);
  std::regex rx_continuation_level_record(nuc_id + continuation_record + " L " + record_data);
  std::regex rx_primary_gamma_record(nuc_id + primary_record + " G " + record_data);
  std::regex rx_continuation_gamma_record(nuc_id + continuation_record + " G " + record_data);
 
  // Open the ENSDF file for parsing
  std::ifstream file_in(filename);

  std::string line; // String to store the current line
                    // of the ENSDF file during parsing

  std::string record; // String that will store the full text
                      // (all lines) of the current ENSDF record

  bool found_decay_scheme = false; // Flag that indicates whether or not
                                   // gamma decay scheme data were found for the
                                   // current nuc_id.

  while (!file_in.eof()) {
    std::getline(file_in, line);
    if (std::regex_match(line, rx_primary_identification_record)) {
      found_decay_scheme = true;
      break;
    }
  }

  //if (!found_decay_scheme) {
  //  std::cout << "Gamma decay scheme data (adopted levels, gammas) for " + nuc_id << std::endl;
  //  std::cout << "could not be found in the file " + filename << std::endl;
  //}
  //else {
  //  std::cout << "Gamma decay scheme data for " + nuc_id << " found. Using ENSDF dataset" << std::endl;
  //  std::cout << line << std::endl;
  //}

  bool no_advance = false; // Flag that prevents advancing through
                           // the ENSDF file when a continuation record
                           // is found

  TEnsdfLevel* p_current_level = nullptr; // Pointer to the current level object
                                          // being filled with gamma ray data 

  while (!file_in.eof()) {
    // Get the next line of the file
    // unless we already did
    if (!no_advance) std::getline(file_in, line);
    no_advance = false;

    // Level Record
    if (std::regex_match(line, rx_primary_level_record)) {
      record = line;
      line = this->process_continuation_records(file_in, record, rx_continuation_level_record);
      no_advance = true;
        
      // Extract the level energy (in keV) as a trimmed string from the ENSDF level record
      std::string level_energy = ensdf_utils::trim_copy(record.substr(9,10)); 

      // Also extract the spin-parity of the level
      std::string spin_parity = ensdf_utils::trim_copy(record.substr(21,18));

      // Add a new level object to this decay scheme object using these data
      this->add_level(TEnsdfLevel(level_energy, spin_parity));

      // Update the current level pointer
      p_current_level = this->get_level(level_energy);
    }

    // Gamma Record
    else if (std::regex_match(line, rx_primary_gamma_record)) {
      record = line;
      line = this->process_continuation_records(file_in,
        record, rx_continuation_gamma_record);
      no_advance = true;

      // Extract the gamma ray's energy and relative (photon)
      // intensity from the ENSDF gamma record
      double gamma_energy = std::stod(record.substr(9,10));
      double gamma_ri = ensdf_utils::str_to_double(record.substr(21,8));

      // If this gamma belongs to a level record, then add its
      // data to the corresponding level object. Gammas that
      // have no assigned level appear before any level records,
      // so p_current_level will be a null pointer for them.
      if (p_current_level != nullptr) {
        p_current_level->add_gamma(TEnsdfGamma(gamma_energy, gamma_ri, p_current_level)); 
      }

    }

    else if (std::regex_match(line, rx_end_record)) {
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
  for(std::vector<TEnsdfLevel*>::iterator j = this->pv_sorted_levels.begin();
    j != this->pv_sorted_levels.end(); ++j)
  {
    
    std::vector<TEnsdfGamma>* p_gammas = (*j)->get_gammas();
    double initial_level_energy = (*j)->get_numerical_energy();

    for(std::vector<TEnsdfGamma>::iterator k = p_gammas->begin();
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
      int e_index = std::distance(this->sorted_level_energies.begin(),
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

std::string TEnsdfDecayScheme::process_continuation_records(std::ifstream &file_in,
  std::string &record, std::regex &rx_cont_record) const {

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

}

void TEnsdfDecayScheme::print_report() {
  // Cycle through each of the levels owned by this decay scheme
  // object in order of increasing energy
  for(std::vector<TEnsdfLevel*>::iterator j = this->pv_sorted_levels.begin();
    j != this->pv_sorted_levels.end(); ++j)
  {

    std::string sp = (*j)->get_spin_parity();
    if (sp.empty()) sp = "UNKNOWN";
    
    std::cout << "Level at " << (*j)->get_string_energy()
      << " keV has spin-parity " << sp << std::endl;
    std::vector<TEnsdfGamma>* p_gammas = (*j)->get_gammas();

    // Cycle through each of the gammas owned by the current level
    // (according to the ENSDF specification, these will already be
    // sorted in order of increasing energy)
    for(std::vector<TEnsdfGamma>::iterator k = p_gammas->begin();
      k != p_gammas->end(); ++k) 
    {
      std::cout << "  has a gamma with energy " << k->get_energy() << " keV";
      std::cout << " (transition to level at "
        << k->get_end_level()->get_string_energy() << " keV)" << std::endl;
      std::cout << "    and relative photon intensity " << k->get_ri() << std::endl; 
    }

  }

}

std::string TEnsdfDecayScheme::get_nuc_id() const {
  return nuc_id;
}

void TEnsdfDecayScheme::set_nuc_id(std::string id) {
  nuc_id = id;
}

bool TEnsdfDecayScheme::compare_level_energies(TEnsdfLevel* first,
  TEnsdfLevel* second)
{
  return first->get_numerical_energy() < second->get_numerical_energy();
}

void TEnsdfDecayScheme::add_level(TEnsdfLevel level) {
  // Add the level to the std::map of level objects. Use the
  // string version of its energy in keV as the key.
  std::string energy_string = level.get_string_energy();
  levels[energy_string] = level;

  // Get a pointer to the just-added level object
  TEnsdfLevel* p_level = &(levels[energy_string]);

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


std::vector<TEnsdfLevel*>* TEnsdfDecayScheme::get_sorted_level_pointers() {
  return &pv_sorted_levels;
}


std::map<std::string, TEnsdfLevel>* TEnsdfDecayScheme::get_levels() {
  return &levels;
}

TEnsdfLevel* TEnsdfDecayScheme::get_level(std::string energy) {
  return &(levels.at(energy));
}


int main() {

  std::string nuc_id = " 40K ";
  std::string filename = "ensdf.040";

  // Create a decay scheme object to store data
  // imported from the ENSDF file
  TEnsdfDecayScheme decay_scheme(nuc_id, filename);

  // Print a report describing the decay scheme
  decay_scheme.print_report();

  std::cout << std::endl << std::endl;

  // Test the decay scheme object by simulating a sample gamma cascade
  decay_scheme.do_cascade(1e10);

}
