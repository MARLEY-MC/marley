#include "marley_utils.hh"
#include "TMarleyStructureDatabase.hh"

TMarleyStructureDatabase::TMarleyStructureDatabase() {
/***
  contbin_width = TMarleyConfigFile::DEFAULT_CONTINUUM_BIN_RESOLUTION;
  contbin_num_subs = TMarleyConfigFile::DEFAULT_CONTINUUM_BIN_SUBINTERVALS;
***/
}

TMarleyStructureDatabase::TMarleyStructureDatabase(
  const TMarleyConfigFile& cf)
{
  add_all_from_config_file(cf);
/***
  contbin_width = cf.get_contbin_width();
  contbin_num_subs = cf.get_contbin_num_subs();
***/
}

void TMarleyStructureDatabase::add_decay_scheme(const std::string nucid,
  TMarleyDecayScheme ds)
{
  // Remove the previous entry (if one exists) filed
  // under the supplied nucid
  decay_scheme_table.erase(nucid);

  // Create a new entry in the database using this nucid
  std::pair<std::unordered_map<std::string,
    TMarleyDecayScheme>::iterator, bool>
    p = decay_scheme_table.emplace(nucid, std::move(ds));

  // Convert an iterator to the added decay scheme
  // into a pointer (this avoids problems that occur
  // when rehashing invalidates iterators to elements
  // of an unordered_map [pointers are not invalidated]).
  TMarleyDecayScheme* ds_p = &(p.first->second);

  // Compute mass number from nucid
  int Z;
  int A = std::stoi(nucid.substr(0,3));

  // Get element symbol from nucid
  std::string element_symbol = nucid.substr(3);
  element_symbol.back() = tolower(element_symbol.back());
  marley_utils::trim_right_inplace(element_symbol);
  // Look up the corresponding atomic number
  std::unordered_map<std::string, int>::const_iterator it
    = marley_utils::atomic_numbers.find(element_symbol);
  if (it != marley_utils::atomic_numbers.end()) Z = it->second;
  else throw std::runtime_error(std::string("Unrecognized")
    + " element symbol '" + element_symbol + "' encountered in"
    + " TMarleyStructureDatabase::add_decay_scheme()");

  // Compute the particle ID for the nucleus based on Z and A
  int pid = marley_utils::get_nucleus_pid(Z, A);

  // Erase the previous entry (if one exists) filed under this
  // pid in the decay scheme table
  pid_decay_scheme_table.erase(pid);

  // Add the decay scheme to the particle-ID-based lookup table
  pid_decay_scheme_table[pid] = ds_p;
}

// Get a pointer to the decay scheme in the database that was filed under the
// given nucid. If the requested decay scheme could not be found, return nullptr.
// Note that the nucid must strictly conform to the ENSDF conventions (3 mass number
// characters [left-padded] followed by 2 upper-case element symbol
// characters [right-padded])
TMarleyDecayScheme* TMarleyStructureDatabase::get_decay_scheme(
  const std::string nucid)
{
  std::unordered_map<std::string, TMarleyDecayScheme>::iterator
    it = decay_scheme_table.find(nucid);
  if (it == decay_scheme_table.end()) return nullptr;
  else return &(it->second);  // Convert the iterator into a pointer to the decay scheme
}

// Add decay schemes to the database from a parsed
// configuration file structure record
void TMarleyStructureDatabase::add_from_record(
  const TMarleyConfigFile::StructureRecord& sr)
{
  // TODO: consider altering the parsing process used here so
  // that all nuclides are loaded from the file in one pass. 
  // Add a decay scheme for each nucid listed in the structure
  // record to the database
  for (const auto& id : sr.nucids) {
    std::cout << "Loading nuclear structure data for "
      << marley_utils::trim_copy(id)
      << " from file " << sr.filename << std::endl;
    add_decay_scheme(id,
      TMarleyDecayScheme(id, sr.filename, sr.format));
  }
}

void TMarleyStructureDatabase::add_all_from_config_file(
  const TMarleyConfigFile& cf)
{
  for (const auto& record : cf.get_structure_records()) {
    add_from_record(record);
  }
}
