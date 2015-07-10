#include "TMarleyStructureDatabase.hh"

void TMarleyStructureDatabase::add_decay_scheme(const std::string nucid,
  TMarleyDecayScheme ds)
{
  // If an old entry using this nucid already exists, erase it.
  if (get_decay_scheme(nucid) != nullptr) {
    decay_scheme_table.erase(nucid);
  }

  // Create a new entry in the database using this nucid
  decay_scheme_table.emplace(nucid, std::move(ds));
}

void TMarleyStructureDatabase::clear() {
  decay_scheme_table.clear();
}

void TMarleyStructureDatabase::remove_decay_scheme(const std::string nucid) {
  // Remove the decay scheme with this nucid if it
  // exists in the database. If it doesn't, do nothing.
  decay_scheme_table.erase(nucid);
}

// Get a pointer to the decay scheme in the database that was filed under the
// given nucid. If the requested decay scheme could not be found, return nullptr.
TMarleyDecayScheme* TMarleyStructureDatabase::get_decay_scheme(const std::string nucid) {
  std::unordered_map<std::string, TMarleyDecayScheme>::iterator
    it = decay_scheme_table.find(nucid);
  if (it == decay_scheme_table.end()) return nullptr;
  else return &(it->second);  // Convert the iterator into a pointer to the decay scheme
}
