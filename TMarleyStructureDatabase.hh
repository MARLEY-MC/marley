#pragma once
#include <string>
#include <unordered_map>

#include "TMarleyDecayScheme.hh"
#include "TMarleyConfigFile.hh"

class TMarleyStructureDatabase {

  public:
    TMarleyStructureDatabase();
    TMarleyStructureDatabase(const TMarleyConfigFile& cf);
    void add_decay_scheme(const std::string nucid,
      TMarleyDecayScheme ds);

    inline void clear() {
      decay_scheme_table.clear();
    }

    inline void remove_decay_scheme(const std::string nucid) {
      // Remove the decay scheme with this nucid if it
      // exists in the database. If it doesn't, do nothing.
      decay_scheme_table.erase(nucid);
    }

    inline TMarleyDecayScheme* get_decay_scheme(const int particle_id) {
      std::unordered_map<int, TMarleyDecayScheme*>::iterator it
        = pid_decay_scheme_table.find(particle_id);
      if (it == pid_decay_scheme_table.end()) return nullptr;
      else return it->second;
    }

    TMarleyDecayScheme* get_decay_scheme(const std::string nucid);

    inline TMarleyDecayScheme* get_decay_scheme(const int Z, const int A)
    {
      int particle_id = marley_utils::get_nucleus_pid(Z, A);
      return get_decay_scheme(particle_id);
    }

    void add_from_record(const TMarleyConfigFile::StructureRecord& sr);
    void add_all_from_config_file(const TMarleyConfigFile& cf);

  private:
    // Lookup table for TMarleyDecayScheme objects.
    // Keys are ENSDF-style nucIDs, values are decay schemes.
    std::unordered_map<std::string, TMarleyDecayScheme> decay_scheme_table;
    // Table for looking up decay schemes by PDG particle ID
    std::unordered_map<int, TMarleyDecayScheme*> pid_decay_scheme_table;
};
