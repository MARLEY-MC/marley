#pragma once
#include <string>
#include <unordered_map>

#include "TMarleyConfigFile.hh"
#include "TMarleyDecayScheme.hh"
#include "TMarleySphericalOpticalModel.hh"

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

    // If the requested optical model object already exists in the lookup table,
    // return it. If not, create it, add it to the table, and then return it.
    inline const TMarleySphericalOpticalModel& get_optical_model(const int Z,
      const int A)
    {
      int pid = marley_utils::get_nucleus_pid(Z, A);

      std::unordered_map<int, TMarleySphericalOpticalModel>::iterator
        it = optical_model_table.find(pid);

      if (it == optical_model_table.end()) {
	// The requested optical model wasn't found, so create it and add it to
	// the table, returning a reference to the stored optical model afterwards.
        return optical_model_table.emplace(pid,
          TMarleySphericalOpticalModel(Z, A)).first->second;
      }
      else return it->second;
    }

  private:
    // Lookup table for TMarleyDecayScheme objects.
    // Keys are ENSDF-style nucIDs, values are decay schemes.
    std::unordered_map<std::string, TMarleyDecayScheme> decay_scheme_table;
    // Table for looking up decay schemes by PDG particle ID
    std::unordered_map<int, TMarleyDecayScheme*> pid_decay_scheme_table;

    // Lookup table for TMarleySphericalOpticalModel objects.
    // Keys are PDG particle IDs, values are optical models.
    std::unordered_map<int, TMarleySphericalOpticalModel> optical_model_table;
};
