#pragma once
#include <string>
#include <unordered_map>

#include "ConfigFile.hh"
#include "DecayScheme.hh"
#include "SphericalOpticalModel.hh"

namespace marley {

  class StructureDatabase {
  
    public:
      StructureDatabase();
      StructureDatabase(const marley::ConfigFile& cf);
      void add_decay_scheme(const std::string nucid,
        marley::DecayScheme ds);
  
      inline void clear() {
        decay_scheme_table.clear();
      }
  
      inline void remove_decay_scheme(const std::string nucid) {
        // Remove the decay scheme with this nucid if it
        // exists in the database. If it doesn't, do nothing.
        decay_scheme_table.erase(nucid);
      }
  
      inline marley::DecayScheme* get_decay_scheme(const int particle_id) {
        std::unordered_map<int, marley::DecayScheme*>::iterator it
          = pid_decay_scheme_table.find(particle_id);
        if (it == pid_decay_scheme_table.end()) return nullptr;
        else return it->second;
      }
  
      marley::DecayScheme* get_decay_scheme(const std::string nucid);
  
      inline marley::DecayScheme* get_decay_scheme(const int Z, const int A)
      {
        int particle_id = marley_utils::get_nucleus_pid(Z, A);
        return get_decay_scheme(particle_id);
      }
  
      void add_from_record(const marley::ConfigFile::StructureRecord& sr);
      void add_all_from_config_file(const marley::ConfigFile& cf);
  
      // If the requested optical model object already exists in the lookup table,
      // return it. If not, create it, add it to the table, and then return it.
      inline /*const*/ marley::SphericalOpticalModel& get_optical_model(const int Z,
        const int A)
      {
        int pid = marley_utils::get_nucleus_pid(Z, A);
  
        std::unordered_map<int, marley::SphericalOpticalModel>::iterator
          it = optical_model_table.find(pid);
  
        if (it == optical_model_table.end()) {
  	// The requested optical model wasn't found, so create it and add it to
  	// the table, returning a reference to the stored optical model afterwards.
          return optical_model_table.emplace(pid,
            marley::SphericalOpticalModel(Z, A)).first->second;
        }
        else return it->second;
      }
  
  /***
      inline double get_contbin_width() const { return contbin_width; }
      inline size_t get_contbin_num_subs() const { return contbin_num_subs; }
  ***/
  
    private:
      // Lookup table for marley::DecayScheme objects.
      // Keys are ENSDF-style nucIDs, values are decay schemes.
      std::unordered_map<std::string, marley::DecayScheme> decay_scheme_table;
      // Table for looking up decay schemes by PDG particle ID
      std::unordered_map<int, marley::DecayScheme*> pid_decay_scheme_table;
  
      // Lookup table for marley::SphericalOpticalModel objects.
      // Keys are PDG particle IDs, values are optical models.
      std::unordered_map<int, marley::SphericalOpticalModel> optical_model_table;
  
  /***
      // Resolution to use when binning the nuclear energy level continuum
      double contbin_width;
      // Number of subintervals to use when computing the partial decay width for
      // a continuum bin
      double contbin_num_subs;
  ***/
  };

}
