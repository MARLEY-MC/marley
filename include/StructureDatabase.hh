#pragma once
#include <memory>
#include <string>
#include <unordered_map>

#include "BackshiftedFermiGasModel.hh"
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
        auto iter = pid_decay_scheme_table.find(particle_id);
        if (iter == pid_decay_scheme_table.end()) return nullptr;
        else return iter->second;
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
      inline marley::SphericalOpticalModel& get_optical_model(int nucleus_pid)
      {
        // TODO: add check for invalid nucleus particle ID value
        auto iter = optical_model_table.find(nucleus_pid);

        if (iter == optical_model_table.end()) {
	  // The requested level density model wasn't found, so create it and add
	  // it to the table, returning a reference to the stored level density
	  // model afterwards.
	  int Z = marley::MassTable::get_particle_Z(nucleus_pid);
	  int A = marley::MassTable::get_particle_A(nucleus_pid);
          return *(optical_model_table.emplace(nucleus_pid,
            std::make_unique<marley::SphericalOpticalModel>(Z, A)).first
            ->second.get());
        }
        else return *(iter->second.get());
      }

      // If the requested optical model object already exists in the lookup table,
      // return it. If not, create it, add it to the table, and then return it.
      inline marley::SphericalOpticalModel& get_optical_model(const int Z,
        const int A)
      {
        int nucleus_pid = marley_utils::get_nucleus_pid(Z, A);
        auto iter = optical_model_table.find(nucleus_pid);

        if (iter == optical_model_table.end()) {
	// The requested level density model wasn't found, so create it and add
	// it to the table, returning a reference to the stored level density
	// model afterwards.
          return *(optical_model_table.emplace(nucleus_pid,
            std::make_unique<marley::SphericalOpticalModel>(Z, A)).first
            ->second.get());
        }
        else return *(iter->second.get());
      }

      // If the requested level density model object already exists in the
      // lookup table, return it. If not, create it, add it to the table, and
      // then return it.
      inline marley::LevelDensityModel& get_level_density_model(const int Z,
        const int A)
      {
        int pid = marley_utils::get_nucleus_pid(Z, A);

        auto iter = level_density_table.find(pid);

        if (iter == level_density_table.end()) {
	// The requested level density model wasn't found, so create it and add
	// it to the table, returning a reference to the stored level density
	// model afterwards.
          return *(level_density_table.emplace(pid,
            std::make_unique<marley::BackshiftedFermiGasModel>(Z, A)).first
            ->second.get());
        }
        else return *(iter->second.get());
      }


    private:
      // Lookup table for marley::DecayScheme objects.
      // Keys are ENSDF-style nucIDs, values are decay schemes.
      std::unordered_map<std::string, marley::DecayScheme> decay_scheme_table;
      // Table for looking up decay schemes by PDG particle ID
      std::unordered_map<int, marley::DecayScheme*> pid_decay_scheme_table;

      // Lookup table for marley::SphericalOpticalModel objects.
      // Keys are PDG particle IDs, values are optical models.
      std::unordered_map<int, std::unique_ptr<marley::SphericalOpticalModel> >
        optical_model_table;

      // Lookup table for marley::LevelDensityModel objects.
      // Keys are PDG particle IDs, values are unique_ptrs to level density
      // models.
      std::unordered_map<int, std::unique_ptr<marley::LevelDensityModel> >
        level_density_table;
  };

}
