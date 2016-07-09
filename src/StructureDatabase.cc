#include "marley_utils.hh"
#include "Logger.hh"
#include "StructureDatabase.hh"

marley::StructureDatabase::StructureDatabase() {}

void marley::StructureDatabase::emplace_decay_scheme(int pdg,
  const std::string& filename, DecayScheme::FileFormat format)
{
  int Z_ds = (pdg % 10000000)/10000;
  int A_ds = (pdg % 10000)/10;

  // Remove the previous entry (if one exists) for the given PDG code
  decay_scheme_table_.erase(pdg);

  // Add the new entry
  decay_scheme_table_.emplace(pdg, std::make_unique<marley::DecayScheme>(
    Z_ds, A_ds, filename, format));
}

marley::DecayScheme* marley::StructureDatabase::get_decay_scheme(
  const int particle_id)
{
  auto iter = decay_scheme_table_.find(particle_id);
  if (iter == decay_scheme_table_.end()) return nullptr;
  else return iter->second.get();
}

marley::DecayScheme* marley::StructureDatabase::get_decay_scheme(const int Z,
  const int A)
{
  int particle_id = marley_utils::get_nucleus_pid(Z, A);
  return get_decay_scheme(particle_id);
}

marley::SphericalOpticalModel& marley::StructureDatabase::get_optical_model(
  int nucleus_pid)
{
  /// @todo add check for invalid nucleus particle ID value
  auto iter = optical_model_table_.find(nucleus_pid);

  if (iter == optical_model_table_.end()) {
    // The requested level density model wasn't found, so create it and add
    // it to the table, returning a reference to the stored level density
    // model afterwards.
    int Z = marley_utils::get_particle_Z(nucleus_pid);
    int A = marley_utils::get_particle_A(nucleus_pid);
    return *(optical_model_table_.emplace(nucleus_pid,
      std::make_unique<marley::SphericalOpticalModel>(Z, A)).first
      ->second.get());
  }
  else return *(iter->second.get());
}

marley::SphericalOpticalModel& marley::StructureDatabase::get_optical_model(
  const int Z, const int A)
{
  int nucleus_pid = marley_utils::get_nucleus_pid(Z, A);
  auto iter = optical_model_table_.find(nucleus_pid);

  if (iter == optical_model_table_.end()) {
  // The requested level density model wasn't found, so create it and add
  // it to the table, returning a reference to the stored level density
  // model afterwards.
    return *(optical_model_table_.emplace(nucleus_pid,
      std::make_unique<marley::SphericalOpticalModel>(Z, A)).first
      ->second.get());
  }
  else return *(iter->second.get());
}

marley::LevelDensityModel& marley::StructureDatabase::get_level_density_model(
  const int Z, const int A)
{
  int pid = marley_utils::get_nucleus_pid(Z, A);

  auto iter = level_density_table_.find(pid);

  if (iter == level_density_table_.end()) {
  // The requested level density model wasn't found, so create it and add
  // it to the table, returning a reference to the stored level density
  // model afterwards.
    return *(level_density_table_.emplace(pid,
      std::make_unique<marley::BackshiftedFermiGasModel>(Z, A)).first
      ->second.get());
  }
  else return *(iter->second.get());
}

void marley::StructureDatabase::remove_decay_scheme(int pdg)
{
  // Remove the decay scheme with this PDG code if it exists in the database.
  // If it doesn't, do nothing.
  decay_scheme_table_.erase(pdg);
}

void marley::StructureDatabase::clear() {
  decay_scheme_table_.clear();
}
