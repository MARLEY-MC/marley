#include <limits>

#include "marley/Generator.hh"
#include "marley/NeutrinoSource.hh"

marley::NeutrinoSource::NeutrinoSource(int particle_id) {
  if (!pdg_is_allowed(particle_id)) throw marley::Error(
    "Creating a neutrino source object that produces"
    " particles with PDG ID number " + std::to_string(particle_id)
    + " is not allowed.");
  else pid_ = particle_id;
}

/// PDG particle IDs of each neutrino that could possibly be produced by a
/// NeutrinoSource object
const std::set<int> marley::NeutrinoSource::pids_ = {
  marley_utils::ELECTRON_NEUTRINO,
  marley_utils::ELECTRON_ANTINEUTRINO,
  marley_utils::MUON_NEUTRINO,
  marley_utils::MUON_ANTINEUTRINO,
  marley_utils::TAU_NEUTRINO,
  marley_utils::TAU_ANTINEUTRINO
};

double marley::NeutrinoSource::sample_incident_neutrino(int& pdg,
  marley::Generator& gen)
{
  static double max = marley_utils::UNKNOWN_MAX;
  pdg = pid_;
  return gen.rejection_sample([this](double E)
    -> double { return this->pdf(E); }, get_Emin(), get_Emax(), max);
}

marley::FermiDiracNeutrinoSource::FermiDiracNeutrinoSource(int particle_id,
  double Emin, double Emax, double temp, double eta)
  : NeutrinoSource(particle_id), Emin_(Emin), Emax_(Emax), temperature_(temp),
  eta_(eta), C_(1.)
{
  // Normalize the source spectrum (not strictly necessary, but having the
  // spectrum approximately normalized makes the default rejection sampling
  // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
  double integral = marley_utils::num_integrate(
    [this](double E) -> double { return this->pdf(E); }, Emin_, Emax_);

  // Update the normalization constant, thereby normalizing this object's
  // pdf in the process.
  C_ /= integral;
}

double marley::FermiDiracNeutrinoSource::pdf(double E) {
  if (E < Emin_ || E > Emax_) return 0.;
  else return (C_ / std::pow(temperature_, 4)) * (std::pow(E, 2)
    / (1 + std::exp((E / temperature_) - eta_)));
}

marley::BetaFitNeutrinoSource::BetaFitNeutrinoSource(int particle_id,
  double Emin, double Emax, double Emean, double beta)
  : NeutrinoSource(particle_id), Emin_(Emin), Emax_(Emax), Emean_(Emean),
  beta_(beta), C_(1.)
{
  // Normalize the source spectrum (not strictly necessary, but having the
  // spectrum approximately normalized makes the default rejection sampling
  // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
  double integral = marley_utils::num_integrate(
    [this](double E) -> double { return this->pdf(E); }, Emin_, Emax_);

  // Update the normalization constant, thereby normalizing this object's
  // pdf in the process.
  C_ /= integral;
}

double marley::BetaFitNeutrinoSource::pdf(double E) {
  if (E < Emin_ || E > Emax_) return 0.;
  else return C_ * std::pow(E / Emean_, beta_ - 1)
    * std::exp(-beta_ * E / Emean_);
}

marley::FunctionNeutrinoSource::FunctionNeutrinoSource(int particle_id,
  double Emin, double Emax, std::function<double(double)> prob_dens_func)
  : NeutrinoSource(particle_id), Emin_(Emin), Emax_(Emax)
{
  // Normalize the supplied spectrum (not strictly necessary, but having the
  // spectrum approximately normalized makes the default rejection sampling
  // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
  double integral = marley_utils::num_integrate(prob_dens_func, Emin, Emax);
  probability_density_ = [prob_dens_func, integral](double E)
    -> double { return prob_dens_func(E) / integral; };
}

marley::DecayAtRestNeutrinoSource::DecayAtRestNeutrinoSource(int particle_id)
  : NeutrinoSource(particle_id)
{
  if (particle_id != marley_utils::ELECTRON_NEUTRINO &&
    particle_id != marley_utils::MUON_ANTINEUTRINO)
  {
    throw marley::Error(std::string("Decay at rest")
      + " neutrino source objects may only produce electron neutrinos"
      + " or muon antineutrinos. PDG particle ID number "
      + std::to_string(particle_id) + " is therefore not allowed.");
  }
}

double marley::DecayAtRestNeutrinoSource::pdf(double E) {
  if (E < Emin_ || E > Emax_) return 0.;
  // Note that both of these source spectra are normalized to 1
  // on the energy interval [0., m_mu_ / 2.]
  /// @todo Refine the approximate decay-at-rest Michel spectra to use more
  /// exact expressions.
  else if (pid_ == marley_utils::ELECTRON_NEUTRINO)
    return 96. * std::pow(E, 2) * m_mu_to_the_minus_four_
      * (m_mu_ - 2*E);
  // Spectrum for muon antineutrinos
  else return 16. * std::pow(E, 2) * m_mu_to_the_minus_four_
    * (3*m_mu_ - 4*E);
}

/// @brief Method called near the end of construction to verify that the
/// newly-created grid source object is valid.
void marley::GridNeutrinoSource::check_for_errors() {
  size_t grid_size = grid_.size();
  if (grid_size < 2) throw marley::Error(std::string("Grid with")
    + " less than 2 gridpoints passed to the constructor of"
    + " marley::GridNeutrinoSource.");

  double sum_of_PDs = 0.;
  for (size_t j = 0; j < grid_size; ++j) {
    auto& pair = grid_.at(j);
    if (pair.first < 0.) throw marley::Error(std::string("All energy")
      + " values used in a marley::GridNeutrinoSource"
      + " object must be nonnegative");

    // Prevent actually sampling an energy value of zero by advancing to
    // the next representable double value.
    else if (pair.first == 0.) pair.first = std::nextafter(0.,
      marley_utils::infinity);
    if (pair.second < 0.) throw marley::Error(std::string("All PDF")
      + " values used in a marley::GridNeutrinoSource object must be"
      + " nonnegative");
    else sum_of_PDs += pair.second;
  }

  if (sum_of_PDs <= 0.) throw marley::Error(std::string("All probability")
    + " density grid point values are zero for the neutrino source");
}
