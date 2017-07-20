#include "marley/marley_utils.hh"
#include "marley/ExitChannel.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/OpticalModel.hh"

void marley::FragmentDiscreteExitChannel::do_decay(double& Ex, int& two_J,
  marley::Parity& Pi, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
{
  Ex = final_level_.energy();
  two_J = final_level_.twoJ();
  Pi = final_level_.parity();
  emitted_particle = marley::Particle(fragment_.get_pid(),
    fragment_.get_mass());
  residual_nucleus = residue_;
}

void marley::GammaDiscreteExitChannel::do_decay(double& Ex, int& two_J,
  marley::Parity& Pi, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
{
  Ex = final_level_.energy();
  two_J = final_level_.twoJ();
  Pi = final_level_.parity();
  emitted_particle = marley::Particle(marley_utils::PHOTON, 0.);
  residual_nucleus = residue_;
}

void marley::FragmentContinuumExitChannel::do_decay(double& Ex,
  int& two_J, marley::Parity& Pi, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& gen)
{
  double Ea;
  Ex = gen.rejection_sample([&Ea, this](double ex)
    -> double { return this->Epdf_(Ea, ex); }, Emin_, Emax_);

  sample_spin_parity(two_J, Pi, gen, Ex, Ea);

  emitted_particle = marley::Particle(fragment_.get_pid(),
    fragment_.get_mass());

  residual_nucleus = gs_residue_;
  double rn_mass = gs_residue_.mass() + Ex;
  residual_nucleus.set_mass(rn_mass);
}

void marley::FragmentContinuumExitChannel::sample_spin_parity(int& twoJ,
  marley::Parity& Pi, marley::Generator& gen, double Exf, double fragment_KE)
{
  // Clear any previous table entries of spin-parities and decay widths
  jpi_widths_table_.clear();

  int two_s = fragment_.get_two_s();
  marley::Parity Pa = fragment_.get_parity();

  marley::OpticalModel& om
    = gen.get_structure_db().get_optical_model(gs_residue_.pdg_code());
  int Zf = om.Z();
  int Af = om.A();

  marley::LevelDensityModel& ldm
    = gen.get_structure_db().get_level_density_model(Zf, Af);

  double continuum_width = 0.;

  // Final nuclear state parity
  marley::Parity Pf;
  // The orbital parity starts as (-1)^0 = 1. Rather than applying parity
  // conservation each time, just find the final state parity Pf for l = 0.
  // Then we can safely flip Pf without further thought for each new l value in
  // the loop.
  if (Pi == Pa) Pf = 1;
  else Pf = -1;
  // For each new iteration, increment l and flip the final-state parity
  for (int l = 0; l <= HauserFeshbachDecay::l_max_; ++l, !Pf) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      for (int twoJf = std::abs(twoJ - two_j);
        twoJf <= twoJ + two_j; twoJf += 2)
      {
        // TODO: since only the spin distribution changes in the twoJf loop,
        // you can optimize this by precomputing most of the level density
        // and the multiplying here by the appropriate spin distribution.
        double width = om.transmission_coefficient(fragment_KE,
          fragment_.get_pid(), two_j, l, two_s)
          * ldm.level_density(Exf, twoJf, Pf);

        // Store the computed decay width for sampling
        continuum_width += width;
        jpi_widths_table_.emplace_back(twoJf, Pf, width);
      }
    }
  }
  // Throw an error if all decays are impossible
  if (continuum_width <= 0.) throw marley::Error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial fragment decay widths are "
    + "zero. cw = " + std::to_string(continuum_width));

  // Sample a final spin and parity
  const auto begin = marley::IteratorToMember<
    std::vector<SpinParityWidth>::const_iterator, SpinParityWidth,
    const double>(jpi_widths_table_.cbegin(),
    &SpinParityWidth::width);
  const auto end = marley::IteratorToMember<
    std::vector<SpinParityWidth>::const_iterator, SpinParityWidth,
    const double>(jpi_widths_table_.cend(),
    &SpinParityWidth::width);
  std::discrete_distribution<size_t> jpi_dist(begin, end);
  size_t jpi_index = gen.sample_from_distribution(jpi_dist);

  // Store the results
  SpinParityWidth Jpi = jpi_widths_table_.at(jpi_index);
  twoJ = Jpi.twoJf;
  Pi = Jpi.Pf;
}

void marley::GammaContinuumExitChannel::do_decay(double& Ex, int& two_J,
  marley::Parity& Pi, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus, marley::Generator& gen)
{
  double Exi = Ex;
  Ex = gen.rejection_sample(Epdf_, Emin_, Emax_);

  int nuc_pid = residual_nucleus.pdg_code();
  int Z = marley_utils::get_particle_Z(nuc_pid);
  int A = marley_utils::get_particle_A(nuc_pid);

  sample_spin_parity(Z, A, two_J, Pi, Exi, Ex, gen);

  emitted_particle = marley::Particle(marley_utils::PHOTON, 0.);
  residual_nucleus = gs_residue_;
  double rn_mass = gs_residue_.mass() + Ex;
  residual_nucleus.set_mass(rn_mass);
}

// Based on some initial parameters (including the nuclear 2J and parity)
// determine the final spin and parity and store it in twoJ and Pi
void marley::GammaContinuumExitChannel::sample_spin_parity(int Z, int A,
  int& twoJ, marley::Parity& Pi, double Exi, double Exf, marley::Generator& gen)
{
  using TrType = marley::GammaStrengthFunctionModel::TransitionType;

  // Clear any previous table entries of spin-parities and decay widths
  jpi_widths_table_.clear();

  double continuum_width = 0.;
  marley::Parity Pf;

  marley::LevelDensityModel& ldm
    = gen.get_structure_db().get_level_density_model(Z, A);

  marley::GammaStrengthFunctionModel& gsfm
    = gen.get_structure_db().get_gamma_strength_function_model(Z, A);

  // Approximate the gamma energy by the energy difference between the two
  // levels
  // TODO: consider adding a nuclear recoil correction here
  double e_gamma = Exi - Exf;
  bool initial_spin_is_zero = twoJ == 0;
  for (int l = 0; l <= marley::HauserFeshbachDecay::l_max_; ++l) {
    int two_l = 2*l;
    int twoJf = twoJ + two_l;
    // Determine the multipolarity being considered in this trip through the
    // loop based on the value of the orbital angular momentum l.
    int mpol = ((l == 0) ? 1 : l);
    // Consider both transition types (equivalently, both final parities). Check
    // for Ji = 0 to Jf = 0 transitions, which are not allowed by angular
    // momentum conservation.
    if (!initial_spin_is_zero || (twoJf > 0)) {

      double tcE = gsfm.transmission_coefficient(TrType::electric, mpol,
        e_gamma);
      double tcM = gsfm.transmission_coefficient(TrType::magnetic, mpol,
        e_gamma);

      if (!initial_spin_is_zero) {

        continuum_width += store_gamma_jpi_width(Exf, twoJf, Pi, tcE, tcM,
          mpol, ldm);

	// Consider the other possible final spin value for this multipolarity
	// if it is allowed (Jf = Ji - l is positive)
        twoJf = twoJ - two_l;
        if (twoJf >= 0) continuum_width += store_gamma_jpi_width(Exf, twoJf,
          Pi, tcE, tcM, mpol, ldm);
      }
      // twoJf > 0 but twoJi == 0
      else continuum_width += store_gamma_jpi_width(Exf, twoJf, Pi, tcE,
        tcM, mpol, ldm);
    }
  }

  // Throw an error if all decays are impossible
  if (continuum_width <= 0.) throw marley::Error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial gamma decay widths are "
    + "zero.");

  // Sample a final spin and parity
  const auto begin = marley::IteratorToMember<
    std::vector<SpinParityWidth>::const_iterator, SpinParityWidth,
    const double>(jpi_widths_table_.cbegin(),
    &SpinParityWidth::width);
  const auto end = marley::IteratorToMember<
    std::vector<SpinParityWidth>::const_iterator, SpinParityWidth,
    const double>(jpi_widths_table_.cend(),
    &SpinParityWidth::width);
  std::discrete_distribution<size_t> jpi_dist(begin, end);
  size_t jpi_index = gen.sample_from_distribution(jpi_dist);

  // Store the results
  SpinParityWidth Jpi = jpi_widths_table_.at(jpi_index);
  twoJ = Jpi.twoJf;
  Pi = Jpi.Pf;
}

// Helper function used when sampling continuum spin-parities for gamma-ray
// transitions
double marley::GammaContinuumExitChannel::store_gamma_jpi_width(double Exf,
  int twoJf, marley::Parity Pi, double tcE, double tcM, int mpol,
  marley::LevelDensityModel& ldm)
{
  marley::Parity Pf;
  double combined_width = 0.;

  // Electric transitions represent a parity flip for odd multipolarities
  if (mpol % 2) Pf = -Pi;
  else Pf = Pi;
  double rho = ldm.level_density(Exf, twoJf, Pf);

  // Compute and store information for the electric transition
  double width = tcE * rho;
  combined_width += width;

  jpi_widths_table_.emplace_back(twoJf, Pf, width);

  // Magnetic transitions have opposite final parity from electric
  // transitions of the same multipolarity
  !Pf;
  // Compute and store information for the magnetic transition
  rho = ldm.level_density(Exf, twoJf, Pf);
  width = tcM * rho;
  combined_width += width;

  jpi_widths_table_.emplace_back(twoJf, Pf, width);

  return combined_width;
}
