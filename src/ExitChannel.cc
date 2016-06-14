#include "ExitChannel.hh"
#include "HauserFeshbachDecay.hh"

//TODO: fix this!
static constexpr int l_max = 2;
static constexpr double DEFAULT_NUMEROV_STEP_SIZE = 0.1;

void marley::FragmentContinuumExitChannel::sample_spin_parity(int& twoJ,
  marley::Parity& Pi, const marley::Fragment& f,
  marley::SphericalOpticalModel& om, marley::Generator& gen,
  double Exf, double Ea)
{
  int fragment_pid = f.get_pid();
  int two_s = f.get_two_s();
  marley::Parity Pa = f.get_parity();

  std::vector<double> widths;
  std::vector<int> twoJfs;
  std::vector<marley::Parity> Pfs;

  double continuum_width = 0;
  int Zf = om.get_Z();
  int Af = om.get_A();

  marley::LevelDensityModel& ldm
    = gen.get_structure_db().get_level_density_model(Zf, Af);

  // Final nuclear state parity
  marley::Parity Pf;
  // The orbital parity starts as (-1)^0 = 1. Rather than applying parity
  // conservation each time, just find the final state parity Pf for l = 0.
  // Then we can safely flip Pf without further thought for each new l value in
  // the loop.
  if (Pi == Pa) Pf = 1;
  else Pf = -1;
  // For each new iteration, increment l and flip the final-state parity
  for (int l = 0; l <= l_max; ++l, !Pf) {
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
        double width = om.transmission_coefficient(Ea, fragment_pid, two_j, l,
          two_s, DEFAULT_NUMEROV_STEP_SIZE) * ldm.level_density(Exf, twoJf, Pf);

        // Store the computed decay width for sampling
        continuum_width += width;
        widths.push_back(width);
        twoJfs.push_back(twoJf);
        Pfs.push_back(Pf);
      }
    }
  }
  // Throw an error if all decays are impossible
  if (continuum_width <= 0.) throw marley::Error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial fragment decay widths are "
    + "zero. cw = " + std::to_string(continuum_width));

  // Sample a final spin and parity
  std::discrete_distribution<size_t> jpi_dist(widths.begin(), widths.end());
  size_t jpi_index = gen.discrete_sample(jpi_dist);

  // Store the results
  twoJ = twoJfs.at(jpi_index);
  Pi = Pfs.at(jpi_index);
}

// Based on some initial parameters (including the nuclear 2J and parity)
// determine the final spin and parity and store it in twoJ and Pi
void marley::GammaContinuumExitChannel::sample_spin_parity(int Z, int A,
  int& twoJ, marley::Parity& Pi, double Exi, double Exf, marley::Generator& gen)
{
  using HFD = marley::HauserFeshbachDecay;
  using TT = HFD::TransitionType;

  std::vector<int> twoJfs;
  std::vector<marley::Parity> Pfs;
  std::vector<double> widths;

  double continuum_width = 0;
  marley::Parity Pf;

  marley::LevelDensityModel& ldm
    = gen.get_structure_db().get_level_density_model(Z, A);
  // Approximate the gamma energy by the energy difference between the two levels
  // TODO: consider adding a nuclear recoil correction here
  double e_gamma = Exi - Exf;
  bool initial_spin_is_zero = twoJ == 0;
  for (int l = 0; l <= l_max; ++l) {
    int two_l = 2*l;
    int twoJf = twoJ + two_l;
    // Determine the multipolarity being considered in this trip through the
    // loop based on the value of the orbital angular momentum l.
    int mpol = ((l == 0) ? 1 : l);
    // Consider both transition types (equivalently, both final parities). Check
    // for Ji = 0 to Jf = 0 transitions, which are not allowed by angular
    // momentum conservation.
    if (!initial_spin_is_zero || (twoJf > 0)) {

      double tcE = HFD::gamma_transmission_coefficient(Z, A, TT::electric,
        mpol, e_gamma);
      double tcM = HFD::gamma_transmission_coefficient(Z, A, TT::magnetic,
        mpol, e_gamma);

      if (!initial_spin_is_zero) {

        continuum_width += store_gamma_pws(Exf, twoJf, Pi, widths, twoJfs,
          Pfs, tcE, tcM, mpol, ldm);

	// Consider the other possible final spin value for this multipolarity
	// if it is allowed (Jf = Ji - l is positive)
        twoJf = twoJ - two_l;
        if (twoJf >= 0) continuum_width += store_gamma_pws(Exf, twoJf,
          Pi, widths, twoJfs, Pfs, tcE, tcM, mpol, ldm);
      }
      // twoJf > 0 but twoJi == 0
      else continuum_width += store_gamma_pws(Exf, twoJf, Pi, widths,
        twoJfs, Pfs, tcE, tcM, mpol, ldm);
    }
  }
  // Throw an error if all decays are impossible
  if (continuum_width <= 0.) throw marley::Error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial gamma decay widths are "
    + "zero.");

  // Sample a final spin and parity
  std::discrete_distribution<size_t> jpi_dist(widths.begin(), widths.end());
  size_t jpi_index = gen.discrete_sample(jpi_dist);

  // Store the results
  twoJ = twoJfs.at(jpi_index);
  Pi = Pfs.at(jpi_index);
}

// Helper function used when sampling continuum spin-parities for gamma-ray
// transitions
double marley::GammaContinuumExitChannel::store_gamma_pws(double Exf, int twoJf,
  marley::Parity Pi, std::vector<double>& widths, std::vector<int>& twoJfs,
  std::vector<marley::Parity>& Pfs, double tcE, double tcM, int mpol,
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

  twoJfs.push_back(twoJf);
  Pfs.push_back(Pf);
  widths.push_back(width);

  // Magnetic transitions have opposite final parity from electric
  // transitions of the same multipolarity
  !Pf;
  // Compute and store information for the magnetic transition
  rho = ldm.level_density(Exf, twoJf, Pf);
  width = tcM * rho;
  combined_width += width;

  twoJfs.push_back(twoJf);
  Pfs.push_back(Pf);
  widths.push_back(width);

  return combined_width;
}
