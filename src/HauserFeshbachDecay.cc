#include <memory>

#include "ExitChannel.hh"
#include "Generator.hh"
#include "Kinematics.hh"
#include "MassTable.hh"
#include "HauserFeshbachDecay.hh"

using TrType = marley::HauserFeshbachDecay::TransitionType;

// Table of nuclear fragments that will be considered when computing
// branching ratios for nuclear de-excitations. Spin-parity values are taken
// from nuclear ground states listed in the 10/2014 release of ENSDF.
const std::vector<marley::Fragment> marley::HauserFeshbachDecay::fragments = {
  marley::Fragment(marley_utils::NEUTRON, 1, marley::Parity(1)),
  marley::Fragment(marley_utils::PROTON, 1, marley::Parity(1)),
  marley::Fragment(marley_utils::DEUTERON, 2, marley::Parity(1)),
  marley::Fragment(marley_utils::TRITON, 1, marley::Parity(1)),
  marley::Fragment(marley_utils::HELION, 1, marley::Parity(1)),
  marley::Fragment(marley_utils::ALPHA, 0, marley::Parity(1)),
};

marley::HauserFeshbachDecay::HauserFeshbachDecay(const marley::Particle&
  compound_nucleus, double Exi, int twoJi, marley::Parity Pi,
  marley::Generator& gen) : compound_nucleus_(compound_nucleus), Exi_(Exi),
  twoJi_(twoJi), Pi_(Pi), gen_(gen)
{
  build_exit_channels();
}

void marley::HauserFeshbachDecay::build_exit_channels()
{
  static const double me
    = marley::MassTable::get_particle_mass(marley_utils::ELECTRON);

  int pid_initial = compound_nucleus_.pdg_code();
  int Zi = marley::MassTable::get_particle_Z(pid_initial);
  int Ai = marley::MassTable::get_particle_A(pid_initial);
  int qi = compound_nucleus_.charge(); // Get net charge of initial ion
  double Mi = compound_nucleus_.mass(); // initial mass
  double Migs = Mi - Exi_; // initial ground-state mass

  total_width_ = 0.; // total compound nucleus decay width

  marley::StructureDatabase& db = gen_.get_structure_db();

  for (const auto& f : fragments) {

    // Get information about the current fragment
    int two_s = f.get_two_s(); // spin
    marley::Parity Pa = f.get_parity();
    int fragment_pid = f.get_pid();
    int Za = f.get_Z(); // atomic number
    double Ma = f.get_mass(); // mass

    // Get information about the final-state nucleus
    int Zf = Zi - f.get_Z(); // atomic number
    int Af = Ai - f.get_A(); // mass number
    int pid_final = marley_utils::get_nucleus_pid(Zf, Af);

    // Approximate the ground state mass of the ion formed when the fragment f
    // is emitted by adding (Za - qi) electron masses to the atomic mass for
    // the final nucleus.
    double Mfgs_ion = marley::MassTable::get_atomic_mass(pid_final)
      + (Za - qi)*me;

    // Get fragment separation energy without a redundant mass table lookup
    double Sa = Mfgs_ion + Ma - Migs;

    // Get discrete level data (if any) and models for the final nucleus
    marley::DecayScheme* ds = db.get_decay_scheme(pid_final);
    marley::SphericalOpticalModel& om = db.get_optical_model(pid_final);

    // Determine the maximum excitation energy available after fragment emission
    // in the final nucleus. Use a simpler functional form as a shortcut if
    // the fragment is neutral [and can therefore be emitted with arbitrarily
    // low kinetic energy due to the absence of a Coulomb barrier])
    double Exf_max;
    double Vc; // Coulomb barrier
    if (f.get_Z() == 0) {
      Vc = 0;
      Exf_max = Exi_ - Sa;
    }
    else {
      Vc = coulomb_barrier(Zf, Af, Za, f.get_A()); // Ea_min = Vc
      Exf_max = std::sqrt(std::pow(Mi - Ma, 2) - 2*Vc*Mi) - Mfgs_ion;
    }

    // The fragment kinetic energy depends on the final level energy. For speed,
    // we precompute here a portion of the fragment energy that doesn't depend on
    // the final level energy.
    double Mconst = (Exi_ - Sa) * (Mi + Mfgs_ion - Ma);

    // Check if emission of this fragment is energetically allowed. If we're
    // below the fragment emission threshold, the partial decay width is
    // approximately zero (neglecting tunneling), so we can skip this fragment
    // rather than slogging through the calculation.  We'll also skip it if
    // we're exactly at threshold to avoid numerical problems.
    // TODO: consider changing this rule to allow tunneling through the Coulomb
    // barrier
    if (Mconst / (2 * Mi) - Vc <= 0.) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0.;

    // If discrete level data is available for the final nuclide, get decay
    // widths for each accessible level
    if (ds) {

      // Get a vector of pointers to levels in the decay scheme. The levels are
      // sorted in order of increasing excitation energy.
      const auto sorted_lps = ds->get_sorted_level_pointers();

      // Use the maximum discrete level energy from the decay scheme object as the
      // lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (sorted_lps.size() > 0) E_c_min = sorted_lps.back()->get_energy();

      // Loop over the final discrete nuclear levels in order of increasing energy
      // until the new level energy exceeds the maximum value. For each
      // energetically allowed level, if a transition to it for a given fragment
      // orbital angular momentum l and total angular momentum j conserves parity,
      // then compute an optical model transmission coefficient and add it to the total.
      for (const auto& level : sorted_lps) {
        double Exf = level->get_energy();
        if (Exf < Exf_max) {
          double discrete_width = 0.;
          double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
          int twoJf = level->get_twoJ();
          marley::Parity Pf = level->get_parity();
          for (int two_j = std::abs(twoJi_ - twoJf); two_j <= twoJi_ + twoJf;
            two_j += 2)
          {
            int j_plus_s = (two_j + two_s) / 2;
            // TODO: consider adding a check that l does not exceed its maximum
            // value l_max that is also used as a cutoff for transmission
            // coefficient calculations in the continuum.
            // For each new iteration, increment l and flip the overall final state parity
            int l = std::abs(two_j - two_s) / 2;
            bool l_is_odd = l % 2;
            marley::Parity P_final_state = Pf * Pa * marley::Parity(!l_is_odd);
            for (; l <= j_plus_s; ++l, !P_final_state)
            {
              // Add to the total transmission coefficient if parity is conserved
              if (Pi_ == P_final_state) discrete_width
                += om.transmission_coefficient(Ea, fragment_pid, two_j, l, two_s,
                DEFAULT_NUMEROV_STEP_SIZE);
            }
          }

          // Store information for this decay channel
	  exit_channels_.push_back(std::make_unique<marley::FragmentDiscreteExitChannel>(
            discrete_width, *level, marley::Particle(pid_final, Mfgs_ion + Exf,
            qi - f.get_Z()), f));
          total_width_ += discrete_width;
        }
        else break;
      }
    }

    // If transitions to the energy continuum are possible, include the
    // continuum in the decay channels
    if (Exf_max > E_c_min) {

      marley::LevelDensityModel& ldm = db.get_level_density_model(Zf, Af);

      // Create a function object for the continuum partial width member
      // function that takes a single argument.
      std::function<double(double)> cpw = [=, &om, &ldm](double Exf) -> double {
        double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
        return fragment_continuum_partial_width(om, ldm, f, Ea, Exf);
      };

      // Numerically integrate using the call wrapper, the integration bounds, and
      // the number of subintervals
      double continuum_width = marley_utils::num_integrate(cpw, E_c_min,
        Exf_max, DEFAULT_CONTINUUM_SUBINTERVALS);

      // Create a normalized probability density function to use for sampling
      // a final excitation energy from the continuum for this exit channel.
      std::function<double(double&, double)> pdf
        = [=, &f, &om, &ldm](double& Ea, double Exf) -> double {
        Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
        return fragment_continuum_partial_width(om, ldm, f, Ea, Exf)
          / continuum_width;
      };

      // Store information for this decay channel
      exit_channels_.push_back(
        std::make_unique<marley::FragmentContinuumExitChannel>(continuum_width,
        E_c_min, Exf_max, pdf, f, marley::Particle(pid_final, Mfgs_ion,
        qi - f.get_Z())));
      total_width_ += continuum_width;
    }
  }

  marley::DecayScheme* ds = db.get_decay_scheme(pid_initial);

  // Let the continuum go down to 0 MeV unless there is a decay scheme object
  // available for the final nuclide (we'll check this in a second).
  double E_c_min = 0.;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if (ds) {

    bool initial_spin_is_zero = twoJi_ == 0;
    // Loop over the final discrete nuclear levels in order of increasing energy
    // until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission coefficient
    // for it
    const auto sorted_lps = ds->get_sorted_level_pointers();

    // Use the maximum discrete level energy from the decay scheme object as the
    // lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if (sorted_lps.size() > 0) E_c_min = sorted_lps.back()->get_energy();

    for (const auto& level_f : sorted_lps) {
      int twoJf = level_f->get_twoJ();
      // 0->0 EM transitions aren't allowed due to angular momentum conservation
      // (photons are spin 1), so if the initial and final spins are zero, skip
      // ahead to the next final level.
      if (initial_spin_is_zero && twoJf == 0) continue;
      double Exf = level_f->get_energy();
      if (Exf < Exi_) {
        // Approximate the gamma energy by the energy difference between the two levels
        // TODO: consider adding a nuclear recoil correction here
        double e_gamma = Exi_ - Exf;
        // Determine the type (electric or magnetic) and multipolarity of the gamma
        // transition between these two levels
        int mpol; // Multipolarity
        auto type = determine_gamma_transition_type(twoJi_, Pi_, level_f, mpol);
        // TODO: allow the user to choose which gamma-ray transition model to use
        // (Weisskopf single-particle estimates, Brink-Axel strength functions, etc.)
        double discrete_width = gamma_transmission_coefficient(Zi, Ai, type,
          mpol, e_gamma);

        // Store information for this decay channel
        exit_channels_.push_back(
          std::make_unique<marley::GammaDiscreteExitChannel>(discrete_width,
          *level_f, marley::Particle(pid_initial, Migs + Exf, qi)));
        total_width_ += discrete_width;
      }
      else break;
    }
  }

  // If gamma transitions to the energy continuum are possible, include them
  // in the possible decay channels
  if (Exi_ > E_c_min) {

    marley::LevelDensityModel& ldm = db.get_level_density_model(Zi, Ai);

    // Numerically integrate the continuum width up to the initial excitation
    // energy
    double continuum_width = marley_utils::num_integrate(
      [&ldm, this](double Exf)
      -> double { return gamma_continuum_partial_width(ldm, Exf); }, E_c_min,
      Exi_, DEFAULT_CONTINUUM_SUBINTERVALS);

    // Normalized probability density used for sampling a final excitation
    // energy in the continuum
    std::function<double(double)> pdf = [=, &ldm](double Exf) -> double {
      return gamma_continuum_partial_width(ldm, Exf) / continuum_width;
    };

    // Store information for this decay channel
    exit_channels_.push_back(
      std::make_unique<marley::GammaContinuumExitChannel>(continuum_width,
      E_c_min, Exi_, pdf, marley::Particle(pid_initial, Migs, qi)));
    total_width_ += continuum_width;
  }
}

bool marley::HauserFeshbachDecay::do_decay(double& Exf, int& twoJf,
  marley::Parity& Pf, marley::Particle& emitted_particle,
  marley::Particle& residual_nucleus)
{
  // Throw an error if all decays are impossible
  if (total_width_ <= 0.) throw marley::Error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial decay widths are zero.");

  // Sample an exit channel
  const auto widths_begin
    = marley::ExitChannel::make_width_iterator(exit_channels_.cbegin());
  const auto widths_end
    = marley::ExitChannel::make_width_iterator(exit_channels_.cend());

  std::discrete_distribution<size_t> exit_channel_dist(widths_begin,
    widths_end);
  size_t exit_channel_index = gen_.discrete_sample(exit_channel_dist);

  marley::ExitChannel* ec = exit_channels_.at(exit_channel_index).get();

  ec->do_decay(Exf, twoJf, Pf, emitted_particle, residual_nucleus, gen_);

  // TODO: consider changing this to a more realistic model
  // instead of isotropic emissions
  double cos_theta_emitted_particle = gen_.uniform_random_double(-1, 1, true);
  double phi_emitted_particle = gen_.uniform_random_double(0,
    marley_utils::two_pi, false);

  // Handle the kinematics calculations for this decay and update the
  // final-state particle objects
  marley::Kinematics::two_body_decay(compound_nucleus_, emitted_particle,
    residual_nucleus, cos_theta_emitted_particle, phi_emitted_particle);

  bool discrete_level = !ec->is_continuum();
  return !discrete_level;
}

TrType marley::HauserFeshbachDecay::determine_gamma_transition_type(int twoJi,
  marley::Parity Pi, marley::Level* level_f, int& l)
{
  int twoJf = level_f->get_twoJ();
  marley::Parity Pf = level_f->get_parity();

  return determine_gamma_transition_type(twoJi, Pi, twoJf, Pf, l);
}

TrType marley::HauserFeshbachDecay::determine_gamma_transition_type(
  marley::Level* level_i, marley::Level* level_f, int& l)
{
  int twoJi = level_i->get_twoJ();
  marley::Parity Pi = level_i->get_parity();

  int twoJf = level_f->get_twoJ();
  marley::Parity Pf = level_f->get_parity();

  return determine_gamma_transition_type(twoJi, Pi, twoJf, Pf, l);
}

// Returns whether a gamma transition from a nuclear state with spin twoJi / 2 and
// parity Pi to a state with spin twoJf / 2 and parity Pf is electric or magnetic.
// Also loads the integer l with the multipolarity of the transition.
TrType marley::HauserFeshbachDecay::determine_gamma_transition_type(int twoJi,
  marley::Parity Pi, int twoJf, marley::Parity Pf, int& l)
{
  // TODO: reconsider how to handle this situation
  if (twoJi == 0 && twoJf == 0) throw marley::Error(
    std::string("0 -> 0 EM transitions are not allowed."));

  int two_delta_J = std::abs(twoJf - twoJi);
  // Odd values of two_delta_J are unphysical because they represent EM
  // transitions where the total angular momentum changes by half (photons are
  // spin 1)
  if (two_delta_J % 2) throw marley::Error(std::string("Unphysical ")
    + "EM transition encountered between nuclear levels with spins 2*Ji = "
    + std::to_string(twoJi) + " and 2*Jf = " + std::to_string(twoJf));

  marley::Parity Pi_times_Pf = Pi * Pf;

  // Determine the multipolarity of this transition. Load l with the result.
  if (two_delta_J == 0) l = 1;
  // We already checked for unphysical odd values of two_delta_J above, so
  // using integer division here will always give the expected result.
  else l = two_delta_J / 2;

  // Determine whether this transition is electric or magnetic based on its
  // multipolarity and whether or not there is a change of parity

  // Pi * Pf = -1 gives electric transitions for odd l
  marley::Parity electric_parity;
  if (l % 2) electric_parity = -1; // l is odd
  // Pi * Pf = +1 gives electric transitions for even l
  else electric_parity = 1; // l is even

  if (Pi_times_Pf == electric_parity) return TransitionType::electric;
  else return TransitionType::magnetic;
}

double marley::HauserFeshbachDecay::gamma_strength_function_coefficient(int Z, int A,
  TransitionType type, int l, double e_gamma)
{
  // TODO: improve this error message
  if (l < 1) throw marley::Error(std::string("Invalid multipolarity")
    + std::to_string(l) + " given for gamma ray strength"
    + " function calculation");

  // The strength, energy, and width of the giant resonance for a transition of
  // type x (E or M) and multipolarity l
  double sigma_xl, e_xl, gamma_xl;

  if (type == TransitionType::electric) {
    if (l == 1) {
      e_xl = 31.2*std::pow(A, -1.0/3.0) + 20.6*std::pow(A, -1.0/6.0);
      gamma_xl = 0.026 * std::pow(e_xl, 1.91);
      sigma_xl = 1.2 * 120 * (A - Z) * Z/(A * marley_utils::pi * gamma_xl)
        * marley_utils::mb;
    }
    if (l > 1) {
      // Values for E2 transitions
      e_xl = 63*std::pow(A, -1.0/3.0);
      gamma_xl = 6.11 - 0.012*A;
      sigma_xl = 0.00014 * std::pow(Z, 2) * e_xl
        / (std::pow(A, 1.0/3.0) * gamma_xl) * marley_utils::mb;
      // If this is an E2 transition, we're done. Otherwise,
      // compute the giant resonance strength iteratively
      for (int i = 2; i < l; ++i) {
        sigma_xl *= 8e-4;
      }
    }
  }
  else if (type == TransitionType::magnetic) {
    // Values for M1 transitions
    // The commented-out version is for RIPL-1
    //double factor_m1 = 1.58e-9*std::pow(A, 0.47);
    // RIPL-2 factor
    const double e_gamma_ref = 7.0; // MeV
    double factor_m1 = gamma_strength_function(Z, A,
      TransitionType::electric, 1, e_gamma_ref)
      / (0.0588 * std::pow(A, 0.878));
    gamma_xl = 4.0;
    e_xl = 41*std::pow(A, -1.0/3.0);
    sigma_xl = (std::pow(std::pow(e_gamma_ref, 2) - std::pow(e_xl, 2), 2)
      + std::pow(e_gamma_ref, 2) * std::pow(gamma_xl, 2))
      * (3 * std::pow(marley_utils::pi, 2) * factor_m1)
      / (e_gamma_ref * std::pow(gamma_xl, 2));
    // If this is an M1 transition, we're done. Otherwise,
    // compute the giant resonance strength iteratively
    for (int i = 1; i < l; ++i) {
      sigma_xl *= 8e-4;
    }
  }
  // TODO: improve this error message
  else throw marley::Error(std::string("Invalid transition type")
    + " given for gamma ray strength function calculation");

  // Now that we have the appropriate giant resonance parameters,
  // calculate the strength function using the Brink-Axel expression.
  // Note that the strength function has units of MeV^(-3)
  double coeff = (sigma_xl * std::pow(gamma_xl, 2)) / ((2*l + 1)
    * std::pow(marley_utils::pi, 2) * (std::pow(std::pow(e_gamma, 2)
    - std::pow(e_xl, 2), 2) + std::pow(e_gamma, 2) * std::pow(gamma_xl, 2)));

  return coeff;
}

// Computes the partial decay width for a gamma transition under the Weisskopf
// single-particle approximation.
double marley::HauserFeshbachDecay::weisskopf_partial_decay_width(int A,
  TransitionType type, int l, double e_gamma)
{
  // Compute double factorial of 2l + 1
  int dfact = 1;
  for (int n = 2*l + 1; n > 0; n -= 2) dfact *= n;

  // Multipolarity factor
  double lambda = (l + 1) / (l * std::pow(dfact, 2))
    * std::pow(3.0 / (l + 3), 2);

  // Estimated nuclear radius (fm)
  double R = marley_utils::r0 * std::pow(A, 1.0/3.0);

  // Electric transition partial decay width
  double el_width = 2 * marley_utils::alpha * lambda
    * std::pow(R * e_gamma / marley_utils::hbar_c, 2*l) * e_gamma;

  if (type == TransitionType::electric) {
    return el_width;
  }

  else if (type == TransitionType::magnetic) {
    double mp = marley::MassTable::get_particle_mass(marley_utils::PROTON);
    return 10 * el_width * std::pow(marley_utils::hbar_c / (mp * R), 2);
  }

  // TODO: improve this error message
  else throw marley::Error(std::string("Invalid transition type")
    + " given for Weisskopf gamma-ray partial width calculation");
}

double marley::HauserFeshbachDecay::gamma_continuum_partial_width(
  marley::LevelDensityModel& ldm, double Exf)
{
  double continuum_width = 0.;
  // Approximate the gamma energy by the energy difference between the two levels
  // TODO: consider adding a nuclear recoil correction here
  double e_gamma = Exi_ - Exf;
  bool initial_spin_is_zero = twoJi_ == 0;
  int cn_pid = compound_nucleus_.pdg_code();
  int Z = marley::MassTable::get_particle_Z(cn_pid);
  int A = marley::MassTable::get_particle_A(cn_pid);
  for (int l = 0; l <= l_max; ++l) {
    int two_l = 2*l;
    int twoJf = twoJi_ + two_l;
    // Determine the multipolarity being considered in this trip through the
    // loop based on the value of the orbital angular momentum l.
    int mpol = ((l == 0) ? 1 : l);
    // Consider both transition types (equivalently, both final parities). Check
    // for Ji = 0 to Jf = 0 transitions, which are not allowed by angular
    // momentum conservation.
    if (!initial_spin_is_zero) {
      continuum_width += gamma_cpw(Z, A, mpol, twoJf, e_gamma, ldm, Exf);
      // Consider the other possible final spin value for this multipolarity if
      // it is allowed (Jf = Ji - l is positive)
      twoJf = twoJi_ - two_l;
      if (twoJf >= 0) continuum_width += gamma_cpw(Z, A, mpol, twoJf, e_gamma,
        ldm, Exf);
    }
    else if (twoJf > 0) continuum_width += gamma_cpw(Z, A, mpol, twoJf, e_gamma,
      ldm, Exf);
  }
  return continuum_width;
}

double marley::HauserFeshbachDecay::fragment_continuum_partial_width(
  marley::SphericalOpticalModel& om, marley::LevelDensityModel& ldm,
  const marley::Fragment& frag, double fragment_KE, double Exf)
{
  int two_s = frag.get_two_s();
  marley::Parity Pa = frag.get_parity();
  double continuum_width = 0.;
  // Final nuclear state parity
  marley::Parity Pf;
  // The orbital parity starts as (-1)^0 = 1. Rather than applying parity
  // conservation each time, just find the final state parity Pf for l = 0.
  // Then we can safely flip Pf without further thought for each new l value in
  // the loop.
  if (Pi_ == Pa) Pf = 1;
  else Pf = -1;
  // For each new iteration, increment l and flip the final-state parity
  for (int l = 0; l <= l_max; ++l, !Pf) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      for (int twoJf = std::abs(twoJi_ - two_j);
        twoJf <= twoJi_ + two_j; twoJf += 2)
      {
        continuum_width += om.transmission_coefficient(fragment_KE,
          frag.get_pid(), two_j, l, two_s, DEFAULT_NUMEROV_STEP_SIZE)
          * ldm.level_density(Exf, twoJf);
          // TODO: since only the spin distribution changes in the twoJf loop,
          // you can optimize this by precomputing most of the level density
          // and the multiplying here by the appropriate spin distribution
      }
    }
  }
  return continuum_width;
}

double marley::HauserFeshbachDecay::get_fragment_emission_threshold(const int Zi,
  const int Ai, const marley::Fragment& f)
{
  int Zf = Zi - f.get_Z();
  int Af = Ai - f.get_A();

  // Get various masses that will be needed for fragment emission kinematics
  // (including nuclear recoil)
  double Migs = marley::MassTable::get_atomic_mass(Zi, Ai);
  double Ma = f.get_mass();
  double Mfgs = marley::MassTable::get_atomic_mass(Zf, Af);
  double me = marley::MassTable::get_particle_mass(marley_utils::ELECTRON);
  // Fragment atomic number
  int Za = f.get_Z();
  // Approximate the ground state rest energy of the negative ion (with
  // charge Za-) formed when the bare fragment f is emitted by adding Za
  // electron masses to the atomic mass for the final nucleus.
  double Mfgs_ion = Mfgs + Za*me;

  // Get Coulomb potential for this fragment
  double Vc;
  if (Za != 0) Vc = coulomb_barrier(Zf, Af, Za, f.get_A());
  else Vc = 0.;

  return Ma - Migs + Vc + marley_utils::real_sqrt(Vc * (2*Ma + Vc)
    + std::pow(Mfgs_ion, 2));
}

double marley::HauserFeshbachDecay::gamma_cpw(int Z, int A, int mpol, int twoJf,
  double e_gamma, marley::LevelDensityModel& ldm, double Exf)
{
  double tcE = gamma_transmission_coefficient(Z, A,
    TransitionType::electric, mpol, e_gamma);
  double tcM = gamma_transmission_coefficient(Z, A,
    TransitionType::magnetic, mpol, e_gamma);

  // Final parity for the electric transition will match
  // the initial parity for even multipolarities
  auto PfE = Pi_;
  bool mpol_is_odd = mpol % 2;
  if (mpol_is_odd) !PfE;

  double rhoE = ldm.level_density(Exf, twoJf, PfE);
  // The magnetic transition corresponds to the opposite parity
  double rhoM = ldm.level_density(Exf, twoJf, -PfE);
  return (tcE * rhoE) + (tcM * rhoM);
}

void marley::HauserFeshbachDecay::print(std::ostream& out) const {

  static constexpr double hbar = 6.58211951e-22; // MeV * s
  int pid_initial = compound_nucleus_.pdg_code();
  int Zi = marley::MassTable::get_particle_Z(pid_initial);
  int Ai = marley::MassTable::get_particle_A(pid_initial);
  marley::LevelDensityModel& ldm
    = gen_.get_structure_db().get_level_density_model(Zi, Ai);
  double two_pi_rho = marley_utils::two_pi
    * ldm.level_density(Exi_, twoJi_, Pi_);

  out << "Compound nucleus " << compound_nucleus_.pdg_code()
    << " with Ex = " << Exi_ << ", spin = " << twoJi_ / 2;
  if (twoJi_ % 2) out << ".5";
  out << ", and parity = " << Pi_ << '\n';
  out << "Total width = " << total_width_ / two_pi_rho << " MeV\n";
  out << "Mean lifetime = " << hbar / (total_width_ / two_pi_rho) << " s\n";
  for (const auto& ec : exit_channels_) {
    double width = ec->width() / two_pi_rho;
    bool continuum = ec->is_continuum();
    bool frag = ec->emits_fragment();
    if (continuum) {
      if (frag) {
        const auto* fcec
          = dynamic_cast<marley::FragmentContinuumExitChannel*>(ec.get());
        if (fcec) out << "  "
          << marley_utils::particle_symbols.at(fcec->get_fragment().get_pid())
          << " emission to the continuum width = " << width << " MeV\n";
      }
      else out << "  gamma-ray emission to the continuum width = "
        << width << " MeV\n";
    }
    else if (frag) {
      const auto* fdec
        = dynamic_cast<marley::FragmentDiscreteExitChannel*>(ec.get());
      if (fdec) out << "  "
        << marley_utils::particle_symbols.at(fdec->get_fragment().get_pid())
        << " emission to level at " << fdec->get_final_level().get_energy()
        << " MeV width = " << width << " MeV\n";
    }
    else {
      const auto* gdec
        = dynamic_cast<marley::GammaDiscreteExitChannel*>(ec.get());
      if (gdec) out << "  gamma-ray emission to level at "
        << gdec->get_final_level().get_energy() << " MeV width = "
        << width << " MeV\n";
    }
  }
}
