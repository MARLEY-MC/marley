#include <memory>

#include "marley/marley_utils.hh"
#include "marley/ExitChannel.hh"
#include "marley/Generator.hh"
#include "marley/marley_kinematics.hh"
#include "marley/MassTable.hh"
#include "marley/HauserFeshbachDecay.hh"

using TrType = marley::GammaStrengthFunctionModel::TransitionType;

namespace {

  // Helper function for
  // marley::HauserFeshbachDecay::gamma_continuum_partial_width(). It doesn't
  // really need access to any class members to do its job, so we've hidden it
  // in an anonymous namespace rather than including it as a private member.
  double gamma_cpw_helper(marley::LevelDensityModel& ldm,
    marley::GammaStrengthFunctionModel& gsfm, marley::Parity Pi, int mpol,
    double e_gamma, double Exf, int twoJf, bool use_TCs)
  {
    double tcE = 1.;
    double tcM = 1.;
    if ( use_TCs ) {
      tcE = gsfm.transmission_coefficient(TrType::electric, mpol, e_gamma);
      tcM = gsfm.transmission_coefficient(TrType::magnetic, mpol, e_gamma);
    }

    // Final parity for the electric transition will match
    // the initial parity for even multipolarities
    auto PfE = Pi;
    bool mpol_is_odd = mpol % 2;
    if (mpol_is_odd) !PfE;

    double rhoE = ldm.level_density(Exf, twoJf, PfE);
    // The magnetic transition corresponds to the opposite parity
    double rhoM = ldm.level_density(Exf, twoJf, -PfE);
    return (tcE * rhoE) + (tcM * rhoM);
  }
}

// Table of nuclear fragments that will be considered when computing
// branching ratios for nuclear de-excitations. Spin-parity values are taken
// from nuclear ground states listed in the 10/2014 release of ENSDF.
const std::vector<marley::Fragment> marley::HauserFeshbachDecay::fragments_ = {
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
  // Remove any pre-existing ExitChannel objects, just in case
  exit_channels_.clear();

  const marley::MassTable& mt = marley::MassTable::Instance();

  static const double me = mt.get_particle_mass(marley_utils::ELECTRON);

  int pid_initial = compound_nucleus_.pdg_code();
  int Zi = marley_utils::get_particle_Z(pid_initial);
  int Ai = marley_utils::get_particle_A(pid_initial);
  int qi = compound_nucleus_.charge(); // Get net charge of initial ion
  double Mi = compound_nucleus_.mass(); // initial mass
  double Migs = Mi - Exi_; // initial ground-state mass

  total_width_ = 0.; // total compound nucleus decay width

  marley::StructureDatabase& db = gen_.get_structure_db();

  for (const auto& f : fragments_) {

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
    double Mfgs_ion = mt.get_atomic_mass(pid_final) + (Za - qi)*me;

    // Get fragment separation energy without a redundant mass table lookup
    double Sa = Mfgs_ion + Ma - Migs;

    // Get discrete level data (if any) and models for the final nucleus
    marley::DecayScheme* ds = db.get_decay_scheme(pid_final);
    marley::OpticalModel& om = db.get_optical_model(pid_final);

    // Determine the maximum excitation energy available after fragment
    // emission in the final nucleus. This is simply the difference
    // between the initial excitation energy and the fragment separation
    // energy.
    double Exf_max = Exi_ - Sa;

    // Check if emission of this fragment is energetically allowed. If we're
    // exactly at threshold, still refuse to emit the fragment to avoid
    // numerical problems.
    if (Exf_max <= 0.) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0.;

    // If discrete level data are available for the final nuclide, get decay
    // widths for each accessible level
    if (ds) {

      // Get a vector of pointers to levels in the decay scheme. The levels are
      // sorted in order of increasing excitation energy.
      const auto& levels = ds->get_levels();

      // Use the maximum discrete level energy from the decay scheme object as
      // the lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (levels.size() > 0) E_c_min = levels.back()->energy();

      // Loop over the final discrete nuclear levels in order of increasing
      // energy until the new level energy exceeds the maximum value. For each
      // energetically allowed level, if a transition to it for a given
      // fragment orbital angular momentum l and total angular momentum j
      // conserves parity, then compute an optical model transmission
      // coefficient and add it to the total.
      for (const auto& level : levels) {
        double Exf = level->energy();
        if (Exf < Exf_max) {
          double discrete_width = 0.;
          double total_KE_CM_frame = Exf_max - Exf;
          int twoJf = level->twoJ();
          marley::Parity Pf = level->parity();
          for (int two_j = std::abs(twoJi_ - twoJf); two_j <= twoJi_ + twoJf;
            two_j += 2)
          {
            int j_plus_s = (two_j + two_s) / 2;
            // TODO: consider adding a check that l does not exceed its maximum
            // value l_max that is also used as a cutoff for transmission
            // coefficient calculations in the continuum.
            // For each new iteration, increment l and flip the overall final
            // state parity
            int l = std::abs(two_j - two_s) / 2;
            bool l_is_odd = l % 2;
            marley::Parity P_final_state = Pf * Pa * marley::Parity(!l_is_odd);
            for (; l <= j_plus_s; ++l, !P_final_state)
            {
              // Add to the total transmission coefficient if parity is
              // conserved
              if (Pi_ == P_final_state) discrete_width
                += om.transmission_coefficient(total_KE_CM_frame, fragment_pid, two_j, l,
                two_s);
            }
          }

          // Store information for this decay channel
	  exit_channels_.push_back(
            std::make_unique<marley::FragmentDiscreteExitChannel>(
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
        double total_KE_CM_frame = Exf_max - Exf;
        return fragment_continuum_partial_width(om, ldm, f, total_KE_CM_frame,
          Exf, true);
      };

      // Numerically integrate using the call wrapper, the integration bounds,
      // and the number of subintervals
      double continuum_width = marley_utils::num_integrate(cpw, E_c_min,
        Exf_max, DEFAULT_CONTINUUM_SUBINTERVALS_);

      // Create a normalized probability density function to use for sampling a
      // final excitation energy from the continuum for this exit channel.
      std::function<double(double&, double, bool)> pdf
        = [=, &f, &om, &ldm](double& total_KE_CM_frame, double Exf,
        bool use_TCs) -> double
      {
        total_KE_CM_frame = Exf_max - Exf;
        return fragment_continuum_partial_width(om, ldm, f, total_KE_CM_frame,
          Exf, use_TCs) / continuum_width;
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

  marley::GammaStrengthFunctionModel&
    gsfm = db.get_gamma_strength_function_model(Zi, Ai);

  // Let the continuum go down to 0 MeV unless there is a decay scheme object
  // available for the final nuclide (we'll check this in a second).
  double E_c_min = 0.;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if (ds) {

    // Loop over the final discrete nuclear levels in order of increasing
    // energy until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission
    // coefficient for it
    const auto& levels = ds->get_levels();

    // Use the maximum discrete level energy from the decay scheme object as
    // the lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if (levels.size() > 0) E_c_min = levels.back()->energy();

    for (const auto& level_f : levels) {
      double Exf = level_f->energy();
      if (Exf < Exi_) {
        double discrete_width = gsfm.transmission_coefficient(Exi_, twoJi_,
          Pi_, *level_f.get());

        // Store information for this decay channel if it is allowed
        if (discrete_width > 0.) {
          exit_channels_.push_back(
            std::make_unique<marley::GammaDiscreteExitChannel>(discrete_width,
            *level_f, marley::Particle(pid_initial, Migs + Exf, qi)));
          total_width_ += discrete_width;
        }
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
      [&ldm, &gsfm, this](double Exf)
      -> double { return gamma_continuum_partial_width(ldm, gsfm, Exf, true); },
      E_c_min, Exi_, DEFAULT_CONTINUUM_SUBINTERVALS_);

    // Normalized probability density used for sampling a final excitation
    // energy in the continuum
    std::function<double(double, bool)> pdf = [=, &ldm, &gsfm](double Exf,
      bool use_TCs) -> double
    { return gamma_continuum_partial_width(ldm, gsfm, Exf, use_TCs)
      / continuum_width; };

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
  size_t exit_channel_index = gen_.sample_from_distribution(exit_channel_dist);

  marley::ExitChannel* ec = exit_channels_.at(exit_channel_index).get();

  ec->do_decay(Exf, twoJf, Pf, emitted_particle, residual_nucleus, gen_);

  // TODO: consider changing this to a more realistic model
  // instead of isotropic emissions
  double cos_theta_emitted_particle = gen_.uniform_random_double(-1, 1, true);
  double phi_emitted_particle = gen_.uniform_random_double(0,
    marley_utils::two_pi, false);

  // Handle the kinematics calculations for this decay and update the
  // final-state particle objects
  marley_kinematics::two_body_decay(compound_nucleus_, emitted_particle,
    residual_nucleus, cos_theta_emitted_particle, phi_emitted_particle);

  bool discrete_level = !ec->is_continuum();
  return !discrete_level;
}

double marley::HauserFeshbachDecay::gamma_continuum_partial_width(
  marley::LevelDensityModel& ldm, marley::GammaStrengthFunctionModel& gsfm,
  double Exf, bool use_TCs)
{
  double continuum_width = 0.;
  // Approximate the gamma energy by the energy difference between the two
  // levels
  // TODO: consider adding a nuclear recoil correction here
  double e_gamma = Exi_ - Exf;
  bool initial_spin_is_zero = twoJi_ == 0;

  for (int l = 0; l <= l_max_; ++l) {
    int two_l = 2*l;
    int twoJf = twoJi_ + two_l;

    // Determine the multipolarity being considered in this trip through the
    // loop based on the value of the orbital angular momentum l.
    int mpol = ((l == 0) ? 1 : l);

    // Consider both transition types (equivalently, both final parities).
    // Check for Ji = 0 to Jf = 0 transitions, which are not allowed by angular
    // momentum conservation.
    if (!initial_spin_is_zero) {
      continuum_width += gamma_cpw_helper(ldm, gsfm, Pi_, mpol, e_gamma,
        Exf, twoJf, use_TCs);
      // Consider the other possible final spin value for this multipolarity if
      // it is allowed (Jf = Ji - l is positive)
      twoJf = twoJi_ - two_l;
      if (twoJf >= 0) continuum_width += gamma_cpw_helper(ldm, gsfm, Pi_, mpol,
        e_gamma, Exf, twoJf, use_TCs);
    }
    else if (twoJf > 0) continuum_width += gamma_cpw_helper(ldm, gsfm, Pi_,
      mpol, e_gamma, Exf, twoJf, use_TCs);
  }
  return continuum_width;
}

double marley::HauserFeshbachDecay::fragment_continuum_partial_width(
  marley::OpticalModel& om, marley::LevelDensityModel& ldm,
  const marley::Fragment& frag, double total_KE_CM_frame, double Exf,
  bool use_TCs)
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
  for (int l = 0; l <= l_max_; ++l, !Pf) {
    int two_l = 2*l;
    for (int two_j = std::abs(two_l - two_s);
      two_j <= two_l + two_s; two_j += 2)
    {
      for (int twoJf = std::abs(twoJi_ - two_j);
        twoJf <= twoJi_ + two_j; twoJf += 2)
      {
        double tc = 1.;
        if ( use_TCs ) {
          tc = om.transmission_coefficient(total_KE_CM_frame,
            frag.get_pid(), two_j, l, two_s);
        }
        continuum_width += tc * ldm.level_density(Exf, twoJf);
        // TODO: since only the spin distribution changes in the twoJf loop,
        // you can optimize this by precomputing most of the level density
        // and the multiplying here by the appropriate spin distribution
      }
    }
  }
  return continuum_width;
}

double marley::HauserFeshbachDecay::get_fragment_emission_threshold(
  const int Zi, const int Ai, const marley::Fragment& f)
{
  const marley::MassTable& mt = marley::MassTable::Instance();

  // Separation energy for the fragment
  double Sa = mt.get_fragment_separation_energy(Zi, Ai, f.get_pid());

  return Sa;
}

void marley::HauserFeshbachDecay::print(std::ostream& out) const {

  // Needed to print results in conventional units
  static constexpr double hbar = 6.58211951e-22; // MeV * s
  int pid_initial = compound_nucleus_.pdg_code();
  int Zi = marley_utils::get_particle_Z(pid_initial);
  int Ai = marley_utils::get_particle_A(pid_initial);
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
        << " emission to level at " << fdec->get_final_level().energy()
        << " MeV width = " << width << " MeV\n";
    }
    else {
      const auto* gdec
        = dynamic_cast<marley::GammaDiscreteExitChannel*>(ec.get());
      if (gdec) out << "  gamma-ray emission to level at "
        << gdec->get_final_level().energy() << " MeV width = "
        << width << " MeV\n";
    }
  }
}
