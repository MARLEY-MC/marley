#include "TMarleyGenerator.hh"
#include "TMarleyKinematics.hh"
#include "TMarleyMassTable.hh"
#include "TMarleyNuclearPhysics.hh"

using TrType = TMarleyNuclearPhysics::TransitionType;

// Table of nuclear fragments that will be considered when computing
// branching ratios for nuclear de-excitations. Spin-parity values are taken
// from nuclear ground states listed in the 10/2014 release of ENSDF.
const std::vector<TMarleyFragment> TMarleyNuclearPhysics::fragments = {
  TMarleyFragment(marley_utils::NEUTRON, 1, 1),
  TMarleyFragment(marley_utils::PROTON, 1, 1),
  TMarleyFragment(marley_utils::DEUTERON, 2, 1),
  TMarleyFragment(marley_utils::TRITON, 1, 1),
  TMarleyFragment(marley_utils::HELION, 1, 1),
  TMarleyFragment(marley_utils::ALPHA, 0, 1),
};

TrType TMarleyNuclearPhysics::determine_gamma_transition_type(int twoJi,
  TMarleyParity Pi, TMarleyLevel* level_f, int& l)
{
  int twoJf = level_f->get_two_J();
  TMarleyParity Pf = level_f->get_parity();

  return determine_gamma_transition_type(twoJi, Pi, twoJf, Pf, l);
}

TrType TMarleyNuclearPhysics::determine_gamma_transition_type(
  TMarleyLevel* level_i, TMarleyLevel* level_f, int& l)
{
  int twoJi = level_i->get_two_J();
  TMarleyParity Pi = level_i->get_parity();

  int twoJf = level_f->get_two_J();
  TMarleyParity Pf = level_f->get_parity();

  return determine_gamma_transition_type(twoJi, Pi, twoJf, Pf, l);
}

// Returns whether a gamma transition from a nuclear state with spin twoJi / 2 and
// parity Pi to a state with spin twoJf / 2 and parity Pf is electric or magnetic.
// Also loads the integer l with the multipolarity of the transition.
TrType TMarleyNuclearPhysics::determine_gamma_transition_type(int twoJi,
  TMarleyParity Pi, int twoJf, TMarleyParity Pf, int& l)
{
  // TODO: reconsider how to handle this situation
  if (twoJi == 0 && twoJf == 0) throw std::runtime_error(
    std::string("0 -> 0 EM transitions are not allowed."));

  int two_delta_J = std::abs(twoJf - twoJi);
  // Odd values of two_delta_J are unphysical because they represent EM
  // transitions where the total angular momentum changes by half (photons are
  // spin 1)
  if (two_delta_J % 2) throw std::runtime_error(std::string("Unphysical ")
    + "EM transition encountered between nuclear levels with spins 2*Ji = "
    + std::to_string(twoJi) + " and 2*Jf = " + std::to_string(twoJf));

  TMarleyParity Pi_times_Pf = Pi * Pf;

  // Determine the multipolarity of this transition. Load l with the result.
  if (two_delta_J == 0) l = 1;
  // We already checked for unphysical odd values of two_delta_J above, so
  // using integer division here will always give the expected result.
  else l = two_delta_J / 2;

  // Determine whether this transition is electric or magnetic based on its
  // multipolarity and whether or not there is a change of parity

  // Pi * Pf = -1 gives electric transitions for odd l
  TMarleyParity electric_parity;
  if (l % 2) electric_parity = -1; // l is odd
  // Pi * Pf = +1 gives electric transitions for even l
  else electric_parity = 1; // l is even

  if (Pi_times_Pf == electric_parity) return TransitionType::electric;
  else return TransitionType::magnetic;
}

double TMarleyNuclearPhysics::gamma_strength_function_coefficient(int Z, int A,
  TransitionType type, int l, double e_gamma)
{
  // TODO: improve this error message
  if (l < 1) throw std::runtime_error(std::string("Invalid multipolarity")
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
  else throw std::runtime_error(std::string("Invalid transition type")
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
double TMarleyNuclearPhysics::weisskopf_partial_decay_width(int A,
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
    double mp = TMarleyMassTable::get_particle_mass(marley_utils::PROTON);
    return 10 * el_width * std::pow(marley_utils::hbar_c / (mp * R), 2);
  }

  // TODO: improve this error message
  else throw std::runtime_error(std::string("Invalid transition type")
    + " given for Weisskopf gamma-ray partial width calculation");
}

double TMarleyNuclearPhysics::gamma_continuum_partial_width(int Z, int A,
  int twoJi, double Exi, double Exf)
{
  double continuum_width = 0;
  // Approximate the gamma energy by the energy difference between the two levels
  // TODO: consider adding a nuclear recoil correction here
  double e_gamma = Exi - Exf;
  bool initial_spin_is_zero = twoJi == 0;
  for (int l = 0; l <= l_max; ++l) {
    int two_l = 2*l;
    int twoJf = twoJi + two_l;
    // Determine the multipolarity being considered in this trip through the
    // loop based on the value of the orbital angular momentum l.
    int mpol = ((l == 0) ? 1 : l);
    // Consider both transition types (equivalently, both final parities). Check
    // for Ji = 0 to Jf = 0 transitions, which are not allowed by angular
    // momentum conservation.
    if (!initial_spin_is_zero) {
      continuum_width += gamma_cpw(Z, A, mpol, twoJf, e_gamma, Exf);
      // Consider the other possible final spin value for this multipolarity if
      // it is allowed (Jf = Ji - l is positive)
      twoJf = twoJi - two_l;
      if (twoJf >= 0) continuum_width += gamma_cpw(Z, A, mpol, twoJf, e_gamma,
        Exf);
    }
    else if (twoJf > 0) continuum_width += gamma_cpw(Z, A, mpol, twoJf, e_gamma,
      Exf);
  }
  return continuum_width;
}

// Computes the partial width (multiplied by 2*pi and the initial nuclear level
// density) predicted by the Hauser-Feshbach statistical model for decay of a
// given nuclear level via gamma-ray emission.
double TMarleyNuclearPhysics::hf_gamma_partial_width(double Ex, int twoJi,
  TMarleyParity Pi, const TMarleyDecayScheme& ds)
{
  int Z = ds.get_Z();
  int A = ds.get_A();
  double discrete_width = 0;
  bool initial_spin_is_zero = twoJi == 0;
  // Loop over the final discrete nuclear levels in order of increasing energy
  // until the new level energy exceeds the maximum value. For each
  // energetically allowed level, compute a gamma ray transmission coefficient
  // for it and add it to the total.
  const std::vector<TMarleyLevel*>* sorted_lps = ds.get_sorted_level_pointers();
  for (const auto& level_f : *sorted_lps) {
    int twoJf = level_f->get_two_J();
    // 0->0 EM transitions aren't allowed due to angular momentum conservation
    // (photons are spin 1), so if the initial and final spins are zero, skip
    // ahead to the next final level.
    if (initial_spin_is_zero && twoJf == 0) continue;
    double Exf = level_f->get_energy();
    if (Exf < Ex) {
      // Approximate the gamma energy by the energy difference between the two levels
      // TODO: consider adding a nuclear recoil correction here
      double e_gamma = Ex - Exf;
      // Determine the type (electric or magnetic) and multipolarity of the gamma
      // transition between these two levels
      int mpol; // Multipolarity
      auto type = determine_gamma_transition_type(twoJi, Pi, level_f, mpol);
      // TODO: allow the user to choose which gamma-ray transition model to use 
      // (Weisskopf single-particle estimates, Brink-Axel strength functions, etc.)
      discrete_width += gamma_transmission_coefficient(Z, A, type, mpol, e_gamma);
    }
    else break;
  }

  // Use the maximum discrete level energy from the decay scheme object as the
  // lower bound for the continuum. If the decay scheme object doesn't have
  // any levels, use the continuum all the way down to zero excitation energy.
  // TODO: consider whether this is the best approach
  double E_c_min;
  if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
  else E_c_min = 0;

  double continuum_width = 0.;

  // If transitions to the energy continuum are possible, add the total
  // continuum contribution to the decay width
  if (Ex > E_c_min) {

    // Create a forwarding call wrapper for the continuum partial width member
    // function that takes a single argument.
    std::function<double(double)> gpw = std::bind(&gamma_continuum_partial_width,
      Z, A, twoJi, Ex, std::placeholders::_1 /*Exf*/);

    // Numerically integrate using the call wrapper, the integration bounds, and
    // the number of subintervals
    continuum_width += marley_utils::num_integrate(gpw, E_c_min, Ex,
      DEFAULT_CONTINUUM_SUBINTERVALS);
  }

  return discrete_width + continuum_width;
}

// Computes the partial width (multiplied by 2*pi and the initial nuclear level
// density) predicted by the Hauser-Feshbach statistical model for decay of a
// given nuclear level via emission of the fragment f.
double TMarleyNuclearPhysics::hf_fragment_partial_width(int Zi, int Ai,
  double Ex, int twoJi, TMarleyParity Pi, const TMarleyFragment& f,
  const TMarleySphericalOpticalModel& om, const TMarleyDecayScheme& ds)
{
  // Get information about the current fragment
  int two_s = f.get_two_s();
  TMarleyParity Pa = f.get_parity();
  int fragment_pid = f.get_pid();
  int Zf = om.get_Z();
  int Af = om.get_A();

  // Get various masses that will be needed for fragment emission kinematics
  // (including nuclear recoil)
  double Migs = TMarleyMassTable::get_atomic_mass(Zi, Ai);
  double Mi = Migs + Ex;
  double Ma = f.get_mass();
  double Mfgs = TMarleyMassTable::get_atomic_mass(Zf, Af);
  double me = TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);
  // Fragment atomic number
  int Za = f.get_Z();
  // Approximate the ground state rest energy of the negative ion (with
  // charge Za-) formed when the bare fragment f is emitted by adding Za
  // electron masses to the atomic mass for the final nucleus.
  double Mfgs_ion = Mfgs + Za*me;
  // Get fragment separation energy without a redundant mass table lookup
  double Sa = Mfgs_ion + Ma - Migs;

  // Determine the maximum excitation energy available after fragment emission
  // in the final nucleus. Use a simpler functional form as a shortcut if
  // the fragment is neutral [and can therefore be emitted with arbitrarily
  // low kinetic energy due to the absence of a Coulomb barrier])
  double Exf_max;
  double Vc; // Coulomb barrier
  if (f.get_Z() == 0) {
    Vc = 0;
    Exf_max = Ex - Sa;
  }
  else {
    int Zf = om.get_Z();
    int Af = om.get_A();
    Vc = coulomb_barrier(Zf, Af, Za, f.get_A()); // Ea_min = Vc
    Exf_max = std::sqrt(std::pow(Mi - Ma, 2) - 2*Vc*Mi) - Mfgs_ion;
  }

  // The fragment kinetic energy depends on the final level energy. For speed,
  // we precompute here a portion of the fragment energy that doesn't depend on
  // the final level energy.
  double Mconst = (Ex - Sa) * (Mi + Mfgs_ion - Ma);

  // Check if emission of this fragment is energetically allowed. If we're below
  // the fragment emission threshold, the partial decay width is zero, so we can
  // return that value immediately rather than slogging through the calculation.
  // We'll also return zero if we're exactly at threshold to avoid numerical problems.
  if (Mconst / (2 * Mi) - Vc <= 0) return 0.;

  // Get a pointer to the vector of sorted pointers to levels in the decay
  // scheme
  const std::vector<TMarleyLevel*>* sorted_lps = ds.get_sorted_level_pointers();

  double discrete_width = fragment_discrete_partial_width(Exf_max, Mconst,
    Mfgs_ion, Mi, twoJi, Pi, fragment_pid, two_s, Pa, om, sorted_lps);

  // Use the maximum discrete level energy from the decay scheme object as the
  // lower bound for the continuum
  // TODO: consider whether this is the best approach
  double E_c_min;
  if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
  else E_c_min = 0;

  double continuum_width = 0.;

  // If transitions to the energy continuum are possible, add the total
  // continuum contribution to the decay width
  if (Exf_max > E_c_min) {

    // Create a forwarding call wrapper for the continuum partial width member
    // function that takes a single argument.
    std::function<double(double)> cpw = std::bind(
      &fragment_continuum_partial_width, Mconst, Mfgs_ion, Mi, twoJi, Pi,
      fragment_pid, two_s, Pa, om, std::placeholders::_1 /*Exf*/);

    // Numerically integrate using the call wrapper, the integration bounds, and
    // the number of subintervals
    continuum_width = marley_utils::num_integrate(cpw, E_c_min, Exf_max,
      DEFAULT_CONTINUUM_SUBINTERVALS);
  }

  return discrete_width + continuum_width;
}

double TMarleyNuclearPhysics::fragment_continuum_partial_width(double Mconst,
  double Mfgs_ion, double Mi, int twoJi, TMarleyParity Pi, int fragment_pid,
  int two_s, TMarleyParity Pa, const TMarleySphericalOpticalModel& om,
  double Exf)
{
  double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
  double continuum_width = 0;
  int Zf = om.get_Z();
  int Af = om.get_A();
  // Final nuclear state parity
  TMarleyParity Pf;
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
      for (int twoJf = std::abs(twoJi - two_j);
        twoJf <= twoJi + two_j; twoJf += 2)
      {
        continuum_width += om.transmission_coefficient(Ea, fragment_pid, two_j, l,
          two_s, DEFAULT_NUMEROV_STEP_SIZE)
          * TMarleyBackshiftedFermiGasModel::level_density(Zf, Af, Exf, twoJf);
          // TODO: since only the spin distribution changes in the twoJf loop,
          // you can optimize this by precomputing most of the level density
          // and the multiplying here by the appropriate spin distribution
      }
    }
  }
  return continuum_width;
}

double TMarleyNuclearPhysics::fragment_discrete_partial_width(double Exf_max,
  double Mconst, double Mfgs_ion, double Mi, int twoJi, TMarleyParity Pi,
  int fragment_pid, int two_s, TMarleyParity Pa,
  const TMarleySphericalOpticalModel& om,
  const std::vector<TMarleyLevel*>* sorted_lps)
{
  double discrete_width = 0;
  // Loop over the final discrete nuclear levels in order of increasing energy
  // until the new level energy exceeds the maximum value. For each
  // energetically allowed level, if a transition to it for a given fragment
  // orbital angular momentum l and total angular momentum j conserves parity,
  // then compute an optical model transmission coefficient and add it to the total.
  for (const auto& level : *sorted_lps) {
    double Exf = level->get_energy();
    if (Exf < Exf_max) {
      double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
      int twoJf = level->get_two_J();
      TMarleyParity Pf = level->get_parity();
      for (int two_j = std::abs(twoJi - twoJf); two_j <= twoJi + twoJf;
        two_j += 2)
      {
        int j_plus_s = (two_j + two_s) / 2;
        // TODO: consider adding a check that l does not exceed its maximum
        // value l_max that is also used as a cutoff for transmission
        // coefficient calculations in the continuum.
        // For each new iteration, increment l and flip the overall final state parity
        int l = std::abs(two_j - two_s) / 2;
        bool l_is_odd = l % 2;
        TMarleyParity P_final_state = Pf * Pa * TMarleyParity(!l_is_odd);
        for (; l <= j_plus_s; ++l, !P_final_state)
        {
          // Add to the total transmission coefficient if parity is conserved
          if (Pi == P_final_state) discrete_width
            += om.transmission_coefficient(Ea, fragment_pid, two_j, l, two_s,
            DEFAULT_NUMEROV_STEP_SIZE);
        }
      }
    }
    else break;
  }
  return discrete_width;
}

// Loads the two product particle objects with the final-state particles from a
// two-body Hauser-Feshbach decay of the initial particle. Also loads Ex, twoJ,
// and Pi (which are originally assumed to contain the initial nuclear
// excitation energy, 2J, and parity) with the final nuclear excitation energy,
// 2J, and parity.  Returns true if the final-state nucleus is in the unbound
// continuum, or false if it is in a bound state level included in a decay
// scheme object.
bool TMarleyNuclearPhysics::hauser_feshbach_decay(int Zi, int Ai,
  const TMarleyParticle& initial_particle, TMarleyParticle& first_product,
  TMarleyParticle& second_product, double& Ex, int& twoJ, TMarleyParity& Pi,
  TMarleyStructureDatabase& db, TMarleyGenerator& gen)
{
  double Mi = initial_particle.get_mass();
  double Migs = Mi - Ex;
  double me = TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);

  double total_width = 0.;

  std::vector<double> widths;
  std::vector<const TMarleyFragment*> frags;
  std::vector<double> final_ion_gs_masses;
  std::vector<TMarleyLevel*> levels;

  std::unordered_map<const TMarleyFragment*, std::function<double(double)> >
    cpws;
  std::unordered_map<const TMarleyFragment*, double> total_cws;
  std::unordered_map<const TMarleyFragment*, double> E_c_mins;
  std::unordered_map<const TMarleyFragment*, double> Exf_maxes;

  for (const auto& f : TMarleyNuclearPhysics::get_fragments()) {
    int Zf = Zi - f.get_Z();
    int Af = Ai - f.get_A();

    TMarleySphericalOpticalModel om(Zf, Af);
    TMarleyDecayScheme* ds = db.get_decay_scheme(Zf, Af);

    // Get information about the current fragment
    int two_s = f.get_two_s();
    TMarleyParity Pa = f.get_parity();
    int fragment_pid = f.get_pid();
  
    // Get various masses that will be needed for fragment emission kinematics
    // (including nuclear recoil)
    double Ma = f.get_mass();
    double Mfgs = TMarleyMassTable::get_atomic_mass(Zf, Af);
    // Fragment atomic number
    int Za = f.get_Z();
    // Approximate the ground state rest energy of the negative ion (with
    // charge Za-) formed when the bare fragment f is emitted by adding Za
    // electron masses to the atomic mass for the final nucleus.
    double Mfgs_ion = Mfgs + Za*me;
    // Get fragment separation energy without a redundant mass table lookup
    double Sa = Mfgs_ion + Ma - Migs;
  
    // Determine the maximum excitation energy available after fragment emission
    // in the final nucleus. Use a simpler functional form as a shortcut if
    // the fragment is neutral [and can therefore be emitted with arbitrarily
    // low kinetic energy due to the absence of a Coulomb barrier])
    double Exf_max;
    double Vc; // Coulomb barrier
    if (f.get_Z() == 0) {
      Vc = 0;
      Exf_max = Ex - Sa;
    }
    else {
      Vc = coulomb_barrier(Zf, Af, Za, f.get_A()); // Ea_min = Vc
      Exf_max = std::sqrt(std::pow(Mi - Ma, 2) - 2*Vc*Mi) - Mfgs_ion;
    }
  
    // The fragment kinetic energy depends on the final level energy. For speed,
    // we precompute here a portion of the fragment energy that doesn't depend on
    // the final level energy.
    double Mconst = (Ex - Sa) * (Mi + Mfgs_ion - Ma);
  
    // Check if emission of this fragment is energetically allowed. If we're below
    // the fragment emission threshold, the partial decay width is zero, so we can
    // skip this fragment rather than slogging through the calculation.
    // We'll also skip it if we're exactly at threshold to avoid numerical problems.
    if (Mconst / (2 * Mi) - Vc <= 0) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0;
  
    // If discrete level data is available for the final nuclide, get decay
    // widths for each accessible level
    if (ds != nullptr) {

      // Get a pointer to the vector of sorted pointers to levels in the decay
      // scheme
      const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

      // Use the maximum discrete level energy from the decay scheme object as the
      // lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
  
      // Loop over the final discrete nuclear levels in order of increasing energy
      // until the new level energy exceeds the maximum value. For each
      // energetically allowed level, if a transition to it for a given fragment
      // orbital angular momentum l and total angular momentum j conserves parity,
      // then compute an optical model transmission coefficient and add it to the total.
      for (const auto& level : *sorted_lps) {
        double Exf = level->get_energy();
        if (Exf < Exf_max) {
          double discrete_width = 0;
          double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
          int twoJf = level->get_two_J();
          TMarleyParity Pf = level->get_parity();
          for (int two_j = std::abs(twoJ - twoJf); two_j <= twoJ + twoJf;
            two_j += 2)
          {
            int j_plus_s = (two_j + two_s) / 2;
            // TODO: consider adding a check that l does not exceed its maximum
            // value l_max that is also used as a cutoff for transmission
            // coefficient calculations in the continuum.
            // For each new iteration, increment l and flip the overall final state parity
            int l = std::abs(two_j - two_s) / 2;
            bool l_is_odd = l % 2;
            TMarleyParity P_final_state = Pf * Pa * TMarleyParity(!l_is_odd);
            for (; l <= j_plus_s; ++l, !P_final_state)
            {
              // Add to the total transmission coefficient if parity is conserved
              if (Pi == P_final_state) discrete_width
                += om.transmission_coefficient(Ea, fragment_pid, two_j, l, two_s,
                DEFAULT_NUMEROV_STEP_SIZE);
            }
          }

          // Store information for this decay channel
          widths.push_back(discrete_width);
          total_width += discrete_width;
          levels.push_back(level);
          frags.push_back(&f);
          final_ion_gs_masses.push_back(Mfgs_ion);
        }
        else break;
      }
    } 

    // If transitions to the energy continuum are possible, include the
    // continuum in the decay channels
    if (Exf_max > E_c_min) {
  
      // Create a forwarding call wrapper for the continuum partial width member
      // function that takes a single argument.
      std::function<double(double)> cpw = std::bind(
        &fragment_continuum_partial_width, Mconst, Mfgs_ion, Mi, twoJ, Pi,
        fragment_pid, two_s, Pa, om, std::placeholders::_1 /*Exf*/);
  
      // Numerically integrate using the call wrapper, the integration bounds, and
      // the number of subintervals
      double continuum_width = marley_utils::num_integrate(cpw, E_c_min,
        Exf_max, DEFAULT_CONTINUUM_SUBINTERVALS);

      // Store information for this decay channel
      widths.push_back(continuum_width);
      total_width += continuum_width;
      levels.push_back(nullptr);
      frags.push_back(&f);
      final_ion_gs_masses.push_back(Mfgs_ion);
      cpws[&f] = cpw;
      total_cws[&f] = continuum_width;
      E_c_mins[&f] = E_c_min;
      Exf_maxes[&f] = Exf_max;
    }
  }

  TMarleyDecayScheme* ds = db.get_decay_scheme(Zi, Ai);

  // Let the continuum go down to 0 MeV unless there is a decay scheme object
  // available for the final nuclide (we'll check this in a second).
  double E_c_min = 0;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if (ds != nullptr) {

    bool initial_spin_is_zero = twoJ == 0;
    // Loop over the final discrete nuclear levels in order of increasing energy
    // until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission coefficient
    // for it
    const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

    // Use the maximum discrete level energy from the decay scheme object as the
    // lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();

    for (const auto& level_f : *sorted_lps) {
      int twoJf = level_f->get_two_J();
      // 0->0 EM transitions aren't allowed due to angular momentum conservation
      // (photons are spin 1), so if the initial and final spins are zero, skip
      // ahead to the next final level.
      if (initial_spin_is_zero && twoJf == 0) continue;
      double Exf = level_f->get_energy();
      if (Exf < Ex) {
        // Approximate the gamma energy by the energy difference between the two levels
        // TODO: consider adding a nuclear recoil correction here
        double e_gamma = Ex - Exf;
        // Determine the type (electric or magnetic) and multipolarity of the gamma
        // transition between these two levels
        int mpol; // Multipolarity
        auto type = determine_gamma_transition_type(twoJ, Pi, level_f, mpol);
        // TODO: allow the user to choose which gamma-ray transition model to use 
        // (Weisskopf single-particle estimates, Brink-Axel strength functions, etc.)
        double discrete_width = gamma_transmission_coefficient(Zi, Ai, type,
          mpol, e_gamma);

        // Store information for this decay channel
        widths.push_back(discrete_width);
        total_width += discrete_width;
        levels.push_back(level_f);
        frags.push_back(nullptr); // A nullptr entry here means a gamma-ray transition
        final_ion_gs_masses.push_back(Migs);
      }
      else break;
    }
  }

  // If gamma transitions to the energy continuum are possible, include them
  // in the possible decay channels
  if (Ex > E_c_min) {

    // Create a forwarding call wrapper for the continuum partial width member
    // function that takes a single argument.
    std::function<double(double)> gpw = std::bind(&gamma_continuum_partial_width,
      Zi, Ai, twoJ, Ex, std::placeholders::_1 /*Exf*/);

    // Numerically integrate using the call wrapper, the integration bounds, and
    // the number of subintervals
    double continuum_width = marley_utils::num_integrate(gpw, E_c_min, Ex,
      DEFAULT_CONTINUUM_SUBINTERVALS);

    // Store information for this decay channel
    widths.push_back(continuum_width);
    total_width += continuum_width;
    levels.push_back(nullptr);
    frags.push_back(nullptr); // A nullptr entry here means a gamma-ray transition
    final_ion_gs_masses.push_back(Migs);
    cpws[nullptr] = gpw;
    total_cws[nullptr] = continuum_width;
    E_c_mins[nullptr] = E_c_min;
    Exf_maxes[nullptr] = Ex;
  }

  // Throw an error if all decays are impossible
  if (total_width <= 0.) throw std::runtime_error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial decay widths are zero.");

  // Sample an exit channel
  std::discrete_distribution<size_t> exit_channel_dist(widths.begin(),
    widths.end());
  size_t exit_channel_index = gen.discrete_sample(exit_channel_dist);

  double Exf, m_frag;
  int pid_frag, pid_final;

  double mfgs = final_ion_gs_masses.at(exit_channel_index);

  const TMarleyLevel* lev = levels.at(exit_channel_index);
  const TMarleyFragment* fr = frags.at(exit_channel_index);
  bool discrete_level = lev != nullptr;
  bool evaporation = fr != nullptr;

  if (discrete_level) {
    Exf = lev->get_energy();
    twoJ = lev->get_two_J();
    Pi = lev->get_parity();
  }
  else {

    double cw = total_cws.at(fr);
    std::function<double(double)> cpw = cpws.at(fr);

    std::function<double(double)> pw = [&cw, &cpw] (double energy)
      -> double { return cpw(energy) / cw; };

    Exf = gen.rejection_sample(pw, E_c_mins.at(fr), Exf_maxes.at(fr));

    // Determine final spin-parity here
    if (evaporation) {
      double Ma = fr->get_mass();
      double Sa = mfgs + Ma - Migs;
      double Mconst = (Ex - Sa) * (Mi + mfgs - Ma);
      double Ea = (Mconst - Exf*(2*mfgs + Exf)) / (2 * Mi);
      int Zf = Zi - fr->get_Z();
      int Af = Ai - fr->get_A();
      TMarleySphericalOpticalModel om(Zf, Af);
      sample_fragment_spin_parity(twoJ, Pi, *fr, om, gen, Exf, Ea);
    }
    else sample_gamma_spin_parity(Zi, Ai, twoJ, Pi, Ex, Exf, gen);
  }

  if (evaporation) {

    pid_frag = fr->get_pid();
    m_frag = fr->get_mass();

    int Zf = Zi - fr->get_Z();
    int Af = Ai - fr->get_A();
    pid_final = marley_utils::get_nucleus_pid(Zf, Af);
  }

  else {
    // Gamma emission
    pid_frag = marley_utils::PHOTON;
    m_frag = TMarleyMassTable::get_particle_mass(pid_frag);

    pid_final = marley_utils::get_nucleus_pid(Zi, Ai);
  }

  first_product = TMarleyParticle(pid_frag, m_frag);
  second_product = TMarleyParticle(pid_final, mfgs + Exf);

  // TODO: consider changing this to a more realistic model
  // instead of isotropic emissions
  double cos_theta_first = gen.uniform_random_double(-1, 1, true);
  double phi_first = gen.uniform_random_double(0, marley_utils::two_pi, false);

  // Handle the kinematics calculations for this decay
  TMarleyKinematics::two_body_decay(initial_particle, first_product,
    second_product, cos_theta_first, phi_first);

  // Update the excitation energy to its final value
  Ex = Exf;

  return !discrete_level;
}

// Based on some initial parameters (including the nuclear 2J and parity)
// determine the final spin and parity and store it in twoJ and Pi
void TMarleyNuclearPhysics::sample_gamma_spin_parity(int Z, int A,
  int& twoJ, TMarleyParity& Pi, double Exi, double Exf, TMarleyGenerator& gen)
{
  std::vector<int> twoJfs;
  std::vector<TMarleyParity> Pfs;
  std::vector<double> widths;

  double continuum_width = 0;
  TMarleyParity Pf;
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

      double tcE = gamma_transmission_coefficient(Z, A,
        TransitionType::electric, mpol, e_gamma);
      double tcM = gamma_transmission_coefficient(Z, A,
        TransitionType::magnetic, mpol, e_gamma);

      if (!initial_spin_is_zero) {

        continuum_width += store_gamma_pws(Z, A, Exf, twoJf, Pi, widths, twoJfs,
          Pfs, tcE, tcM, mpol);

	// Consider the other possible final spin value for this multipolarity
	// if it is allowed (Jf = Ji - l is positive)
        twoJf = twoJ - two_l;
        if (twoJf >= 0) continuum_width += store_gamma_pws(Z, A, Exf, twoJf,
          Pi, widths, twoJfs, Pfs, tcE, tcM, mpol);
      }
      // twoJf > 0 but twoJi == 0
      else continuum_width += store_gamma_pws(Z, A, Exf, twoJf, Pi, widths,
        twoJfs, Pfs, tcE, tcM, mpol);
    }
  }
  // Throw an error if all decays are impossible
  if (continuum_width <= 0.) throw std::runtime_error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial gamma decay widths are "
    + "zero.");

  // Sample a final spin and parity
  std::discrete_distribution<size_t> jpi_dist(widths.begin(), widths.end());
  size_t jpi_index = gen.discrete_sample(jpi_dist);

  // Store the results
  twoJ = twoJfs.at(jpi_index);
  Pi = Pfs.at(jpi_index);
}

void TMarleyNuclearPhysics::sample_fragment_spin_parity(int& twoJ,
  TMarleyParity& Pi, const TMarleyFragment& f,
  const TMarleySphericalOpticalModel& om, TMarleyGenerator& gen,
  double Exf, double Ea)
{
  int fragment_pid = f.get_pid();
  int two_s = f.get_two_s();
  TMarleyParity Pa = f.get_parity();

  std::vector<double> widths;
  std::vector<int> twoJfs;
  std::vector<TMarleyParity> Pfs;

  double continuum_width = 0;
  int Zf = om.get_Z();
  int Af = om.get_A();
  // Final nuclear state parity
  TMarleyParity Pf;
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
        // and the multiplying here by the appropriate spin distribution
        double width = om.transmission_coefficient(Ea, fragment_pid, two_j, l,
          two_s, DEFAULT_NUMEROV_STEP_SIZE)
          * TMarleyBackshiftedFermiGasModel::level_density(Zf, Af, Exf, twoJf);

        // Store the computed decay width for sampling
        continuum_width += width;
        widths.push_back(width);
        twoJfs.push_back(twoJf);
        Pfs.push_back(Pf);
      }
    }
  }
  // Throw an error if all decays are impossible
  if (continuum_width <= 0.) throw std::runtime_error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial fragment decay widths are "
    + "zero. cw = " + std::to_string(continuum_width));

  // Sample a final spin and parity
  std::discrete_distribution<size_t> jpi_dist(widths.begin(), widths.end());
  size_t jpi_index = gen.discrete_sample(jpi_dist);

  // Store the results
  twoJ = twoJfs.at(jpi_index);
  Pi = Pfs.at(jpi_index);

}

void TMarleyNuclearPhysics::hf_test(int Zi, int Ai,
  const TMarleyParticle& initial_particle, double Ex,
  int twoJ, TMarleyParity Pi,
  TMarleyStructureDatabase& db, TMarleyGenerator& gen)
{

  double Mi = initial_particle.get_mass();
  double Migs = Mi - Ex;
  double me = TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);

  double total_width = 0.;

  std::vector<double> widths;
  std::vector<const TMarleyFragment*> frags;
  std::vector<double> final_ion_gs_masses;
  std::vector<TMarleyLevel*> levels;

  std::unordered_map<const TMarleyFragment*, std::function<double(double)> >
    cpws;
  std::unordered_map<const TMarleyFragment*, double> total_cws;
  std::unordered_map<const TMarleyFragment*, double> E_c_mins;
  std::unordered_map<const TMarleyFragment*, double> Exf_maxes;

  for (const auto& f : TMarleyNuclearPhysics::get_fragments()) {
    int Zf = Zi - f.get_Z();
    int Af = Ai - f.get_A();

    TMarleySphericalOpticalModel om(Zf, Af);
    TMarleyDecayScheme* ds = db.get_decay_scheme(Zf, Af);

    // Get information about the current fragment
    int two_s = f.get_two_s();
    TMarleyParity Pa = f.get_parity();
    int fragment_pid = f.get_pid();
  
    // Get various masses that will be needed for fragment emission kinematics
    // (including nuclear recoil)
    double Ma = f.get_mass();
    double Mfgs = TMarleyMassTable::get_atomic_mass(Zf, Af);
    // Fragment atomic number
    int Za = f.get_Z();
    // Approximate the ground state rest energy of the negative ion (with
    // charge Za-) formed when the bare fragment f is emitted by adding Za
    // electron masses to the atomic mass for the final nucleus.
    double Mfgs_ion = Mfgs + Za*me;
    // Get fragment separation energy without a redundant mass table lookup
    double Sa = Mfgs_ion + Ma - Migs;
  
    // Determine the maximum excitation energy available after fragment emission
    // in the final nucleus. Use a simpler functional form as a shortcut if
    // the fragment is neutral [and can therefore be emitted with arbitrarily
    // low kinetic energy due to the absence of a Coulomb barrier])
    double Exf_max;
    double Vc; // Coulomb barrier
    if (f.get_Z() == 0) {
      Vc = 0;
      Exf_max = Ex - Sa;
    }
    else {
      Vc = coulomb_barrier(Zf, Af, Za, f.get_A()); // Ea_min = Vc
      Exf_max = std::sqrt(std::pow(Mi - Ma, 2) - 2*Vc*Mi) - Mfgs_ion;
    }
  
    // The fragment kinetic energy depends on the final level energy. For speed,
    // we precompute here a portion of the fragment energy that doesn't depend on
    // the final level energy.
    double Mconst = (Ex - Sa) * (Mi + Mfgs_ion - Ma);
  
    // Check if emission of this fragment is energetically allowed. If we're below
    // the fragment emission threshold, the partial decay width is zero, so we can
    // skip this fragment rather than slogging through the calculation.
    // We'll also skip it if we're exactly at threshold to avoid numerical problems.
    if (Mconst / (2 * Mi) - Vc <= 0) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0;
  
    // If discrete level data is available for the final nuclide, get decay
    // widths for each accessible level
    if (ds != nullptr) {

      // Get a pointer to the vector of sorted pointers to levels in the decay
      // scheme
      const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

      // Use the maximum discrete level energy from the decay scheme object as the
      // lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
  
      // Loop over the final discrete nuclear levels in order of increasing energy
      // until the new level energy exceeds the maximum value. For each
      // energetically allowed level, if a transition to it for a given fragment
      // orbital angular momentum l and total angular momentum j conserves parity,
      // then compute an optical model transmission coefficient and add it to the total.
      for (const auto& level : *sorted_lps) {
        double Exf = level->get_energy();
        if (Exf < Exf_max) {
          double discrete_width = 0;
          double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
          int twoJf = level->get_two_J();
          TMarleyParity Pf = level->get_parity();
          for (int two_j = std::abs(twoJ - twoJf); two_j <= twoJ + twoJf;
            two_j += 2)
          {
            int j_plus_s = (two_j + two_s) / 2;
            // TODO: consider adding a check that l does not exceed its maximum
            // value l_max that is also used as a cutoff for transmission
            // coefficient calculations in the continuum.
            // For each new iteration, increment l and flip the overall final state parity
            int l = std::abs(two_j - two_s) / 2;
            bool l_is_odd = l % 2;
            TMarleyParity P_final_state = Pf * Pa * TMarleyParity(!l_is_odd);
            for (; l <= j_plus_s; ++l, !P_final_state)
            {
              // Add to the total transmission coefficient if parity is conserved
              if (Pi == P_final_state) discrete_width
                += om.transmission_coefficient(Ea, fragment_pid, two_j, l, two_s,
                DEFAULT_NUMEROV_STEP_SIZE);
            }
          }

          // Store information for this decay channel
          widths.push_back(discrete_width);
          total_width += discrete_width;
          levels.push_back(level);
          frags.push_back(&f);
          final_ion_gs_masses.push_back(Mfgs_ion);
        }
        else break;
      }
    } 

    // If transitions to the energy continuum are possible, include the
    // continuum in the decay channels
    if (Exf_max > E_c_min) {
  
      // Create a forwarding call wrapper for the continuum partial width member
      // function that takes a single argument.
      std::function<double(double)> cpw = std::bind(
        &fragment_continuum_partial_width, Mconst, Mfgs_ion, Mi, twoJ, Pi,
        fragment_pid, two_s, Pa, om, std::placeholders::_1 /*Exf*/);
  
      // Numerically integrate using the call wrapper, the integration bounds, and
      // the number of subintervals
      double continuum_width = marley_utils::num_integrate(cpw, E_c_min,
        Exf_max, DEFAULT_CONTINUUM_SUBINTERVALS);

      // Store information for this decay channel
      widths.push_back(continuum_width);
      total_width += continuum_width;
      levels.push_back(nullptr);
      frags.push_back(&f);
      final_ion_gs_masses.push_back(Mfgs_ion);
      cpws[&f] = cpw;
      total_cws[&f] = continuum_width;
      E_c_mins[&f] = E_c_min;
      Exf_maxes[&f] = Exf_max;
    }
  }

  TMarleyDecayScheme* ds = db.get_decay_scheme(Zi, Ai);

  // Let the continuum go down to 0 MeV unless there is a decay scheme object
  // available for the final nuclide (we'll check this in a second).
  double E_c_min = 0;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if (ds != nullptr) {

    bool initial_spin_is_zero = twoJ == 0;
    // Loop over the final discrete nuclear levels in order of increasing energy
    // until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission coefficient
    // for it
    const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

    // Use the maximum discrete level energy from the decay scheme object as the
    // lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();

    for (const auto& level_f : *sorted_lps) {
      int twoJf = level_f->get_two_J();
      // 0->0 EM transitions aren't allowed due to angular momentum conservation
      // (photons are spin 1), so if the initial and final spins are zero, skip
      // ahead to the next final level.
      if (initial_spin_is_zero && twoJf == 0) continue;
      double Exf = level_f->get_energy();
      if (Exf < Ex) {
        // Approximate the gamma energy by the energy difference between the two levels
        // TODO: consider adding a nuclear recoil correction here
        double e_gamma = Ex - Exf;
        // Determine the type (electric or magnetic) and multipolarity of the gamma
        // transition between these two levels
        int mpol; // Multipolarity
        auto type = determine_gamma_transition_type(twoJ, Pi, level_f, mpol);
        // TODO: allow the user to choose which gamma-ray transition model to use 
        // (Weisskopf single-particle estimates, Brink-Axel strength functions, etc.)
        double discrete_width = gamma_transmission_coefficient(Zi, Ai, type,
          mpol, e_gamma);

        // Store information for this decay channel
        widths.push_back(discrete_width);
        total_width += discrete_width;
        levels.push_back(level_f);
        frags.push_back(nullptr); // A nullptr entry here means a gamma-ray transition
        final_ion_gs_masses.push_back(Migs);
      }
      else break;
    }
  }

  // If gamma transitions to the energy continuum are possible, include them
  // in the possible decay channels
  if (Ex > E_c_min) {

    // Create a forwarding call wrapper for the continuum partial width member
    // function that takes a single argument.
    std::function<double(double)> gpw = std::bind(&gamma_continuum_partial_width,
      Zi, Ai, twoJ, Ex, std::placeholders::_1 /*Exf*/);

    // Numerically integrate using the call wrapper, the integration bounds, and
    // the number of subintervals
    double continuum_width = marley_utils::num_integrate(gpw, E_c_min, Ex,
      DEFAULT_CONTINUUM_SUBINTERVALS);

    // Store information for this decay channel
    widths.push_back(continuum_width);
    total_width += continuum_width;
    levels.push_back(nullptr);
    frags.push_back(nullptr); // A nullptr entry here means a gamma-ray transition
    final_ion_gs_masses.push_back(Migs);
    cpws[nullptr] = gpw;
    total_cws[nullptr] = continuum_width;
    E_c_mins[nullptr] = E_c_min;
    Exf_maxes[nullptr] = Ex;
  }

  // Throw an error if all decays are impossible
  if (total_width <= 0.) throw std::runtime_error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial decay widths are zero.");

  for (size_t i = 0; i < widths.size(); ++i) {
    std::cout << "Exit channel " << i << std::endl;
    std::cout << "  prob = " << widths.at(i) / total_width << std::endl; 

    const TMarleyLevel* lev = levels.at(i);
    const TMarleyFragment* fr = frags.at(i);

    if (fr != nullptr) std::cout << "  " << fr->get_pid() << " emission ";
    else std::cout << "  " << marley_utils::PHOTON << " emission ";
    if (lev != nullptr) std::cout << "to discrete level" << std::endl;
    else std::cout << "to unbound continuum" << std::endl;

    double Exf, m_frag;
    int pid_frag, pid_final;

    double mfgs = final_ion_gs_masses.at(i);

    bool discrete_level = lev != nullptr;
    bool evaporation = fr != nullptr;

    int twoJf = twoJ;
    TMarleyParity Pf = Pi;

    if (discrete_level) {
      Exf = lev->get_energy();
      twoJf = lev->get_two_J();
      Pf = lev->get_parity();
    }
    else {

      double cw = total_cws.at(fr);
      std::function<double(double)> cpw = cpws.at(fr);

      std::function<double(double)> pw = [&cw, &cpw] (double energy)
        -> double { return cpw(energy) / cw; };

      Exf = gen.rejection_sample(pw, E_c_mins.at(fr), Exf_maxes.at(fr));

      // Determine final spin-parity here
      if (evaporation) {
        double Ma = fr->get_mass();
        double Sa = mfgs + Ma - Migs;
        double Mconst = (Ex - Sa) * (Mi + mfgs - Ma);
        double Ea = (Mconst - Exf*(2*mfgs + Exf)) / (2 * Mi);
        int Zf = Zi - fr->get_Z();
        int Af = Ai - fr->get_A();
        TMarleySphericalOpticalModel om(Zf, Af);
        sample_fragment_spin_parity(twoJf, Pf, *fr, om, gen, Exf, Ea);
      }
      else sample_gamma_spin_parity(Zi, Ai, twoJf, Pf, Ex, Exf, gen);
    }

    if (evaporation) {

      pid_frag = fr->get_pid();
      m_frag = fr->get_mass();

      int Zf = Zi - fr->get_Z();
      int Af = Ai - fr->get_A();
      pid_final = marley_utils::get_nucleus_pid(Zf, Af);
    }

    else {
      // Gamma emission
      pid_frag = marley_utils::PHOTON;
      m_frag = TMarleyMassTable::get_particle_mass(pid_frag);

      pid_final = marley_utils::get_nucleus_pid(Zi, Ai);
    }

    TMarleyParticle first_product = TMarleyParticle(pid_frag, m_frag);
    TMarleyParticle second_product = TMarleyParticle(pid_final, mfgs + Exf);

    // TODO: consider changing this to a more realistic model
    // instead of isotropic emissions
    double cos_theta_first = gen.uniform_random_double(-1, 1, true);
    double phi_first = gen.uniform_random_double(0, marley_utils::two_pi, false);

    // Handle the kinematics calculations for this decay
    TMarleyKinematics::two_body_decay(initial_particle, first_product,
      second_product, cos_theta_first, phi_first);

    // Output information about this decay channel
    std::cout << "  Product 1: " << first_product.get_id()
      << " emitted with KE = " << first_product.get_total_energy()
      - first_product.get_mass() << std::endl;
    std::cout << "  Product 2: " << second_product.get_id()
      << " emitted with KE = " << second_product.get_total_energy()
      - second_product.get_mass() << std::endl;
    std::cout << "  Final nuclear Ex = " << Exf << std::endl;
    std::cout << "  Final nuclear 2J = " << twoJf << std::endl;
    std::cout << "  Final nuclear parity = " << Pf << std::endl;
  }

}

void TMarleyNuclearPhysics::hf_test2(int Zi, int Ai,
  const TMarleyParticle& initial_particle, double Ex,
  int twoJ, TMarleyParity Pi,
  TMarleyStructureDatabase& db, TMarleyGenerator& gen,
  size_t num_trials, std::vector<double>& KEs,
  std::vector<int>& f_pids)
{

  KEs.clear();
  f_pids.clear();

  double Mi = initial_particle.get_mass();
  double Migs = Mi - Ex;
  double me = TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);

  double total_width = 0.;

  std::vector<double> widths;
  std::vector<const TMarleyFragment*> frags;
  std::vector<double> final_ion_gs_masses;
  std::vector<TMarleyLevel*> levels;

  std::unordered_map<const TMarleyFragment*, std::function<double(double)> >
    cpws;
  std::unordered_map<const TMarleyFragment*, double> total_cws;
  std::unordered_map<const TMarleyFragment*, double> E_c_mins;
  std::unordered_map<const TMarleyFragment*, double> Exf_maxes;

  for (const auto& f : TMarleyNuclearPhysics::get_fragments()) {
    int Zf = Zi - f.get_Z();
    int Af = Ai - f.get_A();

    TMarleySphericalOpticalModel om(Zf, Af);
    TMarleyDecayScheme* ds = db.get_decay_scheme(Zf, Af);

    // Get information about the current fragment
    int two_s = f.get_two_s();
    TMarleyParity Pa = f.get_parity();
    int fragment_pid = f.get_pid();
  
    // Get various masses that will be needed for fragment emission kinematics
    // (including nuclear recoil)
    double Ma = f.get_mass();
    double Mfgs = TMarleyMassTable::get_atomic_mass(Zf, Af);
    // Fragment atomic number
    int Za = f.get_Z();
    // Approximate the ground state rest energy of the negative ion (with
    // charge Za-) formed when the bare fragment f is emitted by adding Za
    // electron masses to the atomic mass for the final nucleus.
    double Mfgs_ion = Mfgs + Za*me;
    // Get fragment separation energy without a redundant mass table lookup
    double Sa = Mfgs_ion + Ma - Migs;
  
    // Determine the maximum excitation energy available after fragment emission
    // in the final nucleus. Use a simpler functional form as a shortcut if
    // the fragment is neutral [and can therefore be emitted with arbitrarily
    // low kinetic energy due to the absence of a Coulomb barrier])
    double Exf_max;
    double Vc; // Coulomb barrier
    if (f.get_Z() == 0) {
      Vc = 0;
      Exf_max = Ex - Sa;
    }
    else {
      Vc = coulomb_barrier(Zf, Af, Za, f.get_A()); // Ea_min = Vc
      Exf_max = std::sqrt(std::pow(Mi - Ma, 2) - 2*Vc*Mi) - Mfgs_ion;
    }
  
    // The fragment kinetic energy depends on the final level energy. For speed,
    // we precompute here a portion of the fragment energy that doesn't depend on
    // the final level energy.
    double Mconst = (Ex - Sa) * (Mi + Mfgs_ion - Ma);
  
    // Check if emission of this fragment is energetically allowed. If we're below
    // the fragment emission threshold, the partial decay width is zero, so we can
    // skip this fragment rather than slogging through the calculation.
    // We'll also skip it if we're exactly at threshold to avoid numerical problems.
    if (Mconst / (2 * Mi) - Vc <= 0) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0;
  
    // If discrete level data is available for the final nuclide, get decay
    // widths for each accessible level
    if (ds != nullptr) {

      // Get a pointer to the vector of sorted pointers to levels in the decay
      // scheme
      const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

      // Use the maximum discrete level energy from the decay scheme object as the
      // lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
  
      // Loop over the final discrete nuclear levels in order of increasing energy
      // until the new level energy exceeds the maximum value. For each
      // energetically allowed level, if a transition to it for a given fragment
      // orbital angular momentum l and total angular momentum j conserves parity,
      // then compute an optical model transmission coefficient and add it to the total.
      for (const auto& level : *sorted_lps) {
        double Exf = level->get_energy();
        if (Exf < Exf_max) {
          double discrete_width = 0;
          double Ea = (Mconst - Exf*(2*Mfgs_ion + Exf)) / (2 * Mi);
          int twoJf = level->get_two_J();
          TMarleyParity Pf = level->get_parity();
          for (int two_j = std::abs(twoJ - twoJf); two_j <= twoJ + twoJf;
            two_j += 2)
          {
            int j_plus_s = (two_j + two_s) / 2;
            // TODO: consider adding a check that l does not exceed its maximum
            // value l_max that is also used as a cutoff for transmission
            // coefficient calculations in the continuum.
            // For each new iteration, increment l and flip the overall final state parity
            int l = std::abs(two_j - two_s) / 2;
            bool l_is_odd = l % 2;
            TMarleyParity P_final_state = Pf * Pa * TMarleyParity(!l_is_odd);
            for (; l <= j_plus_s; ++l, !P_final_state)
            {
              // Add to the total transmission coefficient if parity is conserved
              if (Pi == P_final_state) discrete_width
                += om.transmission_coefficient(Ea, fragment_pid, two_j, l, two_s,
                DEFAULT_NUMEROV_STEP_SIZE);
            }
          }

          // Store information for this decay channel
          widths.push_back(discrete_width);
          total_width += discrete_width;
          levels.push_back(level);
          frags.push_back(&f);
          final_ion_gs_masses.push_back(Mfgs_ion);
        }
        else break;
      }
    } 

    // If transitions to the energy continuum are possible, include the
    // continuum in the decay channels
    if (Exf_max > E_c_min) {
  
      // Create a forwarding call wrapper for the continuum partial width member
      // function that takes a single argument.
      std::function<double(double)> cpw = std::bind(
        &fragment_continuum_partial_width, Mconst, Mfgs_ion, Mi, twoJ, Pi,
        fragment_pid, two_s, Pa, om, std::placeholders::_1 /*Exf*/);
  
      // Numerically integrate using the call wrapper, the integration bounds, and
      // the number of subintervals
      double continuum_width = marley_utils::num_integrate(cpw, E_c_min,
        Exf_max, DEFAULT_CONTINUUM_SUBINTERVALS);

      // Store information for this decay channel
      widths.push_back(continuum_width);
      total_width += continuum_width;
      levels.push_back(nullptr);
      frags.push_back(&f);
      final_ion_gs_masses.push_back(Mfgs_ion);
      cpws[&f] = cpw;
      total_cws[&f] = continuum_width;
      E_c_mins[&f] = E_c_min;
      Exf_maxes[&f] = Exf_max;
    }
  }

  TMarleyDecayScheme* ds = db.get_decay_scheme(Zi, Ai);

  // Let the continuum go down to 0 MeV unless there is a decay scheme object
  // available for the final nuclide (we'll check this in a second).
  double E_c_min = 0;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if (ds != nullptr) {

    bool initial_spin_is_zero = twoJ == 0;
    // Loop over the final discrete nuclear levels in order of increasing energy
    // until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission coefficient
    // for it
    const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

    // Use the maximum discrete level energy from the decay scheme object as the
    // lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();

    for (const auto& level_f : *sorted_lps) {
      int twoJf = level_f->get_two_J();
      // 0->0 EM transitions aren't allowed due to angular momentum conservation
      // (photons are spin 1), so if the initial and final spins are zero, skip
      // ahead to the next final level.
      if (initial_spin_is_zero && twoJf == 0) continue;
      double Exf = level_f->get_energy();
      if (Exf < Ex) {
        // Approximate the gamma energy by the energy difference between the two levels
        // TODO: consider adding a nuclear recoil correction here
        double e_gamma = Ex - Exf;
        // Determine the type (electric or magnetic) and multipolarity of the gamma
        // transition between these two levels
        int mpol; // Multipolarity
        auto type = determine_gamma_transition_type(twoJ, Pi, level_f, mpol);
        // TODO: allow the user to choose which gamma-ray transition model to use 
        // (Weisskopf single-particle estimates, Brink-Axel strength functions, etc.)
        double discrete_width = gamma_transmission_coefficient(Zi, Ai, type,
          mpol, e_gamma);

        // Store information for this decay channel
        widths.push_back(discrete_width);
        total_width += discrete_width;
        levels.push_back(level_f);
        frags.push_back(nullptr); // A nullptr entry here means a gamma-ray transition
        final_ion_gs_masses.push_back(Migs);
      }
      else break;
    }
  }

  // If gamma transitions to the energy continuum are possible, include them
  // in the possible decay channels
  if (Ex > E_c_min) {

    // Create a forwarding call wrapper for the continuum partial width member
    // function that takes a single argument.
    std::function<double(double)> gpw = std::bind(&gamma_continuum_partial_width,
      Zi, Ai, twoJ, Ex, std::placeholders::_1 /*Exf*/);

    // Numerically integrate using the call wrapper, the integration bounds, and
    // the number of subintervals
    double continuum_width = marley_utils::num_integrate(gpw, E_c_min, Ex,
      DEFAULT_CONTINUUM_SUBINTERVALS);

    // Store information for this decay channel
    widths.push_back(continuum_width);
    total_width += continuum_width;
    levels.push_back(nullptr);
    frags.push_back(nullptr); // A nullptr entry here means a gamma-ray transition
    final_ion_gs_masses.push_back(Migs);
    cpws[nullptr] = gpw;
    total_cws[nullptr] = continuum_width;
    E_c_mins[nullptr] = E_c_min;
    Exf_maxes[nullptr] = Ex;
  }

  // Throw an error if all decays are impossible
  if (total_width <= 0.) throw std::runtime_error(std::string("Cannot ")
    + "continue Hauser-Feshbach decay. All partial decay widths are zero.");

  std::discrete_distribution<size_t> exit_channel_dist(widths.begin(),
    widths.end());

  for (size_t i = 0; i < num_trials; ++i) {

    size_t exit_channel_index = gen.discrete_sample(exit_channel_dist);

    double Exf, m_frag;
    int pid_frag, pid_final;

    double mfgs = final_ion_gs_masses.at(exit_channel_index);

    const TMarleyLevel* lev = levels.at(exit_channel_index);
    const TMarleyFragment* fr = frags.at(exit_channel_index);
    bool discrete_level = lev != nullptr;
    bool evaporation = fr != nullptr;

    if (discrete_level) {
      Exf = lev->get_energy();
      twoJ = lev->get_two_J();
      Pi = lev->get_parity();
    }
    else {

      double cw = total_cws.at(fr);
      std::function<double(double)> cpw = cpws.at(fr);

      std::function<double(double)> pw = [&cw, &cpw] (double energy)
        -> double { return cpw(energy) / cw; };

      Exf = gen.rejection_sample(pw, E_c_mins.at(fr), Exf_maxes.at(fr));

      // Determine final spin-parity here
      if (evaporation) {
        double Ma = fr->get_mass();
        double Sa = mfgs + Ma - Migs;
        double Mconst = (Ex - Sa) * (Mi + mfgs - Ma);
        double Ea = (Mconst - Exf*(2*mfgs + Exf)) / (2 * Mi);
        int Zf = Zi - fr->get_Z();
        int Af = Ai - fr->get_A();
        TMarleySphericalOpticalModel om(Zf, Af);
        sample_fragment_spin_parity(twoJ, Pi, *fr, om, gen, Exf, Ea);
      }
      else sample_gamma_spin_parity(Zi, Ai, twoJ, Pi, Ex, Exf, gen);
    }

    if (evaporation) {

      pid_frag = fr->get_pid();
      m_frag = fr->get_mass();

      int Zf = Zi - fr->get_Z();
      int Af = Ai - fr->get_A();
      pid_final = marley_utils::get_nucleus_pid(Zf, Af);
    }

    else {
      // Gamma emission
      pid_frag = marley_utils::PHOTON;
      m_frag = TMarleyMassTable::get_particle_mass(pid_frag);

      pid_final = marley_utils::get_nucleus_pid(Zi, Ai);
    }

    TMarleyParticle first_product = TMarleyParticle(pid_frag, m_frag);
    TMarleyParticle second_product = TMarleyParticle(pid_final, mfgs + Exf);

    // TODO: consider changing this to a more realistic model
    // instead of isotropic emissions
    double cos_theta_first = gen.uniform_random_double(-1, 1, true);
    double phi_first = gen.uniform_random_double(0, marley_utils::two_pi, false);

    // Handle the kinematics calculations for this decay
    TMarleyKinematics::two_body_decay(initial_particle, first_product,
      second_product, cos_theta_first, phi_first);

    // Store result of this decay
    if (fr != nullptr) f_pids.push_back(fr->get_pid());
    else f_pids.push_back(marley_utils::PHOTON);
    KEs.push_back(first_product.get_total_energy() - first_product.get_mass());

    //std::cout << "Decay #" << i << std::endl;

  }

}

void TMarleyNuclearPhysics::hf_test3(int Zi, int Ai,
  const TMarleyParticle& initial_particle,
  double Ex, int twoJ, TMarleyParity Pi, TMarleyStructureDatabase& db,
  std::unordered_map<const TMarleyFragment*,
    std::function<double(double)> >& funcs,
  std::unordered_map<const TMarleyFragment*, double>& total_c_widths,
  std::unordered_map<const TMarleyFragment*, double>& E_c_mins,
  std::unordered_map<const TMarleyFragment*, double>& Exf_maxes)
{

  funcs.clear();
  total_c_widths.clear();
  E_c_mins.clear();
  Exf_maxes.clear();

  double Mi = initial_particle.get_mass();
  double Migs = Mi - Ex;
  double me = TMarleyMassTable::get_particle_mass(marley_utils::ELECTRON);

  for (const auto& f : TMarleyNuclearPhysics::get_fragments()) {
    int Zf = Zi - f.get_Z();
    int Af = Ai - f.get_A();

    TMarleySphericalOpticalModel om(Zf, Af);
    TMarleyDecayScheme* ds = db.get_decay_scheme(Zf, Af);

    // Get information about the current fragment
    int two_s = f.get_two_s();
    TMarleyParity Pa = f.get_parity();
    int fragment_pid = f.get_pid();
  
    // Get various masses that will be needed for fragment emission kinematics
    // (including nuclear recoil)
    double Ma = f.get_mass();
    double Mfgs = TMarleyMassTable::get_atomic_mass(Zf, Af);
    // Fragment atomic number
    int Za = f.get_Z();
    // Approximate the ground state rest energy of the negative ion (with
    // charge Za-) formed when the bare fragment f is emitted by adding Za
    // electron masses to the atomic mass for the final nucleus.
    double Mfgs_ion = Mfgs + Za*me;
    // Get fragment separation energy without a redundant mass table lookup
    double Sa = Mfgs_ion + Ma - Migs;
  
    // Determine the maximum excitation energy available after fragment emission
    // in the final nucleus. Use a simpler functional form as a shortcut if
    // the fragment is neutral [and can therefore be emitted with arbitrarily
    // low kinetic energy due to the absence of a Coulomb barrier])
    double Exf_max;
    double Vc; // Coulomb barrier
    if (f.get_Z() == 0) {
      Vc = 0;
      Exf_max = Ex - Sa;
    }
    else {
      Vc = coulomb_barrier(Zf, Af, Za, f.get_A()); // Ea_min = Vc
      Exf_max = std::sqrt(std::pow(Mi - Ma, 2) - 2*Vc*Mi) - Mfgs_ion;
    }
  
    // The fragment kinetic energy depends on the final level energy. For speed,
    // we precompute here a portion of the fragment energy that doesn't depend on
    // the final level energy.
    double Mconst = (Ex - Sa) * (Mi + Mfgs_ion - Ma);
  
    // Check if emission of this fragment is energetically allowed. If we're below
    // the fragment emission threshold, the partial decay width is zero, so we can
    // skip this fragment rather than slogging through the calculation.
    // We'll also skip it if we're exactly at threshold to avoid numerical problems.
    if (Mconst / (2 * Mi) - Vc <= 0) continue;

    // Let the continuum go down to 0 MeV unless there is a decay scheme object
    // available for the final nuclide (we'll check this in a second).
    double E_c_min = 0;
  
    // If discrete level data is available for the final nuclide, get decay
    // widths for each accessible level
    if (ds != nullptr) {

      // Get a pointer to the vector of sorted pointers to levels in the decay
      // scheme
      const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

      // Use the maximum discrete level energy from the decay scheme object as the
      // lower bound for the continuum
      // TODO: consider whether this is the best approach
      if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
    } 

    // If transitions to the energy continuum are possible, include the
    // continuum in the decay channels
    if (Exf_max > E_c_min) {
  
      // Create a forwarding call wrapper for the continuum partial width member
      // function that takes a single argument.
      std::function<double(double)> cpw = std::bind(
        &fragment_continuum_partial_width, Mconst, Mfgs_ion, Mi, twoJ, Pi,
        fragment_pid, two_s, Pa, om, std::placeholders::_1 /*Exf*/);
  
      // Numerically integrate using the call wrapper, the integration bounds, and
      // the number of subintervals
      double continuum_width = marley_utils::num_integrate(cpw, E_c_min,
        Exf_max, DEFAULT_CONTINUUM_SUBINTERVALS);

      // Store information for this decay channel
      funcs[&f] = cpw;
      total_c_widths[&f] = continuum_width;
      E_c_mins[&f] = E_c_min;
      Exf_maxes[&f] = Exf_max;
    }
  }

  TMarleyDecayScheme* ds = db.get_decay_scheme(Zi, Ai);

  // Let the continuum go down to 0 MeV unless there is a decay scheme object
  // available for the final nuclide (we'll check this in a second).
  double E_c_min = 0;

  // If discrete level data is available for this nuclide, get gamma decay
  // widths for each accessible level
  if (ds != nullptr) {

    // Loop over the final discrete nuclear levels in order of increasing energy
    // until the new level energy exceeds the maximum value. For each
    // energetically allowed level, compute a gamma ray transmission coefficient
    // for it
    const std::vector<TMarleyLevel*>* sorted_lps = ds->get_sorted_level_pointers();

    // Use the maximum discrete level energy from the decay scheme object as the
    // lower bound for the continuum.
    // TODO: consider whether this is the best approach
    if (sorted_lps->size() > 0) E_c_min = sorted_lps->back()->get_energy();
  }

  // If gamma transitions to the energy continuum are possible, include them
  // in the possible decay channels
  if (Ex > E_c_min) {

    // Create a forwarding call wrapper for the continuum partial width member
    // function that takes a single argument.
    std::function<double(double)> gpw = std::bind(&gamma_continuum_partial_width,
      Zi, Ai, twoJ, Ex, std::placeholders::_1 /*Exf*/);

    // Numerically integrate using the call wrapper, the integration bounds, and
    // the number of subintervals
    double continuum_width = marley_utils::num_integrate(gpw, E_c_min, Ex,
      DEFAULT_CONTINUUM_SUBINTERVALS);

    // Store information for this decay channel
    funcs[nullptr] = gpw;
    total_c_widths[nullptr] = continuum_width;
    E_c_mins[nullptr] = E_c_min;
    Exf_maxes[nullptr] = Ex;
  }

}
