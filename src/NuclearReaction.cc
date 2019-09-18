#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

#include "marley/marley_utils.hh"
#include "marley/Error.hh"
#include "marley/Event.hh"
#include "marley/Generator.hh"
#include "marley/HauserFeshbachDecay.hh"
#include "marley/Level.hh"
#include "marley/Logger.hh"
#include "marley/MatrixElement.hh"
#include "marley/Reaction.hh"

using ME_Type = marley::MatrixElement::TransitionType;

namespace {
  // In cases where no discrete level data are available, a continuum level density
  // is used all the way down to the ground state. To avoid asymptotically approaching
  // Ex = 0 in these cases, the de-excitation cascade will end once the excitation energy of
  // the residual nucleus falls below this (small) value. Excitation energies below this
  // value are considered "close enough" to the ground state for MARLEY not to worry about
  // further de-excitations.
  /// @todo Make this configurable?
  constexpr double CONTINUUM_GS_CUTOFF = 0.001; // MeV
}

marley::NuclearReaction::NuclearReaction(const std::string& filename,
  marley::StructureDatabase& db)
{

  std::regex rx_comment("#.*"); // Matches comment lines

  // Open the reaction data file for parsing
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw marley::Error(std::string("Could not read from the ") +
      "file " + filename);
  }

  std::string line; // String to store the current line
                    // of the reaction data file during parsing

  /// @todo Add more error handling for NuclearReaction::NuclearReaction

  line = marley_utils::get_next_line(file_in, rx_comment, false);

  // Save the reaction description line
  /// @todo Consider automatically generating reaction descriptions from the
  /// particle PDG codes rather than entering the formulae into the reaction
  /// data files by hand.
  description_ = line;

  // Move on to the next non-comment line
  line = marley_utils::get_next_line(file_in, rx_comment, false);

  // Read in the particle IDs
  std::istringstream iss(line);
  iss >> pdg_a_ >> pdg_b_ >> pdg_c_ >> pdg_d_ >> q_d_;

  // For now, set the process type based on whether the projectile
  // and ejectile PDG codes are equal
  /// @todo Note that this will fail to distinguish between NC and EM
  /// reactions when the latter are added to MARLEY. Fix this.
  if ( pdg_a_ == pdg_c_ ) process_type_ = marley::Reaction::ProcessType::NC;
  else process_type_ = marley::Reaction::ProcessType::CC;

  // Get initial and final values of the nuclear
  // charge and mass number from the pids
  Zi_ = (pdg_b_ % 10000000)/10000;
  Ai_ = (pdg_b_ % 10000)/10;
  Zf_ = (pdg_d_ % 10000000)/10000;
  Af_ = (pdg_d_ % 10000)/10;

  const marley::MassTable& mt = marley::MassTable::Instance();

  // Get the particle masses from the mass table
  ma_ = mt.get_particle_mass(pdg_a_);
  mc_ = mt.get_particle_mass(pdg_c_);

  // If the target (particle b) or residue (particle d)
  // has a particle ID greater than 10^9, assume that it
  // is an atom rather than a bare nucleus
  if (pdg_b_ > 1000000000) mb_ = mt.get_atomic_mass(pdg_b_);
  else mb_ = mt.get_particle_mass(pdg_b_);

  if (pdg_d_ > 1000000000) {
    // If particle d is an atom and is ionized as a result of this reaction
    // (e.g., q_d_ != 0), then approximate its ground-state ionized mass by
    // subtracting the appropriate number of electron masses from its atomic
    // (i.e., neutral) ground state mass.
    md_gs_ = mt.get_atomic_mass(pdg_d_)
      - (q_d_ * mt.get_particle_mass(marley_utils::ELECTRON));
  }
  else {
    md_gs_ = mt.get_particle_mass(pdg_d_);
  }

  KEa_threshold_ = (std::pow(mc_ + md_gs_, 2)
    - std::pow(ma_ + mb_, 2))/(2*mb_);

  // Read in all of the level energy (MeV), squared matrix element (B(F) or
  // B(GT) strength), and matrix element type identifier (0 represents
  // B(F), 1 represents B(GT)) triplets.

  // Set the old energy entry to the lowest representable double
  // value. This guarantees that we always read in the first energy
  // value given in the reaction data file
  double old_energy = std::numeric_limits<double>::lowest();
  while (line = marley_utils::get_next_line(file_in, rx_comment, false),
    file_in.good())
  {
    iss.str(line);
    iss.clear();
    /// @todo Consider implementing a sorting procedure rather than strictly
    /// enforcing that energies must be given in ascending order.

    // The order of the entries is important because later uses of the vector of
    // matrix elements assume that they are sorted in order of ascending final
    // level energy.
    double energy, strength;
    int integer_me_type;
    iss >> energy >> strength >> integer_me_type;
    if (old_energy >= energy) throw marley::Error(std::string("Invalid")
      + " reaction dataset. Level energies must be unique and must be"
      + " given in ascending order.");
    // @todo Right now, 0 corresponds to a Fermi transition, and 1 corresponds
    // to a Gamow-Teller transition. As you add new matrix element types,
    // consider changing the convention and its implementation.
    // All of the level pointers owned by the matrix elements will initially be
    // set to nullptr. This may be changed later if discrete level data can be
    // found for the residual nucleus.
    matrix_elements_.emplace_back(energy, strength,
      static_cast<ME_Type>(integer_me_type), nullptr);
    old_energy = energy;
  }

  file_in.close();

  // If discrete level data are available for the residual nucleus, use them to
  // assign values to the level pointers and refine the level energies. If not,
  // keep all of the level pointers nullptr (treat them just like unbound
  // levels) and use the energies given in the reaction dataset.
  marley::DecayScheme* scheme = db.get_decay_scheme(Zf_, Af_);
  if (scheme != nullptr) set_decay_scheme(scheme);

  // Precompute these squared masses for speed
  ma2_ = std::pow(ma_, 2);
  mb2_ = std::pow(mb_, 2);
  mc2_ = std::pow(mc_, 2);
}

// Associate a decay scheme object with this reaction. This will provide
// nuclear structure data for sampling final residue energy levels.
void marley::NuclearReaction::set_decay_scheme(marley::DecayScheme* ds) {

  if (ds == nullptr) throw marley::Error(std::string("Null pointer")
    + " passed to marley::NuclearReaction::set_decay_scheme(marley::Decay"
    + "Scheme* ds)");

  // Check to see if the decay scheme being associated with this
  // reaction is for the correct nuclide. If the PDG code in the decay
  // scheme object does not match the one we'd expect for this reaction's
  // final state nucleus, complain
  int reaction_pdg = marley_utils::get_nucleus_pid(Zf_, Af_);
  int scheme_pdg = marley_utils::get_nucleus_pid(ds->Z(), ds->A());
  if (reaction_pdg != scheme_pdg) throw
    marley::Error(std::string("Nuclear data mismatch: attempted ")
    + "to associate a decay scheme object that has PDG code "
    + std::to_string(scheme_pdg) + " with a reaction object that has PDG code "
    + std::to_string(reaction_pdg));

  // Use the smallest nuclear fragment emission threshold to check for unbound
  // levels. Before looking them up, start by setting the unbound threshold to
  // infinity.
  double unbound_threshold = std::numeric_limits<double>::max();
  MARLEY_LOG_DEBUG() << "unbound_threshold = " << unbound_threshold << '\n';
  for (const auto& f : marley::HauserFeshbachDecay::get_fragments()) {
    double thresh = marley::HauserFeshbachDecay
      ::get_fragment_emission_threshold(Zf_, Af_, f);
    MARLEY_LOG_DEBUG() << f.get_pid() << " emission threshold = " << thresh
      << '\n';
    if (thresh < unbound_threshold) unbound_threshold = thresh;
    MARLEY_LOG_DEBUG() << "unbound_threshold = " << unbound_threshold << '\n';
  }

  // Cycle through each of the level energies given in the reaction dataset.
  for ( auto& mat_el : matrix_elements_ )
  {
    // Get the excitation energy for the level accessed by the transition
    // represented by this matrix element. Use the value from the reaction
    // data file rather than that owned by any previous discrete level
    // assignment. We'll use that value because we need to (re-)assign
    // levels to each matrix element using the DecayScheme ds.
    double en = mat_el.tabulated_level_energy();

    // If the level is above the fragment emission threshold, assign it a null
    // level pointer. Such levels will be handled by a fragment evaporation
    // routine and therefore do not need pointers to level objects describing
    // their de-excitation gammas.
    if ( en > unbound_threshold ) {
      mat_el.set_level(nullptr);
      continue;
    }

    // For each energy, find a pointer to the level with the closest energy
    // owned by the decay scheme object.
    marley::Level* plevel = ds->get_pointer_to_closest_level( en );
    MARLEY_LOG_DEBUG() << "reaction level at " << en
      << " MeV was matched to the decay scheme level at "
      << plevel->energy() << " MeV";

    // Complain if there are duplicates (if there are duplicates, we'll have
    // two different B(F) + B(GT) values for the same level object)
    const auto begin = matrix_elements_.cbegin();
    const auto end = matrix_elements_.cend();
    const auto found = std::find_if(begin, end,
      [plevel](const marley::MatrixElement& me) -> bool
      { return plevel == me.level(); });
    if ( found != end )
    {
      // One of the matrix elements already uses a level pointer equal to plevel
      throw marley::Error("Reaction dataset gives two"
        " level energies that refer to the same DecayScheme level at "
        + std::to_string( plevel->energy() ) + " MeV");
    }

    /// @todo Add check to see if the energy of the chosen level is very
    /// different from the energy given in the reaction dataset. If it is, the
    /// level matchup is likely incorrect.

    // Add the level pointer to the list
    mat_el.set_level( plevel );
  }
}

// Fermi function used to calculate cross-sections
// The form used here is based on http://en.wikipedia.org/wiki/Beta_decay
// but rewritten for convenient use inside this class.
// Input: beta_c (3-momentum magnitude of particle c / total energy of particle
// c), where we assume that the ejectile (particle c) is the light product from
// 2-2 scattering.
double marley::NuclearReaction::fermi_function(double beta_c) const {

  // Don't bother to calculate anything if the light product from
  // this reaction (particle c) is not an electron nor a positron.
  // This situation occurs for neutral current reactions, for example.
  bool electron = (pdg_c_ == marley_utils::ELECTRON);
  if (!electron && pdg_c_ != marley_utils::POSITRON) return 1.;

  // Lorentz factor gamma for particle c
  double gamma_c = std::pow(1. - beta_c*beta_c, -marley_utils::ONE_HALF);

  double s = std::sqrt(1. - std::pow(marley_utils::alpha * Zf_, 2));

  // Estimate the nuclear radius using r_n = r_0*A^(1/3)
  double r_n = marley_utils::r0 * std::pow(Af_, marley_utils::ONE_THIRD);

  // Convert the estimated nuclear radius to natural units (MeV^(-1))
  double rho = r_n / marley_utils::hbar_c;

  double eta = marley_utils::alpha * Zf_ / beta_c;

  // Adjust the value of eta if the light product from this reaction is a
  // positron.
  if (!electron) eta *= -1;

  // Complex variable for the gamma function
  std::complex<double> a(s, eta);
  double b = std::tgamma(1+2*s);

  return 2 * (1 + s) * std::pow(2*beta_c*gamma_c*rho*mc_, 2*s-2)
    * std::exp(marley_utils::pi*eta) * std::norm(marley_utils::gamma(a))
    / std::pow(b, 2);
}

// Return the maximum residue excitation energy E_level that can
// be achieved in the lab frame for a given projectile kinetic energy KEa
// (this corresponds to the final particles all being produced
// at rest in the CM frame). This maximum level energy is used
// to find the allowed levels when creating events.
double marley::NuclearReaction::max_level_energy(double KEa) const {
  // Calculate the total CM frame energy using known quantities
  // from the lab frame
  double E_CM = std::sqrt(std::pow(ma_ + mb_, 2) + 2*mb_*KEa);
  // The maximum level energy is achieved when the final state
  // particles are produced at rest in the CM frame. Subtracting
  // the ground-state rest masses of particles c and d from the
  // total CM energy leaves us with the energy available to create
  // an excited level in the residue (particle d).
  return E_CM - mc_ - md_gs_;
}

double marley::NuclearReaction::threshold_kinetic_energy() const {
  return KEa_threshold_;
}

// Creates an event object by sampling the appropriate quantities and
// performing kinematic calculations
marley::Event marley::NuclearReaction::create_event(int pdg_a, double KEa,
  marley::Generator& gen)
{
  // Check that the projectile supplied to this event is correct. If not, alert
  // the user that this event does not use the requested projectile.
  if (pdg_a != pdg_a_) throw marley::Error(std::string("Could")
    + " not create this event. The requested projectile particle ID, "
    + std::to_string(pdg_a) + ", does not match the projectile"
    + " particle ID, " + std::to_string(pdg_a_) + ", in the reaction dataset.");

  // Sample a final residue energy level. First, check to make sure the given
  // projectile energy is above threshold for this reaction.
  if (KEa < KEa_threshold_) throw std::range_error(std::string("Could")
    + " not create this event. Projectile kinetic energy " + std::to_string(KEa)
    + " MeV is below the threshold value " + std::to_string(KEa_threshold_)
    + " MeV.");

  /// @todo Add more error checks to NuclearReaction::create_event as necessary

  // Create an empty vector of sampling weights (partial total cross
  // sections to each kinematically accessible final level)
  std::vector<double> level_weights;

  // Create a discrete distribution object for level sampling.
  // Its default constructor creates a single weight of 1.
  // We will always explicitly give it weights to use when sampling
  // levels, so we won't worry about its default behavior.
  static std::discrete_distribution<size_t> ldist;

  // Compute the total cross section for a transition to each individual nuclear
  // level, and save the results in the level_weights vector (which will be
  // cleared by summed_xs_helper() before being loaded with the cross sections).
  // The summed_xs_helper() method can also be used for differential
  // (d\sigma/d\cos\theta_c^{CM}) cross sections, so supply a dummy cos_theta_c_cm
  // value and request total cross sections by setting the last argument to false.
  double dummy = 0.;
  double sum_of_xsecs = summed_xs_helper(pdg_a, KEa, dummy,
    &level_weights, false);

  // Note that the elements in matrix_elements_ are given in order of
  // increasing excitation energy (this is currently enforced by the reaction
  // data format and is checked during parsing). This ensures that we can
  // sample a matrix element index from level_weights (which is populated in
  // the same order by summed_xs_helper()) and have it refer to the correct
  // object.

  // Complain if none of the levels we have data for are kinematically allowed
  if ( level_weights.empty() ) {
    throw marley::Error("Could not create this event. The DecayScheme object"
      " associated with this reaction does not contain data for any"
      " kinematically accessible levels for a projectile kinetic energy of "
      + std::to_string(KEa) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) ) + " MeV).");
  }

  // Complain if the total cross section (the sum of all partial level cross
  // sections) is zero or negative (the latter is just to cover all possibilities).
  if ( sum_of_xsecs <= 0. ) {
    throw marley::Error("Could not create this event. All kinematically"
      " accessible levels for a projectile kinetic energy of "
      + std::to_string(KEa) + " MeV (max E_level = "
      + std::to_string( max_level_energy(KEa) )
      + " MeV) have vanishing matrix elements.");
  }

  // Create a list of parameters used to supply the weights to our discrete
  // level sampling distribution
  std::discrete_distribution<size_t>::param_type params(level_weights.begin(),
    level_weights.end());

  // Sample a matrix_element using our discrete distribution and the
  // current set of weights
  size_t me_index = gen.sample_from_distribution(ldist, params);

  const auto& sampled_matrix_el = matrix_elements_.at( me_index );

  // Get a pointer to the selected level (this will be nullptr
  // if the transition is to the unbound continuum)
  const marley::Level* plevel = sampled_matrix_el.level();

  // Get the energy of the selected level.
  double E_level = sampled_matrix_el.level_energy();

  // Update the residue mass based on its excitation energy for the current
  // event
  md_ = md_gs_ + E_level;
  md2_ = std::pow(md_, 2);

  // Compute Mandelstam s, the ejectile's CM frame total energy, the magnitude
  // of the ejectile's CM frame 3-momentum, and the residue's CM frame total
  // energy.
  double s, Ec_cm, pc_cm, Ed_cm;
  two_two_scatter(KEa, s, Ec_cm, pc_cm, Ed_cm);

  // Determine the CM frame velocity of the ejectile
  double beta_c_cm = pc_cm / Ec_cm;

  // Sample a CM frame scattering cosine for the ejectile.
  double cos_theta_c_cm = sample_cos_theta_c_cm(sampled_matrix_el,
    beta_c_cm, gen);

  // Sample a CM frame azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the matrix elements are azimuthally invariant
  double phi_c_cm = gen.uniform_random_double(0., 2.*marley_utils::pi, false);

  // Create the preliminary event object (after 2-2 scattering, but before
  // de-excitation of the residual nucleus)
  marley::Event event = make_event_object(KEa, pc_cm, cos_theta_c_cm, phi_c_cm,
    Ec_cm, Ed_cm, E_level);

  // Get a reference to the residue so that we can handle its de-excitation
  marley::Particle& residue = event.residue();

  bool continuum = (plevel == nullptr);
  int Z = Zf_;
  int A = Af_;
  // Load the initial excitation energy into Ex.  This variable will be changed
  // during every step of the Hauser-Feshbach decay cascade.
  double Ex = E_level;

  if ( continuum ) {

    double cutoff = CONTINUUM_GS_CUTOFF;

    // Load the initial twoJ and parity values into twoJ and P.  These
    // variables will be changed during every step of the Hauser-Feshbach decay
    // cascade.
    // Right now, matrix element types are represented with 0 <-> Fermi, 1 <->
    // Gamow-Teller.  Since the transitions are from the 40Ar ground state (Jpi
    // = 0+), these also correspond to the 40K* state spins.
    /// @todo Come up with a better way of determining the Jpi values that will
    /// work for forbidden transition operators.
    int twoJ;
    if ( sampled_matrix_el.type() == ME_Type::FERMI) twoJ = 0;
    else if ( sampled_matrix_el.type() == ME_Type::GAMOW_TELLER) twoJ = 2;
    else throw marley::Error("Unrecognized matrix element type encountered"
      " during a continuum decay in marley::NuclearReaction::create_event()");

    // Fermi transition gives 0+ -> 0+, GT transition gives 0+ -> 1+
    /// @todo handle target nuclei that are not initially 0+
    /// @todo include possibility of negative parity here.
    marley::Parity P(true); // positive parity
    marley::Particle first, second;

    // The selected level is unbound, so handle its de-excitation using
    // the Hauser-Feshbach statistical model.
    while ( continuum && Ex > cutoff ) {
      marley::HauserFeshbachDecay hfd(residue, Ex, twoJ, P, gen);
      MARLEY_LOG_DEBUG() << hfd;
      continuum = hfd.do_decay(Ex, twoJ, P, first, second);

      MARLEY_LOG_DEBUG() << "Hauser-Feshbach decay to " << first.pdg_code()
        << " and " << second.pdg_code();
      MARLEY_LOG_DEBUG() << second.pdg_code() << " is at Ex = " << Ex << " MeV.";

      residue = second;
      Z = marley_utils::get_particle_Z(residue.pdg_code());
      A = marley_utils::get_particle_A(residue.pdg_code());
      event.add_final_particle(first);
    }
  }

  if ( !continuum ) {
    // Either the selected initial level was bound (so it will only decay via
    // gamma emission) or the Hauser-Feshbach decay process has now accessed a
    // bound level in the residual nucleus. In either case, use gamma-ray decay
    // scheme data to sample the de-excitation gammas and add them to this
    // event's final particle list.
    marley::DecayScheme* dec_scheme
      = gen.get_structure_db().get_decay_scheme(Z, A);
    dec_scheme->do_cascade(*dec_scheme->get_pointer_to_closest_level(Ex),
      event, gen, residue.charge());
  }

  // Return the completed event object
  return event;
}

// Compute the total reaction cross section (summed over all final nuclear levels)
// in units of MeV^(-2) using the center of momentum frame.
double marley::NuclearReaction::total_xs(int pdg_a, double KEa) {
  double dummy_cos_theta = 0.;
  return summed_xs_helper(pdg_a, KEa, dummy_cos_theta, nullptr, false);
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM}
// summed over all final nuclear levels. This is done in units of MeV^(-2)
// using the center of momentum frame.
double marley::NuclearReaction::diff_xs(int pdg_a, double KEa,
  double cos_theta_c_cm)
{
  return summed_xs_helper(pdg_a, KEa, cos_theta_c_cm, nullptr, true);
}

// Compute the differential cross section d\sigma / d\cos\theta_c^{CM} for a
// transition to a particular final nuclear level. This is done in units of
// MeV^(-2) using the center of momentum frame.
double marley::NuclearReaction::diff_xs(const marley::MatrixElement& mat_el,
  double KEa, double cos_theta_c_cm)
{
  // Check that the scattering cosine is within the physically meaningful range
  if ( std::abs(cos_theta_c_cm) > 1. ) return 0.;
  double beta_c_cm;
  double xsec = total_xs(mat_el, KEa, beta_c_cm, true);
  xsec *= mat_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm);
  return xsec;
}

// Helper function for total_xs and diff_xs()
double marley::NuclearReaction::summed_xs_helper(int pdg_a, double KEa,
  double cos_theta_c_cm, std::vector<double>* level_xsecs, bool differential)
{
  // Check that the projectile supplied to this event is correct. If not,
  // return a total cross section of zero since this reaction is not available
  // for the given projectile.
  /// @todo Consider whether you should use an exception if pdg_a != pdg_a_
  if (pdg_a != pdg_a_) return 0.;

  // If we're evaluating a differential cross section, check that the
  // scattering cosine is within the physically meaningful range. If it's
  // not, then just return 0.
  if ( differential && std::abs(cos_theta_c_cm) > 1. ) return 0.;

  // If we've been passed a vector to load with the partial cross sections
  // to each nuclear level, then clear it before storing them
  if ( level_xsecs ) level_xsecs->clear();

  double max_E_level = max_level_energy( KEa );
  double xsec = 0.;
  for ( const auto& mat_el : matrix_elements_ ) {

    // Get the excitation energy for the current level
    double level_energy = mat_el.level_energy();

    // Exit the loop early if you reach a level with an energy that's too high
    if ( level_energy > max_E_level ) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // computing the total xs.
    if ( mat_el.strength() != 0. ) {

      // Set the check_max_E_level flag to false when calculating the total
      // cross section for this level (we've already verified that the
      // current level is kinematically accessible in the check against
      // max_E_level above)
      double beta_c_cm = 0.;
      double partial_xsec = total_xs(mat_el, KEa, beta_c_cm, false);

      // If a differential cross section (d\sigma / d\cos\theta_{CM})
      // is desired, then multiply by the appropriate angular factor
      if ( differential ) {
        partial_xsec *= mat_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm);
      }

      if ( std::isnan(partial_xsec) ) {
        MARLEY_LOG_WARNING() << "Partial cross section for reaction "
          << description_ << " gave NaN result.";
        MARLEY_LOG_DEBUG() << "Parameters were level energy = "
          << mat_el.level_energy() << " MeV, projectile kinetic energy = "
          << KEa << " MeV, and reduced matrix element = " << mat_el.strength();
        MARLEY_LOG_DEBUG() << "The partial cross section to this level"
          << " will be set to zero.";
        partial_xsec = 0.;
      }

      xsec += partial_xsec;

      // Store the partial cross section to the current individual nuclear
      // level if needed (i.e., if level_xsecs is not nullptr)
      if ( level_xsecs ) level_xsecs->push_back( partial_xsec );
    }
  }

  return xsec;
}

// Compute the total reaction cross section (in MeV^(-2)) for a transition to a
// particular nuclear level using the center of momentum frame
double marley::NuclearReaction::total_xs(const marley::MatrixElement& me,
  double KEa, double& beta_c_cm, bool check_max_E_level) const
{
  // Don't bother to compute anything if the matrix element vanishes for this
  // level
  if ( me.strength() == 0. ) return 0.;

  // Also don't proceed further if the reaction is below threshold (equivalently,
  // if the requested level excitation energy E_level exceeds that maximum
  // kinematically-allowed value). To avoid redundant checks of the threshold,
  // skip this check if check_max_E_level is set to false.
  if ( check_max_E_level ) {
    double max_E_level = max_level_energy( KEa );
    if ( me.level_energy() > max_E_level ) return 0.;
  }

  // The final nuclear mass (before nuclear de-excitations) is the sum of the
  // ground state residue mass plus the excitation energy of the accessed level
  double md2 = std::pow(md_gs_ + me.level_energy(), 2);

  // Compute Mandelstam s (the square of the total CM frame energy)
  double s = std::pow(ma_ + mb_, 2) + 2.*mb_*KEa;
  double sqrt_s = std::sqrt(s);

  // Compute CM frame total energies for two of the particles. Also
  // compute the magnitude of the ejectile CM frame momentum.
  double Eb_cm = (s + mb2_ - ma2_) / (2. * sqrt_s);
  double Ec_cm = (s + mc2_ - md2) / (2. * sqrt_s);
  double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc2_);

  // Compute the (dimensionless) speed of the ejectile in the CM frame
  beta_c_cm = pc_cm / Ec_cm;

  // CM frame total energy of the nuclear residue
  double Ed_cm = sqrt_s - Ec_cm;

  // Dot product of the four-momenta of particles c and d
  double pc_dot_pd = Ed_cm*Ec_cm + std::pow(pc_cm, 2);

  // Relative speed of particles c and d, computed with a manifestly
  // Lorentz-invariant expression
  double beta_rel_cd = marley_utils::real_sqrt(
    std::pow(pc_dot_pd, 2) - mc2_*md2) / pc_dot_pd;

  // Common factors for the allowed approximation total cross sections
  // for both CC and NC reactions
  double total_xsec = (marley_utils::GF2 / marley_utils::pi)
    * ( Eb_cm * Ed_cm / s ) * Ec_cm * pc_cm * me.strength();

  if ( process_type_ == ProcessType::CC ) {

    // Calculate a Coulomb correction factor using either a Fermi function
    // or the effective momentum approximation
    double factor_C = coulomb_correction_factor( beta_rel_cd );
    total_xsec *= marley_utils::Vud2 * factor_C;
  }
  else if ( process_type_ == ProcessType::NC )
  {
    // For NC, extra factors are only needed for Fermi transitions (which
    // correspond to CEvNS since they can only access the nuclear ground state)
    if ( me.type() == ME_Type::FERMI ) {
      double Q_w = weak_nuclear_charge();
      total_xsec *= 0.25*std::pow(Q_w, 2);
    }
  }
  else throw marley::Error("Unrecognized process type encountered in"
    " marley::NuclearReaction::total_xs()");

  return total_xsec;
}

// Sample an ejectile scattering cosine in the CM frame.
double marley::NuclearReaction::sample_cos_theta_c_cm(
  const marley::MatrixElement& matrix_el, double beta_c_cm,
  marley::Generator& gen) const
{
  // To avoid wasting time searching for the maximum of these distributions
  // for rejection sampling, set the maximum to the known value before
  // proceeding.
  double max;
  if ( matrix_el.type() == ME_Type::FERMI ) {
    // B(F)
    max = matrix_el.cos_theta_pdf(1., beta_c_cm);
  }
  else if (matrix_el.type() == ME_Type::GAMOW_TELLER) {
    // B(GT)
    max = matrix_el.cos_theta_pdf(-1., beta_c_cm);
  }
  else throw marley::Error("Unrecognized matrix element type "
    + std::to_string(matrix_el.type()) + " encountered while sampling a"
    " CM frame scattering angle");

  // Sample a CM frame scattering cosine using the appropriate distribution for
  // this matrix element.
  return gen.rejection_sample(
    [&matrix_el, &beta_c_cm](double cos_theta_c_cm)
    -> double { return matrix_el.cos_theta_pdf(cos_theta_c_cm, beta_c_cm); },
    -1., 1., max);
}

marley::Event marley::NuclearReaction::make_event_object(double KEa,
  double pc_cm, double cos_theta_c_cm, double phi_c_cm, double Ec_cm,
  double Ed_cm, double E_level)
{
  marley::Event event = marley::Reaction::make_event_object(KEa, pc_cm,
    cos_theta_c_cm, phi_c_cm, Ec_cm, Ed_cm, E_level);

  // Assume that the target is a neutral atom (q_b = 0)
  event.target().set_charge(0);

  // Assign the correct charge to the residue
  event.residue().set_charge(q_d_);

  return event;
}

double marley::NuclearReaction::coulomb_correction_factor(double beta_rel_cd)
  const
{
  // Don't bother to calculate anything if the light product from
  // this reaction (particle c) is not an electron nor a positron.
  // This situation occurs for neutral current reactions, for example.
  /// \todo Revisit this when you have more reaction data files
  bool electron = (pdg_c_ == marley_utils::ELECTRON);
  if (!electron && pdg_c_ != marley_utils::POSITRON) return 1.;

  // Fermi function approach to the Coulomb correction
  double fermi_func = fermi_function(beta_rel_cd);

  // Effective momentum approximation for the Coulomb correction
  bool EMA_ok = false;
  double factor_EMA = ema_factor(beta_rel_cd, EMA_ok);

  // If the effective momentum approximation is invalid because subtracting
  // off the Coulomb potential brings the reaction below threshold, then
  // just use the Fermi function
  if ( !EMA_ok ) return fermi_func;


  // Otherwise, choose the larger of the two factors for antineutrinos, and the
  // smaller of the two factors for neutrinos
  bool is_antineutrino = pdg_a_ < 0;

  double correction_factor = 1.;
  if ( is_antineutrino ) {
    correction_factor = std::max( fermi_func, factor_EMA );
  }
  else {
    correction_factor = std::min( fermi_func, factor_EMA );
  }

  return correction_factor;
}

// Effective momentum approximation for the Coulomb correction factor
double marley::NuclearReaction::ema_factor(double beta_rel_cd, bool& ok) const
{
  // Don't bother to calculate anything if the light product from
  // this reaction (particle c) is not an electron nor a positron.
  // This situation occurs for neutral current reactions, for example.
  bool electron = (pdg_c_ == marley_utils::ELECTRON);
  if (!electron && pdg_c_ != marley_utils::POSITRON) return 1.;

  // Like the Fermi function, this approximation uses a static nuclear Coulomb
  // potential (a sphere at the origin). Typically nuclear recoil is neglected,
  // allowing one to compute the effective lepton momentum in the lab frame. In
  // MARLEY's case, we do this by calculating it in the rest frame of the final
  // nucleus ("FNR" frame). We already have the relative speed of the final
  // nucleus and outgoing lepton, so this is easy.
  double gamma_rel_cd = std::pow(
    1. - std::pow(beta_rel_cd, 2), -marley_utils::ONE_HALF);

  // Check for numerical errors from the square root
  if ( !std::isfinite(gamma_rel_cd) ) {
    MARLEY_LOG_WARNING() << "Invalid beta_rel = " << beta_rel_cd
      << " encountered in marley::NuclearReaction::ema_factor()";
  }

  // Total energy of the outgoing lepton in the FNR frame
  double E_c_FNR = gamma_rel_cd * mc_;

  // Approximate Coulomb potential
  double Vc = -3.*Zf_*marley_utils::alpha / (2. * marley_utils::r0
    * std::pow(Af_, marley_utils::ONE_THIRD));
  if ( !electron ) Vc *= -1;

  // Effective FNR frame total energy
  double E_c_FNR_eff = E_c_FNR - Vc;

  // If subtracting off the Coulomb potential drops the effective energy
  // below the lepton mass, then the expression for the effective momentum
  // will give an imaginary value. Signal this by setting the "ok" flag to
  // false.
  ok = (E_c_FNR_eff >= mc_);

  // Lepton momentum in FNR frame
  double p_c_FNR = marley_utils::real_sqrt( std::pow(E_c_FNR, 2) - mc2_ );

  // Effective momentum in FNR frame
  double p_c_FNR_eff = marley_utils::real_sqrt(
    std::pow(E_c_FNR_eff, 2) - mc2_ );

  // Coulomb correction factor
  double f_EMA2 = std::pow(p_c_FNR_eff / p_c_FNR, 2);

  return f_EMA2;
}

// Factor that appears in the cross section for coherent elastic
// neutrino-nucleus scattering (CEvNS), which corresponds to the Fermi
// component of NC scattering under the allowed approximation
double marley::NuclearReaction::weak_nuclear_charge() const
{
  int Ni = Ai_ - Zi_;
  double Qw = Ni - (1. - 4.*marley_utils::sin2thetaw)*Zi_;
  return Qw;
}
