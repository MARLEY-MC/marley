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
#include "marley/Reaction.hh"

marley::NuclearReaction::NuclearReaction(std::string filename,
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
    int me_type;
    iss >> energy >> strength >> me_type;
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
      static_cast<marley::MatrixElement::TransitionType>(me_type), nullptr);
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
  for (auto& mat_el : matrix_elements_)
  {
    double en = mat_el.level_energy();
    // If the level is above the fragment emission threshold, assign it a null
    // level pointer. Such levels will be handled by a fragment evaporation
    // routine and therefore do not need pointers to level objects describing
    // their de-excitation gammas.
    if (en > unbound_threshold) {
      mat_el.set_level(nullptr);
      continue;
    }

    // For each energy, find a pointer to the level with the closest energy
    // owned by the decay scheme object.
    marley::Level* plevel = ds->get_pointer_to_closest_level(en);
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
    if (found != end)
    {
      // One of the matrix elements already uses a level pointer equal to plevel
      throw marley::Error(std::string("Reaction dataset gives two")
        + " level energies that refer to the same DecayScheme level at "
        + std::to_string(plevel->energy()) + " MeV");
    }

    /// @todo Add check to see if the energy of the chosen level is very
    /// different from the energy given in the reaction dataset. If it is, the
    /// level matchup is likely incorrect.

    // Add the level pointer to the list
    mat_el.set_level(plevel);
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
  double gamma_c = std::pow(1 - beta_c*beta_c, -marley_utils::ONE_HALF);

  double s = std::sqrt(1 - std::pow(marley_utils::alpha * Zf_, 2));

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

  // Calculate the maximum value of E_level that is
  // kinematically allowed
  double max_E_level = max_level_energy(KEa);

  // Create an empty vector of sampling weights
  std::vector<double> level_weights;

  // Create a discrete distribution object for level sampling.
  // Its default constructor creates a single weight of 1.
  // We will always explicitly give it weights to use when sampling
  // levels, so we won't worry about its default behavior.
  static std::discrete_distribution<size_t> ldist;

  // The level pointers owned by the matrix elements are ordered by increasing
  // energy (this is currently enforced by the reaction data format and is
  // checked during parsing). Iterate over the levels, assigning the total cross
  // section for each level as its weight until you reach the end of
  // the matrix elements or a level that is kinematically forbidden. If the
  // matrix element B(F) + B(GT) for a level vanishes, skip computing the total
  // cross section and assign a weight of zero for efficiency.
  bool at_least_one_nonzero_matrix_el = false;

  for (const auto& mat_el : matrix_elements_) {

    double level_energy;

    // Get the excitation energy for the current level
    if (mat_el.level() != nullptr) {
      level_energy = mat_el.level()->energy();
    }
    else {
      // Level is unbound, so just use the energy given in the reaction dataset
      // rather than trying to find its DecayScheme version
      level_energy = mat_el.level_energy();
    }

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // calling total_xs. This will avoid unnecessary numerical integrations.
    double strength = mat_el.strength();
    double xs;

    if (strength == 0.) {
      xs = 0.;
    }
    else {

      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that
      // std::discrete_distribution automatically normalizes the weights, so we
      // don't have to do that ourselves.
      xs = total_xs(level_energy, KEa, strength);
      if (std::isnan(xs)) {
        MARLEY_LOG_WARNING() << "Partial cross section for reaction "
          << description_ << " gave NaN result.";
        MARLEY_LOG_DEBUG() << "Parameters were level energy = " << level_energy
          << " MeV, projectile kinetic energy = " << KEa
          << " MeV, and matrix element = " << strength;
        MARLEY_LOG_DEBUG() << "The partial cross section to this level"
          << " will be set to zero for this event.";
        xs = 0.;
      }
      if (!at_least_one_nonzero_matrix_el && xs != 0)
        at_least_one_nonzero_matrix_el = true;
    }
    level_weights.push_back(xs);
  }

  // Complain if none of the levels we have data for are kinematically allowed
  if (level_weights.empty()) {
    throw marley::Error(std::string("Could not create this event. The ")
      + "DecayScheme object associated with this reaction "
      + "does not contain data for any kinematically accessible levels "
      + "for a projectile kinetic energy of " + std::to_string(KEa)
      + " MeV (max E_level = " + std::to_string(max_E_level)
      + " MeV).");
  }

  // Complain if none of the levels have nonzero weight (this means that
  // all kinematically allowed levels have a vanishing matrix element)
  if (!at_least_one_nonzero_matrix_el) {
    throw marley::Error(std::string("Could not create this event. All ")
      + "kinematically accessible levels for a projectile kinetic energy of "
      + std::to_string(KEa) + " MeV (max E_level = "
      + std::to_string(max_E_level) + " MeV) have vanishing matrix elements.");
  }

  // Create a list of parameters used to supply the weights to our discrete
  // level sampling distribution
  std::discrete_distribution<size_t>::param_type params(level_weights.begin(),
    level_weights.end());

  // Sample a matrix_element using our discrete distribution and the
  // current set of weights
  size_t me_index = gen.sample_discrete(ldist, params);

  const auto& sampled_matrix_el = matrix_elements_.at(me_index);

  // Get a pointer to the selected level
  const marley::Level* plevel = sampled_matrix_el.level();

  // Get the energy of the selected level.
  double E_level;
  if (plevel != nullptr) {
    E_level = plevel->energy();
  }
  else {
    // Level is unbound, so use the level energy found in the reaction dataset
    E_level = sampled_matrix_el.level_energy();
  }

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
  double phi_c_cm = gen.uniform_random_double(0, 2*marley_utils::pi, false);

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

  if (continuum) {
    double cutoff = 0.001; /// @todo Remove hard-coded cutoff value
    // Load the initial twoJ and parity values into twoJ and P.  These
    // variables will be changed during every step of the Hauser-Feshbach decay
    // cascade.
    // Right now, matrix element types are represented with 0 <-> Fermi, 1 <->
    // Gamow-Teller.  Since the transitions are from the 40Ar ground state (Jpi
    // = 0+), these also correspond to the 40K* state spins.
    /// @todo Come up with a better way of determining the Jpi values that will
    /// work for forbidden transition operators.

    marley::MatrixElement::TransitionType me_type = sampled_matrix_el.type();
    int twoJ;
    if (me_type == marley::MatrixElement::TransitionType::FERMI) twoJ = 0;
    else if (me_type == marley::MatrixElement::TransitionType::GAMOW_TELLER) {
      twoJ = 2;
    }
    else throw marley::Error("Unrecognized matrix element type encountered"
      " during a continuum decay in marley::NuclearReaction::create_event()");

    // Fermi transition gives 0+ -> 0+, GT transition gives 0+ -> 1+
    /// @todo include possibility of negative parity here.
    marley::Parity P(true); // positive parity
    marley::Particle first, second;
    // The selected level is unbound, so handle its de-excitation using
    // the Hauser-Feshbach statistical model.
    while (continuum && Ex > cutoff) {
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

  if (!continuum) {
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

// Compute the total reaction cross section (including all final nuclear levels)
// in units of MeV^(-2) using the center of mass frame.
double marley::NuclearReaction::total_xs(int pdg_a, double KEa) {

  // Check that the projectile supplied to this event is correct. If not,
  // return a total cross section of zero since this reaction is not available
  // for the given projectile.
  /// @todo Consider whether you should use an exception if pdg_a != pdg_a_
  if (pdg_a != pdg_a_) return 0.;

  double max_E_level = max_level_energy(KEa);
  double xs = 0;
  for (const auto& mat_el : matrix_elements_) {

    double level_energy;

    // Get the excitation energy for the current level
    if (mat_el.level() != nullptr) {
      level_energy = mat_el.level()->energy();
    }
    else {
      // Level is unbound, so just use the energy given in the reaction dataset
      // rather than trying to find its DecayScheme version
      level_energy = mat_el.level_energy();
    }

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // computing the total xs.
    double matrix_el = mat_el.strength();

    if (matrix_el != 0.) {
      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that
      // std::discrete_distribution automatically normalizes the weights, so we
      // don't have to do that ourselves.
      double md2 = std::pow(md_gs_ + level_energy, 2);

      // Compute Mandelstam s (the square of the total CM frame energy)
      double s = std::pow(ma_ + mb_, 2) + 2*mb_*KEa;
      double sqrt_s = std::sqrt(s);

      // Compute CM frame energies for two of the particles. Also
      // compute the ejectile CM frame momentum.
      double Eb_cm = (s + mb2_ - ma2_) / (2 * sqrt_s);
      double Ec_cm = (s + mc2_ - md2) / (2 * sqrt_s);
      double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc2_);
      double beta_c_cm = pc_cm / Ec_cm;

      xs += fermi_function(beta_c_cm) * matrix_el * pc_cm * Ec_cm
        * Eb_cm * (sqrt_s - Ec_cm) / s;
    }
  }

  // Include the quark mixing matrix element Vud if this is a charged-current
  // reaction. For now, we determine that by asking whether the light product
  // (particle c) is an electron or a positron. If it is neither, we don't
  // multiply by marley_utils::Vud^2.
  if (pdg_c_ == marley_utils::ELECTRON || pdg_c_ == marley_utils::POSITRON) {
    xs *= std::pow(marley_utils::Vud, 2);
  }

  return std::pow(marley_utils::GF, 2) * xs / marley_utils::pi;
}

// Compute the total reaction cross section for a given level (ignoring some
// constants for speed, since this version of the function is intended for use
// in sampling only) using the center of mass frame
double marley::NuclearReaction::total_xs(double E_level, double KEa,
  double matrix_element) const
{
  // Don't bother to compute anything if the matrix element vanishes for this
  // level
  if (matrix_element == 0.) return 0.;
  else {
    double md2 = std::pow(md_gs_ + E_level, 2);

    // Compute Mandelstam s (the square of the total CM frame energy)
    double s = std::pow(ma_ + mb_, 2) + 2*mb_*KEa;
    double sqrt_s = std::sqrt(s);

    // Compute CM frame energies for two of the particles. Also
    // compute the ejectile CM frame momentum.
    double Eb_cm = (s + mb2_ - ma2_) / (2 * sqrt_s);
    double Ec_cm = (s + mc2_ - md2) / (2 * sqrt_s);
    double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc2_);
    double beta_c_cm = pc_cm / Ec_cm;

    return fermi_function(beta_c_cm) * matrix_element * pc_cm * Ec_cm
      * Eb_cm * (sqrt_s - Ec_cm) / s;
  }
}

// Sample an ejectile scattering cosine in the CM frame.
double marley::NuclearReaction::sample_cos_theta_c_cm(
  const marley::MatrixElement& matrix_el, double beta_c_cm,
  marley::Generator& gen) const
{
  using ME_Type = marley::MatrixElement::TransitionType;

  // Choose the correct form factor to use based on the matrix element type.
  // Note that our rejection sampling technique does not require that we
  // normalize the form factors before sampling from them.
  /// @todo Implement a more general way of labeling the matrix elements
  std::function<double(double)> form_factor;

  if (matrix_el.type() == ME_Type::FERMI) {
    // B(F)
    form_factor = [&beta_c_cm](double cos_theta_c_cm)
      -> double { return 1. + beta_c_cm * cos_theta_c_cm; };
  }
  else if (matrix_el.type() == ME_Type::GAMOW_TELLER) {
    // B(GT)
    form_factor = [&beta_c_cm](double cos_theta_c_cm)
      -> double { return (3. - beta_c_cm * cos_theta_c_cm) / 3.; };
  }
  else throw marley::Error("Unrecognized matrix element type "
    + std::to_string(matrix_el.type()) + " encountered while sampling a"
    " CM frame scattering angle");

  // Sample a CM frame scattering cosine using the appropriate form factor for
  // this matrix element.
  return gen.rejection_sample(form_factor, -1, 1);
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
