#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>

#include "marley_utils.hh"
#include "Error.hh"
#include "Event.hh"
#include "Generator.hh"
#include "Level.hh"
#include "Logger.hh"
#include "HauserFeshbachDecay.hh"
#include "Reaction.hh"

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

  // Get the particle masses from the mass table
  ma_ = marley::MassTable::get_particle_mass(pdg_a_);
  mc_ = marley::MassTable::get_particle_mass(pdg_c_);

  // If the target (particle b) or residue (particle d)
  // has a particle ID greater than 10^9, assume that it
  // is an atom rather than a bare nucleus
  if (pdg_b_ > 1000000000) {
    mb_ = marley::MassTable::get_atomic_mass(pdg_b_);
  }
  else {
    mb_ = marley::MassTable::get_particle_mass(pdg_b_);
  }

  if (pdg_d_ > 1000000000) {
    // If particle d is an atom and is ionized as a result of this reaction
    // (e.g., q_d_ != 0), then approximate its ground-state ionized mass by
    // subtracting the appropriate number of electron masses from its atomic
    // (i.e., neutral) ground state mass.
    md_gs_ = marley::MassTable::get_atomic_mass(pdg_d_)
      - (q_d_ * marley::MassTable::get_particle_mass(marley_utils::ELECTRON));
  }
  else {
    md_gs_ = marley::MassTable::get_particle_mass(pdg_d_);
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

    // The order of the entries is important because later uses of the
    // residue_level_energies vector assume that they are sorted in
    // order of ascending energy.
    double energy, strength, strength_id;
    iss >> energy >> strength >> strength_id;
    if (old_energy >= energy) throw marley::Error(std::string("Invalid reaction dataset. ")
      + "Level energies must be unique and must be given in ascending order.");
    residue_level_energies_.push_back(energy);
    residue_level_strengths_.push_back(strength);
    residue_level_strength_ids_.push_back(strength_id);
    // Initialize all of the level pointers to nullptr. This may be changed
    // later if discrete level data can be found for the residual nucleus.
    residue_level_pointers_.push_back(nullptr);
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
    + " passed to marley::NuclearReaction::set_decay_scheme(marley::DecayScheme* ds)");

  // Check to see if the decay scheme being associated with this
  // reaction is for the correct nuclide. If the nuc_id in the decay
  // scheme object does not match the one we'd expect for this reaction's
  // final state nucleus, complain
  std::string reaction_nuc_id = marley_utils::nuc_id(Zf_, Af_);
  std::string scheme_nuc_id = ds->get_nuc_id();
  if (reaction_nuc_id != scheme_nuc_id) throw
    marley::Error(std::string("Nuclear data mismatch: attempted ")
    + "to associate a decay scheme object that has ENSDF nucid "
    + marley_utils::trim_copy(scheme_nuc_id)
    + " with a reaction object that has nucid "
    + marley_utils::trim_copy(reaction_nuc_id));

  // Use the smallest nuclear fragment emission threshold to check for unbound
  // levels. Before looking them up, start by setting the unbound threshold to
  // infinity.
  double unbound_threshold = std::numeric_limits<double>::max();
  LOG_DEBUG() << "unbound_threshold = " << unbound_threshold << std::endl;
  for (const auto& f : marley::HauserFeshbachDecay::get_fragments()) {
    double thresh = marley::HauserFeshbachDecay
      ::get_fragment_emission_threshold(Zf_, Af_, f);
    LOG_DEBUG() << f.get_pid() << " emission threshold = " << thresh << std::endl;
    if (thresh < unbound_threshold) unbound_threshold = thresh;
    LOG_DEBUG() << "unbound_threshold = " << unbound_threshold << std::endl;
  }

  // Cycle through each of the level energies given in the reaction dataset.
  for(size_t j = 0, s = residue_level_energies_.size(); j < s; ++j)
  {
    double en = residue_level_energies_.at(j);
    // If the level is above the fragment emission threshold, assign it a null
    // level pointer. Such levels will be handled by a fragment evaporation routine
    // and therefore do not need pointers to level objects describing their
    // de-excitation gammas.
    if (en > unbound_threshold) {
      residue_level_pointers_.at(j) = nullptr;
      continue;
    }

    // For each energy, find a pointer to the level with the closest energy
    // owned by the decay scheme object.
    marley::Level* plevel = ds->get_pointer_to_closest_level(en);
    LOG_DEBUG() << "reaction level at " << en
      << " MeV was matched to the decay scheme level at "
      << plevel->get_energy() << " MeV";

    // Complain if there are duplicates (if there are duplicates, we'll have
    // two different B(F) + B(GT) values for the same level object)
    if (std::find(residue_level_pointers_.begin(),
      residue_level_pointers_.end(), plevel) != residue_level_pointers_.end())
    {
      // residue_level_pointers already contains plevel
      throw marley::Error(std::string("Reaction dataset gives two")
        + " level energies that refer to the same DecayScheme level at "
        + std::to_string(plevel->get_energy()) + " MeV");
    }

    /// @todo Add check to see if the energy of the chosen level is very
    /// different from the energy given in the reaction dataset. If it is, the
    /// level matchup is likely incorrect.

    // Add the level pointer to the list
    residue_level_pointers_.at(j) = plevel;
  }
}

// Fermi function used to calculate cross-sections
// The form used here is based on http://en.wikipedia.org/wiki/Beta_decay
// but rewritten for convenient use inside this class.
// Input: beta_c (3-momentum magnitude of particle c / total energy of particle c),
// where we assume that the ejectile (particle c) is the light product from
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

// Creates an event object by sampling the appropriate
// quantities and performing kinematic calculations
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

  // The pointers in residue_level_pointers are ordered by increasing energy
  // (this is currently enforced by the reaction data format and is checked
  // during parsing). Iterate over the levels, assigning the total cross
  // section for each level as its weight until you reach the end of
  // residue_level_pointers or a level that is kinematically forbidden. If the
  // matrix element B(F) + B(GT) for a level vanishes, skip computing the total
  // cross section and assign a weight of zero for efficiency.
  bool at_least_one_nonzero_matrix_el = false;

  for(size_t i = 0; i < residue_level_pointers_.size(); ++i) {

    double level_energy;

    // Get the excitation energy for the current level
    if (residue_level_pointers_[i] != nullptr) {
      level_energy = residue_level_pointers_[i]->get_energy();
    }
    else {
      // Level is unbound, so just use the energy given in the reaction dataset
      // rather than trying to find its DecayScheme version
      level_energy = residue_level_energies_[i];
    }

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // calling total_xs. This will avoid unnecessary numerical integrations.
    double matrix_el = residue_level_strengths_[i];
    double xs;

    if (matrix_el == 0.) {
      xs = 0.;
    }
    else {

      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that
      // std::discrete_distribution automatically normalizes the weights, so we
      // don't have to do that ourselves.
      xs = total_xs(level_energy, KEa, matrix_el);
      if (std::isnan(xs)) {
        LOG_WARNING() << "Partial cross section for reaction " << description_
          << " gave NaN result.";
        LOG_DEBUG() << "Parameters were level energy = " << level_energy
          << " MeV, projectile kinetic energy = " << KEa
          << " MeV, and matrix element = " << matrix_el;
        LOG_DEBUG() << "The partial cross section to this level"
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

  // Sample a level index using our discrete distribution and the
  // current set of weights
  size_t l_index = gen.discrete_sample(ldist, params);

  // Get a pointer to the selected level
  marley::Level* plevel = residue_level_pointers_.at(l_index);

  // Get the energy of the selected level.
  double E_level;
  if (plevel != nullptr) {
    E_level = plevel->get_energy();
  }
  else {
    // Level is unbound, so use the level energy found in the reaction dataset
    E_level = residue_level_energies_.at(l_index);
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
  //double matrix_el = residue_level_strengths.at(l_index);
  int m_type = residue_level_strength_ids_.at(l_index);
  double cos_theta_c_cm = sample_cos_theta_c_cm(/*matrix_el,*/ m_type, beta_c_cm,
    gen);

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
    int twoJ = 2*m_type;
    // Fermi transition gives 0+ -> 0+, GT transition gives 0+ -> 1+
    /// @todo include possibility of negative parity here.
    marley::Parity P(true); // positive parity
    marley::Particle first, second;
    // The selected level is unbound, so handle its de-excitation using
    // the Hauser-Feshbach statistical model.
    while (continuum && Ex > cutoff) {
      marley::HauserFeshbachDecay hfd(residue, Ex, twoJ, P, gen);
      LOG_DEBUG() << hfd;
      continuum = hfd.do_decay(Ex, twoJ, P, first, second);

      LOG_DEBUG() << "Hauser-Feshbach decay to " << first.pdg_code()
        << " and " << second.pdg_code();
      LOG_DEBUG() << second.pdg_code() << " is at Ex = " << Ex << " MeV.";

      residue = second;
      Z = marley::MassTable::get_particle_Z(residue.pdg_code());
      A = marley::MassTable::get_particle_A(residue.pdg_code());
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
    dec_scheme->do_cascade(dec_scheme->get_pointer_to_closest_level(Ex),
      &event, gen, residue.charge());
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
  for(size_t i = 0; i < residue_level_pointers_.size(); ++i) {

    double level_energy;

    // Get the excitation energy for the current level
    if (residue_level_pointers_[i] != nullptr) {
      level_energy = residue_level_pointers_[i]->get_energy();
    }
    else {
      // Level is unbound, so just use the energy given in the reaction dataset
      // rather than trying to find its DecayScheme version
      level_energy = residue_level_energies_[i];
    }

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // computing the total xs.
    double matrix_el = residue_level_strengths_[i];

    if (matrix_el != 0.) {
      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that std::discrete_distribution
      // automatically normalizes the weights, so we don't have to do that ourselves.
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
double marley::NuclearReaction::sample_cos_theta_c_cm(/*double matrix_el,*/
  int m_type, double beta_c_cm, marley::Generator& gen) const
{
  // Choose the correct form factor to use based on the matrix element type.
  // Note that our rejection sampling technique does not require that we
  // normalize the form factors before sampling from them.
  /// @todo Implement a more general way of labeling the matrix elements
  std::function<double(double)> form_factor;

  if (m_type == 0) {
    // B(F)
    form_factor = [&beta_c_cm](double cos_theta_c_cm)
      -> double { return 1. + beta_c_cm * cos_theta_c_cm; };
  }
  else if (m_type == 1) {
    // B(GT)
    form_factor = [&beta_c_cm](double cos_theta_c_cm)
      -> double { return (3. - beta_c_cm * cos_theta_c_cm) / 3.; };
  }
  else throw marley::Error(std::string("Unrecognized matrix element")
    + " type " + std::to_string(m_type) + " encountered while sampling a"
    + " CM frame scattering angle");

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
