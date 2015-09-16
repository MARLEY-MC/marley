#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "marley_utils.hh"
#include "TMarleyEvent.hh"
#include "TMarleyGenerator.hh"
#include "TMarleyKinematics.hh"
#include "TMarleyLevel.hh"
#include "TMarleyReaction.hh"

TMarleyReaction::TMarleyReaction(std::string filename, TMarleyDecayScheme* scheme) {

  std::regex rx_comment("#.*"); // Matches comment lines

  // Open the reaction data file for parsing
  std::ifstream file_in(filename);

  // If the file doesn't exist or some other error
  // occurred, complain and give up.
  if (!file_in.good()) {
    throw std::runtime_error(std::string("Could not read from the ") +
      "file " + filename);
  }

  std::string line; // String to store the current line
                    // of the reaction data file during parsing

  // TODO: add more error handling for this function

  line = marley_utils::get_next_line(file_in, rx_comment, false);

  // Read in the particle IDs
  std::istringstream iss(line);
  iss >> pid_a >> pid_b >> pid_c >> pid_d;

  // Get initial and final values of the nuclear
  // charge and mass number from the pids
  Zi = (pid_b % 10000000)/10000;
  Ai = (pid_b % 10000)/10;
  Zf = (pid_d % 10000000)/10000;
  Af = (pid_d % 10000)/10;

  // Get the particle masses from the mass table
  ma = TMarleyMassTable::get_particle_mass(pid_a);
  mc = TMarleyMassTable::get_particle_mass(pid_c);

  // If the target (particle b) or residue (particle d)
  // has a particle ID greater than 10^9, assume that it
  // is an atom rather than a bare nucleus
  if (pid_b > 1000000000) {
    mb = TMarleyMassTable::get_atomic_mass(pid_b);
  }
  else {
    mb = TMarleyMassTable::get_particle_mass(pid_b);
  }

  if (pid_d > 1000000000) {
    md_gs = TMarleyMassTable::get_atomic_mass(pid_d);
  }
  else {
    md_gs = TMarleyMassTable::get_particle_mass(pid_d);
  }

  Ea_threshold = ((mc + md_gs)*(mc + md_gs) - ma*ma - mb*mb)/(2*mb);

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
    // TODO: consider implementing a sorting procedure rather than strictly
    // enforcing that energies must be given in ascending order.

    // The order of the entries is important because later uses of the
    // residue_level_energies vector assume that they are sorted in
    // order of ascending energy.
    double energy, strength, strength_id;
    iss >> energy >> strength >> strength_id;
    if (old_energy >= energy) throw std::runtime_error(std::string("Invalid reaction dataset. ")
      + "Level energies must be unique and must be given in ascending order.");
    residue_level_energies.push_back(energy);
    residue_level_strengths.push_back(strength);
    residue_level_strength_ids.push_back(strength_id);
    old_energy = energy;
  }

  // Compute nuclear fragment evaporation thresholds for the residue
  compute_evaporation_thresholds();

  // If a pointer to a decay scheme object is supplied
  // to the constructor, associate it with this reaction.
  // If not, note that we don't have nuclear structure data
  // associated with this reaction yet.
  if (scheme != nullptr) {
    set_decay_scheme(scheme);
  }
  else {
    ds = nullptr;
  }

  file_in.close();

  // Precompute these squared masses for speed
  ma2 = std::pow(ma, 2);
  mb2 = std::pow(mb, 2);
  mc2 = std::pow(mc, 2);
}

// Associate a decay scheme object with this reaction. This will provide
// nuclear structure data for sampling final residue energy levels.
void TMarleyReaction::set_decay_scheme(TMarleyDecayScheme* scheme) {

  if (scheme == nullptr) throw std::runtime_error(std::string("Null pointer")
    + " passed to TMarleyReaction::set_decay_scheme(TMarleyDecayScheme* scheme)");

  // Check to see if the decay scheme being associated with this
  // reaction is for the correct nuclide. If the nuc_id in the decay
  // scheme object does not match the one we'd expect for this reaction's
  // final state nucleus, complain
  std::string reaction_nuc_id = marley_utils::nuc_id(Zf, Af);
  std::string scheme_nuc_id = scheme->get_nuc_id();
  if (reaction_nuc_id != scheme_nuc_id) throw
    std::runtime_error(std::string("Nuclear data mismatch: attempted ")
    + "to associate a decay scheme object that has ENSDF nucid "
    + marley_utils::trim_copy(scheme_nuc_id)
    + " with a reaction object that has nucid "
    + marley_utils::trim_copy(reaction_nuc_id));

  // The decay scheme object has passed our error checks, so go
  // ahead and associate it with this reaction
  ds = scheme;

  // Clear the list of pointers to the ENSDF levels. This will need to
  // be remade each time we associate nuclear structure data with this
  // reaction object
  residue_level_pointers.clear();

  // Use the smallest nuclear fragment emission threshold to check
  // for unbound levels. Note that the constructor sorts the evaporation
  // threshold objects in order of increasing energy, so all we need to do
  // is to look at the first one.
  double unbound_threshold = evaporation_thresholds.front().get_separation_energy();

  // Cycle through each of the level energies given in the reaction dataset. 
  for(std::vector<double>::iterator it = residue_level_energies.begin();
  it != residue_level_energies.end(); ++it)
  {
    // If the level is above the fragment emission threshold, assign it a null
    // level pointer. Such levels will be handled by a fragment evaporation routine
    // and therefore do not need pointers to level objects describing their
    // de-excitation gammas.
    if ((*it) > unbound_threshold) {
      residue_level_pointers.push_back(nullptr);
      continue;
    }

    // For each one, find a pointer to the level with the closest energy owned
    // by the decay scheme object.
    TMarleyLevel* plevel = ds->get_pointer_to_closest_level(*it);
    //std::cout << "DEBUG: I matched E = " << *it << " MeV to the ENSDF level "
    //  << "with energy " << plevel->get_energy() << " MeV" << std::endl;

    // Complain if there are duplicates (if there are duplicates, we'll have
    // two different B(F) + B(GT) values for the same level object)
    if (std::find(residue_level_pointers.begin(), residue_level_pointers.end(),
      plevel) != residue_level_pointers.end()) {
      // residue_level_pointers already contains plevel
      throw std::runtime_error(std::string("Reaction dataset gives two level energies ")
        + "that refer to the same ENSDF level at "
        + std::to_string(plevel->get_energy()) + " MeV");
    }

    // TODO: add check to see if the energy of the chosen level is very different
    // from the energy given in the reaction dataset. If it is, the level matchup
    // is likely incorrect.

    // Add the ENSDF level pointer to the list
    residue_level_pointers.push_back(plevel);
  }

}

void TMarleyReaction::compute_evaporation_thresholds() {

  // Loop over the particle IDs for each of the possible nuclear evaporation fragments
  for(std::vector<int>::const_iterator it = fragment_pids.begin();
    it != fragment_pids.end(); ++it)
  {
    // Compute the separation energy for the current fragment and add it
    // to the vector of evaporation thresholds
    evaporation_thresholds.push_back(TMarleyEvaporationThreshold(
      TMarleyMassTable::get_fragment_separation_energy(Zf, Af, *it), *it));
  }

  // Sort the vector of evaporation thresholds in order of increasing energy
  std::sort(evaporation_thresholds.begin(), evaporation_thresholds.end(),
    [](const TMarleyEvaporationThreshold &et1, const TMarleyEvaporationThreshold &et2)
    -> bool { return et1.get_separation_energy() < et2.get_separation_energy(); });

  //// DEBUG
  //std::cout << std::endl << std::endl;
  //for(std::vector<TMarleyEvaporationThreshold>::iterator i = evaporation_thresholds.begin();
  //  i != evaporation_thresholds.end(); ++i)
  //{
  //  std::cout << "DEBUG: pid = " << i->get_fragment_pid() << " SE = "
  //    << i->get_separation_energy() << std::endl;
  //}
}

// Fermi function used in calculating cross-sections
// The equation was taken from http://en.wikipedia.org/wiki/Beta_decay
// Input: atomic number Z, mass number A, electron total energy E (MeV), positron or electron?
double TMarleyReaction::fermi_function(int Z, int A, double E, bool electron){
  
  double m = marley_utils::m_e; // electron mass (MeV)
  
  double p = marley_utils::real_sqrt(E*E - m*m); // Electron momentum
  double s = std::sqrt(1 - marley_utils::alpha*marley_utils::alpha*Z*Z);

  // Fitting coefficient for estimating the nuclear radius
  // (taken from Introductory Nuclear Physics by Kenneth S. Krane)
  const double r_0 = 1.2; // fm

  // Estimate the nuclear radius using r_n = r_0*A^(1/3)
  double r_n = r_0*std::pow(A, 1./3);

  // Convert the estimated nuclear radius to natural units (MeV^(-1))
  double rho = r_n/marley_utils::hbar_c;

  double eta = marley_utils::alpha*Z*E/p;
  if (!electron) eta *= -1;

  // Complex variables for the gamma function
  std::complex<double> a(s, eta);
  //std::complex<double> b(1+2*s, 0);
  double b = std::tgamma(1+2*s);

  return 2*(1 + s)*std::pow(2*p*rho, 2*s-2)*std::exp(marley_utils::pi*eta)*std::norm(marley_utils::gamma(a))
    / (b*b);
    // /std::pow(std::abs(marley_utils::gamma(b)), 2);
}


// Input: atomic number Z, electron total energy E (MeV), positron or electron?
double TMarleyReaction::fermi_approx(int Z, double E, bool electron){

  // This is a test for comparison with the "exact" calculation, which involves
  // complex gamma functions. This holds true for Q << mc^2
  
  const double m = marley_utils::m_e; // electron mass (MeV)

  double p = std::sqrt(E*E - m*m);
  
  double eta = marley_utils::alpha*Z*E/p;
  if (!electron) eta *= -1;

  return 2*marley_utils::pi*eta/(1-std::exp(-2*marley_utils::pi*eta));
}

// Return the maximum residue excitation energy E_level that can
// be achieved in the lab frame for a given projectile energy Ea
// (this corresponds to the final particles all being produced
// at rest in the CM frame). This maximum level energy is used
// to find the allowed levels when creating events.
double TMarleyReaction::max_level_energy(double Ea) {
  // Calculate the total CM frame energy using known quantities
  // from the lab frame
  double E_CM = std::sqrt(ma*ma + mb*mb + 2*mb*Ea);
  // The maximum level energy is achieved when the final state
  // particles are produced at rest in the CM frame. Subtracting
  // the ground-state rest masses of particles c and d from the
  // total CM energy leaves us with the energy available to create
  // an excited level in the residue (particle d).
  return E_CM - mc - md_gs;
}

double TMarleyReaction::get_threshold_energy() {
  return Ea_threshold;
}

// Creates an event object by sampling the appropriate
// quantities and performing kinematic calculations
TMarleyEvent TMarleyReaction::create_event(double Ea,
  TMarleyGenerator& gen)
{
  // Sample a final residue energy level. First,
  // check to make sure the given projectile energy
  // is above threshold for this reaction.
  if (Ea < Ea_threshold) {
    throw std::range_error("Could not create this event. Projectile energy "
      + std::to_string(Ea) + " MeV is below the threshold value "
      + std::to_string(Ea_threshold) + " MeV.");
  }

  // Calculate the maximum value of E_level that is
  // kinematically allowed
  double max_E_level = max_level_energy(Ea);

  // Create an empty vector of sampling weights
  std::vector<double> level_weights;

  // Create a discrete distribution object for level sampling.
  // Its default constructor creates a single weight of 1.
  // We will always explicitly give it weights to use when sampling
  // levels, so we won't worry about its default behavior.
  static std::discrete_distribution<unsigned int> ldist;

  // Check that there is a decay scheme
  // object associated with this reaction. Complain
  // if there is a problem.
  if (ds == nullptr) {
    throw std::runtime_error(std::string("Could not create this event. No nuclear ")
      + "structure data is associated with this reaction.");
  } 

  // TODO: add more error checks as necessary
  
  // The pointers in residue_level_pointers are ordered by increasing energy
  // (this is currently enforced by the reaction data format and is checked
  // during parsing). Iterate over the levels, assigning the total cross section
  // for each level as its weight until you reach the end of residue_level_pointers
  // or a level that is kinematically forbidden. If the matrix element B(F) + B(GT)
  // for a level vanishes, skip computing the total cross section and assign a weight
  // of zero for efficiency.
  bool at_least_one_nonzero_matrix_el = false;

  for(unsigned int i = 0; i < residue_level_pointers.size(); ++i) {

    double level_energy;

    // Get the excitation energy for the current level
    if (residue_level_pointers[i] != nullptr) {
      level_energy = residue_level_pointers[i]->get_energy();
    }
    else {
      // Level is unbound, so just use the energy given in the reaction dataset
      // rather than trying to find its ENSDF version
      level_energy = residue_level_energies[i];
    }

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // calling total_xs. This will avoid unnecessary numerical integrations.
    double matrix_el = residue_level_strengths[i];
    double xs;

    if (matrix_el == 0) {
      xs = 0;
    }
    else {
      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that std::discrete_distribution
      // automatically normalizes the weights, so we don't have to do that ourselves.
      xs = total_xs_cm(level_energy, Ea, matrix_el);
      //DEBUG
      if (std::isnan(xs)) {
        std::cout << "DEBUG: this level gave a weight of nan, so I made it zero." << std::endl;
        xs = 0;
      }
      if (!at_least_one_nonzero_matrix_el && xs != 0)
        at_least_one_nonzero_matrix_el = true;
    }
    level_weights.push_back(xs);
  }

  // Complain if none of the levels we have data for are kinematically allowed
  if (level_weights.empty()) {
    throw std::runtime_error(std::string("Could not create this event. The ")
      + "TMarleyDecayScheme object associated with this reaction "
      + "does not contain data for any kinematically accessible levels "
      + "for a projectile energy of " + std::to_string(Ea)
      + " MeV (max E_level = " + std::to_string(max_E_level)
      + " MeV).");
  }

  // Complain if none of the levels have nonzero weight (this means that
  // all kinematically allowed levels have a vanishing matrix element)
  if (!at_least_one_nonzero_matrix_el) {
    throw std::runtime_error(std::string("Could not create this event. All ")
      + "kinematically accessible levels for a projectile energy of "
      + std::to_string(Ea) + " MeV (max E_level = " + std::to_string(max_E_level)
      + " MeV) have vanishing matrix elements.");
  }

  // Create a list of parameters used to supply the weights to our discrete
  // level sampling distribution
  std::discrete_distribution<unsigned int>::param_type params(level_weights.begin(),
    level_weights.end());

  // Sample a level index using our discrete distribution and the
  // current set of weights
  unsigned int l_index = gen.discrete_sample(ldist, params);

  // Get a pointer to the selected level
  TMarleyLevel* plevel = residue_level_pointers.at(l_index);

  // Get the energy of the selected level. This will be
  // needed for sampling a scattering cosine.
  double E_level;
  if (plevel != nullptr) {
    E_level = plevel->get_energy();
  }
  else {
    // Level is unbound, so use the level energy found in the reaction dataset
    E_level = residue_level_energies.at(l_index);
  }

  double md = md_gs + E_level;
  double md2 = std::pow(md, 2);

  // Compute Mandelstam s (the square of the total CM frame energy)
  double s = ma2 + mb2 + 2 * mb * Ea;
  double sqrt_s = std::sqrt(s);

  // Determine the CM frame energy, momentum, and velocity of the ejectile
  double Ec_cm = (s + mc2 - md2) / (2 * sqrt_s);
  double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc2);
  double beta_c_cm = pc_cm / Ec_cm;

  // Sample a CM frame scattering cosine for the ejectile.
  //double matrix_el = residue_level_strengths.at(l_index);
  double m_type = residue_level_strength_ids.at(l_index);
  double cos_theta_c_cm = sample_cos_theta_c_cm(/*matrix_el,*/ m_type, beta_c_cm,
    gen);
  double sin_theta_c_cm = marley_utils::real_sqrt(1
    - std::pow(cos_theta_c_cm, 2));

  // Sample a CM frame azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the matrix elements are azimuthally invariant
  double phi_c_cm = gen.uniform_random_double(0, 2*marley_utils::pi, false);

  // Determine the Cartesian components of the ejectile's CM frame momentum
  double pc_cm_x = sin_theta_c_cm * std::cos(phi_c_cm) * pc_cm;
  double pc_cm_y = sin_theta_c_cm * std::sin(phi_c_cm) * pc_cm;
  double pc_cm_z = cos_theta_c_cm * pc_cm;

  // Determine the residue's CM frame energy. Roundoff errors may cause Ed_cm to
  // dip below md, which is unphysical. Prevent this from occurring by allowing
  // md to be the minimum value of Ed_cm. Also note that, in the CM frame, the
  // residue and ejectile have equal and opposite momenta.
  double Ed_cm = std::max(sqrt_s - Ec_cm, md);

  // Create particle objects representing the projectile and target in the lab frame
  // TODO: edit this to allow for projectile directions other than along the z-axis
  TMarleyParticle projectile(pid_a, Ea, 0, 0,
    marley_utils::real_sqrt(std::pow(Ea, 2) - ma2), ma);

  TMarleyParticle target(pid_b, mb, 0, 0, 0, mb);

  // Create particle objects representing the ejectile and residue in the CM
  // frame.
  TMarleyParticle ejectile(pid_c, Ec_cm, pc_cm_x, pc_cm_y, pc_cm_z, mc);
  TMarleyParticle residue(pid_d, Ed_cm, -pc_cm_x, -pc_cm_y, -pc_cm_z, md);

  // Boost the ejectile and residue into the lab frame.
  // TODO: edit this to allow for projectile directions other than along the z-axis
  double beta_z = Ea / (Ea + mb);
  TMarleyKinematics::lorentz_boost(0, 0, -beta_z, ejectile);
  TMarleyKinematics::lorentz_boost(0, 0, -beta_z, residue);

  // Create the event object and load it with the appropriate information
  TMarleyEvent event(E_level);
  event.set_reaction(this);
  // Add the projectile to this event's initial particle list
  event.add_initial_particle(projectile, TMarleyEvent::ParticleRole::pr_projectile);
  // Add the target to this event's initial particle list
  event.add_initial_particle(target, TMarleyEvent::ParticleRole::pr_target);
  // Add the ejectile to this event's final particle list
  event.add_final_particle(ejectile, TMarleyEvent::ParticleRole::pr_ejectile);
  // Add the residue to this event's final particle list.
  event.add_final_particle(residue, TMarleyEvent::ParticleRole::pr_residue);

  if (plevel != nullptr) {
    // The selected level is bound, so it will only decay via gamma emission.
    // Add the de-excitation gammas to this event's final particle list
    this->ds->do_cascade(plevel, &event, gen);
  }
  else {
    // The selected level is unbound, so handle nucleon emission, etc.
    this->evaporate_particles(E_level, &event, gen);
  }

  // Return the completed event object
  return event;
}

void TMarleyReaction::evaporate_particles(double E_level, TMarleyEvent* p_event,
  TMarleyGenerator& gen)
{
  // Find the highest-energy particle evaporation threshold that is
  // greater than the residue excitation energy E_level
  std::vector<TMarleyEvaporationThreshold>::iterator i_thresh =
    std::upper_bound(evaporation_thresholds.begin(), evaporation_thresholds.end(),
      E_level, [](const double E_level, const TMarleyEvaporationThreshold &et)
      -> bool { return E_level < et.get_separation_energy(); });

  if (i_thresh == evaporation_thresholds.begin()) throw
    std::runtime_error(std::string("Particle evaporation attempted for ")
      + "a residue excitation energy of " + std::to_string(E_level)
      + " MeV, which is less than the particle evaporation threshold of "
      + std::to_string(evaporation_thresholds.front().get_separation_energy())
      + " MeV");

  // Assume that the particle whose evaporation threshold energy is closest to
  // the residue excitation energy (without exceeding it) will be emitted, leaving
  // the residue in its ground state. Get the particle ID and separation energy
  // for the emitted particle.
  --i_thresh;
  int fragment_pid = i_thresh->get_fragment_pid();

  // Get the particle ID for the residue in its final state
  int Zres = Zf - TMarleyMassTable::get_particle_Z(fragment_pid);
  int Ares = Af - TMarleyMassTable::get_particle_A(fragment_pid);
  int pid_res = marley_utils::get_nucleus_pid(Zres, Ares);
  double m_res = TMarleyMassTable::get_atomic_mass(pid_res);

  // Get the mass of the fragment
  double m_fragment = TMarleyMassTable::get_particle_mass(fragment_pid); 

  // Create particle objects to store the recoiling residue and the fragment
  TMarleyParticle new_res(pid_res, m_res);
  TMarleyParticle frag(fragment_pid, m_fragment);

  // TODO: consider changing fragment emission so that it is anisotropic
  // Sample an azimuthal emission angle (phi) uniformly on [0, 2*pi).
  double phi_fragment = gen.uniform_random_double(0, 2*marley_utils::pi, false);
  // Sample a polar emission cosine (cos[theta]) uniformly on [-1, 1].
  double cos_theta_fragment = gen.uniform_random_double(-1, 1, true);

  // Load the particle objects with their correct final energies and momenta
  // based on the two-body decay of the original nucleus
  TMarleyKinematics::two_body_decay(*(p_event->get_residue()), frag, new_res,
    cos_theta_fragment, phi_fragment);

  // Update the residue to adjust for effects of fragment emission
  TMarleyParticle* p_residue = p_event->get_residue();
  *p_residue = new_res;

  // Add the evaporated fragment to this event's list of final particles
  p_event->add_final_particle(frag);

  // Add the evaporated fragment to the residue's children
  p_residue->add_child(&(p_event->get_final_particles()->back())); 
}

// Compute the total reaction cross section (including all final nuclear levels)
// in units of MeV^(-2) using the center of mass frame.
double TMarleyReaction::total_xs_cm(double Ea) {
  double max_E_level = max_level_energy(Ea);
  double xs = 0;
  for(unsigned int i = 0; i < residue_level_pointers.size(); ++i) {

    double level_energy;

    // Get the excitation energy for the current level
    if (residue_level_pointers[i] != nullptr) {
      level_energy = residue_level_pointers[i]->get_energy();
    }
    else {
      // Level is unbound, so just use the energy given in the reaction dataset
      // rather than trying to find its ENSDF version
      level_energy = residue_level_energies[i];
    }

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // Check whether the matrix element (B(F) + B(GT)) is nonvanishing for the
    // current level. If it is, just set the weight equal to zero rather than
    // computing the total xs.
    double matrix_el = residue_level_strengths[i];

    if (matrix_el == 0) {
      xs += 0;
    }
    else {
      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that std::discrete_distribution
      // automatically normalizes the weights, so we don't have to do that ourselves.
      double md2 = std::pow(md_gs + level_energy, 2);

      // Compute Mandelstam s (the square of the total CM frame energy)
      double s = ma2 + mb2 + 2 * mb * Ea;
      double sqrt_s = std::sqrt(s);

      // Compute CM frame energies for two of the particles. Also
      // compute the ejectile CM frame momentum.
      double Eb_cm = (s + mb2 - ma2) / (2 * sqrt_s);
      double Ec_cm = (s + mc2 - md2) / (2 * sqrt_s);
      double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc2);

      xs += fermi_function(Zf, Af, Ec_cm, true) * matrix_el * pc_cm * Ec_cm
        * Eb_cm * (sqrt_s - Ec_cm) / s;
    }
  }
  return std::pow(GF, 2) * std::pow(Vud, 2) * xs / marley_utils::pi;
}

// Compute the total reaction cross section for a given level (ignoring some
// constants for speed, since this version of the function is intended for use
// in sampling only) using the center of mass frame
double TMarleyReaction::total_xs_cm(double E_level, double Ea,
  double matrix_element)
{
  // Don't bother to compute anything if the matrix element vanishes for this level
  if (matrix_element == 0) return 0;
  else {
    double md2 = std::pow(md_gs + E_level, 2);

    // Compute Mandelstam s (the square of the total CM frame energy)
    double s = ma2 + mb2 + 2 * mb * Ea;
    double sqrt_s = std::sqrt(s);

    // Compute CM frame energies for two of the particles. Also
    // compute the ejectile CM frame momentum.
    double Eb_cm = (s + mb2 - ma2) / (2 * sqrt_s);
    double Ec_cm = (s + mc2 - md2) / (2 * sqrt_s);
    double pc_cm = marley_utils::real_sqrt(std::pow(Ec_cm, 2) - mc2);

    return fermi_function(Zf, Af, Ec_cm, true) * matrix_element * pc_cm * Ec_cm
      * Eb_cm * (sqrt_s - Ec_cm) / s;
  }
}

// Sample an ejectile scattering cosine in the CM frame.
double TMarleyReaction::sample_cos_theta_c_cm(/*double matrix_el,*/ int m_type,
  double beta_c_cm, TMarleyGenerator& gen)
{
  // Choose the correct form factor to use based on the matrix element type.
  // Note that our rejection sampling technique does not require that we
  // normalize the form factors before sampling from them.
  // TODO: implement a more general way of labeling these matrix elements
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
  else throw std::runtime_error(std::string("Unrecognized matrix element")
    + " type " + std::to_string(m_type) + " encountered while sampling a"
    + " CM frame scattering angle");

  // Sample a CM frame scattering cosine using the appropriate form factor for
  // this matrix element.
  return gen.rejection_sample(form_factor, -1, 1);
}
