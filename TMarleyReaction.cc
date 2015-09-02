#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "marley_utils.hh"
#include "TMarleyEvent.hh"
#include "TMarleyGenerator.hh"
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

  // Read in all of the level energy (MeV) and squared matrix element (B(F) +
  // B(GT) strength) pairs.

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
    double energy, strength;
    iss >> energy >> strength;
    if (old_energy >= energy) throw std::runtime_error(std::string("Invalid reaction dataset. ")
      + "Level energies must be unique and must be given in ascending order.");
    residue_level_energies.push_back(energy);
    residue_level_strengths.push_back(strength);
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

double TMarleyReaction::ejectile_energy(double E_level, double Ea, double cos_theta_c) {

  // All energies are in MeV
  double md = md_gs + E_level;

  double Etot = Ea + mb;

  // Squared 3-momentum of the projectile
  double p2a = Ea*Ea - ma*ma;

  // Quadratic formula coefficients from 4-momentum
  // conservation solution in the lab frame
  double A = 4*Etot*Etot - 4*p2a*cos_theta_c*cos_theta_c;

  double B = 4*Etot*(md*md + p2a - mc*mc - Etot*Etot);

  double C = std::pow((Etot*Etot + mc*mc - md*md - p2a), 2)
    + 4*mc*mc*p2a*cos_theta_c*cos_theta_c;

  // There are two solutions of the quadratic equation 
  // that correspond to forward and backward scattering
  // solutions. Find both of them and load Ecplus and
  // Ecminus with the results.
  double Ecplus, Ecminus;
  marley_utils::solve_quadratic_equation(A, B, C, Ecplus, Ecminus);

  // Assume we are above the backward threshold, and match the
  // correct solution to the physical ejectile energy.
  // TODO: generalize this to include cases where we're below
  // the backward threshold
  if (cos_theta_c >= 0) {
    return Ecplus; // Forward scattering
  }
  else {
    return Ecminus; // Backward scattering
  }

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
  // All hard-coded energies are in MeV
  double Etot = Ea + mb;

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
  // (this currently enforced by the reaction data format and is checked
  // during parsing). Iterate over the levels, assigning the total cross section
  // for each level as its weight until you reach the end of residue_level_pointers
  // or a level that is kinematically forbidden. If the matrix element B(F) + B(GT)
  // for a level vanishes, skip computing the total cross section and assign a weight
  // of zero for efficiency.
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
      xs = total_xs(level_energy, Ea, matrix_el);
      //DEBUG
      if (std::isnan(xs)) {
        std::cout << "DEBUG: this level gave a weight of nan, so I made it zero." << std::endl;
        xs = 0;
      }
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

  // Create a list of parameters used to supply the weights to our discrete
  // level sampling distribution
  std::discrete_distribution<unsigned int>::param_type params(level_weights.begin(),
    level_weights.end());

  //DEBUG
  //std::vector<double> probs = params.probabilities();
  //for (auto& n : probs)
  //{
  //  int i = &n - &(probs[0]);
  //  std::cout << "DEBUG: level at " << residue_level_pointers.at(i)->get_energy()
  //    << " MeV has probability " << n << " of being selected" << std::endl;
  //}

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

  // Get the matrix element (B(F) + B(GT) value) for the
  // selected level
  double matrix_element = residue_level_strengths.at(l_index);

  // Sample a scattering angle for the ejectile using
  // the differential cross section
  double cos_theta_c = sample_ejectile_scattering_cosine(E_level, Ea,
    matrix_element, gen);
  double theta_c = std::acos(cos_theta_c);

  // Use conservation of 4-momentum to compute the ejectile energy
  // based on the sampled scattering angle
  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);

  // Compute the energy and scattering angle of the residue
  double Ed = Etot - Ec;
  double pc = marley_utils::real_sqrt(Ec*Ec - mc*mc);
  double pd_before_deexcitation = marley_utils::real_sqrt(Ed*Ed - md*md);
  double theta_d = std::asin(pc * std::sin(theta_c) / pd_before_deexcitation);
  double cos_theta_d = std::cos(theta_d);

  // Sample an azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the matrix elements are azimuthally invariant
  double phi = gen.uniform_random_double(0, 2*marley_utils::pi, false);

  // TODO: maybe consider the effect of nuclear recoils that happen when
  // gammas/nucleons are emitted
  // Compute the momentum of the residue when it reaches the ground state
  double Ed_gs = Ed - E_level;
  double pd_gs = marley_utils::real_sqrt(Ed_gs*Ed_gs - md_gs*md_gs);

  // Use the scattering angles to compute Cartesian 3-momentum components
  // for the ejectile and residue
  double pc_x = std::sin(theta_c)*std::cos(phi)*pc;
  double pc_y = std::sin(theta_c)*std::sin(phi)*pc;
  double pc_z = cos_theta_c*pc;
  // The residue scattering angle (theta_d) is measured
  // in the clockwise direction, so the x and y components
  // of the residue 3-momentum pick up minus signs (sin[-x] = -sin[x])
  double pd_x_gs = -std::sin(theta_d)*std::cos(phi)*pd_gs;
  double pd_y_gs = -std::sin(theta_d)*std::sin(phi)*pd_gs;
  double pd_z_gs = cos_theta_d*pd_gs;

  // Print results to std::cout
  //std::cout.precision(15);
  //std::cout << std::scientific;
  //std::cout << "E_level = " << E_level << std::endl;
  //std::cout << "cos_theta_c = " << cos_theta_c << std::endl;
  //std::cout << "Ec = " << Ec << std::endl;
  //std::cout << "Ed = " << Ed << std::endl;
  //std::cout << "cos_theta_d = " << cos_theta_d << std::endl;
  //std::cout << "e- kinetic energy = " << Ec - mc << std::endl;
  //std::cout << "e- mass = " << mc << std::endl;
  //std::cout << "40K kinetic energy = " << Ed - md_gs - E_level << std::endl;
  //std::cout << "Ground state nuclear mass change = " << md_gs - mb << std::endl; 

  // Create the event object and load it with the appropriate information
  TMarleyEvent event(E_level);
  event.set_reaction(this);
  // TODO: edit this to allow for projectile directions other than along the z-axis
  // Add the projectile to this event's initial particle list
  event.add_initial_particle(TMarleyParticle(pid_a, Ea, 0, 0, Ea - ma, ma),
    TMarleyEvent::ParticleRole::pr_projectile);
  // Add the target to this event's initial particle list
  event.add_initial_particle(TMarleyParticle(pid_b, mb, 0, 0, 0, mb),
    TMarleyEvent::ParticleRole::pr_target);
  // Add the ejectile to this event's final particle list
  event.add_final_particle(TMarleyParticle(pid_c, Ec, pc_x, pc_y, pc_z, mc),
    TMarleyEvent::ParticleRole::pr_ejectile);

  if (plevel != nullptr) {
    // The selected level is bound, so it will only decay via gamma emission.
    // Add the residue to this event's final particle list. Don't include
    // its excitation energy since we will soon create the de-excitation gamma rays
    event.add_final_particle(TMarleyParticle(pid_d, Ed_gs, pd_x_gs, pd_y_gs, pd_z_gs, md_gs),
      TMarleyEvent::ParticleRole::pr_residue);

    // Add the de-excitation gammas to this event's final particle list
    this->ds->do_cascade(plevel, &event, gen);
  }
  else {
    // The selected level is unbound, so handle nucleon emission, etc.
    this->evaporate_particles(E_level, &event, Ed_gs, theta_d, phi, gen);
  }

  // Return the completed event object
  return event;
}

void TMarleyReaction::evaporate_particles(double E_level, TMarleyEvent* p_event,
  double Ed_gs, double theta_res, double phi_res, TMarleyGenerator& gen)
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
  double sep_energy = i_thresh->get_separation_energy();

  // Get the particle ID for the residue in its final state
  int Zres = Zf - TMarleyMassTable::get_particle_Z(fragment_pid);
  int Ares = Af - TMarleyMassTable::get_particle_A(fragment_pid);
  int pid_res = marley_utils::get_nucleus_pid(Zres, Ares);
  double m_res = TMarleyMassTable::get_atomic_mass(pid_res);

  // TODO: correct the kinematics for this function to account for nuclear recoil, etc.
  double m_fragment = TMarleyMassTable::get_particle_mass(fragment_pid); 
  double E_fragment = m_fragment + E_level - sep_energy;
  double p_fragment = marley_utils::real_sqrt(E_fragment*E_fragment - m_fragment*m_fragment);

  double E_res = Ed_gs - m_fragment + sep_energy;
  double p_res = marley_utils::real_sqrt(E_res*E_res - m_res*m_res);
  // The residue scattering angle (theta_d) is measured
  // in the clockwise direction, so the x and y components
  // of the residue 3-momentum pick up minus signs (sin[-x] = -sin[x])
  double p_res_x = -std::sin(theta_res)*std::cos(phi_res)*p_res;
  double p_res_y = -std::sin(theta_res)*std::sin(phi_res)*p_res;
  double p_res_z = std::cos(theta_res)*p_res;

  // Add the final-state residue to this event's list of final particles
  p_event->add_final_particle(TMarleyParticle(pid_res, E_res, p_res_x, p_res_y, p_res_z, m_res),
    TMarleyEvent::ParticleRole::pr_residue);

  // Sample an azimuthal emission angle (phi) uniformly on [0, 2*pi).
  double phi_fragment = gen.uniform_random_double(0, 2*marley_utils::pi, false);
  // Sample a polar emission cosine (cos[theta]) uniformly on [-1, 1].
  double cos_theta_fragment = gen.uniform_random_double(-1, 1, true);
  double theta_fragment = std::acos(cos_theta_fragment);
  // Get the Cartesian components of the evaporated fragment's 3-momentum
  double p_fragment_x = std::sin(theta_fragment)*std::cos(phi_fragment)*p_fragment;
  double p_fragment_y = std::sin(theta_fragment)*std::sin(phi_fragment)*p_fragment;
  double p_fragment_z = cos_theta_fragment*p_fragment;

  // Add the evaporated fragment to this event's list of final particles
  p_event->add_final_particle(TMarleyParticle(fragment_pid, E_fragment, p_fragment_x,
    p_fragment_y, p_fragment_z, m_fragment));

  // Add the evaporated fragment to the residue's children
  TMarleyParticle* p_residue = p_event->get_residue();
  p_residue->add_child(&(p_event->get_final_particles()->back())); 
}

// Compute the differential reaction cross section dsigma/dcos_theta
// for a given final residue level, projectile energy, and ejectile
// scattering cosine
double TMarleyReaction::differential_xs(double E_level, double Ea,
  double matrix_element, double cos_theta_c)
{

  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);
  double pc = marley_utils::real_sqrt(Ec*Ec - mc*mc);

  // TODO: adjust this differential cross section expression
  // as needed
  
  // This expression for the differential cross section
  // dsigma/dcos_theta comes from Kuramoto, et al.,
  // Nucl. Phys. A512 (1990) 711-736. It is modified
  // using a nuclear recoil correction factor taken from
  // J. D. Walecka, "Semileptonic Weak Interactions in Nuclei,"
  // In: Muon Physics, Volume II: Weak Interactions.
  // Ed. by V. W. Hughes and C. S. Wu.
 
  // The constants GF and Vud have been removed since,
  // for sampling, they will not be needed (what matters is the
  // relative size of the differential cross section).
  // This might reduce numerical error.
  return (1.0/(2*std::acos(-1)))// * std::pow(GF, 2) * std::pow(Vud, 2)
    *pc*Ec*fermi_function(Zf, Af, Ec, true)*matrix_element
    /(1.0 + (Ea/mb)*(1 - (Ec/pc)*cos_theta_c));
}

// Compute the total reaction cross section (including all final nuclear levels)
// in units of MeV^(-2)
double TMarleyReaction::total_xs(double Ea) {
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
    // calling total_xs. This will avoid unnecessary numerical integrations.
    double matrix_el = residue_level_strengths[i];

    if (matrix_el == 0) {
      xs += 0;
    }
    else {
      // If the matrix element is nonzero, assign a weight to this level equal
      // to the total reaction cross section. Note that std::discrete_distribution
      // automatically normalizes the weights, so we don't have to do that ourselves.
      xs += total_xs(level_energy, Ea, matrix_el);
    }
  }
  return std::pow(GF, 2) * std::pow(Vud, 2) * xs;
}

// Numerically integrate the differential cross section over
// cos(theta_c) = [-1,1] to get the total reaction cross section
// for a given final residue level and projectile energy
double TMarleyReaction::total_xs(double E_level, double Ea, double matrix_element) {
  // TODO: consider switching to a different integration method,
  // using tabulated data, etc.

  // Numerically integrate the differential cross section over
  // the interval cos(theta_c) = [-1,1] using the composite
  // trapezoidal rule over n subintervals
  static const int n = 1e3;

  // Create a forwarding call wrapper for the  differential_xs
  // member function that takes a single argument.
  std::function<double(double)> dxs = std::bind(
    &TMarleyReaction::differential_xs, this, E_level, Ea, matrix_element, std::placeholders::_1);

  // Numerically integrate using the call wrapper, the integration bounds,
  // and the number of subintervals
  return marley_utils::num_integrate(dxs, -1.0, 1.0, n); 
}

// Sample a scattering cosine for the ejectile using the differential
// cross section defined in the member function differential_xs
double TMarleyReaction::sample_ejectile_scattering_cosine(double E_level,
  double Ea, double matrix_element, TMarleyGenerator& gen)
{

  // Get the total cross section for this reaction based
  // on the projectile energy and final residue energy.
  // This will be used as a normalization factor.
  double xstot = total_xs(E_level, Ea, matrix_element);

  // Make a normalized version of the differential cross section
  // to use as our probability density function for sampling.
  // Note that, since we are using a rejection sampling method,
  // normalization is not strictly necessary. That being said, doing so
  // allows us to set our value of epsilon for marley_utils::maximize
  // without having to worry about the absolute scale of the differential cross
  // section (the normalized version has an average value of 0.5,
  // so an epsilon of, say, 1e-8 is perfectly usable).
  std::function<double(double)> ndxs = [this, &xstot, &E_level, &Ea, &matrix_element](double cos_theta_c)
    -> double { return (1.0/xstot)*differential_xs(E_level, Ea, matrix_element, cos_theta_c); };

  // Sample a scattering cosine value using this probability density function
  double cos_theta_c = gen.rejection_sample(ndxs, -1, 1);

  return cos_theta_c;
}
