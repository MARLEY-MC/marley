#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "marley_utils.hh"
#include "TMarleyEvent.hh"
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

  line = get_next_line(file_in, rx_comment, false);

  // Read in the particle IDs
  std::istringstream iss(line);
  iss >> pid_a >> pid_b >> pid_c >> pid_d;

  // Get initial and final values of the nuclear
  // charge and mass number from the pids
  Zi = (pid_b % 10000000)/10000;
  Ai = (pid_b % 10000)/10;
  Zf = (pid_d % 10000000)/10000;
  Af = (pid_b % 10000)/10;

  // Read in the particle masses
  line = get_next_line(file_in, rx_comment, false);
  iss.str(line); 
  iss.clear();
  iss >> ma >> mb >> mc >> md_gs;

  Ea_threshold = ((mc + md_gs)*(mc + md_gs) - ma*ma - mb*mb)/(2*mb);

  // Read in the number of levels that have tabulated
  // B(F) + B(GT) data
  int num_levels;
  line = get_next_line(file_in, rx_comment, false);
  iss.str(line); 
  iss.clear();
  iss >> num_levels;

  // Read in all of the level energies (in MeV)
  int j = 0;
  double entry, old_entry;
  while (j < num_levels && file_in.good()) {
    line = get_next_line(file_in, rx_comment, false);
    iss.str(line);
    iss.clear();
    while (iss >> entry) {
      // TODO: consider implementing a sorting procedure rather than strictly
      // enforcing that energies must be given in ascending order. Note that
      // you will need to alter the order of the B(F) + B(GT) values
      // as well if you sort the energies.
      
      // The order of the entries is important because later uses of the
      // residue_level_energies vector assume that they are sorted in
      // order of ascending energy.
      if (old_entry >= entry) throw std::runtime_error(std::string("Invalid reaction dataset. ")
        + "Level energies must be unique and must be given in ascending order.");
      residue_level_energies.push_back(entry);
      ++j;
      old_entry = entry;
    }
  }

  // Read in all of the level B(F) + B(GT) values
  j = 0;
  while (j < num_levels && file_in.good()) {
    line = get_next_line(file_in, rx_comment, false);
    iss.str(line);
    iss.clear();
    while (iss >> entry) {
      residue_level_strengths.push_back(entry);
      ++j;
    }
  }

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

// Advance to the next line of a file that either matches (match == true)
// or does not match (match == false) a given regular expression
std::string TMarleyReaction::get_next_line(std::ifstream &file_in,
  std::regex &rx, bool match) const {

  std::string line;
  while (!file_in.eof()) {
    // Get the next line of the file
    std::getline(file_in, line);
 
    // Check to see if the new line fulfills the search criteria
    if (std::regex_match(line, rx) == match) {
      // If it does, return it
      return line;
    }
    // If not, keep looking
  }

  // If the end of the file is encountered before a suitable
  // line is found, return an empty string
  return std::string("");
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

  // Cycle through each of the level energies given in the reaction dataset. 
  for(std::vector<double>::iterator it = residue_level_energies.begin();
  it != residue_level_energies.end(); ++it)
  {
    // For each one, find a pointer to the level with the closest energy owned
    // by the decay scheme object.
    TMarleyLevel* plevel = ds->get_pointer_to_closest_level(*it);
    //std::cout << "DEBUG: I matched E = " << *it << " MeV to the ENSDF level "
    //  << "with energy " << plevel->get_numerical_energy() << " MeV" << std::endl;

    // Complain if there are duplicates (if there are duplicates, we'll have
    // two different B(F) + B(GT) values for the same level object)
    if (std::find(residue_level_pointers.begin(), residue_level_pointers.end(),
      plevel) != residue_level_pointers.end()) {
      // residue_level_pointers already contains plevel
      throw std::runtime_error(std::string("Reaction dataset gives two level energies ")
        + "that refer to the same ENSDF level at "
        + std::to_string(plevel->get_numerical_energy()) + " MeV");
    }

    // TODO: add check to see if the energy of the chosen level is very different
    // from the energy given in the reaction dataset. If it is, the level matchup
    // is likely incorrect.

    // Add the ENSDF level pointer to the list
    residue_level_pointers.push_back(plevel);
  }

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
TMarleyEvent TMarleyReaction::create_event(double Ea) {

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

    // Get the excitation energy for the current level
    double level_energy = residue_level_pointers[i]->get_numerical_energy();
    //std::cout << "DEBUG: Found level with energy = " << level_energy << std::endl;

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
  //  std::cout << "DEBUG: level at " << residue_level_pointers.at(i)->get_numerical_energy()
  //    << " MeV has probability " << n << " of being selected" << std::endl;
  //}

  // Sample a level index using our discrete distribution and the
  // current set of weights
  unsigned int l_index = ldist(marley_utils::rand_gen, params);

  // Get a pointer to the selected level
  TMarleyLevel* plevel = residue_level_pointers.at(l_index);

  // Get the energy of the selected level. This will be
  // needed for sampling a scattering cosine.
  double E_level = plevel->get_numerical_energy();
  double md = md_gs + E_level;

  // Get the matrix element (B(F) + B(GT) value) for the
  // selected level
  double matrix_element = residue_level_strengths.at(l_index);

  // Sample a scattering angle for the ejectile using
  // the differential cross section
  double cos_theta_c = sample_ejectile_scattering_cosine(E_level, Ea, matrix_element);
  double theta_c = std::acos(cos_theta_c);

  // Use conservation of 4-momentum to compute the ejectile energy
  // based on the sampled scattering angle
  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);

  // Compute the energy and scattering angle of the residue
  double Ed = Etot - Ec;
  double pc = marley_utils::real_sqrt(Ec*Ec - mc*mc);
  double pd_before_gammas = marley_utils::real_sqrt(Ed*Ed - md*md);
  double theta_d = std::asin(pc * std::sin(theta_c) / pd_before_gammas);
  double cos_theta_d = std::cos(theta_d);

  // Sample an azimuthal scattering angle (phi) uniformly on [0, 2*pi).
  // We can do this because the matrix elements are azimuthally invariant
  double phi = marley_utils::uniform_random_double(0, 2*marley_utils::pi, false);

  // TODO: maybe consider the effect of nuclear recoils that happen when
  // gammas are emitted
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
  std::cout << "E_level = " << E_level << std::endl;
  std::cout << "cos_theta_c = " << cos_theta_c << std::endl;
  std::cout << "Ec = " << Ec << std::endl;
  std::cout << "Ed = " << Ed << std::endl;
  std::cout << "cos_theta_d = " << cos_theta_d << std::endl;
  std::cout << "e- kinetic energy = " << Ec - mc << std::endl;
  std::cout << "e- mass = " << mc << std::endl;
  std::cout << "40K kinetic energy = " << Ed - md_gs - E_level << std::endl;
  std::cout << "Ground state nuclear mass change = " << md_gs - mb << std::endl; 

  // Create the event object and load it with the appropriate information
  TMarleyEvent event(E_level);
  event.set_reaction(this);
  // TODO: edit this to allow for projectile directions other than along the z-axis
  // Add the projectile to this event's initial particle list
  event.add_initial_particle(TMarleyParticle(pid_a, Ea, 0, 0, Ea - ma, ma),
    TMarleyEvent::ParticleRole::projectile);
  // Add the target to this event's initial particle list
  event.add_initial_particle(TMarleyParticle(pid_b, mb, 0, 0, 0, mb),
    TMarleyEvent::ParticleRole::target);
  // Add the ejectile to this event's final particle list
  event.add_final_particle(TMarleyParticle(pid_c, Ec, pc_x, pc_y, pc_z, mc),
    TMarleyEvent::ParticleRole::ejectile);
  // Add the residue to this event's final particle list. Don't include
  // its excitation energy since we will soon create the de-excitation gamma rays
  event.add_final_particle(TMarleyParticle(pid_d, Ed_gs, pd_x_gs, pd_y_gs, pd_z_gs, md_gs),
    TMarleyEvent::ParticleRole::residue);

  // Add the de-excitation gammas to this event's final particle list
  // TODO: consider whether including 
  this->ds->do_cascade(plevel, &event);

  // Return the completed event object
  return event;
}

// Compute the differential reaction cross section dsigma/dcos_theta
// for a given final residue level, projectile energy, and ejectile
// scattering cosine
double TMarleyReaction::differential_xs(double E_level, double Ea,
  double matrix_element, double cos_theta_c)
{

  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);
  double pc = marley_utils::real_sqrt(Ec*Ec - mc*mc);

  // TODO: consider removing the constants GF, Vud, etc. since
  // for sampling they will not be needed (what matters is the
  // relative size of the differential cross section).
  // This may reduce numerical error.
  
  // TODO: adjust this differential cross section expression
  // as needed
  
  // This expression for the differential cross section
  // dsigma/dcos_theta comes from Kuramoto, et al.,
  // Nucl. Phys. A512 (1990) 711-736. It is modified
  // using a nuclear recoil correction factor taken from
  // J. D. Walecka, "Semileptonic Weak Interactions in Nuclei,"
  // In: Muon Physics, Volume II: Weak Interactions.
  // Ed. by V. W. Hughes and C. S. Wu.
  return (1.0/(2*std::acos(-1)))*GF*GF*Vud*Vud
    *pc*Ec*fermi_function(Zf, Af, Ec, true)*matrix_element
    /(1.0 + (Ea/mb)*(1 - (Ec/pc)*cos_theta_c));
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
  double Ea, double matrix_element) {

  // Get the total cross section for this reaction based
  // on the projectile energy and final residue energy.
  // This will be used as a normalization factor.
  double xstot = total_xs(E_level, Ea, matrix_element);

  // Make a normalized version of the differential cross section
  // to use as our probability density function for sampling.
  // Note that, since we are using a rejection sampling method,
  // this is not strictly necessary. That being said, doing so
  // allows us to set our value of epsilon for marley_utils::maximize
  // without having to worry about the absolute scale of the differential cross
  // section (the normalized version has an average value of 0.5,
  // so an epsilon of, say, 1e-8 is perfectly usable).
  std::function<double(double)> ndxs = [this, &xstot, &E_level, &Ea, &matrix_element](double cos_theta_c)
    -> double { return (1.0/xstot)*differential_xs(E_level, Ea, matrix_element, cos_theta_c); };

  // This variable will be loaded with the ejectile cosine that
  // corresponds to the largest differential xs
  // We don't actually use this, but currently it's
  // a required parameter of marley_utils::maximize
  double cmax;

  // Get the maximum value of the normalized differential cross section for
  // the given projectile energy and final residue energy.
  // This is needed to correctly apply rejection sampling.
  double max_ndxs = marley_utils::maximize(ndxs, -1.0, 1.0, 1e-8, cmax);

  double cos_theta_c, height;

  do {
    // Sample cosine value uniformly from [-1,1]
    cos_theta_c = marley_utils::uniform_random_double(-1, 1, true);
    // Sample height uniformly from [0, max_ndxs]
    height = max_ndxs*marley_utils::uniform_random_double(0, 1, true);
  }
  // Keep sampling until you get a height value less than the normalized
  // differential cross section evaluated at cos_theta_c
  while (height > ndxs(cos_theta_c));

  return cos_theta_c;
}
