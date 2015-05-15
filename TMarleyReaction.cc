#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

#include "marley_utils.hh"
#include "TMarleyReaction.hh"
#include "TMarleyLevel.hh"

TMarleyReaction::TMarleyReaction() {
  // TODO: add stuff here
  ds = nullptr;
}

// Associate a decay scheme object with this reaction. This will provide
// nuclear structure data for sampling final residue energy levels.
void TMarleyReaction::set_decay_scheme(TMarleyDecayScheme* scheme) {
  ds = scheme;
}

// TODO: Fix both Fermi functions to use natural units and energies in MeV (this keeps
// everything consistent)

// Fermi function used in calculating cross-sections
// The equation was taken from http://en.wikipedia.org/wiki/Beta_decay
// Input: atomic number Z and electron kinetic energy T (in keV), positron or electron?
double TMarleyReaction::fermi_function(double z, double t, bool electron){
  
  double charge = z;
  double kin_e = t;
  double atomic_mass = 40.;
  
  // Constants
  const double alpha = 7.297e-3; // Fine-structure constant
  const double hbar_c = 1.97e+3; // hbar*c in keV*fm
  const double e = std::exp(1); // Euler's number
  const double pi = 4*std::atan(1);
  const double r_0 = 1.2; // in fm
  const double m = 510.9989; // in keV

  double eta;
  double p = std::sqrt(kin_e*kin_e + 2*m*kin_e); // Electron momentum
  double s = std::sqrt(1 - alpha*alpha*charge*charge);
  double energy = std::sqrt(p*p + m*m);
  double r_n = r_0*std::pow(atomic_mass, 1./3);
  double rho = r_n/hbar_c;
  
  if(electron)
    eta = alpha*charge*energy/p;
  else
    eta = -alpha*charge*energy/p;

  // Complex variables for the gamma function
  std::complex<double> a(s, eta);
  std::complex<double> b(1+2*s, 0);

  return 2*(1 + s)* std::pow(2*p*rho, 2*s-2)*std::pow(e, pi*eta)*std::norm(marley_utils::gamma(a))
    / (std::abs(marley_utils::gamma(b)) * std::abs(marley_utils::gamma(b)));
}

double TMarleyReaction::fermi_approx(int z, double t, bool electron){

  // This is a test for comparison with the "exact" calculation, which involves
  // complex gamma functions. This holds true for Q << mc^2
  
  int charge = z;
  double kin_e = t;

  const double alpha = 7.297e-3; // Fine-structure constant
  const double e = std::exp(1); // Euler's number
  const double pi = 4*std::atan(1);
  const double m = 510.9989; // in keV

  double eta;
  double p = std::sqrt(kin_e*kin_e + 2*kin_e*m);
  double energy = std::sqrt(p*p + m*m);
  
  if(electron)
    eta = alpha*charge*energy/p;
  else
    eta = -alpha*charge*energy/p;

  return 2*pi*eta/(1-std::pow(e, -2*pi*eta));
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

// TODO: change this to return an event object
void TMarleyReaction::create_event(double Ea) {

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

  // Get level information from the decay scheme
  // object associated with this reaction. Complain
  // if problems occur
  if (ds == nullptr) {
    throw std::runtime_error(std::string("Could not create this event. No nuclear ")
      + "structure data is associated with this reaction.");
  } 
  
  std::vector<TMarleyLevel*>* pv_levels = ds->get_sorted_level_pointers();

  if (pv_levels->empty()) {
    throw std::runtime_error(std::string("Could not create this event. The ")
      + "TMarleyDecayScheme object associated with this reaction "
      + "does not contain any level data.");
  }

  // The pointers in pv_levels are ordered by increasing energy.
  // Iterate over the levels, assigning the total cross section
  // for each level as its weight until you reach the end of pv_levels
  // or a level that is kinematically forbidden
  for(std::vector<TMarleyLevel*>::iterator j = pv_levels->begin();
    j != pv_levels->end(); ++j) {

    // Get the excitation energy for the current level
    double level_energy = (*j)->get_numerical_energy();
    std::cout << "DEBUG: Found level with energy = " << level_energy << std::endl;

    // Exit the loop early if you reach a level with an energy that's too high
    if (level_energy > max_E_level) break;

    // TODO: add a check here for whether the matrix element (B(F) + B(GT))
    // is nonvanishing for the current level. If it is, just set the weight
    // equal to zero rather than calling total_xs. This will avoid unnecessary
    // numerical integrations.

    // Assign a weight to this level equal to the total reaction
    // cross section. Note that std::discrete_distribution automatically
    // normalizes the weights, so we don't have to do that ourselves.

    //DEBUG
    double xs = total_xs(level_energy, Ea);
    if (std::isnan(xs)) {
      std::cout << "DEBUG: this level gave a weight of nan, so I made it zero." << std::endl;
      level_weights.push_back(0.0);
    }
    else {
    level_weights.push_back(xs);
    std::cout << "DEBUG: This level was assigned weight = "
      << xs << std::endl;
    }
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
  //for (auto& n : probs)
  //{
  //  int i = &n - &(probs[0]);
  //  std::cout << "DEBUG: level at " << pv_levels->at(i)->get_numerical_energy()
  //    << " MeV has probability " << n << " of being selected" << std::endl;
  //}


  // Sample a level index using our discrete distribution and the
  // current set of weights
  unsigned int l_index = ldist(marley_utils::rand_gen, params);

  // Get a pointer to the selected level
  TMarleyLevel* plevel = pv_levels->at(l_index);

  // Get the energy of the selected level. This will be
  // needed for sampling a scattering cosine.
  double E_level = plevel->get_numerical_energy();
  double md = md_gs + E_level;

  // Sample a scattering angle for the ejectile using
  // the differential cross section
  double cos_theta_c = sample_ejectile_scattering_cosine(E_level, Ea);
  double theta_c = std::acos(cos_theta_c);

  // Use conservation of 4-momentum to compute the ejectile energy
  // based on the sampled scattering angle
  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);

  // Compute the energy and scattering angle of the residue
  double Ed = Etot - Ec;
  double theta_d = std::asin(marley_utils::real_sqrt(Ec*Ec - mc*mc)
    * std::sin(theta_c) / marley_utils::real_sqrt(Ed*Ed - md*md));
  double cos_theta_d = std::cos(theta_d);

  // Print results to std::cout
  //std::cout.precision(15);
  std::cout << std::scientific;
  std::cout << "E_level = " << E_level << std::endl;
  std::cout << "cos_theta_c = " << cos_theta_c << std::endl;
  std::cout << "Ec = " << Ec << std::endl;
  std::cout << "Ed = " << Ed << std::endl;
  std::cout << "cos_theta_d = " << cos_theta_d << std::endl;
  
  // Simulate de-excitation gammas
  this->ds->do_cascade(plevel);
}

// Compute the differential reaction cross section dsigma/dcos_theta
// for a given final residue level, projectile energy, and ejectile
// scattering cosine
double TMarleyReaction::differential_xs(double E_level, double Ea, double cos_theta_c) {

  double Ec = this->ejectile_energy(E_level, Ea, cos_theta_c);
  double pc = marley_utils::real_sqrt(Ec*Ec - mc*mc);

  // TODO: compute the matrix element using B(F) + B(GT)
  // values and/or some other method
  double matrix_element = 1;
  
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
    *pc*Ec*fermi_function(Zf, Ec - mc, true)*matrix_element
    /(1.0 + (Ea/mb)*(1 - (Ec/pc)*cos_theta_c));
}

// Numerically integrate the differential cross section over
// cos(theta_c) = [-1,1] to get the total reaction cross section
// for a given final residue level and projectile energy
double TMarleyReaction::total_xs(double E_level, double Ea) {
  // TODO: consider switching to a different integration method,
  // using tabulated data, etc.

  // Numerically integrate the differential cross section over
  // the interval cos(theta_c) = [-1,1] using the composite
  // trapezoidal rule over n subintervals
  static const int n = 1e3;

  // Create a forwarding call wrapper for the  differential_xs
  // member function that takes a single argument.
  std::function<double(double)> dxs = std::bind(
    &TMarleyReaction::differential_xs, this, E_level, Ea, std::placeholders::_1);

  // Numerically integrate using the call wrapper, the integration bounds,
  // and the number of subintervals
  return marley_utils::num_integrate(dxs, -1.0, 1.0, n); 
}

// Sample a scattering cosine for the ejectile using the differential
// cross section defined in the member function differential_xs
double TMarleyReaction::sample_ejectile_scattering_cosine(double E_level, double Ea) {

  // Get the total cross section for this reaction based
  // on the projectile energy and final residue energy.
  // This will be used as a normalization factor.
  double xstot = total_xs(E_level, Ea);

  // Make a normalized version of the differential cross section
  // to use as our probability density function for sampling.
  // Note that, since we are using a rejection sampling method,
  // this is not strictly necessary. That being said, doing so
  // allows us to set our value of epsilon for marley_utils::maximize
  // without having to worry about the absolute scale of the differential cross
  // section (the normalized version has an average value of 0.5,
  // so an epsilon of, say, 1e-8 is perfectly usable).
  std::function<double(double)> ndxs = [this, &xstot, &E_level, &Ea](double cos_theta_c)
    -> double { return (1.0/xstot)*differential_xs(E_level, Ea, cos_theta_c); };

  // This variable will be loaded with the ejectile cosine that
  // corresponds to the largest differential xs
  // We don't actually use this, but currently it's
  // a required parameter of marley_utils::maximize
  double cmax;

  // Get the maximum value of the normalized differential cross section for
  // the given projectile energy and final residue energy.
  // This is needed to correctly apply rejection sampling.
  double max_ndxs = marley_utils::maximize(ndxs, -1.0, 1.0, 1e-8, cmax);

  // Find the double value that comes immediately after +1.0. This allows
  // us to sample uniformly on [0,1] and [-1,1] rather than [0,1) and [-1,1).
  // This trick comes from
  // http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution/uniform_real_distribution
  static double after_one = std::nextafter(1.0, std::numeric_limits<double>::max());

  // Create a uniform distribution object to perform the rejection sampling.
  // Also create two sets of distribution parameters so that we can sample
  // the horizontal and vertical parts using the same object.
  static std::uniform_real_distribution<double> udist; // Defaults to sampling from [0,1). We will always
                                                       // explicitly supply the upper and lower bounds to
                                                       // this distribution, so we won't worry about the
                                                       // default setting.
  static std::uniform_real_distribution<double>::param_type height_params(0, after_one); // [0,1]
  static std::uniform_real_distribution<double>::param_type cos_params(-1, after_one); // [-1,1]

  double cos_theta_c, height;

  do {
    // Sample cosine value uniformly from [-1,1]
    cos_theta_c = udist(marley_utils::rand_gen, cos_params); 
    // Sample height uniformly from [0, max_ndxs]
    height = max_ndxs*udist(marley_utils::rand_gen, height_params);
  }
  // Keep sampling until you get a height value less than the normalized
  // differential cross section evaluated at cos_theta_c
  while (height > ndxs(cos_theta_c));

  return cos_theta_c;
}
