#include <stdexcept>
#include <string>

#include "TMarleyElectronReaction.hh"
#include "TMarleyGenerator.hh"
#include "TMarleyNuclearReaction.hh"

void TMarleyGenerator::init(const TMarleyConfigFile& cf) {
  // Use the seed from the config file object to prepare the random number
  // generator.  This is an attempt to do a decent job of seeding the random
  // number generator, but optimally accomplishing this can be tricky (see, for
  // example, http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
  seed = cf.get_seed();
  std::seed_seq seed_sequence{seed};
  rand_gen.seed(seed_sequence);

  // Load the nuclear structure data specified in the configuration file
  // into the database.
  structure_db = TMarleyStructureDatabase(cf);

  // Create the reactions. Count them for later reference.
  size_t react_count = 0;
  for (const std::string& filename : cf.get_reaction_filenames()) {
    // TODO: add weighting of the ES reaction by atom fraction
    // TODO: allow the user to specify the atomic target for e- ES reactions
    // TODO: Think about a less crude way of implementing this
    if (filename == std::string("e-ES"))
      reactions.push_back(std::make_unique<TMarleyElectronReaction>(19));
    else {
      std::cout << "Loading reaction data from file " << filename << std::endl;
      reactions.push_back(std::make_unique<TMarleyNuclearReaction>(filename,
        structure_db));
    }
    std::cout << "Added reaction " << reactions.back()->get_description()
      << std::endl;
    ++react_count;
  }

  // Create the neutrino source.
  // TODO: Implement configuration file keywords to adjust neutrino source
  // settings and remove hard-coded stuff here.
  //std::vector<double> Es = { 0., 10., 20. };
  //std::vector<double> densities = { 1., 3., 1. };
  //nu_source = TMarleyGridNeutrinoSource(Es, densities,
  //  marley_utils::ELECTRON_NEUTRINO);
  double m_mu = TMarleyMassTable::get_particle_mass(marley_utils::MUON);
  double m_mu_to_the_minus_four = std::pow(m_mu, -4);

  nu_source = TMarleyFunctionNeutrinoSource(
    [m_mu, m_mu_to_the_minus_four](double E) -> double {
      return 96.*std::pow(E,2)*m_mu_to_the_minus_four*(m_mu - 2*E);
    }, 1., marley_utils::ELECTRON_NEUTRINO, 0., m_mu/2.);

  // Initialize the vector of total cross section values to be all zeros and
  // have as many entries as there are reactions available to this generator.
  total_xs_values.clear();
  total_xs_values.resize(react_count, 0.);

  // Compute the normalization factor for use with the reacting neutrino energy PDF
  std::function<double(double)> unnorm_pdf = std::bind(
    &TMarleyGenerator::unnormalized_Ea_pdf, this, std::placeholders::_1);
  double norm = marley_utils::num_integrate(unnorm_pdf, nu_source.get_Emin(),
    nu_source.get_Emax(), 1e3); // TODO: remove hard-coded value here
  // Create the normalized PDF using the normalization factor
  Ea_pdf = [unnorm_pdf, norm](double Ea)
    -> double { return unnorm_pdf(Ea) / norm; };
}

// Sample a random double uniformly between min and max using the class
// member random number generator rand_gen. The inclusive flag
// determines whether or not max is included in the range. That is,
// when inclusive == false, the sampling is done on the interval [min, max),
// while inclusive == true uses [min, max].
double TMarleyGenerator::uniform_random_double(double min, double max,
  bool inclusive)
{
  // Defaults to sampling from [0,1). We will always
  // explicitly supply the upper and lower bounds to
  // this distribution, so we won't worry about the
  // default setting.
  static std::uniform_real_distribution<double> udist;  

  double max_to_use;

  if (inclusive) { // sample from [min, max]

    // Find the double value that comes immediately after max. This allows
    // us to sample uniformly on [min, max] rather than [min,max). This trick comes from
    // http://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution/uniform_real_distribution
    max_to_use = std::nextafter(max, std::numeric_limits<double>::max());
  }
  else { // sample from [min, max)
    max_to_use = max;
  }

  std::uniform_real_distribution<double>::param_type params(min, max_to_use);

  // Sample a random double from this distribution
  return udist(rand_gen, params);
}

// Sample from a given 1D probability density function f(x)
// on the interval [xmin, xmax] using a simple rejection method:
// (1) find the maximum of the function on [xmin, xmax]
// (2) sample an x value uniformly over the function f(x)'s domain
// (3) sample a y value uniformly over [0, max(f(x))]
// (4) if y <= f(x), accept the sampled x value
// (5) if y > f(x), reject the sampled x value, and return
// to step 2 to try again
// Note that f(x) does not need to be normalized, but its range
// must be nonnegative
double TMarleyGenerator::rejection_sample(std::function<double(double)> f,
  double xmin, double xmax, double max_search_tolerance)
{
  // This variable will be loaded with the value of x
  // that corresponds to the maximum of f(x).
  // We don't actually use this, but currently it's
  // a required parameter of marley_utils::maximize
  double x_at_max;

  // Get the maximum value of f(x). This is needed to
  // correctly apply rejection sampling.
  double fmax = marley_utils::maximize(f, xmin, xmax, max_search_tolerance,
    x_at_max);

  double x, y;

  do {
    // Sample x value uniformly from [xmin, xmax]
    x = uniform_random_double(xmin, xmax, true);
    // Sample y uniformly from [0, fmax]
    y = uniform_random_double(0, fmax, true);
  }
  // Keep sampling until you get a y value less than f(x) 
  // (the probability density function evaluated at the sampled value of x)
  while (y > f(x));

  return x;
}

double TMarleyGenerator::unnormalized_Ea_pdf(double Ea) {
  double pdf = 0.;
  // Sum all of the reaction total cross sections, saving
  // each individual value along the way.
  for (size_t j = 0, s = reactions.size(); j < s; ++j) {
    double tot_xs = reactions.at(j)->total_xs(marley_utils::ELECTRON_NEUTRINO,
      Ea); 
    total_xs_values.at(j) = tot_xs;
    pdf += tot_xs;
  }
  // Multiply the total cross section by the neutrino spectrum
  // from the source object to get the (unnormalized) PDF
  // for sampling reacting neutrino energies.
  pdf *= nu_source.pdf(Ea);
  return pdf;
}

// Loads Ea with the energy of a reacting neutrino, and loads
// r_index with the index of the reaction it undergoes. Input
// values for Ea and r_index are ignored.
void TMarleyGenerator::sample_reaction(double& Ea, size_t& r_index) {
  if (reactions.empty()) throw std::runtime_error(std::string("Cannot sample")
    + " a reaction in TMarleyGenerator::sample_reaction. The vector of"
    + " TMarleyReaction objects owned by this generator is empty.");

  // TODO: protect against nu_source changing E_min or E_max after you compute
  // the normalization factor for Ea_pdf in the TMarleyGenerator constructor
  Ea = rejection_sample(Ea_pdf, nu_source.get_Emin(), nu_source.get_Emax());
  // The total cross section values have already been updated by the final
  // call to unnormalized_Ea_PDF during rejection sampling, so we can
  // now sample a reaction type using our discrete distribution object.
  std::discrete_distribution<size_t>::param_type params(total_xs_values.begin(),
    total_xs_values.end());
  r_index = r_index_dist(rand_gen, params);
}
