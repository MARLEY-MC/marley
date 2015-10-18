#include "TMarleyGenerator.hh"

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

  // Create the reaction.
  // TODO: Here or somewhere, add checks that we have required decay
  // scheme objects in the database (e.g., 40K should be required),
  // or adjust the code to work around missing target nuclide decay schemes.
  for (const std::string& filename : cf.get_reaction_filenames()) {
    // TODO: rewrite the create_event function so that you don't need to construct
    // a reaction object with a decay scheme. Let the generator handle communication
    // between reactions and nuclear data/models.
    std::cout << "Loading reaction data from file " << filename << std::endl;
    reactions.push_back(TMarleyReaction(filename, structure_db));
  }

  // Create the neutrino source.
  // TODO: Implement configuration file keywords to adjust neutrino source
  // settings and remove hard-coded stuff here.
  // 4.36, 50
  nu_source = TMarleyNeutrinoSource(4.36, //DEBUG reactions.front().get_threshold_energy(),
    50., TMarleyNeutrinoSource::NeutrinoType::ElectronNeutrino);
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
