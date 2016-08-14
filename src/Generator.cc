#include <chrono>
#include <cmath>
#include <string>

#include "marley/ConfigurationFile.hh"
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/Logger.hh"
#include "marley/NuclearReaction.hh"


// The default constructor uses the system time as the seed and a
// default-constructured monoenergetic neutrion source. No reactions are
// defined, so the user must call add_reaction() at least once before using a
// default-constructed Generator.
marley::Generator::Generator()
  : seed_(std::chrono::system_clock::now().time_since_epoch().count()),
  source_(new marley::MonoNeutrinoSource),
  structure_db_(new marley::StructureDatabase), dir_vec_{0., 0., 1}
{
  reseed(seed_);

  // Print the MARLEY logo to the logger stream(s) when the first
  // Generator instance is initialized. If other instances are created,
  // don't reprint the logo.
  static bool printed_logo = false;
  if (!printed_logo) {
    MARLEY_LOG_INFO() << '\n' << marley_utils::marley_logo
      << "\nDon't worry about a thing,\n'Cause every little thing"
      << " gonna be all right.\n-- Bob, \"Three Little Birds\"\n\n"
      << "Model of Argon Reaction Low Energy Yields\n"
      << "version " << marley_utils::MARLEY_VERSION << '\n';
    printed_logo = true;
  }

  MARLEY_LOG_INFO() << "Seed for random number generator: " << seed_;
}

marley::Generator::Generator(marley::ConfigurationFile& cf) {
  init(cf);
}

marley::Generator::Generator(const std::string& filename) {
  marley::ConfigurationFile cf(filename);
  init(cf);
}

marley::Event marley::Generator::create_event() {
  double E_nu;
  marley::Reaction& r = sample_reaction(E_nu);
  marley::Event ev = r.create_event(source_->get_pid(), E_nu, *this);
  if (need_to_rotate_events_) rotate_event(ev);
  return ev;
}

void marley::Generator::seed_using_state_string(std::string& state_string) {
  // TODO: add error handling here (check that state_string is valid)
  std::stringstream strstr(state_string);
  strstr >> rand_gen_;
}

void marley::Generator::reseed(uint_fast64_t seed) {
  // This is an attempt to do a decent job of seeding the random number
  // generator, but optimally accomplishing this can be tricky (see, for
  // example, http://www.pcg-random.org/posts/cpp-seeding-surprises.html)
  seed_ = seed;
  std::seed_seq seed_sequence{seed_};
  rand_gen_.seed(seed_sequence);
}

std::string marley::Generator::get_state_string() const {
  std::stringstream ss;
  ss << rand_gen_;
  return ss.str();
}

void marley::Generator::normalize_E_pdf() {

  // Treat monoenergetic sources differently since they can cause
  // problems for the standard numerical integration check
  if (source_->get_Emin() == source_->get_Emax()) {
    double pdf_test = E_pdf(source_->get_Emin());
    if (pdf_test <= 0. || std::isnan(pdf_test)) {
      throw marley::Error(std::string("The total cross section")
        + " for all defined reactions is <= 0 or NaN for the neutrino"
        + " energy defined in a monoenergetic source. Please verify that"
        + " your neutrino source produces particles above threshold for"
        + " at least one reaction.");
    }
  }
  else {
    // Reset the normalization factor to its default of one until
    // we can calculate the new value. This prevents strange things
    // from happening when we lose precision due to an abnormally
    // high or low norm_ value from a previous source or reaction
    // definition.
    norm_ = 1.;

    // Update the normalization factor for use with the reacting neutrino
    // energy probability density function
    norm_ = marley_utils::num_integrate([this](double E)
      -> double { return this->E_pdf(E); }, source_->get_Emin(),
      source_->get_Emax(), 1e3); // TODO: remove hard-coded value here

    if (norm_ <= 0. || std::isnan(norm_)) {
      throw marley::Error(std::string("The integral of the")
        + " cross-section-weighted neutrino flux is <= 0 or NaN. Please verify"
        + " that your neutrino source spectrum produces significant flux"
        + " above the reaction threshold(s).");
    }
  }
}

void marley::Generator::init(marley::ConfigurationFile& cf) {

  // Use the seed from the config file object to prepare the random number
  // generator.
  reseed(cf.get_seed());

  // Transfer ownership of the structure database from the ConfigurationFile
  // object to the Generator. If the ConfigurationFile object does not own a
  // StructureDatabase object (for some strange reason), then create a
  // default-constructed one.
  auto& cf_sdb = cf.get_structure_db();
  if (cf_sdb) structure_db_.reset(cf_sdb.release());
  else {
    structure_db_ = std::make_unique<marley::StructureDatabase>();
    MARLEY_LOG_WARNING() << "A ConfigurationFile that does not"
      << " own a StructureDatabase object was passed to"
      << " marley::Generator::init()";
  }

  // Create the reactions. Count them for later reference.
  size_t react_count = 0;
  for (const std::string& filename : cf.get_reaction_filenames()) {
    MARLEY_LOG_INFO() << "Loading reaction data from file " << filename;
    reactions_.push_back(std::make_unique<marley::NuclearReaction>(filename,
      *structure_db_));
    MARLEY_LOG_INFO() << "Added reaction "
      << reactions_.back()->get_description();
    ++react_count;
  }

  // Transfer ownership of the neutrino source from the ConfigurationFile to
  // the Generator.
  auto& cf_source = cf.get_source();
  if (cf_source) source_.reset(cf_source.release());
  else throw marley::Error(std::string("Cannot finish creating")
    + " the marley::Generator object. The ConfigurationFile passed to"
    + " marley::Generator::init() does not own a NeutrinoSource object.");

  // Initialize the vector of total cross section values to be all zeros and
  // have as many entries as there are reactions available to this generator.
  total_xs_values_.clear();
  total_xs_values_.resize(react_count, 0.);

  // Update the normalization factor for the reacting neutrino energy
  // distribution
  normalize_E_pdf();

  // Print the MARLEY logo to the logger stream(s) when the first
  // Generator instance is initialized. If other instances are created,
  // don't reprint the logo.
  static bool printed_logo = false;
  if (!printed_logo) {
    MARLEY_LOG_INFO() << '\n' << marley_utils::marley_logo
      << "\nDon't worry about a thing,\n'Cause every little thing"
      << " gonna be all right.\n-- Bob, \"Three Little Birds\"\n\n"
      << "Model of Argon Reaction Low Energy Yields\n"
      << "version " << marley_utils::MARLEY_VERSION << '\n';
    printed_logo = true;
  }

  MARLEY_LOG_INFO() << "Seed for random number generator: " << seed_;

  // Initialize the incident neutrino direction settings using the direction
  // given in the configuration file
  std::array<double, 3> nu_dir = cf.get_neutrino_direction();
  set_neutrino_direction(nu_dir);
}

// Sample a random double uniformly between min and max using the class
// member random number generator rand_gen_. The inclusive flag
// determines whether or not max is included in the range. That is,
// when inclusive == false, the sampling is done on the interval [min, max),
// while inclusive == true uses [min, max].
double marley::Generator::uniform_random_double(double min, double max,
  bool inclusive)
{
  // Defaults to sampling from [0,1). We will always
  // explicitly supply the upper and lower bounds to
  // this distribution, so we won't worry about the
  // default setting.
  static std::uniform_real_distribution<double> udist;

  double max_to_use;

  if (inclusive) { // sample from [min, max]

    // Find the double value that comes immediately after max. This allows us
    // to sample uniformly on [min, max] rather than [min,max). This trick
    // comes from http://tinyurl.com/n3ocg3p.
    max_to_use = std::nextafter(max, std::numeric_limits<double>::max());
  }
  else { // sample from [min, max)
    max_to_use = max;
  }

  std::uniform_real_distribution<double>::param_type params(min, max_to_use);

  // Sample a random double from this distribution
  return udist(rand_gen_, params);
}

/// @details The rejection method used by this function consists of the
/// following steps:
/// <ol><li>Find the maximum of the function f(x) on [xmin, xmax]</li>
/// <li>Sample an x value uniformly over f(x)'s domain</li>
/// <li>Sample a y value uniformly over [0, max(f(x))]</li>
/// <li>If y <= f(x), accept the sampled x value</li>
/// <li>If y > f(x), reject the sampled x value, and return
/// to step 2 to try again</li> </ol>
/// Note that f(x) does not need to be normalized, but its range must be
/// nonnegative. In the first step, an iterative method (<a
/// href="http://tinyurl.com/ntqkfck">Brent's method</a>) is
/// used to find the maximum of f(x). The iterations will continue until two
/// successive iterations agree within max_search_tolerance on the location of
/// the maximum of f(x). To avoid problems with functions that yield double
/// values that are small compared to a typical value of max_search_tolerance
/// (say, max_search_tolerance = 1e-8, while many neutrino cross sections of
/// interest for MARLEY are less than 1e-40 cm^2), MARLEY normalizes all
/// probability density functions to unity before using rejection sampling.
/// @todo Check the convergence explanation for the first step.
double marley::Generator::rejection_sample(std::function<double(double)> f,
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

double marley::Generator::E_pdf(double E) {
  double pdf = 0.;
  // Sum all of the reaction total cross sections, saving
  // each individual value along the way.
  for (size_t j = 0, s = reactions_.size(); j < s; ++j) {
    double tot_xs = reactions_.at(j)->total_xs(source_->get_pid(), E);
    total_xs_values_.at(j) = tot_xs;
    pdf += tot_xs;
  }
  // Multiply the total cross section by the neutrino spectrum
  // from the source object to get the (unnormalized) PDF
  // for sampling reacting neutrino energies.
  pdf *= source_->pdf(E);

  //  Divide by the normalization factor (computed when this source
  //  was made available to the Generator) to obtain the normalized PDF.
  return pdf / norm_;
}

marley::Reaction& marley::Generator::sample_reaction(double& E) {
  if (reactions_.empty()) throw marley::Error(std::string("Cannot sample")
    + " a reaction in marley::Generator::sample_reaction(). The vector of"
    + " marley::Reaction objects owned by this generator is empty.");

  // TODO: protect against source_ changing E_min or E_max after you compute
  // the normalization factor norm_ in marley::Generator::init()
  E = rejection_sample([this](double E_nu)
    -> double { return this->E_pdf(E_nu); }, source_->get_Emin(),
    source_->get_Emax());
  // The total cross section values have already been updated by the final call
  // to E_pdg during rejection sampling, so we can now sample a reaction type
  // using our discrete distribution object.
  std::discrete_distribution<size_t>::param_type
    params(total_xs_values_.begin(), total_xs_values_.end());
  size_t r_index = r_index_dist_(rand_gen_, params);
  return *reactions_.at(r_index);
}

const marley::NeutrinoSource& marley::Generator::get_source() {
  if (source_) return *source_;
  else throw marley::Error(std::string("Error")
    + " in marley::Generator::get_source(). The member variable source_ =="
    + " nullptr.");
}

void marley::Generator::set_source(
  std::unique_ptr<marley::NeutrinoSource> source)
{
  // If we're passed a nullptr, then don't bother to do anything
  if (source) {
    // Transfer ownership of the neutrino source the generator, leaving
    // the std::unique_ptr passed to this function null afterwards.
    source_.reset(source.release());

    // Don't bother to renormalize if there are no reactions defined yet
    if (!reactions_.empty()) {
      // Update the neutrino energy probability density function based on the
      // new source spectrum
      normalize_E_pdf();
    }
  }
}

void marley::Generator::add_reaction(std::unique_ptr<marley::Reaction> reaction)
{
  // If we're passed a nullptr, then don't bother to do anything
  if (reaction) {

    // Transfer ownership to a new unique_ptr in the reactions vector, leaving
    // the original empty
    reactions_.push_back(std::move(reaction));

    // Add a new entry in the reaction cross sections vector
    total_xs_values_.push_back(0.);

    // TODO: consider adding a check to see whether source_ is non-null.
    // Right now, this shouldn't be possible, but an explicit check might
    // be good.

    // Update the neutrino energy probability density function by including the
    // cross section for the new reaction
    normalize_E_pdf();
  }
}

void marley::Generator::clear_reactions() {
  reactions_.clear();
  total_xs_values_.clear();
  // Reset the normalization factor to 1. We don't need it until we define
  // one or more new reactions.
  norm_ = 1.;
}

marley::StructureDatabase& marley::Generator::get_structure_db() {
  if (structure_db_) return *structure_db_;
  else throw marley::Error(std::string("Error")
    + " in marley::Generator::get_structure_db(). The member variable"
    + " structure_db_ == nullptr.");
}

void marley::Generator::rotate_event(marley::Event& ev) {

  // Rotate the initial particles
  for (auto* p : ev.get_initial_particles())
    rotation_matrix_.rotate_particle_inplace(*p);

  // Rotate the final particles
  for (auto* p : ev.get_final_particles())
    rotation_matrix_.rotate_particle_inplace(*p);
}

void marley::Generator::set_neutrino_direction(
  const std::array<double, 3> dir_vec)
{

  static constexpr std::array<double, 3> null_three_vector = { 0., 0., 0. };

  if (dir_vec == null_three_vector)
    throw marley::Error(std::string("Null vector")
      + " passed to marley::Generator::set_neutrino_direction()");

  // Get the new incident neutrino direction vector
  dir_vec_ = marley::RotationMatrix::normalize(dir_vec);

  // Print a log message announcing the change of direction
  std::string dir_msg("Incident neutrino direction: (");
  for (size_t i = 0; i < 3; ++i) {
    dir_msg += std::to_string(dir_vec_[i]);
    if (i < 2) dir_msg += ", ";
  }
  MARLEY_LOG_INFO() << dir_msg << ')';

  // Get the default incident neutrino direction
  std::array<double, 3> nu_default_dir = marley::RotationMatrix::normalize(
    marley::ConfigurationFile::get_default_neutrino_direction());

  // Update the rotation matrix and the event rotation flag as needed
  if (dir_vec_ != nu_default_dir)
  {
    need_to_rotate_events_ = true;
    rotation_matrix_ = marley::RotationMatrix(nu_default_dir, dir_vec_);
  }
  else {
    need_to_rotate_events_ = false;
    rotation_matrix_ = marley::RotationMatrix(); // identity matrix
  }
}
