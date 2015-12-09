#pragma once
#include <set>

#include "marley_utils.hh"

class TMarleyGenerator;

class TMarleyNeutrinoSource {
  public:

    inline TMarleyNeutrinoSource(int particle_id, double w = 1.0) {
      if (pids.count(particle_id) == 0) throw std::runtime_error(
        std::string("Creating a neutrino source object that produces")
        + " particles with PDG ID number " + std::to_string(particle_id)
        + " is not allowed.");
      else pid = particle_id;
      if (w < 0.) throw std::runtime_error(
        std::string("Cannot create a neutrino source object with a negative")
        + " weight.");
      else weight = w;
    }

    // Samples a neutrino energy in MeV
    virtual double sample_energy(TMarleyGenerator& gen) = 0;

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual double get_Emax() const = 0;

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual double get_Emin() const = 0;

    // Returns the PDG particle ID for the neutrino type produced
    // by this source
    virtual inline int get_pid() const { return pid; }

    // Returns the weight of this source (relevant for simulations in which
    // multiple sources are used)
    virtual inline int get_weight() const { return weight; }

  protected:
    // Particle ID for the neutrino type produced by this source
    int pid;
    double weight;

  private:
    // Particle IDs of each neutrino that could possibly be produced by a
    // source object
    static const std::set<int> pids;
};

// Monoenergetic neutrino source
class TMarleyMonoNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    inline TMarleyMonoNeutrinoSource(int particle_id, double E)
      : TMarleyNeutrinoSource(particle_id)
    { 
      energy = E;
    }

    virtual inline double sample_energy(TMarleyGenerator& gen) {
      // This line will suppress unused parameter warnings for this function.
      // We don't need the generator here, but we have to include it in the
      // function declraration to be consistent with the other members of the
      // class hierarchy. This trick was taken from the accepted answer at
      // http://stackoverflow.com/questions/3599160/unused-parameter-warnings-in-c-code
      (void)(gen);
      return energy;
    }

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return energy; }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return energy; }

  private:
    double energy;
};

// Supernova cooling neutrino source approximated using a Fermi-Dirac
// distribution
class TMarleyFermiDiracNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    inline TMarleyFermiDiracNeutrinoSource(int particle_id
      = marley_utils::ELECTRON_NEUTRINO,
      double Emin = 0., double Emax = 100., double temp = 3.5,
      double e_t_a = 0.)
      : TMarleyNeutrinoSource(particle_id)
    {
      E_min = Emin;
      E_max = Emax;
      temperature = temp;
      eta = e_t_a;
      C = 1.;
      // Create a call wrapper for the Fermi-Dirac distribution function that
      // allows our numerical integration and rejection sampling routines to
      // manipulate it easily.
      fd_dist = std::bind(&TMarleyFermiDiracNeutrinoSource::fermi_dirac_distribution,
        this, std::placeholders::_1);
      // Normalize the source spectrum (not strictly necessary, but having the
      // spectrum approximately normalized makes the default rejection sampling
      // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
      // TODO: consider removing the hard-coded value here
      double integral = marley_utils::num_integrate(fd_dist, E_min, E_max, 1e4);
      C /= integral;
    }

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return E_max; }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return E_min; }

    virtual double sample_energy(TMarleyGenerator& gen);

    inline double fermi_dirac_distribution(double nu_energy)
    {
      return (C / std::pow(temperature, 3)) * (std::pow(nu_energy, 2)
        / (1 + std::exp(nu_energy / (temperature - eta))));
    }

  private:
    // Temperature for Fermi-Dirac distribution (in MeV)
    double temperature;
    // Temperature offset for Fermi-Dirac distribution
    double eta;
    // Minimum and maximum neutrino energies produced by this source
    double E_min, E_max;
    // Normalization constant (determined during construction)
    double C;
    // Call wrapper for this object's Fermi Dirac distribution function
    std::function<double(double)> fd_dist;
};

//// Represents a histogram probability density function (values are sampled
//// uniformly over a chosen bin)
//class TMarleyHistogram
//{
//  public:
//    inline TMarleyHistogram() {}
//    inline TMarleyHistogram(const std::vector<double>& lows,
//      const std::vector<double>& ws, double xmax)
//    {
//      size_t lows_size = lows.size();
//      size_t ws_size = ws.size();
//      if (lows_size == 0) throw std::runtime_error(std::string("Empty")
//        + " vector of lower bounds passed to the constructor of"
//        + " TMarleyHistogram.");
//      if (ws_size == 0) throw std::runtime_error(std::string("Empty")
//        + " vector of weights passed to the constructor of TMarleyHistogram.");
//      if (lows_size != ws_size) throw std::runtime_error(std::string("The")
//        + " vectors of lower bounds and weights passed to the constructor"
//        + " of TMarleyHistogram have unequal sizes.");
//
//      double last_lower_bound = marley_utils::minus_infinity;
//      for (size_t j = 0; j < lows_size; ++j) {
//        if (ws.at(j) < 0.) throw std::runtime_error(std::string("Negative")
//          + " bin weight passed to the constructor of TMarleyHistogram.");
//        double current_lower_bound = lower_bounds.at(j);
//        if (current_lower_bound <= last_lower_bound) throw std::runtime_error(
//          std::string("Bin lower bounds passed to the constructor of")
//          + " TMarleyHistogram must be monotonically increasing.");
//        last_lower_bound = current_lower_bound;
//      }
//
//      if (last_lower_bound >= x_max) throw std::runtime_error(
//        std::string("The upper bound for the last bin passed")
//        + " to the constructor of TMarleyHistogram must be"
//        + " greater than the lower bound of the last bin.");
//
//      // The parameters passed to the constructor have passed all of our error
//      // checks, so create the histogram object.
//      lower_bounds = lows;
//      weights = ws;
//      x_max = xmax;
//      last_bin_index = weights.size() - 1;
//
//      std::discrete_distribution<size_t>::param_type
//        new_params(weights.begin(), weights.end());
//      bin_dist.param(new_params);
//    }
//
//    inline double get_xmin() const { return lower_bounds.front(); }
//    inline double get_xmax() const { return x_max; }
//
//    // Samples a bin (using the member bin_dist) and then a value
//    // from that bin (assuming a uniform distribution within the
//    // bin, i.e., on the closed interval [x_bin_min, x_bin_max])
//    double sample_value(TMarleyGenerator& gen);
//
//  private:
//    // Discrete distribution used for sampling bins from the histogram
//    std::discrete_distribution<size_t> bin_dist;
//    // Lower bounds for each bin
//    std::vector<double> lower_bounds;
//    // Weights for each bin
//    std::vector<double> weights;
//    // Upper bound for the last bin
//    double x_max;
//    // Index of last bin (stored during construction for speed)
//    size_t last_bin_index;
//};

// Histogram-based neutrino source (energies are sampled uniformly within each
// bin)
class TMarleyHistogramNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    inline TMarleyHistogramNeutrinoSource(int particle_id,
      const std::vector<double>& bin_lows,
      const std::vector<double>& bin_weights)
      : TMarleyNeutrinoSource(particle_id),
      energy_dist(bin_lows.begin(), bin_lows.end(), bin_weights.begin())
    {
      // TODO: Add checks that the bin boundaries and weights make sense
    }

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return energy_dist.max(); }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return energy_dist.min(); }

    virtual double sample_energy(TMarleyGenerator& gen);

  private:
    std::piecewise_constant_distribution<double> energy_dist;
};

// Neutrino source that uses a tabulated continuous probability density
// function to sample neutrino energies
class TMarleyGridNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    inline TMarleyGridNeutrinoSource(int particle_id,
      const std::vector<double>& energies,
      const std::vector<double>& prob_densities)
      : TMarleyNeutrinoSource(particle_id),
      energy_dist(energies.begin(), energies.end(), prob_densities.begin())
    {
      // TODO: Add checks that the energies and probability densities make sense
    }

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return energy_dist.max(); }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return energy_dist.min(); }

    virtual double sample_energy(TMarleyGenerator& gen);

  private:
    std::piecewise_linear_distribution<double> energy_dist;
};

//class TMarleyGridNeutrinoSource : public TMarleyNeutrinoSource {
//  public:
//    inline TMarleyGridNeutrinoSource(int particle_id,
//      const InterpolationGrid& g) : TMarleyNeutrinoSource(particle_id)
//    {
//      // TODO: Check that energy grid values are nondecreasing
//      // TODO: Check for at least one nonzero pair.second value
//      grid = g;
//      size_t grid_size = grid.size();
//      if (grid_size < 2) throw std::runtime_error(std::string("Grid with")
//        + " less than 2 gridpoints passed to the constructor of"
//        + " TMarleyGridNeutrinoSource.");
//      for (size_t j = 0; j < grid_size; ++j) {
//        auto& pair = grid.at(j);
//        if (pair.first < 0.) throw std::runtime_error(std::string("All energy")
//          + " values used in a TMarleyGridNeutrinoSource"
//          + " object must be nonnegative");
//        // Prevent actually sampling an energy value of zero by advancing to
//        // the next representable double value.
//        else if (pair.first == 0.) pair.first = std::next_after(0.,
//          marley_utils::infinity);
//        if (pair.second < 0.) throw std::runtime_error(std::string("All PDF")
//          + " values used in a TMarleyGridNeutrinoSource"
//          + " object must be nonnegative");
//      }
//    }
//  private:
//    InterpolationGrid<double, double> grid;
//};
