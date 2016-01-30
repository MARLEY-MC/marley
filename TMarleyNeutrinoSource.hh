#pragma once
#include <functional>
#include <set>

#include "marley_utils.hh"
#include "InterpolationGrid.hh"

class TMarleyGenerator;

class TMarleyNeutrinoSource {
  public:

    inline TMarleyNeutrinoSource(int particle_id, double w = 1.0) {
      if (!pid_is_allowed(particle_id)) throw std::runtime_error(
        std::string("Creating a neutrino source object that produces")
        + " particles with PDG ID number " + std::to_string(particle_id)
        + " is not allowed.");
      else pid = particle_id;
      if (w < 0.) throw std::runtime_error(
        std::string("Cannot create a neutrino source object with a negative")
        + " weight.");
      else weight = w;
    }

    // Returns true if the particle ID passed to the function is allowed to be
    // used by a neutrino source object, and returns false otherwise.
    static inline bool pid_is_allowed(const int particle_id) {
      return (pids.count(particle_id) > 0);
    }

    // Samples a neutrino energy in MeV
    //virtual double sample_energy(TMarleyGenerator& gen) = 0;

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

    // Probability density function (not necessarily normalized) used by
    // TMarleyGenerator for folding the neutrino spectrum produced by this
    // source with the relevant cross sections.
    virtual double pdf(double E_nu) = 0;

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
    inline TMarleyMonoNeutrinoSource(int particle_id
      = marley_utils::ELECTRON_NEUTRINO, double weight = 1.0,
      double E = 10./* MeV*/) : TMarleyNeutrinoSource(particle_id, weight)
    { 
      energy = E;
    }

    //virtual inline double sample_energy(TMarleyGenerator& gen) {
    //  // This line will suppress unused parameter warnings for this function.
    //  // We don't need the generator here, but we have to include it in the
    //  // function declraration to be consistent with the other members of the
    //  // class hierarchy. This trick was taken from the accepted answer at
    //  // http://stackoverflow.com/questions/3599160/unused-parameter-warnings-in-c-code
    //  (void)(gen);
    //  return energy;
    //}

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return energy; }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return energy; }

    // Monoenergetic spectrum
    virtual inline double pdf(double E_nu) {
      if (energy == E_nu) return 1.;
      else return 0.;
    }

  private:
    double energy;
};

// Supernova cooling neutrino source approximated using a Fermi-Dirac
// distribution
class TMarleyFermiDiracNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    inline TMarleyFermiDiracNeutrinoSource(int particle_id
      = marley_utils::ELECTRON_NEUTRINO, double weight = 1.,
      double Emin = 0., double Emax = 100., double temp = 3.5,
      double e_t_a = 0.)
      : TMarleyNeutrinoSource(particle_id, weight)
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

    //virtual double sample_energy(TMarleyGenerator& gen);

    inline double fermi_dirac_distribution(double nu_energy)
    {
      return (C / std::pow(temperature, 3)) * (std::pow(nu_energy, 2)
        / (1 + std::exp(nu_energy / (temperature - eta))));
    }

    virtual inline double pdf(double E_nu) {
      if (E_nu < E_min || E_nu > E_max) return 0.;
      else return fermi_dirac_distribution(E_nu);
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

// Neutrino source with an arbitrary energy spectrum supplied during
// construction as a std::function<double(double)> object
class TMarleyFunctionNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    inline TMarleyFunctionNeutrinoSource(const std::function<double(double)>&
      prob_dens_func = [](double E) -> double { (void)(E); return 1; },
      double weight = 1., int particle_id = marley_utils::ELECTRON_NEUTRINO,
      double Emin = 0., double Emax = 100.)
      : TMarleyNeutrinoSource(particle_id, weight)
    {
      E_min = Emin;
      E_max = Emax;
      // Normalize the supplied spectrum (not strictly necessary, but having the
      // spectrum approximately normalized makes the default rejection sampling
      // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
      // TODO: consider removing the hard-coded value here
      double integral = marley_utils::num_integrate(prob_dens_func, E_min,
        E_max, 1e4);
      probability_density = [prob_dens_func, integral](double E)
        -> double { return prob_dens_func(E) / integral; };
    }

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return E_max; }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return E_min; }

    virtual inline double pdf(double E_nu) {
      if (E_nu < E_min || E_nu > E_max) return 0.;
      else return probability_density(E_nu);
    }

  private:
    // Minimum and maximum neutrino energies produced by this source
    double E_min, E_max;
    // Call wrapper for this object's probability density function
    std::function<double(double)> probability_density;
};

class TMarleyGridNeutrinoSource : public TMarleyNeutrinoSource {
  public:
    using Grid = InterpolationGrid<double>;
    using InterpolationMethod = Grid::InterpolationMethod;

    inline TMarleyGridNeutrinoSource(int particle_id
      = marley_utils::ELECTRON_NEUTRINO, double weight = 1.,
      InterpolationMethod interp_method = InterpolationMethod::LinearLinear)
      : TMarleyNeutrinoSource(particle_id, weight), grid(interp_method)
    {
      //check_for_errors();
    }

    inline TMarleyGridNeutrinoSource(const Grid& g,
      int particle_id = marley_utils::ELECTRON_NEUTRINO, double weight = 1.)
      : TMarleyNeutrinoSource(particle_id, weight), grid(g)
    {
      check_for_errors();
    }

    inline TMarleyGridNeutrinoSource(const std::vector<double>& Es,
      const std::vector<double>& prob_densities, int particle_id
      = marley_utils::ELECTRON_NEUTRINO, double weight = 1., InterpolationMethod
      interp_method = InterpolationMethod::LinearLinear)
      : TMarleyNeutrinoSource(particle_id, weight), grid(Es, prob_densities,
      interp_method)
    {
      check_for_errors();
    }

    // Method called near the end of construction to verify that the
    // newly-created grid source object is valid.
    inline void check_for_errors() {
      // TODO: Check that energy grid values are nondecreasing
      // TODO: Check for at least one nonzero pair.second value
      size_t grid_size = grid.size();
      if (grid_size < 2) throw std::runtime_error(std::string("Grid with")
        + " less than 2 gridpoints passed to the constructor of"
        + " TMarleyGridNeutrinoSource.");
      for (size_t j = 0; j < grid_size; ++j) {
        auto& pair = grid.at(j);
        if (pair.first < 0.) throw std::runtime_error(std::string("All energy")
          + " values used in a TMarleyGridNeutrinoSource"
          + " object must be nonnegative");
        // Prevent actually sampling an energy value of zero by advancing to
        // the next representable double value.
        else if (pair.first == 0.) pair.first = std::nextafter(0.,
          marley_utils::infinity);
        if (pair.second < 0.) throw std::runtime_error(std::string("All PDF")
          + " values used in a TMarleyGridNeutrinoSource"
          + " object must be nonnegative");
      }
    }

    // Returns the maximum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emax() const { return grid.back().first; }

    // Returns the minimum neutrino energy that can be sampled by this source
    // object
    virtual inline double get_Emin() const { return grid.front().first; }

    virtual inline double pdf(double E_nu) {
      return grid.interpolate(E_nu);
    }

  private:
    Grid grid;
};
