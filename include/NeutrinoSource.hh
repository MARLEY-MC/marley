#pragma once
#include <functional>
#include <set>

#include "marley_utils.hh"
#include "Error.hh"
#include "InterpolationGrid.hh"
#include "MassTable.hh"

namespace marley {

  class Generator;

  class NeutrinoSource {
    public:

      inline NeutrinoSource(int particle_id) {
        if (!pid_is_allowed(particle_id)) throw marley::Error(
          std::string("Creating a neutrino source object that produces")
          + " particles with PDG ID number " + std::to_string(particle_id)
          + " is not allowed.");
        else pid = particle_id;
      }

      // Returns true if the particle ID passed to the function is allowed to be
      // used by a neutrino source object, and returns false otherwise.
      static inline bool pid_is_allowed(const int particle_id) {
        return (pids.count(particle_id) > 0);
      }

      // Samples a neutrino energy in MeV
      //virtual double sample_energy(marley::Generator& gen) = 0;

      // Returns the maximum neutrino energy that can be sampled by this source
      // object
      virtual double get_Emax() const = 0;

      // Returns the minimum neutrino energy that can be sampled by this source
      // object
      virtual double get_Emin() const = 0;

      // Returns the PDG particle ID for the neutrino type produced
      // by this source
      virtual inline int get_pid() const { return pid; }

      // Probability density function (not necessarily normalized) used by
      // marley::Generator for folding the neutrino spectrum produced by this
      // source with the relevant cross sections.
      virtual double pdf(double E_nu) = 0;

    protected:
      // Particle ID for the neutrino type produced by this source
      int pid;

    private:
      // Particle IDs of each neutrino that could possibly be produced by a
      // source object
      static const std::set<int> pids;
  };

  // Monoenergetic neutrino source
  class MonoNeutrinoSource : public NeutrinoSource {
    public:
      inline MonoNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO, double E = 10./* MeV*/)
      : NeutrinoSource(particle_id)
      {
        energy = E;
      }

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
  // distribution (see, for example, Giunti & Kim equation 15.18)
  class FermiDiracNeutrinoSource : public NeutrinoSource {
    public:
      inline FermiDiracNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO, double Emin = 0.,
        double Emax = 100., double temp = 3.5, double e_t_a = 0.)
        : NeutrinoSource(particle_id)
      {
        E_min = Emin;
        E_max = Emax;
        temperature = temp;
        eta = e_t_a;
        C = 1.;

        // Create a call wrapper to allow us to numerically integrate the
        // PDF for this source.
        std::function<double(double)> fd_dist = std::bind(
          &FermiDiracNeutrinoSource::pdf, this, std::placeholders::_1);
        // Normalize the source spectrum (not strictly necessary, but having the
        // spectrum approximately normalized makes the default rejection sampling
        // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
        // TODO: consider removing the hard-coded value here
        double integral = marley_utils::num_integrate(fd_dist, E_min, E_max, 1e4);

        // Update the normalization constant, thereby normalizing this object's
        // pdf in the process.
        C /= integral;
      }

      // Returns the maximum neutrino energy that can be sampled by this source
      // object
      virtual inline double get_Emax() const { return E_max; }

      // Returns the minimum neutrino energy that can be sampled by this source
      // object
      virtual inline double get_Emin() const { return E_min; }

      //virtual double sample_energy(marley::Generator& gen);

      virtual inline double pdf(double E_nu) {
        if (E_nu < E_min || E_nu > E_max) return 0.;
        else return (C / std::pow(temperature, 4)) * (std::pow(E_nu, 2)
          / (1 + std::exp((E_nu / temperature) - eta)));
      }

    private:
      // Temperature for Fermi-Dirac distribution (in MeV)
      double temperature;
      // Dimensionless pinching parameter for Fermi-Dirac distribution
      double eta;
      // Minimum and maximum neutrino energies produced by this source
      double E_min, E_max;
      // Normalization constant (determined during construction)
      double C;
  };

  // Neutrino source with a "beta fit" spectrum (see, for example, equation 7 in
  // arXiv:1511.00806v4). Equation 15.19 in Giunti & Kim also describes
  // this spectrum, but note that our definition of beta is theirs plus
  // one.
  class BetaFitNeutrinoSource : public NeutrinoSource {
    public:
      inline BetaFitNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO,
        double Emin = 0., double Emax = 100., double Emean = 13.,
        double b_e_t_a = 4.5)
        : NeutrinoSource(particle_id)
      {
        E_min = Emin;
        E_max = Emax;
        E_mean = Emean;
        beta = b_e_t_a;
        C = 1.; // Normalization constant (will be updated momentarily)

        // Create a call wrapper to allow us to numerically integrate the
        // PDF for this source.
        std::function<double(double)> spect = std::bind(
          &BetaFitNeutrinoSource::pdf, this, std::placeholders::_1);
        // Normalize the source spectrum (not strictly necessary, but having the
        // spectrum approximately normalized makes the default rejection sampling
        // tolerance of 1e-8 reliable for finding the maximum of the spectrum)
        // TODO: consider removing the hard-coded value here
        double integral = marley_utils::num_integrate(spect, E_min, E_max, 1e4);

        // Update the normalization constant, thereby normalizing this object's
        // pdf in the process.
        C /= integral;
      }

      // Returns the maximum neutrino energy that can be sampled by this source
      // object
      virtual inline double get_Emax() const { return E_max; }

      // Returns the minimum neutrino energy that can be sampled by this source
      // object
      virtual inline double get_Emin() const { return E_min; }

      //virtual double sample_energy(marley::Generator& gen);

      virtual inline double pdf(double E_nu) {
        if (E_nu < E_min || E_nu > E_max) return 0.;
        else return C * std::pow(E_nu / E_mean, beta - 1)
          * std::exp(-beta * E_nu / E_mean);
      }

    private:
      // Mean energy for beta fit distribution (assuming E_min = 0.
      // and E_max = infinity. Truncating the distribution may alter
      // the mean value appreciably)
      double E_mean;
      // Pinching parameter beta
      double beta;
      // Minimum and maximum neutrino energies produced by this source
      double E_min, E_max;
      // Normalization constant (determined during construction)
      double C;
  };

  // Neutrino source with an arbitrary energy spectrum supplied during
  // construction as a std::function<double(double)> object
  class FunctionNeutrinoSource : public NeutrinoSource {
    public:
      inline FunctionNeutrinoSource(const std::function<double(double)>&
        prob_dens_func = [](double E) -> double { (void)(E); return 1; },
        int particle_id = marley_utils::ELECTRON_NEUTRINO,
        double Emin = 0., double Emax = 100.)
        : NeutrinoSource(particle_id)
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

  // Neutrino source with an energy spectrum given by approximate pion-muon
  // decay at rest spectra (see, for example, equation 2 in
  // http://iopscience.iop.org/1742-6596/574/1/012167)
  class DecayAtRestNeutrinoSource : public NeutrinoSource {
    public:
      inline DecayAtRestNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO) : NeutrinoSource(particle_id)
      {

        if (particle_id != marley_utils::ELECTRON_NEUTRINO &&
          particle_id != marley_utils::MUON_ANTINEUTRINO)
        {
          throw marley::Error(std::string("Decay at rest")
            + " neutrino source objects may only produce electron neutrinos"
            + " or muon antineutrinos. PDG ID number "
            + std::to_string(particle_id) + " is therefore not allowed.");
        }
      }

      // Returns the maximum neutrino energy that can be sampled by this source
      // object
      virtual inline double get_Emax() const { return E_max; }

      // Returns the minimum neutrino energy that can be sampled by this source
      // object
      virtual inline double get_Emin() const { return E_min; }

      virtual inline double pdf(double E_nu) {
        if (E_nu < E_min || E_nu > E_max) return 0.;
        /// @todo Refine the decay-at-rest Michel spectra to use
        /// more exact expressions.
        // Note that both of these source spectra are normalized to 1
        // on the energy interval [0., m_mu / 2.]
        else if (pid == marley_utils::ELECTRON_NEUTRINO)
          return 96. * std::pow(E_nu, 2) * m_mu_to_the_minus_four
            * (m_mu - 2*E_nu);
        // Spectrum for muon antineutrinos
        else return 16. * std::pow(E_nu, 2) * m_mu_to_the_minus_four
          * (3*m_mu - 4*E_nu);
      }

    private:
      // Muon mass stuff (m_mu^(-4) pre-computed for speed)
      static constexpr double m_mu = marley_utils::m_mu * marley_utils::micro_amu; // MeV
      static constexpr double m_mu_to_the_minus_four = std::pow(m_mu, -4); // MeV^(-4)
      // Minimum and maximum neutrino energies produced by this source
      static constexpr double E_min = 0.; // MeV
      static constexpr double E_max = m_mu / 2.; // MeV
  };

  class GridNeutrinoSource : public NeutrinoSource {
    public:
      using Grid = InterpolationGrid<double>;
      using InterpolationMethod = Grid::InterpolationMethod;

      inline GridNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO,
        InterpolationMethod interp_method = InterpolationMethod::LinearLinear)
        : NeutrinoSource(particle_id), grid(interp_method)
      {
        //check_for_errors();
      }

      inline GridNeutrinoSource(const Grid& g,
        int particle_id = marley_utils::ELECTRON_NEUTRINO)
        : NeutrinoSource(particle_id), grid(g)
      {
        check_for_errors();
      }

      inline GridNeutrinoSource(const std::vector<double>& Es,
        const std::vector<double>& prob_densities, int particle_id
        = marley_utils::ELECTRON_NEUTRINO, InterpolationMethod
        interp_method = InterpolationMethod::LinearLinear)
        : NeutrinoSource(particle_id), grid(Es, prob_densities,
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
        if (grid_size < 2) throw marley::Error(std::string("Grid with")
          + " less than 2 gridpoints passed to the constructor of"
          + " marley::GridNeutrinoSource.");
        for (size_t j = 0; j < grid_size; ++j) {
          auto& pair = grid.at(j);
          if (pair.first < 0.) throw marley::Error(std::string("All energy")
            + " values used in a marley::GridNeutrinoSource"
            + " object must be nonnegative");
          // Prevent actually sampling an energy value of zero by advancing to
          // the next representable double value.
          else if (pair.first == 0.) pair.first = std::nextafter(0.,
            marley_utils::infinity);
          if (pair.second < 0.) throw marley::Error(std::string("All PDF")
            + " values used in a marley::GridNeutrinoSource"
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

}
