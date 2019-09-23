#pragma once
#include <functional>
#include <set>

#include "marley/marley_utils.hh"
#include "marley/InterpolationGrid.hh"

namespace marley {

  class Generator;

  /// @brief Abstract base class for all objects that describe the incident
  /// neutrino energy distribution
  /// @details Classes derived from NeutrinoSource implement a probability
  /// density function used to sample neutrino energies for each event. Because
  /// MARLEY will weight this probability density by the reaction cross
  /// section(s) during event generation, all NeutrinoSource objects should use
  /// <i>unweighted</i> energy spectra. While normalizing the probability
  /// densities to unity is not strictly required (the cross section weighted
  /// spectra are normalized automatically before sampling), it is encouraged,
  /// and all classes derived from NeutrinoSource currently do so.
  class NeutrinoSource {
    public:

      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      NeutrinoSource(int particle_id);

      virtual ~NeutrinoSource() = default;

      /// @brief Get the maximum neutrino energy (MeV) that can be sampled by
      /// this source
      virtual double get_Emax() const = 0;

      /// @brief Get the minimum neutrino energy (MeV) that can be sampled by
      /// this source
      virtual double get_Emin() const = 0;

      /// @brief Get the PDG particle ID for the neutrino type produced by this
      /// source
      inline virtual int get_pid() const;

      /// @brief Probability density function describing the incident neutrino
      /// energy distribution
      /// @details The neutrino spectrum produced by this source will be folded
      /// with the relevant cross sections by a Generator object during event
      /// creation
      /// @param E neutrino energy (MeV)
      /// @return Probability density (MeV<sup> -1</sup>)
      virtual double pdf(double E) const = 0;

      /// Returns true if the Particle Data Group code passed to the function
      /// is allowed to be used by a neutrino source object, and returns false
      /// otherwise.
      static inline bool pdg_is_allowed(const int pdg);

      /// @brief Samples an incident neutrino energy and loads pdg with
      /// the PDG code of the appropriate neutrino type
      /// @details This function uses a probability density function that
      /// represents the <i>incident</i> neutrino spectrum, i.e., it is
      /// not weighted by the reaction cross section(s). To sample a
      /// neutrino from the <i>reacting</i> (i.e., cross section weighted)
      /// spectrum, use Generator::sample_reaction().
      /// @param[out] pdg PDG code of the sampled neutrino
      /// @param[in] gen Generator to use for random sampling
      /// @return The energy of the sampled neutrino (MeV)
      virtual double sample_incident_neutrino(int& pdg,
        marley::Generator& gen) const;

    protected:

      int pid_; ///< PDG particle ID for the neutrinos produced by this source

    private:
      /// PDG particle IDs for each neutrino that could possibly be produced by
      /// a NeutrinoSource object. Attempting to create a NeutrinoSource object
      /// that produces a particle that does not appear in this set will throw
      /// an exception.
      static const std::set<int> pids_;
  };

  /// @brief Monoenergetic neutrino source
  class MonoNeutrinoSource : public NeutrinoSource {
    public:
      /// @param particle_id neutrino PDG particle ID
      /// @param E neutrino energy (MeV)
      inline MonoNeutrinoSource(int particle_id =
        marley_utils::ELECTRON_NEUTRINO, double E = 10.);

      inline virtual double get_Emax() const override;

      inline virtual double get_Emin() const override;

      inline virtual double pdf(double E) const override;

    protected:
      double energy_; ///< neutrino energy (MeV)
  };

  /// @brief Supernova cooling neutrino source approximated using a Fermi-Dirac
  /// distribution
  /// @details Neutrino energies from this source are sampled from a
  /// Fermi-Dirac distribution with temperature @f$T@f$ (MeV) and pinching
  /// parameter @f$\eta@f$. The probability density function is given by
  /// @f{align*}{ P(E) &= \frac{C\,E^2} {T^4\left[1+\exp\left(\frac{E}{T} -
  /// \eta\right)\right]} & \text{E}_\text{min} \leq E \leq
  /// \text{E}_\text{max} @f} where @f$C@f$ is a normalization constant.
  class FermiDiracNeutrinoSource : public NeutrinoSource {
    public:

      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      /// @param Emin minimum allowed neutrino energy (MeV)
      /// @param Emax maximum allowed neutrino energy (MeV)
      /// @param temp temperature (MeV)
      /// @param eta pinching parameter
      FermiDiracNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO, double Emin = 0.,
        double Emax = 50., double temp = 3.5, double eta = 0.);

      inline virtual double get_Emin() const override;

      inline virtual double get_Emax() const override;

      virtual double pdf(double E) const override;

    protected:

      double Emin_; ///< minimum neutrino energy (MeV)
      double Emax_; ///< maximum neutrino energy (MeV)
      double temperature_;  ///< temperature (MeV)
      double eta_; ///< dimensionless pinching parameter
      double C_; ///< normalization constant (MeV<sup>2</sup>)
  };

  /// @brief "Beta-fit" neutrino source
  /// @details Neutrino energies from this source are sampled from a "beta-fit"
  /// spectrum (see, for example, equation 7 in this <a
  /// href="http://arxiv.org/abs/1511.00806">preprint</a>) with mean
  /// neutrino energy @f$E_\text{mean}@f$ and fit parameter
  /// @f$\beta@f$.  The probability density function is given by @f{align*}{
  /// P(E) &= C\left(E/\,E_\text{mean}\right)^{\beta - 1}
  /// \exp\left(-\beta\,E/\,E_\text{mean}\right) & \text{E}_\text{min} \leq E
  /// \leq \text{E}_\text{max} @f} where @f$C@f$ is a normalization constant.
  class BetaFitNeutrinoSource : public NeutrinoSource {
    public:
      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      /// @param Emin minimum allowed neutrino energy (MeV)
      /// @param Emax maximum allowed neutrino energy (MeV)
      /// @param Emean mean neutrino energy (MeV)
      /// @param beta dimensionless fit parameter
      BetaFitNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO, double Emin = 0.,
        double Emax = 50., double Emean = 13., double beta = 4.5);

      inline virtual double get_Emin() const override;
      inline virtual double get_Emax() const override;

      virtual double pdf(double E) const override;

    protected:
      double Emin_; ///< minimum neutrino energy (MeV)
      double Emax_; ///< maximum neutrino energy (MeV)
      /// @brief mean neutrino energy
      /// @note This is exactly the mean neutrino energy only if @f$
      /// E_\text{min} = 0 @f$ and @f$ Emax = \infty @f$. Truncating the
      /// distribution may alter the mean value appreciably.
      double Emean_;
      double beta_; ///< pinching parameter
      double C_; ///< dimensionless normalization constant
  };

  /// @brief Neutrino source with an arbitrary energy spectrum described by a
  /// std::function<double(double)> object
  class FunctionNeutrinoSource : public NeutrinoSource {
    public:
      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      /// @param Emin minimum allowed neutrino energy (MeV)
      /// @param Emax maximum allowed neutrino energy (MeV)
      /// @param prob_dens_func a std::function<double(double)> object to be
      /// used as a probability density function
      FunctionNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO, double Emin = 0., double Emax = 50.,
        std::function<double(double)> prob_dens_func
        = [](double) -> double { return 1.; });

      inline virtual double get_Emax() const override;

      inline virtual double get_Emin() const override;

      inline virtual double pdf(double E) const override;

    private:
      double Emin_; ///< minimum neutrino energy (MeV)
      double Emax_; ///< maximum neutrino energy (MeV)
      /// @brief user-supplied probability density function
      std::function<double(double)> probability_density_;
  };

  /// @brief Muon decay-at-rest neutrino source
  /// @details Neutrino energies from this source are sampled from the
  /// appropriate Michel spectrum for muon decay-at-rest.  For a @f$\nu_e@f$
  /// source, the probability density function is given by @f{align*}{ P(E) &=
  /// 96E^2m_\mu^{-4}(m_\mu - 2E) & 0 < E < m_\mu/2 @f} where @f$m_\mu@f$ is
  /// the muon mass. For a @f$\bar{\nu}_\mu@f$ source, the probability density
  /// function is given by @f{align*}{ P(E) &= 16E^2m_\mu^{-4}(3m_\mu - 4E) & 0
  /// < E < m_\mu/2. @f}
  class DecayAtRestNeutrinoSource : public NeutrinoSource {
    public:
      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      DecayAtRestNeutrinoSource(int particle_id
        = marley_utils::ELECTRON_NEUTRINO);

      inline virtual double get_Emax() const override;

      inline virtual double get_Emin() const override;

      virtual double pdf(double E) const override;

    private:
      // Muon mass stuff (m_mu^(-4) pre-computed for speed)
      static constexpr double m_mu_ = marley_utils::m_mu
        * marley_utils::micro_amu; // MeV
      static constexpr double m_mu_to_the_minus_four_
        = 1. / (m_mu_ * m_mu_ * m_mu_ * m_mu_);
      static constexpr double Emin_ = 0.; // MeV
      static constexpr double Emax_ = m_mu_ / 2.; // MeV
  };

  /// @brief Neutrino source that uses a tabulated energy spectrum
  class GridNeutrinoSource : public NeutrinoSource {
    public:
      using Grid = InterpolationGrid<double>;
      using Method = Grid::InterpolationMethod;

      /// @param g InterpolationGrid<double> object that
      /// describes this source's probability density function
      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      inline GridNeutrinoSource(const Grid& g, int particle_id
        = marley_utils::ELECTRON_NEUTRINO);

      /// @param Es vector of neutrino energy gridpoints (MeV)
      /// @param PDs vector of probability densities (MeV<sup> -1</sup>)
      /// @param particle_id PDG particle ID for the neutrinos produced by this
      /// source
      /// @param method specifier indicating which interpolation rule should be
      /// used between grid points
      inline GridNeutrinoSource(const std::vector<double>& Es,
        const std::vector<double>& PDs, int particle_id
        = marley_utils::ELECTRON_NEUTRINO, Method method
        = Method::LinearLinear);

      inline virtual double get_Emax() const override;

      inline virtual double get_Emin() const override;

      inline virtual double pdf(double E) const override;

    protected:
      Grid grid_;

    private:
      /// Method called near the end of construction to verify that the
      /// newly-created grid source object is valid.
      void check_for_errors();
  };

  // Inline function definitions
  inline int NeutrinoSource::get_pid() const { return pid_; }
  inline bool NeutrinoSource::pdg_is_allowed(const int pdg)
    { return (pids_.count(pdg) > 0); }

  inline MonoNeutrinoSource::MonoNeutrinoSource(int particle_id, double E)
    : NeutrinoSource(particle_id), energy_(E) {}
  inline double MonoNeutrinoSource::get_Emax() const { return energy_; }
  inline double MonoNeutrinoSource::get_Emin() const { return energy_; }
  inline double MonoNeutrinoSource::pdf(double E) const
    { if (energy_ == E) return 1.; else return 0.; }

  inline double FermiDiracNeutrinoSource::get_Emax() const { return Emax_; }
  inline double FermiDiracNeutrinoSource::get_Emin() const { return Emin_; }

  inline double BetaFitNeutrinoSource::get_Emax() const { return Emax_; }
  inline double BetaFitNeutrinoSource::get_Emin() const { return Emin_; }

  inline double FunctionNeutrinoSource::get_Emax() const { return Emax_; }
  inline double FunctionNeutrinoSource::get_Emin() const { return Emin_; }
  inline double FunctionNeutrinoSource::pdf(double E) const {
    if (E < Emin_ || E > Emax_) return 0.;
    else return probability_density_(E);
  }

  inline double DecayAtRestNeutrinoSource::get_Emax() const { return Emax_; }
  inline double DecayAtRestNeutrinoSource::get_Emin() const { return Emin_; }

  inline double GridNeutrinoSource::get_Emax() const
    { return grid_.back().first; }
  inline double GridNeutrinoSource::get_Emin() const
    { return grid_.front().first; }
  inline double GridNeutrinoSource::pdf(double E) const
    { return grid_.interpolate(E); }

  inline GridNeutrinoSource::GridNeutrinoSource(const Grid& g, int particle_id)
    : NeutrinoSource(particle_id), grid_(g) { check_for_errors(); }

  inline GridNeutrinoSource::GridNeutrinoSource(const std::vector<double>& Es,
    const std::vector<double>& prob_densities, int particle_id, Method method)
    : NeutrinoSource(particle_id), grid_(Es, prob_densities, method)
    { check_for_errors(); }
}
