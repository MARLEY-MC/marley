#pragma once
#include <vector>

#include "marley/Fragment.hh"
#include "marley/Generator.hh"
#include "marley/IteratorToPointerMember.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/Parity.hh"

namespace marley {

  /// @brief Abstract base class for compound nucleus de-excitation channels
  class ExitChannel {

    public:

      /// @param width Partial decay width (MeV)
      ExitChannel(double width) : width_(width) {}

      virtual ~ExitChannel() = default;

      /// @brief Returns true if this channel accesses the particle-unbound
      /// continuum of nuclear levels or false otherwise
      virtual bool is_continuum() const = 0;

      /// @brief Returns true if this decay channel involves fragment emission
      /// or false if it involves gamma-ray emission.
      virtual bool emits_fragment() const = 0;

      /// @brief Simulates a nuclear decay into this channel
      /// @details The excitation energy, spin, and parity values are loaded
      /// with their final values as this function returns.
      /// @param[in,out] Ex The nuclear excitation energy
      /// @param[in,out] two_J Two times the nuclear spin
      /// @param[in,out] Pi The nuclear parity
      /// @param[out] emitted_particle Particle emitted in the de-excitation
      /// @param[out] residual_nucleus Final-state nucleus after particle
      /// emission
      /// @param gen Generator to use for random sampling
      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen) = 0;

      /// @brief Convert an iterator that points to an ExitChannel object into
      /// an iterator to the ExitChannel's width_ member variable.
      /// @details This is used to load a std::discrete_distribution with decay
      /// widths for sampling without redundant storage.
      template<typename It> static inline
        marley::IteratorToPointerMember<It, marley::ExitChannel, double>
        make_width_iterator(It it)
      {
        return marley::IteratorToPointerMember<It, marley::ExitChannel,
          double>(it, &marley::ExitChannel::width_);
      }

      /// @brief Get the partial decay width to this channel
      inline double width() const { return width_; }

    protected:

      /// @brief Partial decay width (MeV)
      double width_;
  };

  /// @brief Abstract base class for ExitChannel objects that lead to
  /// discrete nuclear levels in the final state
  class DiscreteExitChannel : public ExitChannel {
    public:

      /// @param width Partial decay width (MeV)
      /// @param[in] flev Reference to the final-state nuclear level
      /// @param residue Particle object to use as the final-state nucleus
      DiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue) : ExitChannel(width), final_level_(flev),
        residue_(residue) {}

      inline virtual bool is_continuum() const final override { return false; }

      /// @brief Get a const reference to the final-state nuclear level
      inline const marley::Level& get_final_level() const
        { return final_level_; }
      /// @brief Get a non-const reference to the final-state nuclear level
      inline marley::Level& get_final_level() { return final_level_; }

    protected:
      /// @brief Reference to the final-state nuclear level
      marley::Level& final_level_;
      /// @brief Residual nucleus Particle object
      marley::Particle residue_;
  };

  /// @brief %Fragment emission ExitChannel that leads to a discrete nuclear
  /// level in the final state
  class FragmentDiscreteExitChannel : public DiscreteExitChannel {

    public:

      /// @param width Partial decay width (MeV)
      /// @param[in] flev Reference to the final-state nuclear level
      /// @param residue Particle object to use as the final-state nucleus
      /// @param frag Reference to the emitted Fragment
      FragmentDiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue, const marley::Fragment& frag)
        : DiscreteExitChannel(width, flev, residue), fragment_(frag) {}

      inline virtual bool emits_fragment() const final override { return true; }

      /// @brief Get a reference to the emitted Fragment
      inline const marley::Fragment& get_fragment() const { return fragment_; }

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
        override;

    protected:

      /// @brief Fragment emitted by this exit channel
      const marley::Fragment& fragment_;
  };

  /// @brief %Gamma emission exit channel that leads to a discrete nuclear
  /// level in the final state
  class GammaDiscreteExitChannel : public DiscreteExitChannel {

    public:

      /// @param width Partial decay width (MeV)
      /// @param[in] flev Reference to the final-state nuclear level
      /// @param residue Particle object to use as the final-state nucleus
      GammaDiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue) : DiscreteExitChannel(width, flev, residue) {}

      inline virtual bool emits_fragment() const final override
        { return false; }

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
        override;
  };


  /// @brief Abstract base class for ExitChannel objects that lead to the
  /// unbound continuum in the final state
  class ContinuumExitChannel : public ExitChannel
  {
    public:

      /// @param width Partial decay width (MeV)
      /// @param Emin Minimum accessible nuclear excitation energy (MeV)
      /// @param Emax Minimum accessible nuclear excitation energy (MeV)
      /// @param gs_residue Residual nucleus Particle object whose mass
      /// corresponds to the ground state
      /// @note Calls to do_decay() for this channel will use the mass of the
      /// gs_residue Particle and the sampled excitation energy to determine
      /// the mass of the final-state nucleus
      ContinuumExitChannel(double width, double Emin, double Emax,
        marley::Particle gs_residue) : marley::ExitChannel(width),
        Emin_(Emin), Emax_(Emax), gs_residue_(gs_residue) {}

      inline virtual bool is_continuum() const final override { return true; }

    protected:

      /// @brief A spin-parity value with its corresponding partial decay width
      /// @details This struct is used to sample final-state nuclear
      /// spin-parities in classes derived from ContinuumExitChannel
      struct SpinParityWidth {

        /// @param twoJ Two times the nuclear spin
        /// @param p Nuclear parity
        /// @param w Partial decay width (MeV) for the given spin-parity
        SpinParityWidth(int twoJ, marley::Parity p, double w)
          : twoJf(twoJ), Pf(p), width(w) {}

        int twoJf; ///< Final nuclear spin
        marley::Parity Pf; ///< Final nuclear parity
        double width; ///< Partial decay width (MeV)
      };

      double Emin_; ///< Minimum accessible nuclear excitation energy (MeV)
      double Emax_; ///< Maximum accessible nuclear excitation energy (MeV)
      marley::Particle gs_residue_; ///< Ground-state residual nucleus

      /// @brief Table of possible final-state spin-parities together
      /// with their partial decay widths
      std::vector<SpinParityWidth> jpi_widths_table_;
  };

  /// @brief %Fragment emission ExitChannel that leads to the unbound continuum
  /// in the final state
  class FragmentContinuumExitChannel : public ContinuumExitChannel
  {
    public:

      /// @param width Partial decay width (MeV)
      /// @param Emin Minimum accessible nuclear excitation energy (MeV)
      /// @param Emax Minimum accessible nuclear excitation energy (MeV)
      /// @param Epdf std::function describing the final-state nuclear
      /// excitation energy distribution
      /// @param frag Reference to the emitted Fragment
      /// @param gs_residue Residual nucleus Particle object whose mass
      /// corresponds to the ground state
      /// @note Calls to do_decay() for this channel will use the mass of the
      /// gs_residue Particle and the sampled excitation energy to determine
      /// the mass of the final-state nucleus
      FragmentContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double&, double, bool)> Epdf,
        const marley::Fragment& frag, marley::Particle gs_residue)
        : marley::ContinuumExitChannel(width, Emin, Emax, gs_residue),
        Epdf_(Epdf), fragment_(frag) {}

      inline virtual bool emits_fragment() const final override { return true; }

      /// @brief Get a reference to the emitted Fragment
      inline const marley::Fragment& get_fragment() const { return fragment_; }

      /// @brief Sample a final-state spin and parity for the residual nucleus
      /// @param[out] twoJ Two times the final-state nuclear spin
      /// @param[out] Pi The final-state nuclear parity
      /// @param gen Reference to the Generator object to use for random
      /// sampling
      /// @param Exf The final-state nuclear excitation energy
      /// @param Ea The final-state fragment kinetic energy
      void sample_spin_parity(int& twoJ, marley::Parity& Pi,
        marley::Generator& gen, double Exf, double Ea);

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override;

    protected:
      /// @brief Probability density function describing the distribution
      /// of final-state nuclear excitation energies within the continuum.
      /// @details The first argument is a reference to a double that will
      /// be loaded with the final-state fragment kinetic energy (MeV).
      /// The second is the final-state nuclear excitation energy (MeV).
      /// The return value is a probability density (MeV<sup> -1</sup>)
      /// for sampling the final-state nuclear excitation energy.
      std::function<double(double&, double, bool)> Epdf_;
      const marley::Fragment& fragment_; ///< Emitted Fragment
  };

  /// @brief %Gamma emission exit channel that leads to the unbound continuum
  /// in the final state
  class GammaContinuumExitChannel : public ContinuumExitChannel
  {
    public:

      /// @param width Partial decay width (MeV)
      /// @param Emin Minimum accessible nuclear excitation energy (MeV)
      /// @param Emax Minimum accessible nuclear excitation energy (MeV)
      /// @param Epdf std::function describing the final-state nuclear
      /// excitation energy distribution
      /// @param gs_residue Residual nucleus Particle object whose mass
      /// corresponds to the ground state
      /// @note Calls to do_decay() for this channel will use the mass of the
      /// gs_residue Particle and the sampled excitation energy to determine
      /// the mass of the final-state nucleus
      GammaContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double, bool)> Epdf, marley::Particle gs_residue)
        : marley::ContinuumExitChannel(width, Emin, Emax, gs_residue),
        Epdf_(Epdf) {}

      inline virtual bool emits_fragment() const final override
        { return false; }

      /// @brief Sample a final-state spin and parity for the residual nucleus
      /// @param Z Atomic number
      /// @param A Mass number
      /// @param[out] twoJ Two times the final-state nuclear spin
      /// @param[out] Pi The final-state nuclear parity
      /// @param Exi The initial nuclear excitation energy
      /// @param Exf The final nuclear excitation energy
      /// @param gen Reference to the Generator object to use for random
      /// sampling
      void sample_spin_parity(int Z, int A, int& twoJ, marley::Parity& Pi,
        double Exi, double Exf, marley::Generator& gen);

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override;

    protected:

      /// @brief Probability density function describing the distribution
      /// of final-state nuclear excitation energies within the continuum.
      /// @details The argument is the final-state nuclear excitation energy
      /// (MeV). The return value is a probability density (MeV<sup> -1</sup>)
      /// for sampling the final-state nuclear excitation energy.
      std::function<double(double, bool)> Epdf_;

    private:

      /// @brief Helper function for building the table of final spin-parities
      /// and widths
      /// @param Exf Final nuclear excitation energy
      /// @param twoJf Two times the final nuclear spin
      /// @param Pi Final nuclear parity
      /// @param tcE Electric &gamma;-ray transition transmission coefficient
      /// @param tcM Magnetic &gamma;-ray transition transmission coefficient
      /// @param mpol Multipolarity of the transition
      /// @param ldm Reference to a LevelDensityModel object representing
      /// the density of nuclear levels in the continuum
      double store_gamma_jpi_width(double Exf, int twoJf, marley::Parity Pi,
        double tcE, double tcM, int mpol, marley::LevelDensityModel& ldm);
  };
}
