#pragma once
#include <vector>

#include "Fragment.hh"
#include "Generator.hh"
#include "IteratorToPointerMember.hh"
#include "Level.hh"
#include "MassTable.hh"
#include "Parity.hh"

namespace marley {

  class ExitChannel {

    public:

      ExitChannel(double width) : width_(width) {}

      // Whether this decay channel accesses the unbound continuum
      // (true) or a discrete level (false)
      virtual bool is_continuum() const = 0;

      // Whether or not this decay channel involves fragment emission (true)
      // or gamma-ray emission (false)
      virtual bool emits_fragment() const = 0;

      // Determine final values for several quantities of interest for the
      // residual nuclide after the decay occurs.  Input values are the initial
      // state's excitation energy Ex, 2*(total angular momentum) two_J, parity
      // Pi, and two marley::Particles. These variables are then loaded with
      // the corresponding values for the final state, with the first particle
      // representing the emitted particle and the second representing the
      // residual nucleus.
      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen) = 0;

      // Convert an iterator that points to this marley::ExitChannel object
      // into an iterator to the exit channel's decay width member variable.
      // This is used to load a std::discrete_distribution with decay widths
      // for sampling without redundant storage.
      template<typename It> static inline
        marley::IteratorToPointerMember<It, marley::ExitChannel, double>
        make_width_iterator(It it)
      {
        return marley::IteratorToPointerMember<It, marley::ExitChannel,
          double>(it, &marley::ExitChannel::width_);
      }

      inline double get_width() const { return width_; }

    protected:
      // Decay width to this exit channel
      double width_;
  };

  // Exit channel that leads to a discrete nuclear level in the final state
  class DiscreteExitChannel : public ExitChannel {
    public:
      DiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue) : ExitChannel(width), final_level_(flev),
        residue_(residue) {}

      inline virtual bool is_continuum() const final override { return false; }

      inline marley::Level& get_final_level() { return final_level_; }
      inline const marley::Level& get_final_level() const { return final_level_; }

    protected:
      marley::Level& final_level_; // Pointer to final level in a nuclear decay scheme
      marley::Particle residue_; // particle object representing the residual nucleus
      //int l;      // Orbital angular momentum of the outgoing fragment
  };

  // Fragment emission exit channel that leads to a discrete nuclear level in
  // the final state
  class FragmentDiscreteExitChannel : public DiscreteExitChannel {
    public:
      FragmentDiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue, const marley::Fragment& frag)
        : DiscreteExitChannel(width, flev, residue), fragment_(frag) {}

      inline virtual bool emits_fragment() const final override { return true; }

      inline const marley::Fragment& get_fragment() const { return fragment_; }

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
        override;

    protected:
      const marley::Fragment& fragment_; // emitted fragment in this exit channel
  };

  // Gamma emission exit channel that leads to a discrete nuclear level in
  // the final state
  class GammaDiscreteExitChannel : public DiscreteExitChannel {
    public:
      GammaDiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue) : DiscreteExitChannel(width, flev, residue) {}

      inline virtual bool emits_fragment() const final override
        { return false; }

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
        override;
  };


  // Exit channel that leads to the unbound continuum in the final state
  class ContinuumExitChannel : public ExitChannel
  {
    public:

      // Constructor takes a marley::Particle representing the residual
      // nucleus in its ground state. Its mass will be changed to match
      // the sampled final excitation energy each time.
      ContinuumExitChannel(double width, double Emin, double Emax,
        marley::Particle gs_residue) : marley::ExitChannel(width),
        Emin_(Emin), Emax_(Emax), gs_residue_(gs_residue) {}

      inline virtual bool is_continuum() const final override { return true; }

    protected:

      // member struct used for sampling final-state spins and parities
      struct SpinParityWidth {
        SpinParityWidth(int twoJ, marley::Parity p, double w)
          : twoJf(twoJ), Pf(p), width(w) {}
        int twoJf; // final nuclear spin
        marley::Parity Pf; // final nuclear parity
        double width; // partial width for this spin-parity combination
      };

      //virtual void sample_spin_parity(double Ea, int& twoJ,
      //  marley::Parity& P) = 0;

      double Emin_;
      double Emax_;
      marley::Particle gs_residue_;
      std::vector<SpinParityWidth> jpi_widths_table_;
  };

  // Fragment emission exit channel that leads to the unbound continuum in the
  // final state
  class FragmentContinuumExitChannel : public ContinuumExitChannel
  {
    public:
      FragmentContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double&, double)> Epdf,
        const marley::Fragment& frag, marley::Particle gs_residue)
        : marley::ContinuumExitChannel(width, Emin, Emax, gs_residue),
        Epdf_(Epdf), fragment_(frag) {}

      inline virtual bool emits_fragment() const final override { return true; }

      inline const marley::Fragment& get_fragment() const { return fragment_; }

      void sample_spin_parity(int& twoJ, marley::Parity& Pi,
        marley::Generator& gen, double Exf, double Ea);

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override;

    protected:
      std::function<double(double&, double)> Epdf_;
      const marley::Fragment& fragment_; // fragment emitted in this exit channel
  };

  // Gamma emission exit channel that leads to the unbound continuum in the
  // final state
  class GammaContinuumExitChannel : public ContinuumExitChannel
  {
    public:
      GammaContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double)> Epdf, marley::Particle gs_residue)
        : marley::ContinuumExitChannel(width, Emin, Emax, gs_residue),
        Epdf_(Epdf) {}

      inline virtual bool emits_fragment() const final override
        { return false; }

      void sample_spin_parity(int Z, int A, int& twoJ, marley::Parity& Pi,
        double Exi, double Exf, marley::Generator& gen);

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override;

    protected:
      std::function<double(double)> Epdf_;

      // Helper function for building the table of gamma-ray spin-parities and
      // widths
      double store_gamma_jpi_width(double Exf, int twoJf, marley::Parity Pi,
        double tcE, double tcM, int mpol, marley::LevelDensityModel& ldm);
  };
}
