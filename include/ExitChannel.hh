#pragma once

#include "Fragment.hh"
#include "Generator.hh"
#include "Level.hh"
#include "Parity.hh"
#include "IteratorToPointerMember.hh"

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

      inline virtual bool is_continuum() const override { return false; }

      inline marley::Level& get_final_level() { return final_level_; }

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

      inline virtual bool emits_fragment() const override { return true; }

      inline virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*gen unused*/)
        override
      {
        Ex = final_level_.get_energy();
        two_J = final_level_.get_two_J();
        Pi = final_level_.get_parity();
        emitted_particle = marley::Particle(fragment_.get_pid(),
          fragment_.get_mass());
        residual_nucleus = residue_;
      }

    protected:
      const marley::Fragment& fragment_; // emitted fragment in this exit channel
  };

  // Gamma emission exit channel that leads to a discrete nuclear level in
  // the final state
  class GammaDiscreteExitChannel : public DiscreteExitChannel {
    public:
      GammaDiscreteExitChannel(double width, marley::Level& flev,
        marley::Particle residue) : DiscreteExitChannel(width, flev, residue) {}

      inline virtual bool emits_fragment() const override { return false; }

      inline virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*gen unused*/)
        override
      {
        Ex = final_level_.get_energy();
        two_J = final_level_.get_two_J();
        Pi = final_level_.get_parity();
        emitted_particle = marley::Particle(marley_utils::PHOTON, 0.);
        residual_nucleus = residue_;
      }
  };


  // Exit channel that leads to the unbound continuum in the final state
  class ContinuumExitChannel : public ExitChannel
  {
    public:
      // Constructor takes a marley::Particle representing the residual
      // nucleus in its ground state. Its mass will be changed to match
      // the sampled final excitation energy each time.
      ContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double)> Epdf, marley::Particle gs_residue)
        : marley::ExitChannel(width), Emin_(Emin), Emax_(Emax),
        Epdf_(Epdf), gs_residue_(gs_residue) {}

      inline virtual bool is_continuum() const override { return true; }

      //virtual void get_post_decay_parameters(double& Ex, int& two_J,
      //  marley::Parity& Pi, marley::Generator& gen) = 0;

      // Create and return a marley::Particle object representing the
      // emitted fragment or gamma ray from this exit channel
      //virtual marley::Particle emitted_particle() = 0;

      // Create and return a marley::Particle object representing the
      // residual nucleus after particle emission
      //virtual marley::Particle residual_particle() = 0;

    protected:
      double Emin_;
      double Emax_;
      std::function<double(double)> Epdf_;
      marley::Particle gs_residue_;
  };

  // Fragment emission exit channel that leads to the unbound continuum in the
  // final state
  class FragmentContinuumExitChannel : public ContinuumExitChannel
  {
    public:
      FragmentContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double)> Epdf, const marley::Fragment& frag,
        marley::Particle gs_residue)
        : marley::ContinuumExitChannel(width, Emin, Emax, Epdf, gs_residue),
        fragment_(frag) {}

      inline virtual bool emits_fragment() const override { return true; }

      inline const marley::Fragment& get_fragment() const { return fragment_; }

      inline virtual void do_decay(double& Ex, int& /*two_J*/,
        marley::Parity& /*Pi*/, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override
      {
        Ex = gen.rejection_sample(Epdf_, Emin_, Emax_);

        // TODO: fix these
        //two_J = final_level_.get_two_J();
        //Pi = final_level_.get_parity();
        emitted_particle = marley::Particle(fragment_.get_pid(),
          fragment_.get_mass());

        residual_nucleus = gs_residue_;
        double rn_mass = gs_residue_.get_mass() + Ex;
        residual_nucleus.set_mass(rn_mass);
      }

    protected:
      const marley::Fragment& fragment_; // fragment emitted in this exit channel
  };

  // Gamma emission exit channel that leads to the unbound continuum in the
  // final state
  class GammaContinuumExitChannel : public ContinuumExitChannel
  {
    public:
      GammaContinuumExitChannel(double width, double Emin, double Emax,
        std::function<double(double)> Epdf, marley::Particle gs_residue)
        : marley::ContinuumExitChannel(width, Emin, Emax, Epdf, gs_residue) {}

      inline virtual bool emits_fragment() const override { return false; }

      inline virtual void do_decay(double& Ex, int& /*two_J*/,
        marley::Parity& /*Pi*/, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override
      {
        Ex = gen.rejection_sample(Epdf_, Emin_, Emax_);
        // TODO: fix these
        //two_J = final_level_.get_two_J();
        //Pi = final_level_.get_parity();
        emitted_particle = marley::Particle(marley_utils::PHOTON, 0.);
        residual_nucleus = gs_residue_;
        double rn_mass = gs_residue_.get_mass() + Ex;
        residual_nucleus.set_mass(rn_mass);
      }
  };
}
