#pragma once

#include "Fragment.hh"
#include "Generator.hh"
#include "IteratorToPointerMember.hh"
#include "Level.hh"
#include "MassTable.hh"
#include "NuclearPhysics.hh"
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
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
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
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
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
        marley::Particle gs_residue) : marley::ExitChannel(width),
        Emin_(Emin), Emax_(Emax), gs_residue_(gs_residue) {}

      inline virtual bool is_continuum() const override { return true; }

    protected:

      // member struct used for sampling final-state spins and parities
      //struct SpinParityWidth {
      //  int twoJf; // final nuclear spin
      //  marley::Parity Pf; // final nuclear parity
      //  double width; // partial width for this spin-parity combination
      //};

      //virtual void sample_spin_parity(double Ea, int& twoJ,
      //  marley::Parity& P) = 0;

      double Emin_;
      double Emax_;
      marley::Particle gs_residue_;
      //std::vector<SpinParityWidth> jpi_widths_table_;
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

      inline virtual bool emits_fragment() const override { return true; }

      inline const marley::Fragment& get_fragment() const { return fragment_; }

      inline virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override
      {
        double Ea;
        Ex = gen.rejection_sample([&Ea, this](double ex)
          -> double { return this->Epdf_(Ea, ex); }, Emin_, Emax_);

        marley::NuclearPhysics::sample_fragment_spin_parity(two_J, Pi,
          fragment_, gen.get_structure_db().get_optical_model(
          gs_residue_.get_id()), gen, Ex, Ea);

        emitted_particle = marley::Particle(fragment_.get_pid(),
          fragment_.get_mass());

        residual_nucleus = gs_residue_;
        double rn_mass = gs_residue_.get_mass() + Ex;
        residual_nucleus.set_mass(rn_mass);
      }

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

      inline virtual bool emits_fragment() const override { return false; }

      inline virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen)
        override
      {
        double Exi = Ex;
        Ex = gen.rejection_sample(Epdf_, Emin_, Emax_);

        int nuc_pid = residual_nucleus.get_id();
        int Z = marley::MassTable::get_particle_Z(nuc_pid);
        int A = marley::MassTable::get_particle_A(nuc_pid);

        marley::NuclearPhysics::sample_gamma_spin_parity(Z, A, two_J, Pi, Exi,
          Ex, gen);

        emitted_particle = marley::Particle(marley_utils::PHOTON, 0.);
        residual_nucleus = gs_residue_;
        double rn_mass = gs_residue_.get_mass() + Ex;
        residual_nucleus.set_mass(rn_mass);
      }

    protected:
      std::function<double(double)> Epdf_;
  };
}
