/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see ${MARLEY}/COPYING or
// visit http://opensource.org/licenses/GPL-3.0

#pragma once

// standard library includes
#include <memory>
#include <vector>

// MARLEY includes
#include "marley/ChebyshevInterpolatingFunction.hh"
#include "marley/Fragment.hh"
#include "marley/Generator.hh"
#include "marley/IteratorToPointerMember.hh"
#include "marley/Level.hh"
#include "marley/MassTable.hh"
#include "marley/Parity.hh"
#include "marley/marley_utils.hh"

namespace marley {

  /// @brief Abstract base class for compound nucleus de-excitation channels
  class ExitChannel {

    public:

      /// @param pdgi PDG code for the initial nucleus
      /// @param qi Net charge (in units of the elementary charge) of the
      /// initial atom or ion
      /// @param Exi Initial nuclear excitation energy @f$ E_x @f$ (MeV)
      /// @param twoJi Two times the initial nuclear spin @f$ J @f$
      /// @param Pi Initial nuclear parity @f$ \Pi @f$
      /// @param rho_i Initial nuclear level density @f$ \rho(E_x, J, \Pi) @f$
      /// in the vicinity of the initial level (MeV<sup> -1</sup>). This
      /// value will be used to compute an overall normalization factor
      /// for decay widths in this exit channel.
      /// @param sdb Reference to the StructureDatabase that will be used
      /// in decay width calculations by this ExitChannel object
      ExitChannel(int pdgi, int qi, double Exi, int twoJi, marley::Parity Pi,
        double rho_i, marley::StructureDatabase& sdb ) : pdgi_( pdgi ),
        qi_( qi ), Exi_( Exi ), twoJi_( twoJi ), Pi_( Pi ), sdb_( sdb )
      {
        one_over_two_pi_rho_i_ = std::pow( 2. * marley_utils::pi * rho_i, -1 );
      }

      virtual ~ExitChannel() = default;

      /// @brief Returns true if this channel accesses the particle-unbound
      /// continuum of nuclear levels or false otherwise
      virtual bool is_continuum() const = 0;

      /// @brief Returns true if this channel involves fragment emission
      /// or false if it involves gamma-ray emission.
      virtual bool emits_fragment() const = 0;

      /// @brief Simulates a nuclear decay into this channel
      /// @param[out] Exf The final nuclear excitation energy
      /// @param[out] two_Jf Two times the final nuclear spin
      /// @param[out] Pf The final nuclear parity
      /// @param[out] emitted_particle Particle emitted in the de-excitation
      /// @param[out] residual_nucleus Final-state nucleus after particle
      /// emission
      /// @param gen Generator to use for random sampling
      virtual void do_decay(double& Exf, int& two_Jf,
        marley::Parity& Pf, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& gen) const = 0;

      /// @brief Convert an iterator that points to an ExitChannel object into
      /// an iterator to the ExitChannel's width_ member variable.
      /// @details This is used to load a std::discrete_distribution with decay
      /// widths for sampling without redundant storage.
      template<typename It> static inline
        marley::IteratorToPointerMember<It, double> make_width_iterator(It it)
      {
        return marley::IteratorToPointerMember<It,
          double>(it, &marley::ExitChannel::width_);
      }

      /// @brief Get the total decay width into this channel (MeV)
      inline double width() const { return width_; }

      /// @brief Returns the PDG code for the particle (gamma-ray or nuclear
      /// fragment) emitted by decays into this ExitChannel
      virtual int emitted_particle_pdg() const = 0;

      /// @brief Returns the PDG code for the final nucleus
      virtual int final_nucleus_pdg() const = 0;

    protected:

      /// Helper function that initializes the width_ member variable upon
      /// construction
      virtual void compute_total_width() = 0;

      /// @brief Helper function that prepares Particle objects representing
      /// the products of the two-body decay
      /// @param[out] emitted_particle Particle emitted in the binary decay
      /// @param[out] residual_nucleus Particle representing the daughter
      /// nucleus
      /// @param[in] Exf Excitation energy of the daughter nucleus
      virtual void prepare_products(marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, double Exf) const;

      /// PDG code for the initial nucleus
      int pdgi_;

      /// Net charge (in units of the elementary charge) of the initial atom
      /// or ion
      int qi_;

      /// Initial nuclear excitation energy @f$ E_x @f$ (MeV)
      double Exi_;

      /// Two times the initial nuclear spin @f$ J @f$
      int twoJi_;

      /// Initial nuclear parity @f$ \Pi @f$
      marley::Parity Pi_;

      /// Normalization factor needed so that the decay widths have correct units
      /// @details This factor is equal to @f$ \left[ 2\pi\rho(E_x, J, \Pi)
      /// \right]^{-1} @f$ where @f$\rho@f$ is the initial nuclear level density
      /// (MeV<sup> -1</sup>) in the vicinity of the initial nuclear level. This
      /// level has excitation energy @f$ E_x @f$, total spin @f$ J @f$, and
      /// parity @f$ \Pi @f$.
      double one_over_two_pi_rho_i_;

      /// Total decay width into this channel (MeV)
      double width_;

      /// StructureDatabase to use in calculations
      marley::StructureDatabase& sdb_;
  };

  /// @brief Abstract base class for ExitChannel objects that lead to
  /// discrete nuclear levels in the final state
  class DiscreteExitChannel : virtual public ExitChannel {
    public:

      /// @param[in] flev Reference to the final-state nuclear level
      DiscreteExitChannel( const marley::Level& flev ) : final_level_( flev ) {}

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
        const final override;

      inline virtual bool is_continuum() const final override { return false; }

      /// @brief Get a const reference to the final-state nuclear level
      inline const marley::Level& get_final_level() const
        { return final_level_; }

    protected:

      /// @brief Reference to the final-state nuclear level
      const marley::Level& final_level_;
  };

  /// @brief Abstract base class for ExitChannel objects that represent
  /// emission of a nuclear fragment
  class FragmentExitChannel : virtual public ExitChannel {
    public:

      /// @param fragment Fragment emitted in this exit channel
      FragmentExitChannel(const marley::Fragment& fragment)
        : fragment_pdg_( fragment.get_pid() ) {}

      virtual int emitted_particle_pdg() const final override
        { return fragment_pdg_; }

      virtual bool emits_fragment() const final override
        { return true; }

      virtual int final_nucleus_pdg() const final override;

    protected:

      /// @brief Helper function that returns that maximum possible excitation
      /// energy for the daughter nucleus after emission of the fragment
      double max_Exf() const;

      /// @brief PDG code identifying the emitted fragment
      int fragment_pdg_;
  };

  /// @brief Abstract base class for ExitChannel objects that represent
  /// emission of a gamma-ray
  class GammaExitChannel : virtual public ExitChannel {
    public:

      GammaExitChannel() {}

      virtual int emitted_particle_pdg() const final override
        { return marley_utils::PHOTON; }

      virtual bool emits_fragment() const final override
        { return false; }

      inline virtual int final_nucleus_pdg() const final override
        { return pdgi_; }

    protected:

      // Returns the gamma-ray energy corresponding to a particular
      // final nuclear excitation energy
      // @param Exf Final nuclear excitation energy (MeV)
      // @return Energy of the gamma-ray emitted in this exit channel (MeV)
      double gamma_energy( double Exf ) const;

      marley::GammaStrengthFunctionModel::TransitionType get_transition_type(
        int mpol, marley::Parity Pf ) const;
  };

  /// @brief Abstract base class for ExitChannel objects that lead to the
  /// unbound continuum in the final state
  class ContinuumExitChannel : virtual public ExitChannel
  {
    public:

      /// @param Ec_min Minimum accessible nuclear excitation energy in the
      /// continuum. Below this value, only discrete nuclear levels are
      /// assumed to be present.
      ContinuumExitChannel(double Ec_min) : E_c_min_( Ec_min ) {}

      /// Helper function that initializes the width_ member variable upon
      /// construction
      virtual void compute_total_width() final override;

      virtual void do_decay(double& Ex, int& two_J,
        marley::Parity& Pi, marley::Particle& emitted_particle,
        marley::Particle& residual_nucleus, marley::Generator& /*unused*/)
        const final override;

      virtual double differential_width( double Exf,
        bool store_jpi_widths = false ) const = 0;

      inline virtual bool is_continuum() const final override { return true; }

      /// @brief Sets the flag that will skip sampling of a final-state
      /// nuclear spin-parity value in do_decay()
      /// @details The skipping functionality should only be used for testing
      /// purposes!
      inline void set_skip_jpi_sampling(bool skip_it) const
        { skip_jpi_sampling_ = skip_it; }

      /// @brief A spin-parity value with its corresponding partial decay width
      /// @details This struct is used to sample final-state nuclear
      /// spin-parities in classes derived from ContinuumExitChannel
      struct SpinParityWidth {

        /// @param twoJ Two times the nuclear spin
        /// @param p Nuclear parity
        /// @param w Partial decay width (MeV) for the given spin-parity
        SpinParityWidth(int twoJ, marley::Parity p, double w)
          : twoJf(twoJ), Pf(p), diff_width(w) {}

        int twoJf; ///< Final nuclear spin
        marley::Parity Pf; ///< Final nuclear parity
        double diff_width; ///< Partial differential decay width (MeV)
      };

      /// Minimum accessible nuclear excitation energy (MeV) in the continuum
      double E_c_min_;

      /// @brief Table of possible final-state spin-parities together
      /// with their partial differential decay widths
      mutable std::vector<SpinParityWidth> jpi_widths_table_;

      /// @brief Flag that allows skipping the sampling of a final
      /// nuclear spin-parity (useful only for testing purposes)
      mutable bool skip_jpi_sampling_ = false;

      /// @brief Chebyshev polynomial interpolant to the cumulative
      /// density function for the final-state nuclear excitation energy
      /// @details This pointer will be initialized lazily during the
      /// first call to do_decay()
      mutable std::unique_ptr<marley::ChebyshevInterpolatingFunction> Exf_cdf_;

      double sample_Exf(marley::Generator& gen) const;

      void sample_spin_parity(double Exf, int& two_Jf, marley::Parity& Pf,
        marley::Generator& gen) const;

      /// @brief Returns the maximum accessible excitation energy to be
      /// used when integrating over the continuum
      virtual double E_c_max() const = 0;
  };

  /// @brief %Fragment emission ExitChannel that leads to a discrete nuclear
  /// level in the final state
  class FragmentDiscreteExitChannel : public DiscreteExitChannel,
    public FragmentExitChannel
  {
    public:

      /// @copydoc marley::ExitChannel::ExitChannel()
      /// @copydoc marley::DiscreteExitChannel( const marley::Level& )
      /// @copydoc marley::FragmentExitChannel( const marley::Fragment& )
      FragmentDiscreteExitChannel(int pdgi, int qi, double Exi, int twoJi,
        marley::Parity Pi, double rho_i, marley::StructureDatabase& sdb,
        const marley::Level& flev, const marley::Fragment& frag)
        : ExitChannel( pdgi, qi, Exi, twoJi, Pi, rho_i, sdb ),
        DiscreteExitChannel( flev ), FragmentExitChannel( frag )
      {
        this->compute_total_width();
      }

      virtual void compute_total_width() final override;
  };

  /// @brief %Gamma emission exit channel that leads to a discrete nuclear
  /// level in the final state
  class GammaDiscreteExitChannel : public DiscreteExitChannel,
    public GammaExitChannel
  {
    public:

      /// @copydoc marley::ExitChannel::ExitChannel()
      /// @copydoc marley::DiscreteExitChannel( const marley::Level& )
      /// @copydoc marley::GammaExitChannel()
      GammaDiscreteExitChannel(int pdgi, int qi, double Exi, int twoJi,
        marley::Parity Pi, double rho_i, marley::StructureDatabase& sdb,
        const marley::Level& flev) : ExitChannel( pdgi, qi, Exi, twoJi, Pi,
        rho_i, sdb ), DiscreteExitChannel( flev ), GammaExitChannel()
      {
        this->compute_total_width();
      }

      virtual void compute_total_width() final override;
  };


  /// @brief %Fragment emission ExitChannel that leads to the unbound continuum
  /// in the final state
  class FragmentContinuumExitChannel : public ContinuumExitChannel,
    public FragmentExitChannel
  {
    public:

      /// @copydoc marley::ExitChannel::ExitChannel()
      /// @copydoc marley::ContinuumExitChannel( double )
      /// @copydoc marley::FragmentExitChannel( const marley::Fragment& )
      FragmentContinuumExitChannel(int pdgi, int qi, double Exi, int twoJi,
        marley::Parity Pi, double rho_i, marley::StructureDatabase& sdb,
        double Ec_min, const marley::Fragment& frag)
        : ExitChannel( pdgi, qi, Exi, twoJi, Pi, rho_i, sdb ),
        ContinuumExitChannel( Ec_min ), FragmentExitChannel( frag )
      {
        this->compute_total_width();
      }

      virtual double differential_width( double Exf,
        bool store_jpi_widths = false ) const final override;

      inline virtual double E_c_max() const final override
        { return this->max_Exf(); }
  };

  /// @brief %Gamma emission exit channel that leads to the unbound continuum
  /// in the final state
  class GammaContinuumExitChannel : public ContinuumExitChannel,
    public GammaExitChannel
  {
    public:

      /// @copydoc marley::ExitChannel::ExitChannel()
      /// @copydoc marley::ContinuumExitChannel( double )
      /// @copydoc marley::GammaExitChannel()
      GammaContinuumExitChannel(int pdgi, int qi, double Exi, int twoJi,
        marley::Parity Pi, double rho_i, marley::StructureDatabase& sdb,
        double Ec_min) : ExitChannel( pdgi, qi, Exi, twoJi, Pi, rho_i, sdb ),
        ContinuumExitChannel( Ec_min ), GammaExitChannel()
      {
        this->compute_total_width();
      }

      virtual double differential_width( double Exf,
        bool store_jpi_widths = false ) const final override;

      inline virtual double E_c_max() const final override
        { return Exi_; }
  };
}
