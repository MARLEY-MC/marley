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
#include <cmath>
#include <memory>
#include <ostream>
#include <unordered_map>

#include "marley/marley_utils.hh"
#include "marley/ExitChannel.hh"
#include "marley/Fragment.hh"
#include "marley/GammaStrengthFunctionModel.hh"
#include "marley/Level.hh"
#include "marley/LevelDensityModel.hh"
#include "marley/MassTable.hh"
#include "marley/OpticalModel.hh"
#include "marley/Parity.hh"
#include "marley/StructureDatabase.hh"

namespace marley {

  /// @brief Monte Carlo implementation of the Hauser-Feshbach statistical
  /// model for decays of highly-excited nuclei
  class HauserFeshbachDecay {

    public:

      /// @param compound_nucleus Particle object that represents
      /// the excited nucleus
      /// @param Exi Initial excitation energy (MeV)
      /// @param twoJi Two times the initial nuclear spin
      /// @param Pi Initial nuclear parity
      /// @param gen Reference to the Generator to use for random sampling
      HauserFeshbachDecay(const marley::Particle& compound_nucleus,
        double Exi, int twoJi, marley::Parity Pi, marley::Generator& gen);

      /// @brief Simulates a decay of the compound nucleus
      /// @param[out] Exf Final nuclear excitation energy (MeV)
      /// @param[out] twoJf Two times the final nuclear spin
      /// @param[out] Pf Final nuclear parity
      /// @param[out] emitted_particle Particle object representing the nuclear
      /// fragment or &gamma;-ray emitted during the decay
      /// @param[out] residual_nucleus Particle object representing the
      /// final-state nucleus
      bool do_decay(double& Exf, int& twoJf, marley::Parity& Pf,
        marley::Particle& emitted_particle, marley::Particle& residual_nucleus);

      /// @brief Get a const reference to the table of Fragment objects
      /// that will be considered when simulating compound nucleus decays
      inline static const std::vector<marley::Fragment>& get_fragments();

      /// @brief Maximum value of the orbital angular momentum to use when
      /// considering compound nucleus decays to the continuum of nuclear
      /// levels
      /// @todo Make this a user-controlled value specified in the configuration
      /// file
      static constexpr int l_max_ = 5;

      /// @brief Print information about the possible decay channels to a
      /// std::ostream
      void print(std::ostream& out) const;

      /// @brief Get a non-const reference to the owned vector of ExitChannel
      /// pointers
      inline std::vector<std::unique_ptr<marley::ExitChannel> >&
        exit_channels();

      /// @brief Get a const reference to the owned vector of ExitChannel
      /// pointers
      inline const std::vector<std::unique_ptr<marley::ExitChannel> >&
        exit_channels() const;

      /// @brief Helper function for do_decay(). Samples an ExitChannel
      /// using the partial decay widths as weights
      const std::unique_ptr<marley::ExitChannel>& sample_exit_channel() const;

    private:

      /// @brief Helper function called by the constructor. Loads
      /// exit_channels_ with ExitChannel objects representing all of the
      /// possible decay modes
      void build_exit_channels();

      /// @brief Helper function used when computing the partial decay width
      /// for emission of a nuclear fragment to the continuum of nuclear levels
      double fragment_continuum_partial_width(marley::OpticalModel& om,
        marley::LevelDensityModel& ldm, const marley::Fragment& frag,
        double fragment_KE, double Exf);

      /// @brief Helper function used when computing the partial decay
      /// width for gamma emission to the continuum of nuclear levels
      double gamma_continuum_partial_width(marley::LevelDensityModel& ldm,
        marley::GammaStrengthFunctionModel& gsfm, double Exf);

      /// @brief Mass of a charged pion (MeV)
      /// @hideinitializer
      static constexpr double mpiplus_ = 139.57018;

      /// @brief Squared charged pion Compton wavelength (fm)
      /// @hideinitializer
      static constexpr double lambda_piplus2_ = (marley_utils::hbar_c
        / mpiplus_) * (marley_utils::hbar_c / mpiplus_);

      /// @brief Particle object that represents the compound nucleus before it
      /// decays
      const marley::Particle& compound_nucleus_;
      const double Exi_; ///< Initial nuclear excitation energy
      const int twoJi_; ///< Two times the initial nuclear spin
      const marley::Parity Pi_; ///< Two times the initial nuclear parity

      /// @brief Generator to use for obtaining discrete level data/nuclear
      /// models and simulating statistical decays
      marley::Generator& gen_;

      /// @brief Total decay width for the compound nucleus (times an overall
      /// factor of @f$2\pi\rho(Ex_\text{i},J_\text{i},\Pi_\text{i})@f$ which
      /// is omitted for speed)
      double total_width_ = 0.;

      /// @brief Table of nuclear fragments that will be considered when
      /// computing branching ratios for nuclear de-excitations
      /// @hideinitializer
      static const std::vector<marley::Fragment> fragments_;

      /// @brief Table of exit channels used for sampling decays
      std::vector<std::unique_ptr<marley::ExitChannel> > exit_channels_;
  };

  // Inline function definitions
  inline const std::vector<marley::Fragment>&
    HauserFeshbachDecay::get_fragments() { return fragments_; }

  inline std::vector<std::unique_ptr<marley::ExitChannel> >&
    HauserFeshbachDecay::exit_channels() { return exit_channels_; }

  inline const std::vector<std::unique_ptr<marley::ExitChannel> >&
    HauserFeshbachDecay::exit_channels() const { return exit_channels_; }
}

/// @brief Operator for printing HauserFeshbachDecay objects to a std::ostream
inline std::ostream& operator<<(std::ostream& out,
  const marley::HauserFeshbachDecay& hfd)
{
  hfd.print(out);
  return out;
}
