/// @file
/// @copyright Copyright (C) 2016-2020 Steven Gardiner
/// @license GNU General Public License, version 3
//
// This file is part of MARLEY (Model of Argon Reaction Low Energy Yields)
//
// MARLEY is free software: you can redistribute it and/or modify it under the
// terms of version 3 of the GNU General Public License as published by the
// Free Software Foundation.
//
// For the full text of the license please see COPYING or
// visit http://opensource.org/licenses/GPL-3.0
//
// Please respect the MCnet academic usage guidelines. See GUIDELINES
// or visit https://www.montecarlonet.org/GUIDELINES for details.

#pragma once
#include <memory>
#include <ostream>

#include "marley/ExitChannel.hh"
#include "marley/Parity.hh"

namespace marley {

  // Forward-declare the StructureDatabase class
  class StructureDatabase;

  /// @brief Monte Carlo implementation of the Hauser-Feshbach statistical
  /// model for decays of highly-excited nuclei
  class HauserFeshbachDecay {

    public:

      /// @param compound_nucleus Particle object that represents
      /// the excited nucleus
      /// @param Exi Initial excitation energy (MeV)
      /// @param twoJi Two times the initial nuclear spin
      /// @param Pi Initial nuclear parity
      /// @param sdb Reference to the StructureDatabase to use in calculations
      HauserFeshbachDecay( const marley::Particle& compound_nucleus,
        double Exi, int twoJi, marley::Parity Pi,
        marley::StructureDatabase& sdb );

      /// @brief Simulates a decay of the compound nucleus
      /// @param[out] Exf Final nuclear excitation energy (MeV)
      /// @param[out] twoJf Two times the final nuclear spin
      /// @param[out] Pf Final nuclear parity
      /// @param[out] emitted_particle Particle object representing the nuclear
      /// fragment or &gamma;-ray emitted during the decay
      /// @param[out] residual_nucleus Particle object representing the
      /// final-state nucleus
      /// param gen Generator to use for random sampling
      bool do_decay( double& Exf, int& twoJf, marley::Parity& Pf,
        marley::Particle& emitted_particle, marley::Particle& residual_nucleus,
        marley::Generator& gen );

      /// @brief Print information about the possible decay channels to a
      /// std::ostream
      void print( std::ostream& out ) const;

      /// @brief Get a non-const reference to the owned vector of ExitChannel
      /// pointers
      inline std::vector< std::unique_ptr<marley::ExitChannel> >&
        exit_channels();

      /// @brief Get a const reference to the owned vector of ExitChannel
      /// pointers
      inline const std::vector< std::unique_ptr<marley::ExitChannel> >&
        exit_channels() const;

      /// @brief Helper function for do_decay(). Samples an ExitChannel
      /// using the partial decay widths as weights
      const std::unique_ptr<marley::ExitChannel>& sample_exit_channel(
        marley::Generator& gen) const;

    private:

      /// @brief Helper function called by the constructor. Loads
      /// exit_channels_ with ExitChannel objects representing all of the
      /// possible decay modes
      void build_exit_channels( marley::StructureDatabase& sdb );

      /// @brief Particle object that represents the compound nucleus before it
      /// decays
      const marley::Particle compound_nucleus_;
      const double Exi_; ///< Initial nuclear excitation energy
      const int twoJi_; ///< Two times the initial nuclear spin
      const marley::Parity Pi_; ///< Two times the initial nuclear parity

      /// @brief Total decay width (MeV) for the compound nucleus
      double total_width_ = 0.;

      /// @brief Table of ExitChannel objects used for sampling decays
      std::vector< std::unique_ptr<marley::ExitChannel> > exit_channels_;
  };

  // Inline function definitions
  inline std::vector<std::unique_ptr<marley::ExitChannel> >&
    HauserFeshbachDecay::exit_channels() { return exit_channels_; }

  inline const std::vector<std::unique_ptr<marley::ExitChannel> >&
    HauserFeshbachDecay::exit_channels() const { return exit_channels_; }
}

/// @brief Operator for printing a HauserFeshbachDecay object to a std::ostream
inline std::ostream& operator<<( std::ostream& out,
  const marley::HauserFeshbachDecay& hfd )
{
  hfd.print(out);
  return out;
}
