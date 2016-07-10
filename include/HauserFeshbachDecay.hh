#pragma once
#include <cmath>
#include <ostream>
#include <unordered_map>

#include "ExitChannel.hh"
#include "Fragment.hh"
#include "GammaStrengthFunctionModel.hh"
#include "Level.hh"
#include "LevelDensityModel.hh"
#include "MassTable.hh"
#include "Parity.hh"
#include "SphericalOpticalModel.hh"
#include "StructureDatabase.hh"
#include "marley_utils.hh"

namespace marley {

  class HauserFeshbachDecay {
    public:

      HauserFeshbachDecay(const marley::Particle& compound_nucleus,
        double Exi, int twoJi, marley::Parity Pi, marley::Generator& gen);

      void build_exit_channels();

      bool do_decay(double& Exf, int& twoJf, marley::Parity& Pf,
        marley::Particle& emitted_particle, marley::Particle& residual_nucleus);

      // Approximates the electric potential energy (in MeV) of two nuclei
      // (with atomic numbers Z, z and mass numbers A, a) that are just
      // touching by modeling them as uniformly charged spheres with radii
      // given by R = (1.2 fm) * A^(1/3).
      static inline double coulomb_barrier(int Z, int A, int z, int a) {
        return Z * z * marley_utils::e2
          / (marley_utils::r0 * (std::pow(A, 1.0/3.0) + std::pow(a, 1.0/3.0)));
      }

      // Table of nuclear fragments that will be considered when computing
      // branching ratios for nuclear de-excitations.
      inline static const std::vector<marley::Fragment>& get_fragments() {
        return fragments;
      }

      static double get_fragment_emission_threshold(const int Zi, const int Ai,
        const marley::Fragment& f);

      // Maximum value of the orbital angular momentum to use in
      // Hauser-Feshbach calculations
      // TODO: make this a user-controlled value specified in the configuration
      // file
      static constexpr int l_max = 2;

      // Default step size for computing optical model transmission
      // coefficients via the Numerov method
      // TODO: make this a user-controlled value specified in the configuration
      // file
      static constexpr double DEFAULT_NUMEROV_STEP_SIZE = 0.1; // fm

      void print(std::ostream& out) const;

    protected:

      // Mass of a charged pion
      static constexpr double mpiplus = 139.57018; // MeV
      // Squared pion Compton wavelength
      static constexpr double lambda_piplus2 = (marley_utils::hbar_c
        / mpiplus) * (marley_utils::hbar_c / mpiplus); // fm

      // Table of nuclear fragments that will be considered when computing
      // branching ratios for nuclear de-excitations.
      static const std::vector<marley::Fragment> fragments;

      // Default number of subintervals to use when integrating fragment and
      // gamma decay widths over the energy continuum.
      // TODO: remove this when you don't need to use it in old testing code
      // anymore
      static constexpr int DEFAULT_CONTINUUM_SUBINTERVALS = 50;

      // Helper functions used internally while computing decay widths using
      // the Hauser-Feshbach statistical model
      double fragment_continuum_partial_width(marley::SphericalOpticalModel& om,
        marley::LevelDensityModel& ldm, const marley::Fragment& frag,
        double fragment_KE, double Exf);

      double gamma_continuum_partial_width(marley::LevelDensityModel& ldm,
        marley::GammaStrengthFunctionModel& gsfm, double Exf);

      // Particle object representing the compound nucleus before it decays
      const marley::Particle& compound_nucleus_;
      const double Exi_; // initial nuclear excitation energy
      const int twoJi_; // initial spin
      const marley::Parity Pi_; // initial parity

      // Generator to use for obtaining discrete level data/nuclear models and
      // simulating statistical decays
      marley::Generator& gen_;

      // Total decay width for this compound nucleus (times an overall factor
      // of 2*pi*rho[Exi_,twoJi_,Pi_] which is omitted for speed)
      double total_width_ = 0.;

      // Table of exit channels used for sampling decays
      std::vector<std::unique_ptr<marley::ExitChannel> > exit_channels_;
  };

}

// Operator for printing HauserFeshbachDecay objects to a std::ostream
inline std::ostream& operator<<(std::ostream& out,
  const marley::HauserFeshbachDecay& hfd)
{
  hfd.print(out);
  return out;
}
