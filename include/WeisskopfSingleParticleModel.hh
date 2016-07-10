#pragma once
#include "GammaStrengthFunctionModel.hh"
#include "LevelDensityModel.hh"
#include "StructureDatabase.hh"

namespace marley {

  // Model for strength functions to use when modeling gamma-ray competition
  // during Hauser-Fesbach decay processes
  class WeisskopfSingleParticleModel : public GammaStrengthFunctionModel
  {
    public:

      /// @param Z Atomic number
      /// @param A Mass number
      /// @param D0 Level spacing parameter (MeV)
      WeisskopfSingleParticleModel(int Z, int A, double D0 = 1.);

      virtual double strength_function(TransitionType type, int l,
        double e_gamma) override;

      virtual double transmission_coefficient(TransitionType type, int l,
        double e_gamma) override;

    private:

      double D0_; // Level spacing parameter (MeV)

      double partial_decay_width(TransitionType type, int l, double e_gamma);
  };

}
