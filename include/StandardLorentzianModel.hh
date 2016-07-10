#pragma once
#include "GammaStrengthFunctionModel.hh"

namespace marley {

  // Model for strength functions to use when modeling gamma-ray competition
  // during Hauser-Fesbach decay processes
  class StandardLorentzianModel : public GammaStrengthFunctionModel
  {

    public:

      StandardLorentzianModel(int Z, int A);

      virtual double strength_function(TransitionType type, int l,
        double e_gamma) override;

      virtual double transmission_coefficient(TransitionType type, int l,
        double e_gamma) override;

    private:

      double strength_function_coefficient(TransitionType type,
        int l, double e_gamma);

      /// @todo Consider other more elegant ways of storing these parameters
      double e_E1_; ///< E1 giant resonance energy (MeV)
      double sigma_E1_; ///< E1 giant resonance strength (mb)
      double gamma_E1_; ///< E1 giant resonance width (MeV)

      double e_E2_; ///< E2 giant resonance energy (MeV)
      double sigma_E2_; ///< E2 giant resonance strength (mb)
      double gamma_E2_; ///< E2 giant resonance width (MeV)

      double e_M1_; ///< M1 giant resonance energy (MeV)
      double sigma_M1_; ///< M1 giant resonance strength (mb)
      double gamma_M1_; ///< M1 giant resonance width (MeV)
  };

}
