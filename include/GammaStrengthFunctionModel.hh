#pragma once

namespace marley {

  class Level;
  class Parity;

  // Abstract base class
  // Model for strength functions to use when modeling gamma-ray competition
  // during Hauser-Fesbach decay processes
  class GammaStrengthFunctionModel {

    public:

      GammaStrengthFunctionModel(int Z, int A);

      enum class TransitionType { electric, magnetic, unphysical };

      static TransitionType determine_transition_type(int twoJi,
        marley::Parity Pi, int twoJf, marley::Parity Pf, int& l);

      static TransitionType determine_transition_type(int twoJi,
        marley::Parity Pi, marley::Level& level_f, int& l);

      static TransitionType determine_transition_type(marley::Level& level_i,
        marley::Level& level_f, int& l);

      virtual double strength_function(TransitionType type, int l,
        double e_gamma) = 0;

      virtual double transmission_coefficient(TransitionType type, int l,
        double e_gamma) = 0;

      double transmission_coefficient(double Exi, int twoJ,
        marley::Parity Pi, marley::Level& level_f);

    protected:

      static void check_multipolarity(int l);

      int Z_; ///< Atomic number
      int A_; ///< Mass number
  };

}
