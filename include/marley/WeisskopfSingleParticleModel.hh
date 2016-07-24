#pragma once
#include "marley/GammaStrengthFunctionModel.hh"
#include "marley/LevelDensityModel.hh"
#include "marley/StructureDatabase.hh"

namespace marley {

  /// @brief Implements the Weisskopf single-particle estimates of the
  /// gamma-ray strength functions
  /// @details Under this model, the electric (@f$f_{\text{E}\ell}@f$) and
  /// magnetic (@f$f_{\text{M}\ell}@f$) gamma-ray strength functions are
  /// independent of the gamma energy and are given by
  /// @f[ f_{\text{E}\ell} = \frac{2\alpha\Lambda}{\text{D}_0}
  /// \left(\frac{R}{\hbar c}\right)^{\!2\ell} @f]
  /// and
  /// @f[ f_{\text{M}\ell} = 10\left(\frac{\hbar c}
  /// {R\,m_\text{p}}\right)^{\!2}
  /// f_{\text{E}\ell}, @f]
  /// where @f$\alpha@f$ is the fine-structure constant, @f$\text{D}_0@f$
  /// is the level spacing parameter (MeV), @f$R = (1.2\text{ fm})A^{1/3}@f$
  /// is the approximate nuclear radius, @f$m_\text{p}@f$ is the proton mass
  /// (MeV), and @f$\Lambda@f$ is a function of the multipolarity @f$\ell@f$
  /// given by
  /// @f[ \Lambda \equiv \left(\frac{3}{\ell + 3}\right)^{\!2}
  /// \left(\frac{\ell + 1}{\ell\left[(2\ell + 1)!!\right]^2}
  /// \right).@f]
  /// These estimates are typically only good to an order of
  /// magnitude, so using a more sophisticated model, e.g., the
  /// StandardLorentzianModel, is strongly recommended.
  class WeisskopfSingleParticleModel : public GammaStrengthFunctionModel
  {
    public:

      /// @param Z Atomic number of the desired nuclide
      /// @param A Mass number of the desired nuclide
      /// @param D0 %Level spacing parameter (MeV)
      WeisskopfSingleParticleModel(int Z, int A, double D0 = 1.);

      /// @copydoc marley::GammaStrengthFunctionModel::strength_function
      /// @note As described above, the Weisskopf estimates of the gamma-ray
      /// strength functions are independent of the gamma energy, and so the
      /// parameter e_gamma is ignored by this function.
      virtual double strength_function(TransitionType type, int l,
        double e_gamma) override;

      virtual double transmission_coefficient(TransitionType type, int l,
        double e_gamma) override;

    private:

      double D0_; // Level spacing parameter (MeV)

      double partial_decay_width(TransitionType type, int l, double e_gamma);
  };

}
