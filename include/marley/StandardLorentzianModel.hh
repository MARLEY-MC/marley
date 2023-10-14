/// @file
/// @copyright Copyright (C) 2016-2023 Steven Gardiner
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
#include "marley/GammaStrengthFunctionModel.hh"

namespace marley {

  /// @brief Implements Brink-Axel strength functions based on the <a
  /// href="https://www-nds.iaea.org/RIPL-3/">Reference Input Parameter
  /// Library</a>'s Standard Lorentzian (SLO) Model
  /// @details <p>Under this model, the gamma-ray strength functions are
  /// taken to have a Lorentzian shape with an energy-independent width.
  /// If @f$E_{\text{X}\ell}@f$, @f$\Gamma_{\text{X}\ell}@f$, and
  /// @f$\sigma_{\text{X}\ell}@f$ are respectively the energy, width, and peak
  /// cross section of the @f$\text{X}\ell@f$ giant resonance, then the
  /// standard Lorentzian model gamma-ray strength function is given by @f[
  /// f_{\text{X}\ell}(E_\gamma) = \frac{\sigma_{\text{X}\ell}}
  /// {(2\ell+1)\pi^2(\hbar c)^2}\left[\frac{\Gamma_{\text{X}\ell}^2
  /// E_\gamma^{3-2\ell}}{\left(E_\gamma^2 - E_{\text{X}\ell}^2\right)^2 +
  /// E_\gamma^2\Gamma_{\text{X}\ell}^2}\right] @f] where @f$E_\gamma@f$
  /// is the gamma-ray energy and the type of transition @f$\text{X}@f$
  /// is either electric @f$\text{(E)}@f$ or magnetic @f$\text{(M)}@f$.</p>
  /// <p>The giant resonance parameters @f$E_{\text{X}\ell}@f$,
  /// @f$\Gamma_{\text{X}\ell}@f$, and @f$\sigma_{\text{X}\ell}@f$
  /// used by StandardLorentzianModel are the same as those used by
  /// default in the <a href="http://talys.eu">TALYS</a> nuclear reaction code,
  /// <a href="http://www.talys.eu/fileadmin/talys/user/docs/talys1.6.pdf">
  /// version 1.6</a>. More details about these parameters are given in the
  /// table below.
  /// <table><tr><th>Transition</th><th>Parameters</th><th>Units</th>
  /// <th>Source</th></tr>
  /// <tr><td rowspan="3">Electric dipole (E1)</td><td>@f$E_{\text{E}1}
  /// = 31.2A^{-1/3} + 20.6A^{-1/6}@f$</td><td>MeV</td>
  /// <td rowspan="3">Empirical fit for spherical nuclei from
  /// the <a href="https://www-nds.iaea.org/RIPL-2/handbook/ripl2.pdf">
  /// RIPL-2 handbook</a>, p. 129 </td></tr>
  /// <tr><td>@f$\Gamma_{\text{E}1} = 0.026{E_{\text{E}1}}^{1.91}@f$
  /// </td><td>MeV</td></tr>
  /// <tr><td>@f$\displaystyle\sigma_{\text{E}1} = 1.2\left(\frac{120NZ}{\pi A
  /// \,\Gamma_{\text{E}1}}\right)@f$
  /// </td><td>mb</td></tr>
  /// <tr><td rowspan="3">Electric quadrupole (E2)</td><td>@f$E_{\text{E}2}
  /// = 63A^{-1/3}@f$</td><td>MeV</td>
  /// <td rowspan="3"> Global fit given by Kopecky in the
  /// <a href="https://www-nds.iaea.org/ripl/readme/ripl_handbook.ps">RIPL-1
  /// handbook</a>, p. 103</td></tr>
  /// <tr><td>@f$\Gamma_{\text{E}2} = 6.11 - 0.012A@f$
  /// </td><td>MeV</td></tr>
  /// <tr><td>@f$\displaystyle\sigma_{\text{E}2}
  /// = \frac{0.00014Z^2 E_{\text{E}2}}{A^{1/3}\Gamma_{\text{E}2}}@f$
  /// </td><td>mb</td></tr>
  /// <tr><td rowspan="3">Magnetic dipole (M1)</td><td>@f$E_{\text{M}1}
  /// = 41A^{-1/3}@f$</td><td>MeV</td>
  /// <td rowspan="3"> Global SLO model fit given in the
  /// <a href="https://www-nds.iaea.org/RIPL-2/handbook/ripl2.pdf">RIPL-2
  /// handbook</a>, p. 132
  /// </td></tr>
  /// <tr><td>@f$\Gamma_{\text{M}1} = 4@f$
  /// </td><td>MeV</td></tr>
  /// <tr><td><p>@f$\displaystyle\sigma_{\text{M}1}
  /// = 3\,\pi^2\hbar^2c^2\left[\frac{\left(B_\text{n}^2
  /// - E_{\text{M}1}^2\right)^2 + B_\text{n}^2\,
  /// \Gamma_{\text{M}1}^2}{B_\text{n}\,\Gamma_{\text{M}1}^2}
  /// \right]\Bigg[\frac{f_{\text{E}1}(B_\text{n})}{0.0588A^{0.878}}\Bigg]
  /// @f$</p><p>where @f$B_\text{n}@f$ = 7 MeV and @f$f_{\text{E}1}@f$
  /// is calculated using the E1 parameters above.</p>
  /// </td><td>mb</td></tr>
  /// <tr><td rowspan="3">Other electric transitions (E3+)</td>
  /// <td>@f$E_{\text{E}\ell}
  /// = E_{\text{E}2}@f$</td><td>MeV</td>
  /// <td rowspan="3">Default approximation used by the
  /// <a href="http://talys.eu">TALYS</a> nuclear code,
  /// <a href="http://www.talys.eu/fileadmin/talys/user/docs/talys1.6.pdf">
  /// version 1.6</a></td></tr>
  /// <tr><td>@f$\Gamma_{\text{E}\ell} = \Gamma_{\text{E}2}@f$
  /// </td><td>MeV</td></tr>
  /// <tr><td><p>@f$\displaystyle\sigma_{\text{E}\ell} =
  /// (0.0008)^{\ell - 2}\,\sigma_{\text{E}2}@f$
  /// </td><td>mb</td></tr>
  /// <tr><td rowspan="3">Other magnetic transitions (M2+)</td>
  /// <td>@f$M_{\text{M}\ell}
  /// = M_{\text{M}1}@f$</td><td>MeV</td>
  /// <td rowspan="3">Default approximation used by the
  /// <a href="http://talys.eu">TALYS</a> nuclear code,
  /// <a href="http://www.talys.eu/fileadmin/talys/user/docs/talys1.6.pdf">
  /// version 1.6</a></td></tr>
  /// <tr><td>@f$\Gamma_{\text{M}\ell} = \Gamma_{\text{M}1}@f$
  /// </td><td>MeV</td></tr>
  /// <tr><td><p>@f$\displaystyle\sigma_{\text{M}\ell} =
  /// (0.0008)^{\ell - 1}\,\sigma_{\text{M}1}@f$
  /// </td><td>mb</td></tr>
  /// </table>
  /// </p>
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
