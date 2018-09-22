// Functions to calculate the Coulomb wavefunctions using the GNU Scientific
// Library
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_coulomb.h"

#include "marley/coulomb_wavefunctions.hh"

std::complex<double> coulomb_H_plus(int l, double eta, double rho) {
  // H+ = G + i F
  gsl_sf_result F, Fp, G, Gp;
  double exp_F, exp_G;
  int code = gsl_sf_coulomb_wave_FG_e(eta, rho, static_cast<double>(l), 0,
    &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  if (code == GSL_EOVRFLW) {
    double Fl = F.val * std::exp(exp_F);
    double Gl = G.val * std::exp(exp_G);
    return std::complex<double>(Gl, Fl);
  }

  return std::complex<double>(G.val, F.val);
}
