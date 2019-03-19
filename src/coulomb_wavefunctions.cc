// Functions to calculate the Coulomb wavefunctions using the GNU Scientific
// Library
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_coulomb.h"

#include "marley/coulomb_wavefunctions.hh"
#include "marley/Logger.hh"

std::complex<double> coulomb_H_plus(int l, double eta, double rho) {

  // Enable the MARLEY gsl error handler upon using this function
  // for the first time
  static bool error_handler_set = false;
  if ( !error_handler_set ) {
    gsl_set_error_handler( &marley_gsl_error_handler );
    error_handler_set = true;
  }

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

/// Custom GSL error handler used by MARLEY when computing the Coulomb wavefunctions
void marley_gsl_error_handler(const char* reason, const char* file, int line,
  int gsl_errno)
{
  // Overflows are expected sometimes, and they are handled explicitly by coulomb_H_plus()
  if ( gsl_errno == GSL_EOVRFLW ) return;

  MARLEY_LOG_ERROR() << "GSL error: " << reason << " in file " << file
    << " at line " << line << " with error code " << gsl_errno;
}
