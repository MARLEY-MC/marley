// Unified GSL interface header for MARLEY.
// When GSL is present on the system this header simply pulls in the real
// GSL headers that MARLEY needs.  When the built-in GSL fallback is in use
// it provides the necessary declarations directly.

#pragma once

#ifdef MARLEY_FOUND_GSL

// ---- System GSL: pull in only the headers actually used by MARLEY ----
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_coulomb.h"
#include "gsl/gsl_cdf.h"

#else

// ---- Built-in GSL fallback declarations ----

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum {
  GSL_SUCCESS  = 0,
  GSL_EOVRFLW  = 16
};

typedef struct gsl_sf_result_struct {
  double val;
  double err;
} gsl_sf_result;

typedef void gsl_error_handler_t( const char* reason, const char* file,
  int line, int gsl_errno );

gsl_error_handler_t* gsl_set_error_handler( gsl_error_handler_t* new_handler );

int gsl_sf_coulomb_wave_FG_e( const double eta, const double x,
  const double lam_F, const int k_lam_G, gsl_sf_result* F, gsl_sf_result* Fp,
  gsl_sf_result* G, gsl_sf_result* Gp, double* exp_F, double* exp_G );

int gsl_sf_coulomb_wave_FG_array( double lam_min, int kmax, double eta,
  double x, double* fc_array, double* gc_array, double* F_exponent,
  double* G_exponent );

double gsl_cdf_chisq_Q( const double x, const double nu );

__END_DECLS

#endif /* MARLEY_FOUND_GSL */
