/*
 * Built-in GSL fallback implementation for MARLEY.
 *
 * This single translation unit contains the minimal amalgamated subset of GSL
 * needed by MARLEY when gsl-config is unavailable.
 *
 * MARLEY call sites:
 *   gsl_set_error_handler         -- src/coulomb_wavefunctions.cc
 *   gsl_sf_coulomb_wave_FG_e      -- src/coulomb_wavefunctions.cc
 *   gsl_cdf_chisq_Q               -- include/marley/tests/Histogram.hh
 *   gsl_sf_coulomb_wave_FG_array  -- retained for possible future use
 *
 * Provenance:
 *   Upstream mirror: https://github.com/sjgardiner/gsl
 *   Target tag:      release-2-8 (GSL 2.8)
 *   Commit used:     17f5f257099a6a9ba539d950e95a2110279c3193 (master)
 *
 * Differences versus upstream:
 *   - Required content is inlined into this single source file.
 *   - gsl_sf_gamma_inc_e is omitted (would require specfunc/expint.c).
 *   - gsl_cdf_gamma_P and gsl_cdf_chisq_P are omitted; only _Q variants are
 *     needed by MARLEY.
 *   - gsl_cdf_gamma_Q and gsl_cdf_chisq_Q deliberately reject non-physical
 *     domains (a <= 0, b <= 0, nu <= 0) through the GSL error path.
 *
 * Files amalgamated (dependency order):
 *   err/error.c, err/stream.c
 *   sys/fdiv.c, sys/infnan.c, sys/coerce.c
 *   sys/minmax.c, sys/prec.c
 *   complex/inline.c, complex/math.c
 *   specfunc/log.c, specfunc/trig.c, specfunc/exp.c, specfunc/gamma.c
 *   specfunc/elementary.c, specfunc/zeta.c, specfunc/psi.c
 *   specfunc/airy.c, specfunc/coulomb.c
 *   specfunc/cheb_eval.c
 *   specfunc/cheb_eval_mode.c
 *   specfunc/erfc.c (partial: gsl_sf_erfc_e only)
 *   specfunc/gamma_inc.c (partial: gsl_sf_gamma_inc_P, gsl_sf_gamma_inc_Q only)
 *   cdf/gamma.c (partial: gsl_cdf_gamma_Q only)
 *   cdf/chisq.c (partial: gsl_cdf_chisq_Q only)
 *
 * Licensing:
 *   Upstream GSL code is licensed under GPLv3 (and parts LGPLv3).
 *   See upstream GSL licensing files for complete terms and notices.
 */

// Ignore eveything here unless we actually need to build the built-in
// fallback for the GSL library
#ifndef MARLEY_FOUND_GSL

/* begin inlined header: specfunc/error.h */
#ifndef MARLEY_BUILTIN_GSL_SPECFUNC_ERROR_H
#define MARLEY_BUILTIN_GSL_SPECFUNC_ERROR_H
#define OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)
#define INTERNAL_OVERFLOW_ERROR(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; return GSL_EOVRFLW; } while(0)
#define INTERNAL_UNDERFLOW_ERROR(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; return GSL_EUNDRFLW; } while(0)
#define DOMAIN_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("domain error", GSL_EDOM); } while(0)
#define DOMAIN_ERROR_MSG(msg, result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ((msg), GSL_EDOM); } while(0)
#define DOMAIN_ERROR_E10(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; (result)->e10 = 0 ; GSL_ERROR ("domain error", GSL_EDOM); } while(0)
#define OVERFLOW_ERROR_E10(result) do { (result)->val = GSL_POSINF; (result)->err = GSL_POSINF; (result)->e10 = 0; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define UNDERFLOW_ERROR_E10(result) do { (result)->val = 0.0; (result)->err = GSL_DBL_MIN; (result)->e10 = 0; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)
#define OVERFLOW_ERROR_2(r1,r2) do { (r1)->val = GSL_POSINF; (r1)->err = GSL_POSINF; (r2)->val = GSL_POSINF ; (r2)->err=GSL_POSINF; GSL_ERROR ("overflow", GSL_EOVRFLW); } while(0)
#define UNDERFLOW_ERROR_2(r1,r2) do { (r1)->val = 0.0; (r1)->err = GSL_DBL_MIN; (r2)->val = 0.0 ; (r2)->err = GSL_DBL_MIN; GSL_ERROR ("underflow", GSL_EUNDRFLW); } while(0)
#define DOMAIN_ERROR_2(r1,r2) do { (r1)->val = GSL_NAN; (r1)->err = GSL_NAN;  (r2)->val = GSL_NAN; (r2)->err = GSL_NAN;  GSL_ERROR ("domain error", GSL_EDOM); } while(0)
#define MAXITER_ERROR(result) do { (result)->val = GSL_NAN; (result)->err = GSL_NAN; GSL_ERROR ("too many iterations error", GSL_EMAXITER); } while(0)
#endif
/* end inlined header: specfunc/error.h */

/* begin inlined header: specfunc/check.h */
#ifndef MARLEY_BUILTIN_GSL_SPECFUNC_CHECK_H
#define MARLEY_BUILTIN_GSL_SPECFUNC_CHECK_H
#define CHECK_UNDERFLOW(r) if (fabs((r)->val) < GSL_DBL_MIN) GSL_ERROR("underflow", GSL_EUNDRFLW);
#endif
/* end inlined header: specfunc/check.h */

/* begin inlined header: specfunc/eval.h */
#ifndef MARLEY_BUILTIN_GSL_SPECFUNC_EVAL_H
#define MARLEY_BUILTIN_GSL_SPECFUNC_EVAL_H
#define EVAL_RESULT(fn)    gsl_sf_result result;    int status = fn;    if (status != GSL_SUCCESS) {      GSL_ERROR_VAL(#fn, status, result.val);    } ;    return result.val;
#define EVAL_DOUBLE(fn)    int status = fn;    if (status != GSL_SUCCESS) {      GSL_ERROR_VAL(#fn, status, result);    } ;    return result;
#endif
/* end inlined header: specfunc/eval.h */

/* begin inlined header: specfunc/chebyshev.h */
#ifndef MARLEY_BUILTIN_GSL_CHEBYSHEV_H
#define MARLEY_BUILTIN_GSL_CHEBYSHEV_H
struct cheb_series_struct {
  double * c;
  int order;
  double a;
  double b;
  int order_sp;
};
typedef struct cheb_series_struct cheb_series;
#endif
/* end inlined header: specfunc/chebyshev.h */
/* Error handling ----------------------------------------------------------*/
/* begin inlined source: err/error.c */
/* err/error.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
gsl_error_handler_t * gsl_error_handler = NULL;

void
gsl_error (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_error_handler)
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  gsl_stream_printf ("ERROR", file, line, reason);

  fflush (stdout);
  fprintf (stderr, "Default GSL error handler invoked.\n");
  fflush (stderr);

  abort ();
}

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler)
{
  gsl_error_handler_t * previous_handler = gsl_error_handler;
  gsl_error_handler = new_handler;
  return previous_handler;
}



/* end inlined source: err/error.c */
/* begin inlined source: err/stream.c */
/* err/stream.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
FILE * gsl_stream = NULL ;
static void (*gsl_stream_handler)(const char * label, const char * file,
                                   int line, const char * reason) = NULL;

void
gsl_stream_printf (const char *label, const char *file, int line,
                   const char *reason)
{
  if (gsl_stream == NULL)
    {
      gsl_stream = stderr;
    }
  if (gsl_stream_handler)
    {
      (*gsl_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);

}

/* end inlined source: err/stream.c */
/* System utilities --------------------------------------------------------*/
/* begin inlined source: sys/fdiv.c */
/* sys/fdiv.c
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
double
gsl_fdiv (const double x, const double y)
{
  return x / y;
}
/* end inlined source: sys/fdiv.c */
/* begin inlined source: sys/infnan.c */
/* sys/infnan.c
 *
 * Copyright (C) 2001, 2004, 2007, 2010 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <math.h>

#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
double gsl_nan (void)
{
  return gsl_fdiv (0.0, 0.0);
}

double gsl_posinf (void)
{
  return gsl_fdiv (+1.0, 0.0);
}

double gsl_neginf (void)
{
  return gsl_fdiv (-1.0, 0.0);
}


int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

#if defined(_MSC_VER) /* Microsoft Visual C++ */
#include <float.h>
int
gsl_isnan (const double x)
{
  return _isnan(x);
}

int
gsl_isinf (const double x)
{
  int fpc = _fpclass(x);

  if (fpc == _FPCLASS_PINF)
    return +1;
  else if (fpc == _FPCLASS_NINF)
    return -1;
  else
    return 0;
}

int
gsl_finite (const double x)
{
  return _finite(x);
}
#else

# if HAVE_DECL_ISFINITE
int
gsl_finite (const double x)
{
  return isfinite(x);
}
# elif HAVE_DECL_FINITE
int
gsl_finite (const double x)
{
  return finite(x);
}
# elif HAVE_IEEE_COMPARISONS
int
gsl_finite (const double x)
{
  const double y = x - x;
  int status = (y == y);
  return status;
}
# else
# error "cannot define gsl_finite without HAVE_DECL_FINITE or HAVE_IEEE_COMPARISONS"
# endif

# if HAVE_DECL_ISNAN
int
gsl_isnan (const double x)
{
  return isnan(x);
}
#elif HAVE_IEEE_COMPARISONS
int
gsl_isnan (const double x)
{
  int status = (x != x);
  return status;
}
# else
# error "cannot define gsl_isnan without HAVE_DECL_ISNAN or HAVE_IEEE_COMPARISONS"
# endif

# if HAVE_DECL_ISINF
int
gsl_isinf (const double x)
{
  /* isinf(3): In glibc 2.01 and earlier, isinf() returns a
     non-zero value (actually: 1) if x is an infinity (positive or
     negative).  (This is all that C99 requires.) */

  if (isinf(x))
    {
      return (x > 0) ? 1 : -1;
    }
  else
    {
      return 0;
    }
}
# else

int
gsl_isinf (const double x)
{
  if (! gsl_finite(x) && ! gsl_isnan(x))
    {
      return (x > 0 ? +1 : -1);
    }
  else
    {
      return 0;
    }
}

# endif
#endif

/* end inlined source: sys/infnan.c */
/* begin inlined source: sys/coerce.c */
/* sys/coerce.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
double
gsl_coerce_double (const double x)
{
  volatile double y;
  y = x;
  return y;
}

float
gsl_coerce_float (const float x)
{
  volatile float y;
  y = x;
  return y;
}

/* The following function is not needed, but is included for completeness */

long double
gsl_coerce_long_double (const long double x)
{
  volatile long double y;
  y = x;
  return y;
}

/* end inlined source: sys/coerce.c */
/* sys/minmax.c and sys/prec.c use the COMPILE_INLINE_STATIC / build.h
 * mechanism to emit linkable definitions for the inline helpers declared
 * in gsl_minmax.h / gsl_precision.h.  build.h now uses a #ifndef guard
 * so subsequent inclusions from sys/prec.c and complex/inline.c below are
 * silently skipped, while the idempotent macro setup from the first
 * inclusion remains in force. */
/* begin inlined source: sys/minmax.c */
/* sys/minmax.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* Compile all the inline functions */

#define COMPILE_INLINE_STATIC
/* begin inlined header: build.h */
/* Compile subsequent inline functions as static functions */

/* Allow multiple inclusion in a single translation unit (needed for the
 * amalgamated gsl_marley.c build where sys/minmax.c, sys/prec.c, and
 * complex/inline.c each set COMPILE_INLINE_STATIC and include this header).
 * The macro setup below is idempotent, so later inclusions are harmless. */
#ifndef __GSL_BUILD_H__
#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else
#error must be called with COMPILE_INLINE_STATIC
#endif

#endif /* __GSL_BUILD_H__ */

/* end inlined header: build.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
/* Define some static functions which are always available */

double gsl_max (double a, double b)
{
  return GSL_MAX (a, b);
}

double gsl_min (double a, double b)
{
  return GSL_MIN (a, b);
}

/* end inlined source: sys/minmax.c */
/* begin inlined source: sys/prec.c */
/* sys/prec.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* Compile all the inline functions */

#define COMPILE_INLINE_STATIC
/* begin inlined header: build.h */
/* Compile subsequent inline functions as static functions */

/* Allow multiple inclusion in a single translation unit (needed for the
 * amalgamated gsl_marley.c build where sys/minmax.c, sys/prec.c, and
 * complex/inline.c each set COMPILE_INLINE_STATIC and include this header).
 * The macro setup below is idempotent, so later inclusions are harmless. */
#ifndef __GSL_BUILD_H__
#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else
#error must be called with COMPILE_INLINE_STATIC
#endif

#endif /* __GSL_BUILD_H__ */

/* end inlined header: build.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_mode.h */
/* gsl_mode.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_MODE_H__
#define __GSL_MODE_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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


/* Some functions can take a mode argument. This
 * is a rough method to do things like control
 * the precision of the algorithm. This mainly
 * occurs in special functions, but we figured
 * it was ok to have a general facility.
 *
 * The mode type is 32-bit field. Most of
 * the fields are currently unused. Users
 * '|' various predefined constants to get
 * a desired mode.
 */
typedef unsigned int gsl_mode_t;


/* Here are the predefined constants.
 * Note that the precision constants
 * are special because they are used
 * to index arrays, so do not change
 * them. The precision information is
 * in the low order 3 bits of gsl_mode_t
 * (the third bit is currently unused).
 */

/* Note that "0" is double precision,
 * so that you get that by default if
 * you forget a flag.
 */
#define GSL_PREC_DOUBLE  0
#define GSL_PREC_SINGLE  1
#define GSL_PREC_APPROX  2

#ifdef HAVE_INLINE
INLINE_FUN unsigned int GSL_MODE_PREC(gsl_mode_t mt);

INLINE_FUN unsigned int
GSL_MODE_PREC(gsl_mode_t mt)
{ return  (mt & (unsigned int)7); }
#else  /* HAVE_INLINE */
#define GSL_MODE_PREC(mt) ((mt) & (unsigned int)7)
#endif /* HAVE_INLINE */


/* Here are some predefined generic modes.
 */
#define GSL_MODE_DEFAULT  0


__END_DECLS

#endif /* __GSL_MODE_H__ */

/* end inlined header: gsl/gsl_mode.h */
const double gsl_prec_eps[_GSL_PREC_T_NUM] = {
  GSL_DBL_EPSILON,
  GSL_FLT_EPSILON,
  GSL_SFLT_EPSILON
};

const double gsl_prec_sqrt_eps[_GSL_PREC_T_NUM] = {
  GSL_SQRT_DBL_EPSILON,
  GSL_SQRT_FLT_EPSILON,
  GSL_SQRT_SFLT_EPSILON
};

const double gsl_prec_root3_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT3_DBL_EPSILON,
  GSL_ROOT3_FLT_EPSILON,
  GSL_ROOT3_SFLT_EPSILON
};

const double gsl_prec_root4_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT4_DBL_EPSILON,
  GSL_ROOT4_FLT_EPSILON,
  GSL_ROOT4_SFLT_EPSILON
};

const double gsl_prec_root5_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT5_DBL_EPSILON,
  GSL_ROOT5_FLT_EPSILON,
  GSL_ROOT5_SFLT_EPSILON
};

const double gsl_prec_root6_eps[_GSL_PREC_T_NUM] = {
  GSL_ROOT6_DBL_EPSILON,
  GSL_ROOT6_FLT_EPSILON,
  GSL_ROOT6_SFLT_EPSILON
};

/* end inlined source: sys/prec.c */

/* Complex arithmetic ------------------------------------------------------*/
/* complex/inline.c also uses build.h (skipped on second inclusion) and
 * gsl_complex_math.h (skipped by its own include guard); it exists solely
 * to emit linkable versions of the complex inline helpers, which are not
 * needed in this single-TU build. Including it is harmless. */
/* begin inlined source: complex/inline.c */
/* complex/inline.c
 *
 * Copyright (C) 2008 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Compile all the inline functions */

#define COMPILE_INLINE_STATIC
/* begin inlined header: build.h */
/* Compile subsequent inline functions as static functions */

/* Allow multiple inclusion in a single translation unit (needed for the
 * amalgamated gsl_marley.c build where sys/minmax.c, sys/prec.c, and
 * complex/inline.c each set COMPILE_INLINE_STATIC and include this header).
 * The macro setup below is idempotent, so later inclusions are harmless. */
#ifndef __GSL_BUILD_H__
#define __GSL_BUILD_H__

#ifdef COMPILE_INLINE_STATIC
#ifndef HIDE_INLINE_STATIC  /* skip if inline functions are hidden */

#undef __GSL_INLINE_H__
#define __GSL_INLINE_H__  /* first, ignore the gsl_inline.h header file */

#undef INLINE_DECL
#define INLINE_DECL       /* disable inline in declarations */

#undef INLINE_FUN
#define INLINE_FUN        /* disable inline in definitions */

#ifndef HAVE_INLINE       /* enable compilation of definitions in .h files */
#define HAVE_INLINE
#endif

/* Compile range checking code for inline functions used in the library */
#undef GSL_RANGE_CHECK
#define GSL_RANGE_CHECK 1

/* Use the global variable gsl_check_range to enable/disable range checking at
   runtime */
#undef GSL_RANGE_COND
#define GSL_RANGE_COND(x) (gsl_check_range && (x))

#endif
#else
#error must be called with COMPILE_INLINE_STATIC
#endif

#endif /* __GSL_BUILD_H__ */

/* end inlined header: build.h */
/* begin inlined header: gsl/gsl_complex_math.h */
/* complex/gsl_complex_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Jorma Olavi T�htinen, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_MATH_H__
#define __GSL_COMPLEX_MATH_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_complex.h */
/* complex/gsl_complex.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Copyright (C) 2020, 2021 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_H__
#define __GSL_COMPLEX_H__

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


/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

typedef const double *       gsl_const_complex_packed ;
typedef const float *        gsl_const_complex_packed_float  ;
typedef const long double *  gsl_const_complex_packed_long_double ;


/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;

typedef const double *       gsl_const_complex_packed_array ;
typedef const float *        gsl_const_complex_packed_array_float  ;
typedef const long double *  gsl_const_complex_packed_array_long_double ;


/* Yes... this seems weird. Trust us. The point is just that
   sometimes you want to make it obvious that something is
   an output value. The fact that it lacks a 'const' may not
   be enough of a clue for people in some contexts.
 */
typedef double *       gsl_complex_packed_ptr ;
typedef float *        gsl_complex_packed_float_ptr  ;
typedef long double *  gsl_complex_packed_long_double_ptr ;

typedef const double *       gsl_const_complex_packed_ptr ;
typedef const float *        gsl_const_complex_packed_float_ptr  ;
typedef const long double *  gsl_const_complex_packed_long_double_ptr ;

/*
 * If <complex.h> is included, use the C99 complex type.  Otherwise
 * define a type bit-compatible with C99 complex. The GSL_REAL and GSL_IMAG
 * macros require C11 functionality also (_Generic)
 */

/* older gcc compilers claim to be C11 compliant but do not support _Generic */
#if defined(__GNUC__) && (__GNUC__ < 7)
#  define GSL_COMPLEX_LEGACY 1
#endif

#if !defined(GSL_COMPLEX_LEGACY) &&          \
     defined(_Complex_I) &&                  \
     defined(complex) &&                     \
     defined(I) &&                           \
     defined(__STDC__) && (__STDC__ == 1) && \
     defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) /* C11 */

#  define GSL_COMPLEX_DEFINE(R, C) typedef R _Complex C ;

#  define GSL_COMPLEX_P(zp)        (&(zp))
#  define GSL_COMPLEX_EQ(z1,z2)    ((z1) == (z2))
#  define GSL_SET_COMPLEX(zp,x,y)  (*(zp) = (x) + I * (y))

#  define GSL_REAL(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[0])

#  define GSL_IMAG(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[1])

#  define GSL_COMPLEX_P_REAL(zp)   GSL_REAL(*(zp))
#  define GSL_COMPLEX_P_IMAG(zp)   GSL_IMAG(*(zp))
#  define GSL_SET_REAL(zp,x)       do { GSL_REAL(*(zp)) = (x); } while(0)
#  define GSL_SET_IMAG(zp,x)       do { GSL_IMAG(*(zp)) = (x); } while(0)

#else /* legacy complex definitions */

/*
 * According to the C17 standard, 6.2.5 paragraph 13:
 *
 * "Each complex type has the same representation and alignment requirements
 * as an array type containing exactly two elements of the corresponding real
 * type; the first element is equal to the real part, and the second element to
 * the imaginary part, of the complex number."
 */

/*#  define GSL_COMPLEX_DEFINE(R, C) typedef R C[2]*/
#  define GSL_COMPLEX_DEFINE(R, C) typedef struct { R dat[2]; } C ;

#  define GSL_REAL(z)              ((z).dat[0])
#  define GSL_IMAG(z)              ((z).dat[1])
#  define GSL_COMPLEX_P(zp)        ((zp)->dat)
#  define GSL_COMPLEX_P_REAL(zp)   ((zp)->dat[0])
#  define GSL_COMPLEX_P_IMAG(zp)   ((zp)->dat[1])
#  define GSL_COMPLEX_EQ(z1,z2)    (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

#  define GSL_SET_COMPLEX(zp,x,y)  do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#  define GSL_SET_REAL(zp,x)       do {(zp)->dat[0]=(x);} while(0)
#  define GSL_SET_IMAG(zp,y)       do {(zp)->dat[1]=(y);} while(0)

#endif

GSL_COMPLEX_DEFINE(double, gsl_complex)
GSL_COMPLEX_DEFINE(long double, gsl_complex_long_double)
GSL_COMPLEX_DEFINE(float, gsl_complex_float)

#define GSL_SET_COMPLEX_PACKED(zp,n,x,y) do {*((zp)+2*(n))=(x); *((zp)+(2*(n)+1))=(y);} while(0)

__END_DECLS

#endif /* __GSL_COMPLEX_H__ */

/* end inlined header: gsl/gsl_complex.h */
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS           /* empty */
#define __END_DECLS             /* empty */
#endif

__BEGIN_DECLS

/* Complex numbers */

gsl_complex gsl_complex_polar (double r, double theta); /* r= r e^(i theta) */

INLINE_DECL gsl_complex gsl_complex_rect (double x, double y);  /* r= real+i*imag */

#ifdef HAVE_INLINE
INLINE_FUN gsl_complex
gsl_complex_rect (double x, double y)
{                               /* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
#endif

#define GSL_COMPLEX_ONE (gsl_complex_rect(1.0,0.0))
#define GSL_COMPLEX_ZERO (gsl_complex_rect(0.0,0.0))
#define GSL_COMPLEX_NEGONE (gsl_complex_rect(-1.0,0.0))

/* Properties of complex numbers */

double gsl_complex_arg (gsl_complex z); /* return arg(z), -pi< arg(z) <=+pi */
double gsl_complex_abs (gsl_complex z);   /* return |z|   */
double gsl_complex_abs2 (gsl_complex z);  /* return |z|^2 */
double gsl_complex_logabs (gsl_complex z); /* return log|z| */

/* Complex arithmetic operators */

gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b);  /* r=a+b */
gsl_complex gsl_complex_sub (gsl_complex a, gsl_complex b);  /* r=a-b */
gsl_complex gsl_complex_mul (gsl_complex a, gsl_complex b);  /* r=a*b */
gsl_complex gsl_complex_div (gsl_complex a, gsl_complex b);  /* r=a/b */

gsl_complex gsl_complex_add_real (gsl_complex a, double x);  /* r=a+x */
gsl_complex gsl_complex_sub_real (gsl_complex a, double x);  /* r=a-x */
gsl_complex gsl_complex_mul_real (gsl_complex a, double x);  /* r=a*x */
gsl_complex gsl_complex_div_real (gsl_complex a, double x);  /* r=a/x */

gsl_complex gsl_complex_add_imag (gsl_complex a, double y);  /* r=a+iy */
gsl_complex gsl_complex_sub_imag (gsl_complex a, double y);  /* r=a-iy */
gsl_complex gsl_complex_mul_imag (gsl_complex a, double y);  /* r=a*iy */
gsl_complex gsl_complex_div_imag (gsl_complex a, double y);  /* r=a/iy */

gsl_complex gsl_complex_conjugate (gsl_complex z);  /* r=conj(z) */
gsl_complex gsl_complex_inverse (gsl_complex a);    /* r=1/a */
gsl_complex gsl_complex_negative (gsl_complex a);    /* r=-a */

/* Elementary Complex Functions */

gsl_complex gsl_complex_sqrt (gsl_complex z);  /* r=sqrt(z) */
gsl_complex gsl_complex_sqrt_real (double x);  /* r=sqrt(x) (x<0 ok) */

gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b);  /* r=a^b */
gsl_complex gsl_complex_pow_real (gsl_complex a, double b);  /* r=a^b */

gsl_complex gsl_complex_exp (gsl_complex a);    /* r=exp(a) */
gsl_complex gsl_complex_log (gsl_complex a);    /* r=log(a) (base e) */
gsl_complex gsl_complex_log10 (gsl_complex a);  /* r=log10(a) (base 10) */
gsl_complex gsl_complex_log_b (gsl_complex a, gsl_complex b);   /* r=log_b(a) (base=b) */

/* Complex Trigonometric Functions */

gsl_complex gsl_complex_sin (gsl_complex a);  /* r=sin(a) */
gsl_complex gsl_complex_cos (gsl_complex a);  /* r=cos(a) */
gsl_complex gsl_complex_sec (gsl_complex a);  /* r=sec(a) */
gsl_complex gsl_complex_csc (gsl_complex a);  /* r=csc(a) */
gsl_complex gsl_complex_tan (gsl_complex a);  /* r=tan(a) */
gsl_complex gsl_complex_cot (gsl_complex a);  /* r=cot(a) */

/* Inverse Complex Trigonometric Functions */

gsl_complex gsl_complex_arcsin (gsl_complex a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arcsin_real (double a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arccos (gsl_complex a);  /* r=arccos(a) */
gsl_complex gsl_complex_arccos_real (double a);  /* r=arccos(a) */
gsl_complex gsl_complex_arcsec (gsl_complex a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arcsec_real (double a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arccsc (gsl_complex a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arccsc_real (double a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arctan (gsl_complex a);  /* r=arctan(a) */
gsl_complex gsl_complex_arccot (gsl_complex a);  /* r=arccot(a) */

/* Complex Hyperbolic Functions */

gsl_complex gsl_complex_sinh (gsl_complex a);  /* r=sinh(a) */
gsl_complex gsl_complex_cosh (gsl_complex a);  /* r=coshh(a) */
gsl_complex gsl_complex_sech (gsl_complex a);  /* r=sech(a) */
gsl_complex gsl_complex_csch (gsl_complex a);  /* r=csch(a) */
gsl_complex gsl_complex_tanh (gsl_complex a);  /* r=tanh(a) */
gsl_complex gsl_complex_coth (gsl_complex a);  /* r=coth(a) */

/* Inverse Complex Hyperbolic Functions */

gsl_complex gsl_complex_arcsinh (gsl_complex a);  /* r=arcsinh(a) */
gsl_complex gsl_complex_arccosh (gsl_complex a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arccosh_real (double a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arcsech (gsl_complex a);  /* r=arcsech(a) */
gsl_complex gsl_complex_arccsch (gsl_complex a);  /* r=arccsch(a) */
gsl_complex gsl_complex_arctanh (gsl_complex a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arctanh_real (double a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arccoth (gsl_complex a);  /* r=arccoth(a) */

__END_DECLS

#endif /* __GSL_COMPLEX_MATH_H__ */

/* end inlined header: gsl/gsl_complex_math.h */
/* end inlined source: complex/inline.c */
/* begin inlined source: complex/math.c */
/* complex/math.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jorma Olavi T�htinen, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Basic complex arithmetic functions

 * Original version by Jorma Olavi T�htinen <jotahtin@cc.hut.fi>
 *
 * Modified for GSL by Brian Gough, 3/2000
 */

/* The following references describe the methods used in these
 * functions,
 *
 *   T. E. Hull and Thomas F. Fairgrieve and Ping Tak Peter Tang,
 *   "Implementing Complex Elementary Functions Using Exception
 *   Handling", ACM Transactions on Mathematical Software, Volume 20
 *   (1994), pp 215-244, Corrigenda, p553
 *
 *   Hull et al, "Implementing the complex arcsin and arccosine
 *   functions using exception handling", ACM Transactions on
 *   Mathematical Software, Volume 23 (1997) pp 299-335
 *
 *   Abramowitz and Stegun, Handbook of Mathematical Functions, "Inverse
 *   Circular Functions in Terms of Real and Imaginary Parts", Formulas
 *   4.4.37, 4.4.38, 4.4.39
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <math.h>
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_complex.h */
/* complex/gsl_complex.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Copyright (C) 2020, 2021 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_H__
#define __GSL_COMPLEX_H__

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


/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

typedef const double *       gsl_const_complex_packed ;
typedef const float *        gsl_const_complex_packed_float  ;
typedef const long double *  gsl_const_complex_packed_long_double ;


/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;

typedef const double *       gsl_const_complex_packed_array ;
typedef const float *        gsl_const_complex_packed_array_float  ;
typedef const long double *  gsl_const_complex_packed_array_long_double ;


/* Yes... this seems weird. Trust us. The point is just that
   sometimes you want to make it obvious that something is
   an output value. The fact that it lacks a 'const' may not
   be enough of a clue for people in some contexts.
 */
typedef double *       gsl_complex_packed_ptr ;
typedef float *        gsl_complex_packed_float_ptr  ;
typedef long double *  gsl_complex_packed_long_double_ptr ;

typedef const double *       gsl_const_complex_packed_ptr ;
typedef const float *        gsl_const_complex_packed_float_ptr  ;
typedef const long double *  gsl_const_complex_packed_long_double_ptr ;

/*
 * If <complex.h> is included, use the C99 complex type.  Otherwise
 * define a type bit-compatible with C99 complex. The GSL_REAL and GSL_IMAG
 * macros require C11 functionality also (_Generic)
 */

/* older gcc compilers claim to be C11 compliant but do not support _Generic */
#if defined(__GNUC__) && (__GNUC__ < 7)
#  define GSL_COMPLEX_LEGACY 1
#endif

#if !defined(GSL_COMPLEX_LEGACY) &&          \
     defined(_Complex_I) &&                  \
     defined(complex) &&                     \
     defined(I) &&                           \
     defined(__STDC__) && (__STDC__ == 1) && \
     defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) /* C11 */

#  define GSL_COMPLEX_DEFINE(R, C) typedef R _Complex C ;

#  define GSL_COMPLEX_P(zp)        (&(zp))
#  define GSL_COMPLEX_EQ(z1,z2)    ((z1) == (z2))
#  define GSL_SET_COMPLEX(zp,x,y)  (*(zp) = (x) + I * (y))

#  define GSL_REAL(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[0])

#  define GSL_IMAG(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[1])

#  define GSL_COMPLEX_P_REAL(zp)   GSL_REAL(*(zp))
#  define GSL_COMPLEX_P_IMAG(zp)   GSL_IMAG(*(zp))
#  define GSL_SET_REAL(zp,x)       do { GSL_REAL(*(zp)) = (x); } while(0)
#  define GSL_SET_IMAG(zp,x)       do { GSL_IMAG(*(zp)) = (x); } while(0)

#else /* legacy complex definitions */

/*
 * According to the C17 standard, 6.2.5 paragraph 13:
 *
 * "Each complex type has the same representation and alignment requirements
 * as an array type containing exactly two elements of the corresponding real
 * type; the first element is equal to the real part, and the second element to
 * the imaginary part, of the complex number."
 */

/*#  define GSL_COMPLEX_DEFINE(R, C) typedef R C[2]*/
#  define GSL_COMPLEX_DEFINE(R, C) typedef struct { R dat[2]; } C ;

#  define GSL_REAL(z)              ((z).dat[0])
#  define GSL_IMAG(z)              ((z).dat[1])
#  define GSL_COMPLEX_P(zp)        ((zp)->dat)
#  define GSL_COMPLEX_P_REAL(zp)   ((zp)->dat[0])
#  define GSL_COMPLEX_P_IMAG(zp)   ((zp)->dat[1])
#  define GSL_COMPLEX_EQ(z1,z2)    (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

#  define GSL_SET_COMPLEX(zp,x,y)  do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#  define GSL_SET_REAL(zp,x)       do {(zp)->dat[0]=(x);} while(0)
#  define GSL_SET_IMAG(zp,y)       do {(zp)->dat[1]=(y);} while(0)

#endif

GSL_COMPLEX_DEFINE(double, gsl_complex)
GSL_COMPLEX_DEFINE(long double, gsl_complex_long_double)
GSL_COMPLEX_DEFINE(float, gsl_complex_float)

#define GSL_SET_COMPLEX_PACKED(zp,n,x,y) do {*((zp)+2*(n))=(x); *((zp)+(2*(n)+1))=(y);} while(0)

__END_DECLS

#endif /* __GSL_COMPLEX_H__ */

/* end inlined header: gsl/gsl_complex.h */
/* begin inlined header: gsl/gsl_complex_math.h */
/* complex/gsl_complex_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Jorma Olavi T�htinen, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_MATH_H__
#define __GSL_COMPLEX_MATH_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_complex.h */
/* complex/gsl_complex.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Copyright (C) 2020, 2021 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_H__
#define __GSL_COMPLEX_H__

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


/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

typedef const double *       gsl_const_complex_packed ;
typedef const float *        gsl_const_complex_packed_float  ;
typedef const long double *  gsl_const_complex_packed_long_double ;


/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;

typedef const double *       gsl_const_complex_packed_array ;
typedef const float *        gsl_const_complex_packed_array_float  ;
typedef const long double *  gsl_const_complex_packed_array_long_double ;


/* Yes... this seems weird. Trust us. The point is just that
   sometimes you want to make it obvious that something is
   an output value. The fact that it lacks a 'const' may not
   be enough of a clue for people in some contexts.
 */
typedef double *       gsl_complex_packed_ptr ;
typedef float *        gsl_complex_packed_float_ptr  ;
typedef long double *  gsl_complex_packed_long_double_ptr ;

typedef const double *       gsl_const_complex_packed_ptr ;
typedef const float *        gsl_const_complex_packed_float_ptr  ;
typedef const long double *  gsl_const_complex_packed_long_double_ptr ;

/*
 * If <complex.h> is included, use the C99 complex type.  Otherwise
 * define a type bit-compatible with C99 complex. The GSL_REAL and GSL_IMAG
 * macros require C11 functionality also (_Generic)
 */

/* older gcc compilers claim to be C11 compliant but do not support _Generic */
#if defined(__GNUC__) && (__GNUC__ < 7)
#  define GSL_COMPLEX_LEGACY 1
#endif

#if !defined(GSL_COMPLEX_LEGACY) &&          \
     defined(_Complex_I) &&                  \
     defined(complex) &&                     \
     defined(I) &&                           \
     defined(__STDC__) && (__STDC__ == 1) && \
     defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) /* C11 */

#  define GSL_COMPLEX_DEFINE(R, C) typedef R _Complex C ;

#  define GSL_COMPLEX_P(zp)        (&(zp))
#  define GSL_COMPLEX_EQ(z1,z2)    ((z1) == (z2))
#  define GSL_SET_COMPLEX(zp,x,y)  (*(zp) = (x) + I * (y))

#  define GSL_REAL(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[0])

#  define GSL_IMAG(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[1])

#  define GSL_COMPLEX_P_REAL(zp)   GSL_REAL(*(zp))
#  define GSL_COMPLEX_P_IMAG(zp)   GSL_IMAG(*(zp))
#  define GSL_SET_REAL(zp,x)       do { GSL_REAL(*(zp)) = (x); } while(0)
#  define GSL_SET_IMAG(zp,x)       do { GSL_IMAG(*(zp)) = (x); } while(0)

#else /* legacy complex definitions */

/*
 * According to the C17 standard, 6.2.5 paragraph 13:
 *
 * "Each complex type has the same representation and alignment requirements
 * as an array type containing exactly two elements of the corresponding real
 * type; the first element is equal to the real part, and the second element to
 * the imaginary part, of the complex number."
 */

/*#  define GSL_COMPLEX_DEFINE(R, C) typedef R C[2]*/
#  define GSL_COMPLEX_DEFINE(R, C) typedef struct { R dat[2]; } C ;

#  define GSL_REAL(z)              ((z).dat[0])
#  define GSL_IMAG(z)              ((z).dat[1])
#  define GSL_COMPLEX_P(zp)        ((zp)->dat)
#  define GSL_COMPLEX_P_REAL(zp)   ((zp)->dat[0])
#  define GSL_COMPLEX_P_IMAG(zp)   ((zp)->dat[1])
#  define GSL_COMPLEX_EQ(z1,z2)    (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

#  define GSL_SET_COMPLEX(zp,x,y)  do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#  define GSL_SET_REAL(zp,x)       do {(zp)->dat[0]=(x);} while(0)
#  define GSL_SET_IMAG(zp,y)       do {(zp)->dat[1]=(y);} while(0)

#endif

GSL_COMPLEX_DEFINE(double, gsl_complex)
GSL_COMPLEX_DEFINE(long double, gsl_complex_long_double)
GSL_COMPLEX_DEFINE(float, gsl_complex_float)

#define GSL_SET_COMPLEX_PACKED(zp,n,x,y) do {*((zp)+2*(n))=(x); *((zp)+(2*(n)+1))=(y);} while(0)

__END_DECLS

#endif /* __GSL_COMPLEX_H__ */

/* end inlined header: gsl/gsl_complex.h */
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS           /* empty */
#define __END_DECLS             /* empty */
#endif

__BEGIN_DECLS

/* Complex numbers */

gsl_complex gsl_complex_polar (double r, double theta); /* r= r e^(i theta) */

INLINE_DECL gsl_complex gsl_complex_rect (double x, double y);  /* r= real+i*imag */

#ifdef HAVE_INLINE
INLINE_FUN gsl_complex
gsl_complex_rect (double x, double y)
{                               /* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
#endif

#define GSL_COMPLEX_ONE (gsl_complex_rect(1.0,0.0))
#define GSL_COMPLEX_ZERO (gsl_complex_rect(0.0,0.0))
#define GSL_COMPLEX_NEGONE (gsl_complex_rect(-1.0,0.0))

/* Properties of complex numbers */

double gsl_complex_arg (gsl_complex z); /* return arg(z), -pi< arg(z) <=+pi */
double gsl_complex_abs (gsl_complex z);   /* return |z|   */
double gsl_complex_abs2 (gsl_complex z);  /* return |z|^2 */
double gsl_complex_logabs (gsl_complex z); /* return log|z| */

/* Complex arithmetic operators */

gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b);  /* r=a+b */
gsl_complex gsl_complex_sub (gsl_complex a, gsl_complex b);  /* r=a-b */
gsl_complex gsl_complex_mul (gsl_complex a, gsl_complex b);  /* r=a*b */
gsl_complex gsl_complex_div (gsl_complex a, gsl_complex b);  /* r=a/b */

gsl_complex gsl_complex_add_real (gsl_complex a, double x);  /* r=a+x */
gsl_complex gsl_complex_sub_real (gsl_complex a, double x);  /* r=a-x */
gsl_complex gsl_complex_mul_real (gsl_complex a, double x);  /* r=a*x */
gsl_complex gsl_complex_div_real (gsl_complex a, double x);  /* r=a/x */

gsl_complex gsl_complex_add_imag (gsl_complex a, double y);  /* r=a+iy */
gsl_complex gsl_complex_sub_imag (gsl_complex a, double y);  /* r=a-iy */
gsl_complex gsl_complex_mul_imag (gsl_complex a, double y);  /* r=a*iy */
gsl_complex gsl_complex_div_imag (gsl_complex a, double y);  /* r=a/iy */

gsl_complex gsl_complex_conjugate (gsl_complex z);  /* r=conj(z) */
gsl_complex gsl_complex_inverse (gsl_complex a);    /* r=1/a */
gsl_complex gsl_complex_negative (gsl_complex a);    /* r=-a */

/* Elementary Complex Functions */

gsl_complex gsl_complex_sqrt (gsl_complex z);  /* r=sqrt(z) */
gsl_complex gsl_complex_sqrt_real (double x);  /* r=sqrt(x) (x<0 ok) */

gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b);  /* r=a^b */
gsl_complex gsl_complex_pow_real (gsl_complex a, double b);  /* r=a^b */

gsl_complex gsl_complex_exp (gsl_complex a);    /* r=exp(a) */
gsl_complex gsl_complex_log (gsl_complex a);    /* r=log(a) (base e) */
gsl_complex gsl_complex_log10 (gsl_complex a);  /* r=log10(a) (base 10) */
gsl_complex gsl_complex_log_b (gsl_complex a, gsl_complex b);   /* r=log_b(a) (base=b) */

/* Complex Trigonometric Functions */

gsl_complex gsl_complex_sin (gsl_complex a);  /* r=sin(a) */
gsl_complex gsl_complex_cos (gsl_complex a);  /* r=cos(a) */
gsl_complex gsl_complex_sec (gsl_complex a);  /* r=sec(a) */
gsl_complex gsl_complex_csc (gsl_complex a);  /* r=csc(a) */
gsl_complex gsl_complex_tan (gsl_complex a);  /* r=tan(a) */
gsl_complex gsl_complex_cot (gsl_complex a);  /* r=cot(a) */

/* Inverse Complex Trigonometric Functions */

gsl_complex gsl_complex_arcsin (gsl_complex a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arcsin_real (double a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arccos (gsl_complex a);  /* r=arccos(a) */
gsl_complex gsl_complex_arccos_real (double a);  /* r=arccos(a) */
gsl_complex gsl_complex_arcsec (gsl_complex a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arcsec_real (double a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arccsc (gsl_complex a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arccsc_real (double a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arctan (gsl_complex a);  /* r=arctan(a) */
gsl_complex gsl_complex_arccot (gsl_complex a);  /* r=arccot(a) */

/* Complex Hyperbolic Functions */

gsl_complex gsl_complex_sinh (gsl_complex a);  /* r=sinh(a) */
gsl_complex gsl_complex_cosh (gsl_complex a);  /* r=coshh(a) */
gsl_complex gsl_complex_sech (gsl_complex a);  /* r=sech(a) */
gsl_complex gsl_complex_csch (gsl_complex a);  /* r=csch(a) */
gsl_complex gsl_complex_tanh (gsl_complex a);  /* r=tanh(a) */
gsl_complex gsl_complex_coth (gsl_complex a);  /* r=coth(a) */

/* Inverse Complex Hyperbolic Functions */

gsl_complex gsl_complex_arcsinh (gsl_complex a);  /* r=arcsinh(a) */
gsl_complex gsl_complex_arccosh (gsl_complex a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arccosh_real (double a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arcsech (gsl_complex a);  /* r=arcsech(a) */
gsl_complex gsl_complex_arccsch (gsl_complex a);  /* r=arccsch(a) */
gsl_complex gsl_complex_arctanh (gsl_complex a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arctanh_real (double a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arccoth (gsl_complex a);  /* r=arccoth(a) */

__END_DECLS

#endif /* __GSL_COMPLEX_MATH_H__ */

/* end inlined header: gsl/gsl_complex_math.h */
/**********************************************************************
 * Complex numbers
 **********************************************************************/

gsl_complex
gsl_complex_polar (double r, double theta)
{                               /* return z = r exp(i theta) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, r * cos (theta), r * sin (theta));
  return z;
}

/**********************************************************************
 * Properties of complex numbers
 **********************************************************************/

double
gsl_complex_arg (gsl_complex z)
{                               /* return arg(z),  -pi < arg(z) <= +pi */
  double x = GSL_REAL (z);
  double y = GSL_IMAG (z);

  if (x == 0.0 && y == 0.0)
    {
      return 0;
    }

  return atan2 (y, x);
}

double
gsl_complex_abs (gsl_complex z)
{                               /* return |z| */
  return hypot (GSL_REAL (z), GSL_IMAG (z));
}

double
gsl_complex_abs2 (gsl_complex z)
{                               /* return |z|^2 */
  double x = GSL_REAL (z);
  double y = GSL_IMAG (z);

  return (x * x + y * y);
}

double
gsl_complex_logabs (gsl_complex z)
{                               /* return log|z| */
  double xabs = fabs (GSL_REAL (z));
  double yabs = fabs (GSL_IMAG (z));
  double max, u;

  if (xabs >= yabs)
    {
      max = xabs;
      u = yabs / xabs;
    }
  else
    {
      max = yabs;
      u = xabs / yabs;
    }

  /* Handle underflow when u is close to 0 */

  return log (max) + 0.5 * log1p (u * u);
}


/***********************************************************************
 * Complex arithmetic operators
 ***********************************************************************/

gsl_complex
gsl_complex_add (gsl_complex a, gsl_complex b)
{                               /* z=a+b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar + br, ai + bi);
  return z;
}

gsl_complex
gsl_complex_add_real (gsl_complex a, double x)
{                               /* z=a+x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) + x, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_add_imag (gsl_complex a, double y)
{                               /* z=a+iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), GSL_IMAG (a) + y);
  return z;
}


gsl_complex
gsl_complex_sub (gsl_complex a, gsl_complex b)
{                               /* z=a-b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar - br, ai - bi);
  return z;
}

gsl_complex
gsl_complex_sub_real (gsl_complex a, double x)
{                               /* z=a-x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) - x, GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_sub_imag (gsl_complex a, double y)
{                               /* z=a-iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), GSL_IMAG (a) - y);
  return z;
}

gsl_complex
gsl_complex_mul (gsl_complex a, gsl_complex b)
{                               /* z=a*b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar * br - ai * bi, ar * bi + ai * br);
  return z;
}

gsl_complex
gsl_complex_mul_real (gsl_complex a, double x)
{                               /* z=a*x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x * GSL_REAL (a), x * GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_mul_imag (gsl_complex a, double y)
{                               /* z=a*iy */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -y * GSL_IMAG (a), y * GSL_REAL (a));
  return z;
}

gsl_complex
gsl_complex_div (gsl_complex a, gsl_complex b)
{                               /* z=a/b */
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  double s = 1.0 / gsl_complex_abs (b);

  double sbr = s * br;
  double sbi = s * bi;

  double zr = (ar * sbr + ai * sbi) * s;
  double zi = (ai * sbr - ar * sbi) * s;

  gsl_complex z;
  GSL_SET_COMPLEX (&z, zr, zi);
  return z;
}

gsl_complex
gsl_complex_div_real (gsl_complex a, double x)
{                               /* z=a/x */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a) / x, GSL_IMAG (a) / x);
  return z;
}

gsl_complex
gsl_complex_div_imag (gsl_complex a, double y)
{                               /* z=a/(iy) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_IMAG (a) / y,  - GSL_REAL (a) / y);
  return z;
}

gsl_complex
gsl_complex_conjugate (gsl_complex a)
{                               /* z=conj(a) */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_negative (gsl_complex a)
{                               /* z=-a */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, -GSL_REAL (a), -GSL_IMAG (a));
  return z;
}

gsl_complex
gsl_complex_inverse (gsl_complex a)
{                               /* z=1/a */
  double s = 1.0 / gsl_complex_abs (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, (GSL_REAL (a) * s) * s, -(GSL_IMAG (a) * s) * s);
  return z;
}

/**********************************************************************
 * Elementary complex functions
 **********************************************************************/

gsl_complex
gsl_complex_sqrt (gsl_complex a)
{                               /* z=sqrt(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, 0, 0);
    }
  else
    {
      double x = fabs (GSL_REAL (a));
      double y = fabs (GSL_IMAG (a));
      double w;

      if (x >= y)
        {
          double t = y / x;
          w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + t * t)));
        }
      else
        {
          double t = x / y;
          w = sqrt (y) * sqrt (0.5 * (t + sqrt (1.0 + t * t)));
        }

      if (GSL_REAL (a) >= 0.0)
        {
          double ai = GSL_IMAG (a);
          GSL_SET_COMPLEX (&z, w, ai / (2.0 * w));
        }
      else
        {
          double ai = GSL_IMAG (a);
          double vi = (ai >= 0) ? w : -w;
          GSL_SET_COMPLEX (&z, ai / (2.0 * vi), vi);
        }
    }

  return z;
}

gsl_complex
gsl_complex_sqrt_real (double x)
{                               /* z=sqrt(x) */
  gsl_complex z;

  if (x >= 0)
    {
      GSL_SET_COMPLEX (&z, sqrt (x), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, 0.0, sqrt (-x));
    }

  return z;
}

gsl_complex
gsl_complex_exp (gsl_complex a)
{                               /* z=exp(a) */
  double rho = exp (GSL_REAL (a));
  double theta = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, rho * cos (theta), rho * sin (theta));
  return z;
}

gsl_complex
gsl_complex_pow (gsl_complex a, gsl_complex b)
{                               /* z=a^b */
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0.0)
    {
      if (GSL_REAL (b) == 0 && GSL_IMAG (b) == 0.0)
        {
          GSL_SET_COMPLEX (&z, 1.0, 0.0);
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0.0, 0.0);
        }
    }
  else if (GSL_REAL (b) == 1.0 && GSL_IMAG (b) == 0.0)
    {
      return a;
    }
  else if (GSL_REAL (b) == -1.0 && GSL_IMAG (b) == 0.0)
    {
      return gsl_complex_inverse (a);
    }
  else
    {
      double logr = gsl_complex_logabs (a);
      double theta = gsl_complex_arg (a);

      double br = GSL_REAL (b), bi = GSL_IMAG (b);

      double rho = exp (logr * br - bi * theta);
      double beta = theta * br + bi * logr;

      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_pow_real (gsl_complex a, double b)
{                               /* z=a^b */
  gsl_complex z;

  if (GSL_REAL (a) == 0 && GSL_IMAG (a) == 0)
    {
      if (b == 0)
        {
          GSL_SET_COMPLEX (&z, 1, 0);
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0, 0);
        }
    }
  else
    {
      double logr = gsl_complex_logabs (a);
      double theta = gsl_complex_arg (a);
      double rho = exp (logr * b);
      double beta = theta * b;
      GSL_SET_COMPLEX (&z, rho * cos (beta), rho * sin (beta));
    }

  return z;
}

gsl_complex
gsl_complex_log (gsl_complex a)
{                               /* z=log(a) */
  double logr = gsl_complex_logabs (a);
  double theta = gsl_complex_arg (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, logr, theta);
  return z;
}

gsl_complex
gsl_complex_log10 (gsl_complex a)
{                               /* z = log10(a) */
  return gsl_complex_mul_real (gsl_complex_log (a), 1 / log (10.));
}

gsl_complex
gsl_complex_log_b (gsl_complex a, gsl_complex b)
{
  return gsl_complex_div (gsl_complex_log (a), gsl_complex_log (b));
}

/***********************************************************************
 * Complex trigonometric functions
 ***********************************************************************/

gsl_complex
gsl_complex_sin (gsl_complex a)
{                               /* z = sin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (I == 0.0)
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      GSL_SET_COMPLEX (&z, sin (R), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, sin (R) * cosh (I), cos (R) * sinh (I));
    }

  return z;
}

gsl_complex
gsl_complex_cos (gsl_complex a)
{                               /* z = cos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (I == 0.0)
    {
      /* avoid returing negative zero (-0.0) for the imaginary part  */

      GSL_SET_COMPLEX (&z, cos (R), 0.0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, cos (R) * cosh (I), sin (R) * sinh (-I));
    }

  return z;
}

gsl_complex
gsl_complex_tan (gsl_complex a)
{                               /* z = tan(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (fabs (I) < 1)
    {
      double D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);

      GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) / D, 0.5 * sinh (2 * I) / D);
    }
  else
    {
      double D = pow (cos (R), 2.0) + pow (sinh (I), 2.0);
      double F = 1 + pow(cos (R)/sinh (I), 2.0);

      GSL_SET_COMPLEX (&z, 0.5 * sin (2 * R) / D, 1 / (tanh (I) * F));
    }

  return z;
}

gsl_complex
gsl_complex_sec (gsl_complex a)
{                               /* z = sec(a) */
  gsl_complex z = gsl_complex_cos (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_csc (gsl_complex a)
{                               /* z = csc(a) */
  gsl_complex z = gsl_complex_sin (a);
  return gsl_complex_inverse(z);
}


gsl_complex
gsl_complex_cot (gsl_complex a)
{                               /* z = cot(a) */
  gsl_complex z = gsl_complex_tan (a);
  return gsl_complex_inverse (z);
}

/**********************************************************************
 * Inverse Complex Trigonometric Functions
 **********************************************************************/

gsl_complex
gsl_complex_arcsin (gsl_complex a)
{                               /* z = arcsin(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      z = gsl_complex_arcsin_real (R);
    }
  else
    {
      double x = fabs (R), y = fabs (I);
      double r = hypot (x + 1, y), s = hypot (x - 1, y);
      double A = 0.5 * (r + s);
      double B = x / A;
      double y2 = y * y;

      double real, imag;

      const double A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = asin (B);
        }
      else
        {
          if (x <= 1)
            {
              double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (x / sqrt (D));
            }
          else
            {
              double Apx = A + x;
              double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan (x / (y * sqrt (D)));
            }
        }

      if (A <= A_crossover)
        {
          double Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      GSL_SET_COMPLEX (&z, (R >= 0) ? real : -real, (I >= 0) ? imag : -imag);
    }

  return z;
}

gsl_complex
gsl_complex_arcsin_real (double a)
{                               /* z = arcsin(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (a), 0.0);
    }
  else
    {
      if (a < 0.0)
        {
          GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, M_PI_2, -acosh (a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccos (gsl_complex a)
{                               /* z = arccos(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      z = gsl_complex_arccos_real (R);
    }
  else
    {
      double x = fabs (R), y = fabs (I);
      double r = hypot (x + 1, y), s = hypot (x - 1, y);
      double A = 0.5 * (r + s);
      double B = x / A;
      double y2 = y * y;

      double real, imag;

      const double A_crossover = 1.5, B_crossover = 0.6417;

      if (B <= B_crossover)
        {
          real = acos (B);
        }
      else
        {
          if (x <= 1)
            {
              double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
              real = atan (sqrt (D) / x);
            }
          else
            {
              double Apx = A + x;
              double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
              real = atan ((y * sqrt (D)) / x);
            }
        }

      if (A <= A_crossover)
        {
          double Am1;

          if (x < 1)
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
            }
          else
            {
              Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
            }

          imag = log1p (Am1 + sqrt (Am1 * (A + 1)));
        }
      else
        {
          imag = log (A + sqrt (A * A - 1));
        }

      GSL_SET_COMPLEX (&z, (R >= 0) ? real : M_PI - real, (I >= 0) ? -imag : imag);
    }

  return z;
}

gsl_complex
gsl_complex_arccos_real (double a)
{                               /* z = arccos(a) */
  gsl_complex z;

  if (fabs (a) <= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (a), 0);
    }
  else
    {
      if (a < 0.0)
        {
          GSL_SET_COMPLEX (&z, M_PI, -acosh (-a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0, acosh (a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arctan (gsl_complex a)
{                               /* z = arctan(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);
  gsl_complex z;

  if (I == 0)
    {
      GSL_SET_COMPLEX (&z, atan (R), 0);
    }
  else
    {
      /* FIXME: This is a naive implementation which does not fully
         take into account cancellation errors, overflow, underflow
         etc.  It would benefit from the Hull et al treatment. */

      double r = hypot (R, I);

      double imag;

      double u = 2 * I / (1 + r * r);

      /* FIXME: the following cross-over should be optimized but 0.1
         seems to work ok */

      if (fabs (u) < 0.1)
        {
          imag = 0.25 * (log1p (u) - log1p (-u));
        }
      else
        {
          double A = hypot (R, I + 1);
          double B = hypot (R, I - 1);
          imag = 0.5 * log (A / B);
        }

      if (R == 0)
        {
          if (I > 1)
            {
              GSL_SET_COMPLEX (&z, M_PI_2, imag);
            }
          else if (I < -1)
            {
              GSL_SET_COMPLEX (&z, -M_PI_2, imag);
            }
          else
            {
              GSL_SET_COMPLEX (&z, 0, imag);
            };
        }
      else
        {
          GSL_SET_COMPLEX (&z, 0.5 * atan2 (2 * R, ((1 + r) * (1 - r))), imag);
        }
    }

  return z;
}

gsl_complex
gsl_complex_arcsec (gsl_complex a)
{                               /* z = arcsec(a) */
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arccos (z);
}

gsl_complex
gsl_complex_arcsec_real (double a)
{                               /* z = arcsec(a) */
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, acos (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
        {
          GSL_SET_COMPLEX (&z, 0, acosh (1 / a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, M_PI, -acosh (-1 / a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccsc (gsl_complex a)
{                               /* z = arccsc(a) */
  gsl_complex z = gsl_complex_inverse (a);
  return gsl_complex_arcsin (z);
}

gsl_complex
gsl_complex_arccsc_real (double a)
{                               /* z = arccsc(a) */
  gsl_complex z;

  if (a <= -1.0 || a >= 1.0)
    {
      GSL_SET_COMPLEX (&z, asin (1 / a), 0.0);
    }
  else
    {
      if (a >= 0.0)
        {
          GSL_SET_COMPLEX (&z, M_PI_2, -acosh (1 / a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, -M_PI_2, acosh (-1 / a));
        }
    }

  return z;
}

gsl_complex
gsl_complex_arccot (gsl_complex a)
{                               /* z = arccot(a) */
  gsl_complex z;

  if (GSL_REAL (a) == 0.0 && GSL_IMAG (a) == 0.0)
    {
      GSL_SET_COMPLEX (&z, M_PI_2, 0);
    }
  else
    {
      z = gsl_complex_inverse (a);
      z = gsl_complex_arctan (z);
    }

  return z;
}

/**********************************************************************
 * Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex
gsl_complex_sinh (gsl_complex a)
{                               /* z = sinh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, sinh (R) * cos (I), cosh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_cosh (gsl_complex a)
{                               /* z = cosh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, cosh (R) * cos (I), sinh (R) * sin (I));
  return z;
}

gsl_complex
gsl_complex_tanh (gsl_complex a)
{                               /* z = tanh(a) */
  double R = GSL_REAL (a), I = GSL_IMAG (a);

  gsl_complex z;

  if (fabs(R) < 1.0)
    {
      double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);

      GSL_SET_COMPLEX (&z, sinh (R) * cosh (R) / D, 0.5 * sin (2 * I) / D);
    }
  else
    {
      double D = pow (cos (I), 2.0) + pow (sinh (R), 2.0);
      double F = 1 + pow (cos (I) / sinh (R), 2.0);

      GSL_SET_COMPLEX (&z, 1.0 / (tanh (R) * F), 0.5 * sin (2 * I) / D);
    }

  return z;
}

gsl_complex
gsl_complex_sech (gsl_complex a)
{                               /* z = sech(a) */
  gsl_complex z = gsl_complex_cosh (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_csch (gsl_complex a)
{                               /* z = csch(a) */
  gsl_complex z = gsl_complex_sinh (a);
  return gsl_complex_inverse (z);
}

gsl_complex
gsl_complex_coth (gsl_complex a)
{                               /* z = coth(a) */
  gsl_complex z = gsl_complex_tanh (a);
  return gsl_complex_inverse (z);
}

/**********************************************************************
 * Inverse Complex Hyperbolic Functions
 **********************************************************************/

gsl_complex
gsl_complex_arcsinh (gsl_complex a)
{                               /* z = arcsinh(a) */
  gsl_complex z = gsl_complex_mul_imag(a, 1.0);
  z = gsl_complex_arcsin (z);
  z = gsl_complex_mul_imag (z, -1.0);
  return z;
}

gsl_complex
gsl_complex_arccosh (gsl_complex a)
{                               /* z = arccosh(a) */
  gsl_complex z = gsl_complex_arccos (a);
  z = gsl_complex_mul_imag (z, GSL_IMAG(z) > 0 ? -1.0 : 1.0);
  return z;
}

gsl_complex
gsl_complex_arccosh_real (double a)
{                               /* z = arccosh(a) */
  gsl_complex z;

  if (a >= 1)
    {
      GSL_SET_COMPLEX (&z, acosh (a), 0);
    }
  else
    {
      if (a >= -1.0)
        {
          GSL_SET_COMPLEX (&z, 0, acos (a));
        }
      else
        {
          GSL_SET_COMPLEX (&z, acosh (-a), M_PI);
        }
    }

  return z;
}

gsl_complex
gsl_complex_arctanh (gsl_complex a)
{                               /* z = arctanh(a) */
  if (GSL_IMAG (a) == 0.0)
    {
      return gsl_complex_arctanh_real (GSL_REAL (a));
    }
  else
    {
      gsl_complex z = gsl_complex_mul_imag(a, 1.0);
      z = gsl_complex_arctan (z);
      z = gsl_complex_mul_imag (z, -1.0);
      return z;
    }
}

gsl_complex
gsl_complex_arctanh_real (double a)
{                               /* z = arctanh(a) */
  gsl_complex z;

  if (a > -1.0 && a < 1.0)
    {
      GSL_SET_COMPLEX (&z, atanh (a), 0);
    }
  else
    {
      GSL_SET_COMPLEX (&z, atanh (1 / a), (a < 0) ? M_PI_2 : -M_PI_2);
    }

  return z;
}

gsl_complex
gsl_complex_arcsech (gsl_complex a)
{                               /* z = arcsech(a); */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arccosh (t);
}

gsl_complex
gsl_complex_arccsch (gsl_complex a)
{                               /* z = arccsch(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arcsinh (t);
}

gsl_complex
gsl_complex_arccoth (gsl_complex a)
{                               /* z = arccoth(a) */
  gsl_complex t = gsl_complex_inverse (a);
  return gsl_complex_arctanh (t);
}
/* end inlined source: complex/math.c */

/* Special functions -------------------------------------------------------
 * Inlined in bottom-up dependency order so that each file's callee
 * definitions precede the callers. */
/* begin inlined source: specfunc/log.c */
/* specfunc/log.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_log.h */
/* specfunc/gsl_sf_log.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_LOG_H__
#define __GSL_SF_LOG_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_e(const double x, gsl_sf_result * result);
double gsl_sf_log(const double x);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_abs_e(const double x, gsl_sf_result * result);
double gsl_sf_log_abs(const double x);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx(const double x);


/* Log(1 + x) - x
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx_mx(const double x);

__END_DECLS

#endif /* __GSL_SF_LOG_H__ */

/* end inlined header: gsl/gsl_sf_log.h */
/* inlined: error.h */

/* inlined: chebyshev.h */
/* begin inlined source: specfunc/cheb_eval.c */

#ifndef MARLEY_BUILTIN_GSL_CHEB_EVAL_C
#define MARLEY_BUILTIN_GSL_CHEB_EVAL_C

static inline int
cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  {
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

#endif /* MARLEY_BUILTIN_GSL_CHEB_EVAL_C */
/* end inlined source: specfunc/cheb_eval.c */

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* Chebyshev expansion for log(1 + x(t))/x(t)
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static double lopx_data[21] = {
  2.16647910664395270521272590407,
 -0.28565398551049742084877469679,
  0.01517767255690553732382488171,
 -0.00200215904941415466274422081,
  0.00019211375164056698287947962,
 -0.00002553258886105542567601400,
  2.9004512660400621301999384544e-06,
 -3.8873813517057343800270917900e-07,
  4.7743678729400456026672697926e-08,
 -6.4501969776090319441714445454e-09,
  8.2751976628812389601561347296e-10,
 -1.1260499376492049411710290413e-10,
  1.4844576692270934446023686322e-11,
 -2.0328515972462118942821556033e-12,
  2.7291231220549214896095654769e-13,
 -3.7581977830387938294437434651e-14,
  5.1107345870861673561462339876e-15,
 -7.0722150011433276578323272272e-16,
  9.7089758328248469219003866867e-17,
 -1.3492637457521938883731579510e-17,
  1.8657327910677296608121390705e-18
};
static cheb_series lopx_cs = {
  lopx_data,
  20,
  -1, 1,
  10
};

/* Chebyshev expansion for (log(1 + x(t)) - x(t))/x(t)^2
 *
 * x(t) = (4t-1)/(2(4-t))
 * t(x) = (8x+1)/(2(x+2))
 * -1/2 < x < 1/2
 * -1 < t < 1
 */
static double lopxmx_data[20] = {
 -1.12100231323744103373737274541,
  0.19553462773379386241549597019,
 -0.01467470453808083971825344956,
  0.00166678250474365477643629067,
 -0.00018543356147700369785746902,
  0.00002280154021771635036301071,
 -2.8031253116633521699214134172e-06,
  3.5936568872522162983669541401e-07,
 -4.6241857041062060284381167925e-08,
  6.0822637459403991012451054971e-09,
 -8.0339824424815790302621320732e-10,
  1.0751718277499375044851551587e-10,
 -1.4445310914224613448759230882e-11,
  1.9573912180610336168921438426e-12,
 -2.6614436796793061741564104510e-13,
  3.6402634315269586532158344584e-14,
 -4.9937495922755006545809120531e-15,
  6.8802890218846809524646902703e-16,
 -9.5034129794804273611403251480e-17,
  1.3170135013050997157326965813e-17
};
static cheb_series lopxmx_cs = {
  lopxmx_data,
  19,
  -1, 1,
  9
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_log_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_log_abs_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    result->val = log(fabs(x));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta)
{
  /* CHECK_POINTER(lnr) */
  /* CHECK_POINTER(theta) */

  if(zr != 0.0 || zi != 0.0) {
    const double ax = fabs(zr);
    const double ay = fabs(zi);
    const double min = GSL_MIN(ax, ay);
    const double max = GSL_MAX(ax, ay);
    lnr->val = log(max) + 0.5 * log(1.0 + (min/max)*(min/max));
    lnr->err = 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
    theta->val = atan2(zi, zr);
    theta->err = GSL_DBL_EPSILON * fabs(lnr->val);
    return GSL_SUCCESS;
  }
  else {
    DOMAIN_ERROR_2(lnr, theta);
  }
}


int
gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < GSL_ROOT6_DBL_EPSILON) {
    const double c1 = -0.5;
    const double c2 =  1.0/3.0;
    const double c3 = -1.0/4.0;
    const double c4 =  1.0/5.0;
    const double c5 = -1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c7 = -1.0/8.0;
    const double c8 =  1.0/9.0;
    const double c9 = -1.0/10.0;
    const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
    result->val = x * (1.0 + x*(c1 + x*(c2 + x*(c3 + x*(c4 + x*t)))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    double t = 0.5*(8.0*x + 1.0)/(x+2.0);
    gsl_sf_result c;
    cheb_eval_e(&lopx_cs, t, &c);
    result->val = x * c.val;
    result->err = fabs(x * c.err);
    return GSL_SUCCESS;
  }
  else {
    result->val = log(1.0 + x);
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= -1.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < GSL_ROOT5_DBL_EPSILON) {
    const double c1 = -0.5;
    const double c2 =  1.0/3.0;
    const double c3 = -1.0/4.0;
    const double c4 =  1.0/5.0;
    const double c5 = -1.0/6.0;
    const double c6 =  1.0/7.0;
    const double c7 = -1.0/8.0;
    const double c8 =  1.0/9.0;
    const double c9 = -1.0/10.0;
    const double t  =  c5 + x*(c6 + x*(c7 + x*(c8 + x*c9)));
    result->val = x*x * (c1 + x*(c2 + x*(c3 + x*(c4 + x*t))));
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(fabs(x) < 0.5) {
    double t = 0.5*(8.0*x + 1.0)/(x+2.0);
    gsl_sf_result c;
    cheb_eval_e(&lopxmx_cs, t, &c);
    result->val = x*x * c.val;
    result->err = x*x * c.err;
    return GSL_SUCCESS;
  }
  else {
    const double lterm = log(1.0 + x);
    result->val = lterm - x;
    result->err = GSL_DBL_EPSILON * (fabs(lterm) + fabs(x));
    return GSL_SUCCESS;
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_log(const double x)
{
  EVAL_RESULT(gsl_sf_log_e(x, &result));
}

double gsl_sf_log_abs(const double x)
{
  EVAL_RESULT(gsl_sf_log_abs_e(x, &result));
}

double gsl_sf_log_1plusx(const double x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_e(x, &result));
}

double gsl_sf_log_1plusx_mx(const double x)
{
  EVAL_RESULT(gsl_sf_log_1plusx_mx_e(x, &result));
}
/* end inlined source: specfunc/log.c */
/* begin inlined source: specfunc/trig.c */
/* specfunc/trig.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_log.h */
/* specfunc/gsl_sf_log.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_LOG_H__
#define __GSL_SF_LOG_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_e(const double x, gsl_sf_result * result);
double gsl_sf_log(const double x);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_abs_e(const double x, gsl_sf_result * result);
double gsl_sf_log_abs(const double x);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx(const double x);


/* Log(1 + x) - x
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx_mx(const double x);

__END_DECLS

#endif /* __GSL_SF_LOG_H__ */

/* end inlined header: gsl/gsl_sf_log.h */
/* begin inlined header: gsl/gsl_sf_trig.h */
/* specfunc/gsl_sf_trig.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_TRIG_H__
#define __GSL_SF_TRIG_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Sin(x) with GSL semantics. This is actually important
 * because we want to control the error estimate, and trying
 * to guess the error for the standard library implementation
 * every time it is used would be a little goofy.
 */
int gsl_sf_sin_e(double x, gsl_sf_result * result);
double gsl_sf_sin(const double x);


/* Cos(x) with GSL semantics.
 */
int gsl_sf_cos_e(double x, gsl_sf_result * result);
double gsl_sf_cos(const double x);


/* Hypot(x,y) with GSL semantics.
 */
int gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_hypot(const double x, const double y);


/* Sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_e(const double zr, const double zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* Cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_e(const double zr, const double zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* Log(Sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_e(const double zr, const double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* Sinc(x) = sin(pi x) / (pi x)
 *
 * exceptions: none
 */
int gsl_sf_sinc_e(double x, gsl_sf_result * result);
double gsl_sf_sinc(const double x);


/* Log(Sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnsinh_e(const double x, gsl_sf_result * result);
double gsl_sf_lnsinh(const double x);


/* Log(Cosh(x))
 *
 * exceptions: none
 */
int gsl_sf_lncosh_e(const double x, gsl_sf_result * result);
double gsl_sf_lncosh(const double x);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_polar_to_rect(const double r, const double theta, gsl_sf_result * x, gsl_sf_result * y);

/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_rect_to_polar(const double x, const double y, gsl_sf_result * r, gsl_sf_result * theta);

/* Sin(x) for quantity with an associated error.
 */
int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result);


/* Force an angle to lie in the range (-pi,pi].
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_symm_e(double * theta);
double gsl_sf_angle_restrict_symm(const double theta);


/* Force an angle to lie in the range [0, 2pi)
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_pos_e(double * theta);
double gsl_sf_angle_restrict_pos(const double theta);


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result);

int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_TRIG_H__ */

/* end inlined header: gsl/gsl_sf_trig.h */
/* inlined: error.h */

/* inlined: chebyshev.h */
/* inlined once already: specfunc/cheb_eval.c */

/* sinh(x) series
 * double-precision for |x| < 1.0
 */
inline
static
int
sinh_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 1.0/6.0;
  const double c1 = 1.0/120.0;
  const double c2 = 1.0/5040.0;
  const double c3 = 1.0/362880.0;
  const double c4 = 1.0/39916800.0;
  const double c5 = 1.0/6227020800.0;
  const double c6 = 1.0/1307674368000.0;
  const double c7 = 1.0/355687428096000.0;
  *result = x*(1.0 + y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))))));
  return GSL_SUCCESS;
}


/* cosh(x)-1 series
 * double-precision for |x| < 1.0
 */
inline
static
int
cosh_m1_series(const double x, double * result)
{
  const double y = x*x;
  const double c0 = 0.5;
  const double c1 = 1.0/24.0;
  const double c2 = 1.0/720.0;
  const double c3 = 1.0/40320.0;
  const double c4 = 1.0/3628800.0;
  const double c5 = 1.0/479001600.0;
  const double c6 = 1.0/87178291200.0;
  const double c7 = 1.0/20922789888000.0;
  const double c8 = 1.0/6402373705728000.0;
  *result = y*(c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*c8))))))));
  return GSL_SUCCESS;
}


/* Chebyshev expansion for f(t) = sinc((t+1)/2), -1 < t < 1
 */
static double sinc_data[17] = {
  1.133648177811747875422,
 -0.532677564732557348781,
 -0.068293048346633177859,
  0.033403684226353715020,
  0.001485679893925747818,
 -0.000734421305768455295,
 -0.000016837282388837229,
  0.000008359950146618018,
  0.000000117382095601192,
 -0.000000058413665922724,
 -0.000000000554763755743,
  0.000000000276434190426,
  0.000000000001895374892,
 -0.000000000000945237101,
 -0.000000000000004900690,
  0.000000000000002445383,
  0.000000000000000009925
};
static cheb_series sinc_cs = {
  sinc_data,
  16,
  -1, 1,
  10
};


/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (sin(x)/x - 1)/(x*x)
 */
static double sin_data[12] = {
  -0.3295190160663511504173,
   0.0025374284671667991990,
   0.0006261928782647355874,
  -4.6495547521854042157541e-06,
  -5.6917531549379706526677e-07,
   3.7283335140973803627866e-09,
   3.0267376484747473727186e-10,
  -1.7400875016436622322022e-12,
  -1.0554678305790849834462e-13,
   5.3701981409132410797062e-16,
   2.5984137983099020336115e-17,
  -1.1821555255364833468288e-19
};
static cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  11
};

/* Chebyshev expansion for f(t) = g((t+1)Pi/8), -1<t<1
 * g(x) = (2(cos(x) - 1)/(x^2) + 1) / x^2
 */
static double cos_data[11] = {
  0.165391825637921473505668118136,
 -0.00084852883845000173671196530195,
 -0.000210086507222940730213625768083,
  1.16582269619760204299639757584e-6,
  1.43319375856259870334412701165e-7,
 -7.4770883429007141617951330184e-10,
 -6.0969994944584252706997438007e-11,
  2.90748249201909353949854872638e-13,
  1.77126739876261435667156490461e-14,
 -7.6896421502815579078577263149e-17,
 -3.7363121133079412079201377318e-18
};
static cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  10
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

/* I would have prefered just using the library sin() function.
 * But after some experimentation I decided that there was
 * no good way to understand the error; library sin() is just a black box.
 * So we have to roll our own.
 */
int
gsl_sf_sin_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double sgn_x = GSL_SIGN(x);
    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = x * (1.0 - x2/6.0);
      result->err = fabs(x*x2*x2 / 100.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = sgn_x;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_cos_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = 1.0 - 0.5*x2;
      result->err = fabs(x2*x2/12.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = 1.0;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }

      if(octant > 1) {
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}


int
gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0 && y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double a = fabs(x);
    const double b = fabs(y);
    const double min = GSL_MIN_DBL(a,b);
    const double max = GSL_MAX_DBL(a,b);
    const double rat = min/max;
    const double root_term = sqrt(1.0 + rat*rat);

    if(max < GSL_DBL_MAX/root_term) {
      result->val = max * root_term;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
}


int
gsl_sf_complex_sin_e(const double zr, const double zi,
                        gsl_sf_result * szr, gsl_sf_result * szi)
{
  /* CHECK_POINTER(szr) */
  /* CHECK_POINTER(szi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    szr->val = sin(zr)*(ch_m1 + 1.0);
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    szr->val = sin(zr)*ch;
    szi->val = cos(zr)*sh;
    szr->err = 2.0 * GSL_DBL_EPSILON * fabs(szr->val);
    szi->err = 2.0 * GSL_DBL_EPSILON * fabs(szi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(szr, szi);
  }
}


int
gsl_sf_complex_cos_e(const double zr, const double zi,
                        gsl_sf_result * czr, gsl_sf_result * czi)
{
  /* CHECK_POINTER(czr) */
  /* CHECK_POINTER(czi) */

  if(fabs(zi) < 1.0) {
    double ch_m1, sh;
    sinh_series(zi, &sh);
    cosh_m1_series(zi, &ch_m1);
    czr->val =  cos(zr)*(ch_m1 + 1.0);
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else if(fabs(zi) < GSL_LOG_DBL_MAX) {
    double ex = exp(zi);
    double ch = 0.5*(ex+1.0/ex);
    double sh = 0.5*(ex-1.0/ex);
    czr->val =  cos(zr)*ch;
    czi->val = -sin(zr)*sh;
    czr->err = 2.0 * GSL_DBL_EPSILON * fabs(czr->val);
    czi->err = 2.0 * GSL_DBL_EPSILON * fabs(czi->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR_2(czr,czi);
  }
}


int
gsl_sf_complex_logsin_e(const double zr, const double zi,
                           gsl_sf_result * lszr, gsl_sf_result * lszi)
{
  /* CHECK_POINTER(lszr) */
  /* CHECK_POINTER(lszi) */

  if(zi > 60.0) {
    lszr->val = -M_LN2 + zi;
    lszi->val =  0.5*M_PI - zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else if(zi < -60.0) {
    lszr->val = -M_LN2 - zi;
    lszi->val = -0.5*M_PI + zr;
    lszr->err = 2.0 * GSL_DBL_EPSILON * fabs(lszr->val);
    lszi->err = 2.0 * GSL_DBL_EPSILON * fabs(lszi->val);
  }
  else {
    gsl_sf_result sin_r, sin_i;
    int status;
    gsl_sf_complex_sin_e(zr, zi, &sin_r, &sin_i); /* ok by construction */
    status = gsl_sf_complex_log_e(sin_r.val, sin_i.val, lszr, lszi);
    if(status == GSL_EDOM) {
      DOMAIN_ERROR_2(lszr, lszi);
    }
  }
  return gsl_sf_angle_restrict_symm_e(&(lszi->val));
}


int
gsl_sf_lnsinh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 1.0) {
    double eps;
    sinh_series(x, &eps);
    result->val = log(eps);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = x + log(0.5*(1.0 - exp(-2.0*x)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lncosh_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x) < 1.0) {
    double eps;
    cosh_m1_series(x, &eps);
    return gsl_sf_log_1plusx_e(eps, result);
  }
  else if(fabs(x) < -0.5*GSL_LOG_DBL_EPSILON) {
    result->val = fabs(x) + log(0.5*(1.0 + exp(-2.0*fabs(x))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = -M_LN2 + fabs(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/*
inline int gsl_sf_sincos_e(const double theta, double * s, double * c)
{
  double tan_half = tan(0.5 * theta);
  double den = 1. + tan_half*tan_half;
  double cos_theta = (1.0 - tan_half*tan_half) / den;
  double sin_theta = 2.0 * tan_half / den;
}
*/

int
gsl_sf_polar_to_rect(const double r, const double theta,
                          gsl_sf_result * x, gsl_sf_result * y)
{
  double t   = theta;
  int status = gsl_sf_angle_restrict_symm_e(&t);
  double c = cos(t);
  double s = sin(t);
  x->val = r * cos(t);
  y->val = r * sin(t);
  x->err  = r * fabs(s * GSL_DBL_EPSILON * t);
  x->err += 2.0 * GSL_DBL_EPSILON * fabs(x->val);
  y->err  = r * fabs(c * GSL_DBL_EPSILON * t);
  y->err += 2.0 * GSL_DBL_EPSILON * fabs(y->val);
  return status;
}


int
gsl_sf_rect_to_polar(const double x, const double y,
                          gsl_sf_result * r, gsl_sf_result * theta)
{
  int stat_h = gsl_sf_hypot_e(x, y, r);
  if(r->val > 0.0) {
    theta->val = atan2(y, x);
    theta->err = 2.0 * GSL_DBL_EPSILON * fabs(theta->val);
    return stat_h;
  }
  else {
    DOMAIN_ERROR(theta);
  }
}


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.8539812564849853515625e-01;
  const double P2 = 4 * 3.7748947079307981766760e-08;
  const double P3 = 4 * 2.6951514290790594840552e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double y = GSL_SIGN(theta) * 2 * floor(fabs(theta)/TwoPi);
  double r = ((theta - y*P1) - y*P2) - y*P3;

  if(r >  M_PI) { r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
  else if (r < -M_PI) r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */

  result->val = r;

  if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
    result->val = GSL_NAN;
    result->err = GSL_NAN;
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val - theta);
    return GSL_SUCCESS;
  }
  else {
    double delta = fabs(result->val - theta);
    result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
    return GSL_SUCCESS;
  }
}


int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result)
{
  /* synthetic extended precision constants */
  const double P1 = 4 * 7.85398125648498535156e-01;
  const double P2 = 4 * 3.77489470793079817668e-08;
  const double P3 = 4 * 2.69515142907905952645e-15;
  const double TwoPi = 2*(P1 + P2 + P3);

  const double y = 2*floor(theta/TwoPi);

  double r = ((theta - y*P1) - y*P2) - y*P3;

  if(r > TwoPi) {r = (((r-2*P1)-2*P2)-2*P3); }  /* r-TwoPi */
  else if (r < 0) { /* may happen due to FP rounding */
    r = (((r+2*P1)+2*P2)+2*P3); /* r+TwoPi */
  }

  result->val = r;

  if(fabs(theta) > 0.0625/GSL_DBL_EPSILON) {
    result->val = GSL_NAN;
    result->err = fabs(result->val);
    GSL_ERROR ("error", GSL_ELOSS);
  }
  else if(fabs(theta) > 0.0625/GSL_SQRT_DBL_EPSILON) {
    result->err = GSL_DBL_EPSILON * fabs(result->val - theta);
    return GSL_SUCCESS;
  }
  else {
    double delta = fabs(result->val - theta);
    result->err = 2.0 * GSL_DBL_EPSILON * ((delta < M_PI) ? delta : M_PI);
    return GSL_SUCCESS;
  }
}


int gsl_sf_angle_restrict_symm_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_symm_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_angle_restrict_pos_e(double * theta)
{
  gsl_sf_result r;
  int stat = gsl_sf_angle_restrict_pos_err_e(*theta, &r);
  *theta = r.val;
  return stat;
}


int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_s = gsl_sf_sin_e(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}


int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_c = gsl_sf_cos_e(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}


#if 0
int
gsl_sf_sin_pi_x_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(-100.0 < x && x < 100.0) {
    result->val = sin(M_PI * x) / (M_PI * x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double N = floor(x + 0.5);
    const double f = x - N;

    if(N < INT_MAX && N > INT_MIN) {
      /* Make it an integer if we can. Saves another
       * call to floor().
       */
      const int intN    = (int)N;
      const double sign = ( GSL_IS_ODD(intN) ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI * f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }
    else if(N > 2.0/GSL_DBL_EPSILON || N < -2.0/GSL_DBL_EPSILON) {
      /* All integer-valued floating point numbers
       * bigger than 2/eps=2^53 are actually even.
       */
      result->val = 0.0;
      result->err = 0.0;
    }
    else {
      const double resN = N - 2.0*floor(0.5*N); /* 0 for even N, 1 for odd N */
      const double sign = ( fabs(resN) > 0.5 ? -1.0 : 1.0 );
      result->val = sign * sin(M_PI*f);
      result->err = GSL_DBL_EPSILON * fabs(result->val);
    }

    return GSL_SUCCESS;
  }
}
#endif


int gsl_sf_sinc_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double ax = fabs(x);

    if(ax < 0.8) {
      /* Do not go to the limit of the fit since
       * there is a zero there and the Chebyshev
       * accuracy will go to zero.
       */
      return cheb_eval_e(&sinc_cs, 2.0*ax-1.0, result);
    }
    else if(ax < 100.0) {
      /* Small arguments are no problem.
       * We trust the library sin() to
       * roughly machine precision.
       */
      result->val = sin(M_PI * ax)/(M_PI * ax);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Large arguments must be handled separately.
       */
      const double r = M_PI*ax;
      gsl_sf_result s;
      int stat_s = gsl_sf_sin_e(r, &s);
      result->val = s.val/r;
      result->err = s.err/r + 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_s;
    }
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_sin(const double x)
{
  EVAL_RESULT(gsl_sf_sin_e(x, &result));
}

double gsl_sf_cos(const double x)
{
  EVAL_RESULT(gsl_sf_cos_e(x, &result));
}

double gsl_sf_hypot(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_hypot_e(x, y, &result));
}

double gsl_sf_lnsinh(const double x)
{
  EVAL_RESULT(gsl_sf_lnsinh_e(x, &result));
}

double gsl_sf_lncosh(const double x)
{
  EVAL_RESULT(gsl_sf_lncosh_e(x, &result));
}

double gsl_sf_angle_restrict_symm(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_symm_e(&result));
}

double gsl_sf_angle_restrict_pos(const double theta)
{
  double result = theta;
  EVAL_DOUBLE(gsl_sf_angle_restrict_pos_e(&result));
}

#if 0
double gsl_sf_sin_pi_x(const double x)
{
  EVAL_RESULT(gsl_sf_sin_pi_x_e(x, &result));
}
#endif

double gsl_sf_sinc(const double x)
{
  EVAL_RESULT(gsl_sf_sinc_e(x, &result));
}
/* end inlined source: specfunc/trig.c */
/* begin inlined source: specfunc/exp.c */
/* specfunc/exp.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
#include <stdlib.h>
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_gamma.h */
/* specfunc/gsl_sf_gamma.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_GAMMA_H__
#define __GSL_SF_GAMMA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
double gsl_sf_lngamma(const double x);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
double gsl_sf_gamma(const double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
double gsl_sf_gammastar(const double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
double gsl_sf_gammainv(const double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_taylorcoeff(const int n, const double x);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_fact(const unsigned int n);


/* n!! = n(n-2)(n-4) ...
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_doublefact(const unsigned int n);


/* log(n!)
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lnfact(const unsigned int n);


/* log(n!!)
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lndoublefact(const unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochhammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_lnpoch(const double a, const double x);


/* Logarithm of Pochhammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_e(const double a, const double x, gsl_sf_result * result, double * sgn);


/* Pochhammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_poch(const double a, const double x);


/* Relative Pochhammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_pochrel(const double a, const double x);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a >= 0, x >= 0
 *   Q(a,0) := 1
 *   Q(0,x) := 0, x != 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_Q(const double a, const double x);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_P(const double a, const double x);


/* Non-normalized Incomplete Gamma Function
 *
 * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * x >= 0.0
 *   Gamma(a, 0) := Gamma(a)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc(const double a, const double x);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_lnbeta(const double a, const double b);

int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_beta(const double a, const double b);


/* Normalized Incomplete Beta Function
 * B_x(a,b)/B(a,b)
 *
 * a > 0, b > 0, 0 <= x <= 1
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_beta_inc(const double a, const double b, const double x);


/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297

__END_DECLS

#endif /* __GSL_SF_GAMMA_H__ */

/* end inlined header: gsl/gsl_sf_gamma.h */
/* begin inlined header: gsl/gsl_sf_exp.h */
/* specfunc/gsl_sf_exp.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_EXP_H__
#define __GSL_SF_EXP_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
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

/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e(const double x, gsl_sf_result * result);
double gsl_sf_exp(const double x);


/* Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_exp_mult(const double x, const double y);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
double gsl_sf_expm1(const double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
double gsl_sf_exprel(const double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
double gsl_sf_exprel_2(const double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_exprel_n(const int n, const double x);

int gsl_sf_exprel_n_CF_e(const double n, const double x, gsl_sf_result * result);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);

/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);

__END_DECLS

#endif /* __GSL_SF_EXP_H__ */

/* end inlined header: gsl/gsl_sf_exp.h */
/* inlined: error.h */

/* Evaluate the continued fraction for exprel.
 * [Abramowitz+Stegun, 4.2.41]
 */
static
int
exprel_n_CF(const double N, const double x, gsl_sf_result * result)
{
  const double RECUR_BIG = GSL_SQRT_DBL_MAX;
  const int maxiter = 5000;
  int n = 1;
  double Anm2 = 1.0;
  double Bnm2 = 0.0;
  double Anm1 = 0.0;
  double Bnm1 = 1.0;
  double a1 = 1.0;
  double b1 = 1.0;
  double a2 = -x;
  double b2 = N+1;
  double an, bn;

  double fn;

  double An = b1*Anm1 + a1*Anm2;   /* A1 */
  double Bn = b1*Bnm1 + a1*Bnm2;   /* B1 */

  /* One explicit step, before we get to the main pattern. */
  n++;
  Anm2 = Anm1;
  Bnm2 = Bnm1;
  Anm1 = An;
  Bnm1 = Bn;
  An = b2*Anm1 + a2*Anm2;   /* A2 */
  Bn = b2*Bnm1 + a2*Bnm2;   /* B2 */

  fn = An/Bn;

  while(n < maxiter) {
    double old_fn;
    double del;
    n++;
    Anm2 = Anm1;
    Bnm2 = Bnm1;
    Anm1 = An;
    Bnm1 = Bn;
    an = ( GSL_IS_ODD(n) ? ((n-1)/2)*x : -(N+(n/2)-1)*x );
    bn = N + n - 1;
    An = bn*Anm1 + an*Anm2;
    Bn = bn*Bnm1 + an*Bnm2;

    if(fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG) {
      An /= RECUR_BIG;
      Bn /= RECUR_BIG;
      Anm1 /= RECUR_BIG;
      Bnm1 /= RECUR_BIG;
      Anm2 /= RECUR_BIG;
      Bnm2 /= RECUR_BIG;
    }

    old_fn = fn;
    fn = An/Bn;
    del = old_fn/fn;

    if(fabs(del - 1.0) < 2.0*GSL_DBL_EPSILON) break;
  }

  result->val = fn;
  result->err = 4.0*(n+1.0)*GSL_DBL_EPSILON*fabs(fn);

  if(n == maxiter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_exp_e(const double x, gsl_sf_result * result)
{
  if(x > GSL_LOG_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else if(x < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else {
    result->val = exp(x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}

int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result)
{
  if(x > INT_MAX-1) {
    OVERFLOW_ERROR_E10(result);
  }
  else if(x < INT_MIN+1) {
    UNDERFLOW_ERROR_E10(result);
  }
  else {
    const int N = (x > GSL_LOG_DBL_MAX || x < GSL_LOG_DBL_MIN) ? (int) floor(x/M_LN10) : 0;
    result->val = exp(x-N*M_LN10);
    result->err = 2.0 * (fabs(x)+1.0) * GSL_DBL_EPSILON * fabs(result->val);
    result->e10 = N;
    return GSL_SUCCESS;
  }
}


int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    const double ex = exp(x);
    result->val = y * ex;
    result->err = (2.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double lnr = x + ly;

    if(lnr > GSL_LOG_DBL_MAX - 0.01) {
      OVERFLOW_ERROR(result);
    }
    else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
      UNDERFLOW_ERROR(result);
    }
    else {
      const double sy   = GSL_SIGN(y);
      const double M    = floor(x);
      const double N    = floor(ly);
      const double a    = x  - M;
      const double b    = ly - N;
      const double berr = 2.0 * GSL_DBL_EPSILON * (fabs(ly) + fabs(N));
      result->val  = sy * exp(M+N) * exp(a+b);
      result->err  = berr * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * (M + N + 1.0) * fabs(result->val);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    const double ex = exp(x);
    result->val = y * ex;
    result->err = (2.0 + fabs(x)) * GSL_DBL_EPSILON * fabs(result->val);
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double l10_val = (x + ly)/M_LN10;

    if(l10_val > INT_MAX-1) {
      OVERFLOW_ERROR_E10(result);
    }
    else if(l10_val < INT_MIN+1) {
      UNDERFLOW_ERROR_E10(result);
    }
    else {
      const double sy  = GSL_SIGN(y);
      const int    N   = (int) floor(l10_val);
      const double arg_val = (l10_val - N) * M_LN10;
      const double arg_err = 2.0 * GSL_DBL_EPSILON * (fabs(x) + fabs(ly) + M_LN10*abs(N));

      result->val  = sy * exp(arg_val);
      result->err  = arg_err * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      result->e10 = N;

      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_exp_mult_err_e(const double x, const double dx,
                             const double y, const double dy,
                             gsl_sf_result * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = fabs(dy * exp(x));
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    double ex = exp(x);
    result->val  = y * ex;
    result->err  = ex * (fabs(dy) + fabs(y*dx));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double lnr = x + ly;

    if(lnr > GSL_LOG_DBL_MAX - 0.01) {
      OVERFLOW_ERROR(result);
    }
    else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
      UNDERFLOW_ERROR(result);
    }
    else {
      const double sy  = GSL_SIGN(y);
      const double M   = floor(x);
      const double N   = floor(ly);
      const double a   = x  - M;
      const double b   = ly - N;
      const double eMN = exp(M+N);
      const double eab = exp(a+b);
      result->val  = sy * eMN * eab;
      result->err  = eMN * eab * 2.0*GSL_DBL_EPSILON;
      result->err += eMN * eab * fabs(dy/y);
      result->err += eMN * eab * fabs(dx);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_exp_mult_err_e10_e(const double x, const double dx,
                             const double y, const double dy,
                             gsl_sf_result_e10 * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = fabs(dy * exp(x));
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    const double ex = exp(x);
    result->val  = y * ex;
    result->err  = ex * (fabs(dy) + fabs(y*dx));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    result->e10 = 0;
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double l10_val = (x + ly)/M_LN10;

    if(l10_val > INT_MAX-1) {
      OVERFLOW_ERROR_E10(result);
    }
    else if(l10_val < INT_MIN+1) {
      UNDERFLOW_ERROR_E10(result);
    }
    else {
      const double sy  = GSL_SIGN(y);
      const int    N   = (int) floor(l10_val);
      const double arg_val = (l10_val - N) * M_LN10;
      const double arg_err = dy/fabs(y) + dx + 2.0*GSL_DBL_EPSILON*fabs(arg_val);

      result->val  = sy * exp(arg_val);
      result->err  = arg_err * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      result->e10 = N;

      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_expm1_e(const double x, gsl_sf_result * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    result->val = -1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    result->val = exp(x) - 1.0;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    result->val = x * (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_LOG_DBL_MAX) {
    result->val = exp(x) - 1.0;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int gsl_sf_exprel_e(const double x, gsl_sf_result * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    result->val = -1.0/x;
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    result->val = (exp(x) - 1.0)/x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    result->val = (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_LOG_DBL_MAX) {
    result->val = (exp(x) - 1.0)/x;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int gsl_sf_exprel_2_e(double x, gsl_sf_result * result)
{
  const double cut = 0.002;

  if(x < GSL_LOG_DBL_MIN) {
    result->val = -2.0/x*(1.0 + 1.0/x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < -cut) {
    result->val = 2.0*(exp(x) - 1.0 - x)/(x*x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < cut) {
    result->val = (1.0 + 1.0/3.0*x*(1.0 + 0.25*x*(1.0 + 0.2*x*(1.0 + 1.0/6.0*x))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_LOG_DBL_MAX) {
    result->val = 2.0*(exp(x) - 1.0 - x)/(x*x);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int
gsl_sf_exprel_n_CF_e(const double N, const double x, gsl_sf_result * result)
{
  return exprel_n_CF(N, x, result);
}

int
gsl_sf_exprel_n_e(const int N, const double x, gsl_sf_result * result)
{
  if(N < 0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(fabs(x) < GSL_ROOT3_DBL_EPSILON * N) {
    result->val = 1.0 + x/(N+1) * (1.0 + x/(N+2));
    result->err = 2.0 * GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(N == 0) {
    return gsl_sf_exp_e(x, result);
  }
  else if(N == 1) {
    return gsl_sf_exprel_e(x, result);
  }
  else if(N == 2) {
    return gsl_sf_exprel_2_e(x, result);
  }
  else {
    if(x > N && (-x + N*(1.0 + log(x/N)) < GSL_LOG_DBL_EPSILON)) {
      /* x is much larger than n.
       * Ignore polynomial part, so
       * exprel_N(x) ~= e^x N!/x^N
       */
      gsl_sf_result lnf_N;
      double lnr_val;
      double lnr_err;
      double lnterm;
      gsl_sf_lnfact_e(N, &lnf_N);
      lnterm = N*log(x);
      lnr_val  = x + lnf_N.val - lnterm;
      lnr_err  = GSL_DBL_EPSILON * (fabs(x) + fabs(lnf_N.val) + fabs(lnterm));
      lnr_err += lnf_N.err;
      return gsl_sf_exp_err_e(lnr_val, lnr_err, result);
    }
    else if(x > N) {
      /* Write the identity
       *   exprel_n(x) = e^x n! / x^n (1 - Gamma[n,x]/Gamma[n])
       * then use the asymptotic expansion
       * Gamma[n,x] ~ x^(n-1) e^(-x) (1 + (n-1)/x + (n-1)(n-2)/x^2 + ...)
       */
      double ln_x = log(x);
      gsl_sf_result lnf_N;
      double lg_N;
      double lnpre_val;
      double lnpre_err;
      gsl_sf_lnfact_e(N, &lnf_N);    /* log(N!)       */
      lg_N  = lnf_N.val - log(N);       /* log(Gamma(N)) */
      lnpre_val  = x + lnf_N.val - N*ln_x;
      lnpre_err  = GSL_DBL_EPSILON * (fabs(x) + fabs(lnf_N.val) + fabs(N*ln_x));
      lnpre_err += lnf_N.err;
      if(lnpre_val < GSL_LOG_DBL_MAX - 5.0) {
        int stat_eG;
        gsl_sf_result bigG_ratio;
        gsl_sf_result pre;
        int stat_ex = gsl_sf_exp_err_e(lnpre_val, lnpre_err, &pre);
        double ln_bigG_ratio_pre = -x + (N-1)*ln_x - lg_N;
        double bigGsum = 1.0;
        double term = 1.0;
        int k;
        for(k=1; k<N; k++) {
          term *= (N-k)/x;
          bigGsum += term;
        }
        stat_eG = gsl_sf_exp_mult_e(ln_bigG_ratio_pre, bigGsum, &bigG_ratio);
        if(stat_eG == GSL_SUCCESS) {
          result->val  = pre.val * (1.0 - bigG_ratio.val);
          result->err  = pre.val * (2.0*GSL_DBL_EPSILON + bigG_ratio.err);
          result->err += pre.err * fabs(1.0 - bigG_ratio.val);
          result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
          return stat_ex;
        }
        else {
          result->val = 0.0;
          result->err = 0.0;
          return stat_eG;
        }
      }
      else {
        OVERFLOW_ERROR(result);
      }
    }
    else if(x > -10.0*N) {
      return exprel_n_CF(N, x, result);
    }
    else {
      /* x -> -Inf asymptotic:
       * exprel_n(x) ~ e^x n!/x^n - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
       *             ~ - n/x (1 + (n-1)/x + (n-1)(n-2)/x + ...)
       */
      double sum  = 1.0;
      double term = 1.0;
      int k;
      for(k=1; k<N; k++) {
        term *= (N-k)/x;
        sum  += term;
      }
      result->val = -N/x * sum;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result)
{
  const double adx = fabs(dx);

  /* CHECK_POINTER(result) */

  if(x + adx > GSL_LOG_DBL_MAX) {
    OVERFLOW_ERROR(result);
  }
  else if(x - adx < GSL_LOG_DBL_MIN) {
    UNDERFLOW_ERROR(result);
  }
  else {
    const double ex  = exp(x);
    const double edx = exp(adx);
    result->val  = ex;
    result->err  = ex * GSL_MAX_DBL(GSL_DBL_EPSILON, edx - 1.0/edx);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int
gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result)
{
  const double adx = fabs(dx);

  /* CHECK_POINTER(result) */

  if(x + adx > INT_MAX - 1) {
    OVERFLOW_ERROR_E10(result);
  }
  else if(x - adx < INT_MIN + 1) {
    UNDERFLOW_ERROR_E10(result);
  }
  else {
    const int    N  = (int)floor(x/M_LN10);
    const double ex = exp(x-N*M_LN10);
    result->val = ex;
    result->err = ex * (2.0 * GSL_DBL_EPSILON * (fabs(x) + 1.0) + adx);
    result->e10 = N;
    return GSL_SUCCESS;
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_exp(const double x)
{
  EVAL_RESULT(gsl_sf_exp_e(x, &result));
}

double gsl_sf_exp_mult(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_exp_mult_e(x, y, &result));
}

double gsl_sf_expm1(const double x)
{
  EVAL_RESULT(gsl_sf_expm1_e(x, &result));
}

double gsl_sf_exprel(const double x)
{
  EVAL_RESULT(gsl_sf_exprel_e(x, &result));
}

double gsl_sf_exprel_2(const double x)
{
  EVAL_RESULT(gsl_sf_exprel_2_e(x, &result));
}

double gsl_sf_exprel_n(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_exprel_n_e(n, x, &result));
}
/* end inlined source: specfunc/exp.c */
/* begin inlined source: specfunc/gamma.c */
/* specfunc/gamma.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_exp.h */
/* specfunc/gsl_sf_exp.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_EXP_H__
#define __GSL_SF_EXP_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
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

/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e(const double x, gsl_sf_result * result);
double gsl_sf_exp(const double x);


/* Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_exp_mult(const double x, const double y);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
double gsl_sf_expm1(const double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
double gsl_sf_exprel(const double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
double gsl_sf_exprel_2(const double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_exprel_n(const int n, const double x);

int gsl_sf_exprel_n_CF_e(const double n, const double x, gsl_sf_result * result);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);

/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);

__END_DECLS

#endif /* __GSL_SF_EXP_H__ */

/* end inlined header: gsl/gsl_sf_exp.h */
/* begin inlined header: gsl/gsl_sf_log.h */
/* specfunc/gsl_sf_log.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_LOG_H__
#define __GSL_SF_LOG_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Provide a logarithm function with GSL semantics.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_e(const double x, gsl_sf_result * result);
double gsl_sf_log(const double x);


/* Log(|x|)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_abs_e(const double x, gsl_sf_result * result);
double gsl_sf_log_abs(const double x);


/* Complex Logarithm
 *   exp(lnr + I theta) = zr + I zi
 * Returns argument in [-pi,pi].
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_log_e(const double zr, const double zi, gsl_sf_result * lnr, gsl_sf_result * theta);


/* Log(1 + x)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx(const double x);


/* Log(1 + x) - x
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_log_1plusx_mx_e(const double x, gsl_sf_result * result);
double gsl_sf_log_1plusx_mx(const double x);

__END_DECLS

#endif /* __GSL_SF_LOG_H__ */

/* end inlined header: gsl/gsl_sf_log.h */
/* begin inlined header: gsl/gsl_sf_psi.h */
/* specfunc/gsl_sf_psi.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_PSI_H__
#define __GSL_SF_PSI_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Poly-Gamma Functions
 *
 * psi(m,x) := (d/dx)^m psi(0,x) = (d/dx)^{m+1} log(gamma(x))
 */


/* Di-Gamma Function  psi(n) = psi(0,n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_int(const int n);


/* Di-Gamma Function psi(x) = psi(0, x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi(const double x);


/* Di-Gamma Function Re[psi(1 + I y)]
 *
 * exceptions: none
 */
int     gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result);
double  gsl_sf_psi_1piy(const double y);


/* Di-Gamma Function psi(z) for general complex argument z = x + iy
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  );


/* Tri-Gamma Function psi^(1)(n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_1_int(const int n);


/* Tri-Gamma Function psi^(1)(x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_1_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi_1(const double x);


/* Poly-Gamma Function psi^(n)(x)
 *
 * n >= 0, x > 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result);
double  gsl_sf_psi_n(const int n, const double x);


__END_DECLS

#endif /* __GSL_SF_PSI_H__ */

/* end inlined header: gsl/gsl_sf_psi.h */
/* begin inlined header: gsl/gsl_sf_trig.h */
/* specfunc/gsl_sf_trig.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_TRIG_H__
#define __GSL_SF_TRIG_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Sin(x) with GSL semantics. This is actually important
 * because we want to control the error estimate, and trying
 * to guess the error for the standard library implementation
 * every time it is used would be a little goofy.
 */
int gsl_sf_sin_e(double x, gsl_sf_result * result);
double gsl_sf_sin(const double x);


/* Cos(x) with GSL semantics.
 */
int gsl_sf_cos_e(double x, gsl_sf_result * result);
double gsl_sf_cos(const double x);


/* Hypot(x,y) with GSL semantics.
 */
int gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_hypot(const double x, const double y);


/* Sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_e(const double zr, const double zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* Cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_e(const double zr, const double zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* Log(Sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_e(const double zr, const double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* Sinc(x) = sin(pi x) / (pi x)
 *
 * exceptions: none
 */
int gsl_sf_sinc_e(double x, gsl_sf_result * result);
double gsl_sf_sinc(const double x);


/* Log(Sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnsinh_e(const double x, gsl_sf_result * result);
double gsl_sf_lnsinh(const double x);


/* Log(Cosh(x))
 *
 * exceptions: none
 */
int gsl_sf_lncosh_e(const double x, gsl_sf_result * result);
double gsl_sf_lncosh(const double x);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_polar_to_rect(const double r, const double theta, gsl_sf_result * x, gsl_sf_result * y);

/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_rect_to_polar(const double x, const double y, gsl_sf_result * r, gsl_sf_result * theta);

/* Sin(x) for quantity with an associated error.
 */
int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result);


/* Force an angle to lie in the range (-pi,pi].
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_symm_e(double * theta);
double gsl_sf_angle_restrict_symm(const double theta);


/* Force an angle to lie in the range [0, 2pi)
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_pos_e(double * theta);
double gsl_sf_angle_restrict_pos(const double theta);


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result);

int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_TRIG_H__ */

/* end inlined header: gsl/gsl_sf_trig.h */
/* begin inlined header: gsl/gsl_sf_gamma.h */
/* specfunc/gsl_sf_gamma.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_GAMMA_H__
#define __GSL_SF_GAMMA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
double gsl_sf_lngamma(const double x);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
double gsl_sf_gamma(const double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
double gsl_sf_gammastar(const double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
double gsl_sf_gammainv(const double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_taylorcoeff(const int n, const double x);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_fact(const unsigned int n);


/* n!! = n(n-2)(n-4) ...
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_doublefact(const unsigned int n);


/* log(n!)
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lnfact(const unsigned int n);


/* log(n!!)
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lndoublefact(const unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochhammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_lnpoch(const double a, const double x);


/* Logarithm of Pochhammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_e(const double a, const double x, gsl_sf_result * result, double * sgn);


/* Pochhammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_poch(const double a, const double x);


/* Relative Pochhammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_pochrel(const double a, const double x);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a >= 0, x >= 0
 *   Q(a,0) := 1
 *   Q(0,x) := 0, x != 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_Q(const double a, const double x);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_P(const double a, const double x);


/* Non-normalized Incomplete Gamma Function
 *
 * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * x >= 0.0
 *   Gamma(a, 0) := Gamma(a)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc(const double a, const double x);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_lnbeta(const double a, const double b);

int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_beta(const double a, const double b);


/* Normalized Incomplete Beta Function
 * B_x(a,b)/B(a,b)
 *
 * a > 0, b > 0, 0 <= x <= 1
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_beta_inc(const double a, const double b, const double x);


/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297

__END_DECLS

#endif /* __GSL_SF_GAMMA_H__ */

/* end inlined header: gsl/gsl_sf_gamma.h */
/* inlined: error.h */
/* inlined: check.h */

/* inlined: chebyshev.h */
/* inlined once already: specfunc/cheb_eval.c */

#define LogRootTwoPi_  0.9189385332046727418


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

static struct {int n; double f; long i; } fact_table[GSL_SF_FACT_NMAX + 1] = {
    { 0,  1.0,     1L     },
    { 1,  1.0,     1L     },
    { 2,  2.0,     2L     },
    { 3,  6.0,     6L     },
    { 4,  24.0,    24L    },
    { 5,  120.0,   120L   },
    { 6,  720.0,   720L   },
    { 7,  5040.0,  5040L  },
    { 8,  40320.0, 40320L },

    { 9,  362880.0,     362880L    },
    { 10, 3628800.0,    3628800L   },
    { 11, 39916800.0,   39916800L  },
    { 12, 479001600.0,  479001600L },

    { 13, 6227020800.0,                               0 },
    { 14, 87178291200.0,                              0 },
    { 15, 1307674368000.0,                            0 },
    { 16, 20922789888000.0,                           0 },
    { 17, 355687428096000.0,                          0 },
    { 18, 6402373705728000.0,                         0 },
    { 19, 121645100408832000.0,                       0 },
    { 20, 2432902008176640000.0,                      0 },
    { 21, 51090942171709440000.0,                     0 },
    { 22, 1124000727777607680000.0,                   0 },
    { 23, 25852016738884976640000.0,                  0 },
    { 24, 620448401733239439360000.0,                 0 },
    { 25, 15511210043330985984000000.0,               0 },
    { 26, 403291461126605635584000000.0,              0 },
    { 27, 10888869450418352160768000000.0,            0 },
    { 28, 304888344611713860501504000000.0,           0 },
    { 29, 8841761993739701954543616000000.0,          0 },
    { 30, 265252859812191058636308480000000.0,        0 },
    { 31, 8222838654177922817725562880000000.0,       0 },
    { 32, 263130836933693530167218012160000000.0,     0 },
    { 33, 8683317618811886495518194401280000000.0,    0 },
    { 34, 2.95232799039604140847618609644e38,  0 },
    { 35, 1.03331479663861449296666513375e40,  0 },
    { 36, 3.71993326789901217467999448151e41,  0 },
    { 37, 1.37637530912263450463159795816e43,  0 },
    { 38, 5.23022617466601111760007224100e44,  0 },
    { 39, 2.03978820811974433586402817399e46,  0 },
    { 40, 8.15915283247897734345611269600e47,  0 },
    { 41, 3.34525266131638071081700620534e49,  0 },
    { 42, 1.40500611775287989854314260624e51,  0 },
    { 43, 6.04152630633738356373551320685e52,  0 },
    { 44, 2.65827157478844876804362581101e54,  0 },
    { 45, 1.19622220865480194561963161496e56,  0 },
    { 46, 5.50262215981208894985030542880e57,  0 },
    { 47, 2.58623241511168180642964355154e59,  0 },
    { 48, 1.24139155925360726708622890474e61,  0 },
    { 49, 6.08281864034267560872252163321e62,  0 },
    { 50, 3.04140932017133780436126081661e64,  0 },
    { 51, 1.55111875328738228022424301647e66,  0 },
    { 52, 8.06581751709438785716606368564e67,  0 },
    { 53, 4.27488328406002556429801375339e69,  0 },
    { 54, 2.30843697339241380472092742683e71,  0 },
    { 55, 1.26964033536582759259651008476e73,  0 },
    { 56, 7.10998587804863451854045647464e74,  0 },
    { 57, 4.05269195048772167556806019054e76,  0 },
    { 58, 2.35056133128287857182947491052e78,  0 },
    { 59, 1.38683118545689835737939019720e80,  0 },
    { 60, 8.32098711274139014427634118320e81,  0 },
    { 61, 5.07580213877224798800856812177e83,  0 },
    { 62, 3.14699732603879375256531223550e85,  0 },
    { 63, 1.982608315404440064116146708360e87,  0 },
    { 64, 1.268869321858841641034333893350e89,  0 },
    { 65, 8.247650592082470666723170306800e90,  0 },
    { 66, 5.443449390774430640037292402480e92,  0 },
    { 67, 3.647111091818868528824985909660e94,  0 },
    { 68, 2.480035542436830599600990418570e96,  0 },
    { 69, 1.711224524281413113724683388810e98,  0 },
    { 70, 1.197857166996989179607278372170e100,  0 },
    { 71, 8.504785885678623175211676442400e101,  0 },
    { 72, 6.123445837688608686152407038530e103,  0 },
    { 73, 4.470115461512684340891257138130e105,  0 },
    { 74, 3.307885441519386412259530282210e107,  0 },
    { 75, 2.480914081139539809194647711660e109,  0 },
    { 76, 1.885494701666050254987932260860e111,  0 },
    { 77, 1.451830920282858696340707840860e113,  0 },
    { 78, 1.132428117820629783145752115870e115,  0 },
    { 79, 8.946182130782975286851441715400e116,  0 },
    { 80, 7.156945704626380229481153372320e118,  0 },
    { 81, 5.797126020747367985879734231580e120,  0 },
    { 82, 4.753643337012841748421382069890e122,  0 },
    { 83, 3.945523969720658651189747118010e124,  0 },
    { 84, 3.314240134565353266999387579130e126,  0 },
    { 85, 2.817104114380550276949479442260e128,  0 },
    { 86, 2.422709538367273238176552320340e130,  0 },
    { 87, 2.107757298379527717213600518700e132,  0 },
    { 88, 1.854826422573984391147968456460e134,  0 },
    { 89, 1.650795516090846108121691926250e136,  0 },
    { 90, 1.485715964481761497309522733620e138,  0 },
    { 91, 1.352001527678402962551665687590e140,  0 },
    { 92, 1.243841405464130725547532432590e142,  0 },
    { 93, 1.156772507081641574759205162310e144,  0 },
    { 94, 1.087366156656743080273652852570e146,  0 },
    { 95, 1.032997848823905926259970209940e148,  0 },
    { 96, 9.916779348709496892095714015400e149,  0 },
    { 97, 9.619275968248211985332842594960e151,  0 },
    { 98, 9.426890448883247745626185743100e153,  0 },
    { 99, 9.332621544394415268169923885600e155,  0 },
    { 100, 9.33262154439441526816992388563e157,  0 },
    { 101, 9.42594775983835942085162312450e159,  0 },
    { 102, 9.61446671503512660926865558700e161,  0 },
    { 103, 9.90290071648618040754671525458e163,  0 },
    { 104, 1.02990167451456276238485838648e166,  0 },
    { 105, 1.08139675824029090050410130580e168,  0 },
    { 106, 1.146280563734708354534347384148e170,  0 },
    { 107, 1.226520203196137939351751701040e172,  0 },
    { 108, 1.324641819451828974499891837120e174,  0 },
    { 109, 1.443859583202493582204882102460e176,  0 },
    { 110, 1.588245541522742940425370312710e178,  0 },
    { 111, 1.762952551090244663872161047110e180,  0 },
    { 112, 1.974506857221074023536820372760e182,  0 },
    { 113, 2.231192748659813646596607021220e184,  0 },
    { 114, 2.543559733472187557120132004190e186,  0 },
    { 115, 2.925093693493015690688151804820e188,  0 },
    { 116, 3.393108684451898201198256093590e190,  0 },
    { 117, 3.96993716080872089540195962950e192,  0 },
    { 118, 4.68452584975429065657431236281e194,  0 },
    { 119, 5.57458576120760588132343171174e196,  0 },
    { 120, 6.68950291344912705758811805409e198,  0 },
    { 121, 8.09429852527344373968162284545e200,  0 },
    { 122, 9.87504420083360136241157987140e202,  0 },
    { 123, 1.21463043670253296757662432419e205,  0 },
    { 124, 1.50614174151114087979501416199e207,  0 },
    { 125, 1.88267717688892609974376770249e209,  0 },
    { 126, 2.37217324288004688567714730514e211,  0 },
    { 127, 3.01266001845765954480997707753e213,  0 },
    { 128, 3.85620482362580421735677065923e215,  0 },
    { 129, 4.97450422247728744039023415041e217,  0 },
    { 130, 6.46685548922047367250730439554e219,  0 },
    { 131, 8.47158069087882051098456875820e221,  0 },
    { 132, 1.11824865119600430744996307608e224,  0 },
    { 133, 1.48727070609068572890845089118e226,  0 },
    { 134, 1.99294274616151887673732419418e228,  0 },
    { 135, 2.69047270731805048359538766215e230,  0 },
    { 136, 3.65904288195254865768972722052e232,  0 },
    { 137, 5.01288874827499166103492629211e234,  0 },
    { 138, 6.91778647261948849222819828311e236,  0 },
    { 139, 9.61572319694108900419719561353e238,  0 },
    { 140, 1.34620124757175246058760738589e241,  0 },
    { 141, 1.89814375907617096942852641411e243,  0 },
    { 142, 2.69536413788816277658850750804e245,  0 },
    { 143, 3.85437071718007277052156573649e247,  0 },
    { 144, 5.55029383273930478955105466055e249,  0 },
    { 145, 8.04792605747199194484902925780e251,  0 },
    { 146, 1.17499720439091082394795827164e254,  0 },
    { 147, 1.72724589045463891120349865931e256,  0 },
    { 148, 2.55632391787286558858117801578e258,  0 },
    { 149, 3.80892263763056972698595524351e260,  0 },
    { 150, 5.71338395644585459047893286526e262,  0 },
    { 151, 8.62720977423324043162318862650e264,  0 },
    { 152, 1.31133588568345254560672467123e267,  0 },
    { 153, 2.00634390509568239477828874699e269,  0 },
    { 154, 3.08976961384735088795856467036e271,  0 },
    { 155, 4.78914290146339387633577523906e273,  0 },
    { 156, 7.47106292628289444708380937294e275,  0 },
    { 157, 1.17295687942641442819215807155e278,  0 },
    { 158, 1.85327186949373479654360975305e280,  0 },
    { 159, 2.94670227249503832650433950735e282,  0 },
    { 160, 4.71472363599206132240694321176e284,  0 },
    { 161, 7.59070505394721872907517857094e286,  0 },
    { 162, 1.22969421873944943411017892849e289,  0 },
    { 163, 2.00440157654530257759959165344e291,  0 },
    { 164, 3.28721858553429622726333031164e293,  0 },
    { 165, 5.42391066613158877498449501421e295,  0 },
    { 166, 9.00369170577843736647426172359e297,  0 },
    { 167, 1.50361651486499904020120170784e300,  0 },
    { 168, 2.52607574497319838753801886917e302,  0 },
    { 169, 4.26906800900470527493925188890e304,  0 },
    { 170, 7.25741561530799896739672821113e306,  0 },

    /*
    { 171, 1.24101807021766782342484052410e309,  0 },
    { 172, 2.13455108077438865629072570146e311,  0 },
    { 173, 3.69277336973969237538295546352e313,  0 },
    { 174, 6.42542566334706473316634250653e315,  0 },
    { 175, 1.12444949108573632830410993864e318,  0 },
    { 176, 1.97903110431089593781523349201e320,  0 },
    { 177, 3.50288505463028580993296328086e322,  0 },
    { 178, 6.23513539724190874168067463993e324,  0 },
    { 179, 1.11608923610630166476084076055e327,  0 },
    { 180, 2.00896062499134299656951336898e329,  0 },
    { 181, 3.63621873123433082379081919786e331,  0 },
    { 182, 6.61791809084648209929929094011e333,  0 },
    { 183, 1.21107901062490622417177024204e336,  0 },
    { 184, 2.22838537954982745247605724535e338,  0 },
    { 185, 4.12251295216718078708070590390e340,  0 },
    { 186, 7.66787409103095626397011298130e342,  0 },
    { 187, 1.43389245502278882136241112750e345,  0 },
    { 188, 2.69571781544284298416133291969e347,  0 },
    { 189, 5.09490667118697324006491921822e349,  0 },
    { 190, 9.68032267525524915612334651460e351,  0 },
    { 191, 1.84894163097375258881955918429e354,  0 },
    { 192, 3.54996793146960497053355363384e356,  0 },
    { 193, 6.85143810773633759312975851330e358,  0 },
    { 194, 1.32917899290084949306717315158e361,  0 },
    { 195, 2.59189903615665651148098764559e363,  0 },
    { 196, 5.08012211086704676250273578535e365,  0 },
    { 197, 1.00078405584080821221303894971e368,  0 },
    { 198, 1.98155243056480026018181712043e370,  0 },
    { 199, 3.94328933682395251776181606966e372,  0 },
    { 200, 7.88657867364790503552363213932e374,  0 }
    */
};

static struct {int n; double f; long i; } doub_fact_table[GSL_SF_DOUBLEFACT_NMAX + 1] = {
  { 0,  1.000000000000000000000000000,    1L    },
  { 1,  1.000000000000000000000000000,    1L    },
  { 2,  2.000000000000000000000000000,    2L    },
  { 3,  3.000000000000000000000000000,    3L    },
  { 4,  8.000000000000000000000000000,    8L    },
  { 5,  15.00000000000000000000000000,    15L   },
  { 6,  48.00000000000000000000000000,    48L   },
  { 7,  105.0000000000000000000000000,    105L  },
  { 8,  384.0000000000000000000000000,    384L  },
  { 9,  945.0000000000000000000000000,    945L  },
  { 10, 3840.000000000000000000000000,    3840L   },
  { 11, 10395.00000000000000000000000,    10395L  },
  { 12, 46080.00000000000000000000000,       46080L       },
  { 13, 135135.0000000000000000000000,       135135L      },
  { 14, 645120.00000000000000000000000,      645120L      },
  { 15, 2.02702500000000000000000000000e6,   2027025L     },
  { 16, 1.03219200000000000000000000000e7,   10321920L    },
  { 17, 3.4459425000000000000000000000e7,    34459425L    },
  { 18, 1.85794560000000000000000000000e8,   185794560L   },
  { 19, 6.5472907500000000000000000000e8,            0 },
  { 20, 3.7158912000000000000000000000e9,            0 },
  { 21, 1.37493105750000000000000000000e10,          0 },
  { 22, 8.1749606400000000000000000000e10,           0 },
  { 23, 3.1623414322500000000000000000e11,           0 },
  { 24, 1.96199055360000000000000000000e12,          0 },
  { 25, 7.9058535806250000000000000000e12,           0 },
  { 26, 5.1011754393600000000000000000e13,           0 },
  { 27, 2.13458046676875000000000000000e14,          0 },
  { 28, 1.42832912302080000000000000000e15,          0 },
  { 29, 6.1902833536293750000000000000e15,           0 },
  { 30, 4.2849873690624000000000000000e16,           0 },
  { 31, 1.91898783962510625000000000000e17,          0 },
  { 32, 1.37119595809996800000000000000e18,          0 },
  { 33, 6.3326598707628506250000000000e18,           0 },
  { 34, 4.6620662575398912000000000000e19,           0 },
  { 35, 2.21643095476699771875000000000e20,          0 },
  { 36, 1.67834385271436083200000000000e21,          0 },
  { 37, 8.2007945326378915593750000000e21,           0 },
  { 38, 6.3777066403145711616000000000e22,           0 },
  { 39, 3.1983098677287777081562500000e23,           0 },
  { 40, 2.55108265612582846464000000000e24,          0 },
  { 41, 1.31130704576879886034406250000e25,          0 },
  { 42, 1.07145471557284795514880000000e26,          0 },
  { 43, 5.6386202968058350994794687500e26,           0 },
  { 44, 4.7144007485205310026547200000e27,           0 },
  { 45, 2.53737913356262579476576093750e28,          0 },
  { 46, 2.16862434431944426122117120000e29,          0 },
  { 47, 1.19256819277443412353990764062e30,          0 },
  { 48, 1.04093968527333324538616217600e31,          0 },
  { 49, 5.8435841445947272053455474391e31,           0 },
  { 50, 5.2046984263666662269308108800e32,           0 },
  { 51, 2.98022791374331087472622919392e33,          0 },
  { 52, 2.70644318171066643800402165760e34,          0 },
  { 53, 1.57952079428395476360490147278e35,          0 },
  { 54, 1.46147931812375987652217169510e36,          0 },
  { 55, 8.6873643685617511998269581003e36,           0 },
  { 56, 8.1842841814930553085241614926e37,           0 },
  { 57, 4.9517976900801981839013661172e38,           0 },
  { 58, 4.7468848252659720789440136657e39,           0 },
  { 59, 2.92156063714731692850180600912e40,       0 },
  { 60, 2.84813089515958324736640819942e41,       0 },
  { 61, 1.78215198865986332638610166557e42,       0 },
  { 62, 1.76584115499894161336717308364e43,       0 },
  { 63, 1.12275575285571389562324404931e44,       0 },
  { 64, 1.13013833919932263255499077353e45,       0 },
  { 65, 7.2979123935621403215510863205e45,        0 },
  { 66, 7.4589130387155293748629391053e46,        0 },
  { 67, 4.8896013036866340154392278347e47,        0 },
  { 68, 5.0720608663265599749067985916e48,        0 },
  { 69, 3.3738248995437774706530672060e49,        0 },
  { 70, 3.5504426064285919824347590141e50,        0 },
  { 71, 2.39541567867608200416367771623e51,       0 },
  { 72, 2.55631867662858622735302649017e52,       0 },
  { 73, 1.74865344543353986303948473285e53,       0 },
  { 74, 1.89167582070515380824123960272e54,       0 },
  { 75, 1.31149008407515489727961354964e55,       0 },
  { 76, 1.43767362373591689426334209807e56,       0 },
  { 77, 1.00984736473786927090530243322e57,       0 },
  { 78, 1.12138542651401517752540683649e58,       0 },
  { 79, 7.9777941814291672401518892225e58,        0 },
  { 80, 8.9710834121121214202032546920e59,        0 },
  { 81, 6.4620132869576254645230302702e60,        0 },
  { 82, 7.3562883979319395645666688474e61,        0 },
  { 83, 5.3634710281748291355541151243e62,        0 },
  { 84, 6.1792822542628292342360018318e63,        0 },
  { 85, 4.5589503739486047652209978556e64,        0 },
  { 86, 5.3141827386660331414429615754e65,        0 },
  { 87, 3.9662868253352861457422681344e66,        0 },
  { 88, 4.6764808100261091644698061863e67,        0 },
  { 89, 3.5299952745484046697106186396e68,        0 },
  { 90, 4.2088327290234982480228255677e69,        0 },
  { 91, 3.2122956998390482494366629620e70,        0 },
  { 92, 3.8721261107016183881809995223e71,        0 },
  { 93, 2.98743500085031487197609655470e72,       0 },
  { 94, 3.6397985440595212848901395509e73,        0 },
  { 95, 2.83806325080779912837729172696e74,       0 },
  { 96, 3.4942066022971404334945339689e75,        0 },
  { 97, 2.75292135328356515452597297515e76,       0 },
  { 98, 3.4243224702511976248246432895e77,        0 },
  { 99, 2.72539213975072950298071324540e78,       0 },
  { 100, 3.4243224702511976248246432895e79,       0 },
  { 101, 2.75264606114823679801052037785e80,      0 },
  { 102, 3.4928089196562215773211361553e81,       0 },
  { 103, 2.83522544298268390195083598919e82,      0 },
  { 104, 3.6325212764424704404139816015e83,       0 },
  { 105, 2.97698671513181809704837778865e84,      0 },
  { 106, 3.8504725530290186668388204976e85,       0 },
  { 107, 3.1853757851910453638417642339e86,       0 },
  { 108, 4.1585103572713401601859261374e87,       0 },
  { 109, 3.4720596058582394465875230149e88,       0 },
  { 110, 4.5743613929984741762045187512e89,       0 },
  { 111, 3.8539861625026457857121505465e90,       0 },
  { 112, 5.1232847601582910773490610013e91,       0 },
  { 113, 4.3550043636279897378547301176e92,       0 },
  { 114, 5.8405446265804518281779295415e93,       0 },
  { 115, 5.0082550181721881985329396352e94,       0 },
  { 116, 6.7750317668333241206863982681e95,       0 },
  { 117, 5.8596583712614601922835393732e96,       0 },
  { 118, 7.9945374848633224624099499564e97,       0 },
  { 119, 6.9729934618011376288174118541e98,       0 },
  { 120, 9.5934449818359869548919399477e99,       0 },
  { 121, 8.4373220887793765308690683435e100,      0 },
  { 122, 1.17040028778399040849681667362e102,       0 },
  { 123, 1.03779061691986331329689540625e103,       0 },
  { 124, 1.45129635685214810653605267528e104,       0 },
  { 125, 1.29723827114982914162111925781e105,       0 },
  { 126, 1.82863340963370661423542637086e106,       0 },
  { 127, 1.64749260436028300985882145742e107,       0 },
  { 128, 2.34065076433114446622134575470e108,       0 },
  { 129, 2.12526545962476508271787968008e109,       0 },
  { 130, 3.04284599363048780608774948111e110,       0 },
  { 131, 2.78409775210844225836042238090e111,       0 },
  { 132, 4.0165567115922439040358293151e112,        0 },
  { 133, 3.7028500103042282036193617666e113,        0 },
  { 134, 5.3821859935336068314080112822e114,        0 },
  { 135, 4.9988475139107080748861383849e115,        0 },
  { 136, 7.3197729512057052907148953438e116,        0 },
  { 137, 6.8484210940576700625940095873e117,        0 },
  { 138, 1.01012866726638733011865555744e119,       0 },
  { 139, 9.5193053207401613870056733264e119,        0 },
  { 140, 1.41418013417294226216611778042e121,       0 },
  { 141, 1.34222205022436275556779993902e122,       0 },
  { 142, 2.00813579052557801227588724819e123,       0 },
  { 143, 1.91937753182083874046195391280e124,       0 },
  { 144, 2.89171553835683233767727763739e125,       0 },
  { 145, 2.78309742114021617366983317355e126,       0 },
  { 146, 4.2219046860009752130088253506e127,        0 },
  { 147, 4.0911532090761177752946547651e128,        0 },
  { 148, 6.2484189352814433152530615189e129,        0 },
  { 149, 6.0958182815234154851890356000e130,        0 },
  { 150, 9.3726284029221649728795922783e131,        0 },
  { 151, 9.2046856051003573826354437561e132,        0 },
  { 152, 1.42463951724416907587769802630e134,       0 },
  { 153, 1.40831689758035467954322289468e135,       0 },
  { 154, 2.19394485655602037685165496051e136,       0 },
  { 155, 2.18289119124954975329199548675e137,       0 },
  { 156, 3.4225539762273917878885817384e138,        0 },
  { 157, 3.4271391702617931126684329142e139,        0 },
  { 158, 5.4076352824392790248639591467e140,        0 },
  { 159, 5.4491512807162510491428083336e141,        0 },
  { 160, 8.6522164519028464397823346347e142,        0 },
  { 161, 8.7731335619531641891199214170e143,        0 },
  { 162, 1.40165906520826112324473821082e145,       0 },
  { 163, 1.43002077059836576282654719098e146,       0 },
  { 164, 2.29872086694154824212137066574e147,       0 },
  { 165, 2.35953427148730350866380286512e148,       0 },
  { 166, 3.8158766391229700819214753051e149,        0 },
  { 167, 3.9404222333837968594685507847e150,        0 },
  { 168, 6.4106727537265897376280785126e151,        0 },
  { 169, 6.6593135744186166925018508262e152,        0 },
  { 170, 1.08981436813352025539677334714e154,       0 },
  { 171, 1.13874262122558345441781649128e155,       0 },
  { 172, 1.87448071318965483928245015709e156,       0 },
  { 173, 1.97002473472025937614282252992e157,       0 },
  { 174, 3.2615964409499994203514632733e158,        0 },
  { 175, 3.4475432857604539082499394274e159,        0 },
  { 176, 5.7404097360719989798185753611e160,        0 },
  { 177, 6.1021516157960034176023927864e161,        0 },
  { 178, 1.02179293302081581840770641427e163,       0 },
  { 179, 1.09228513922748461175082830877e164,       0 },
  { 180, 1.83922727943746847313387154568e165,       0 },
  { 181, 1.97703610200174714726899923887e166,       0 },
  { 182, 3.3473936485761926211036462131e167,        0 },
  { 183, 3.6179760666631972795022686071e168,        0 },
  { 184, 6.1592043133801944228307090322e169,        0 },
  { 185, 6.6932557233269149670791969232e170,        0 },
  { 186, 1.14561200228871616264651187999e172,       0 },
  { 187, 1.25163882026213309884380982464e173,       0 },
  { 188, 2.15375056430278638577544233437e174,       0 },
  { 189, 2.36559737029543155681480056857e175,       0 },
  { 190, 4.0921260721752941329733404353e176,        0 },
  { 191, 4.5182909772642742735162690860e177,        0 },
  { 192, 7.8568820585765647353088136358e178,        0 },
  { 193, 8.7203015861200493478863993359e179,        0 },
  { 194, 1.52423511936385355864990984535e181,       0 },
  { 195, 1.70045880929340962283784787050e182,       0 },
  { 196, 2.98750083395315297495382329688e183,       0 },
  { 197, 3.3499038543080169569905603049e184,        0 },
  { 198, 5.9152516512272428904085701278e185,        0 },
  { 199, 6.6663086700729537444112150067e186,        0 },
  { 200, 1.18305033024544857808171402556e188,       0 },
  { 201, 1.33992804268466370262665421635e189,       0 },
  { 202, 2.38976166709580612772506233164e190,       0 },
  { 203, 2.72005392664986731633210805920e191,       0 },
  { 204, 4.8751138008754445005591271565e192,        0 },
  { 205, 5.5761105496322279984808215214e193,        0 },
  { 206, 1.00427344298034156711518019425e195,       0 },
  { 207, 1.15425488377387119568553005492e196,       0 },
  { 208, 2.08888876139911045959957480403e197,       0 },
  { 209, 2.41239270708739079898275781478e198,       0 },
  { 210, 4.3866663989381319651591070885e199,        0 },
  { 211, 5.0901486119543945858536189892e200,        0 },
  { 212, 9.2997327657488397661373070276e201,        0 },
  { 213, 1.08420165434628604678682084470e203,       0 },
  { 214, 1.99014281187025170995338370390e204,       0 },
  { 215, 2.33103355684451500059166481610e205,       0 },
  { 216, 4.2987084736397436934993088004e206,        0 },
  { 217, 5.0583428183525975512839126509e207,        0 },
  { 218, 9.3711844725346412518284931849e208,        0 },
  { 219, 1.10777707721921886373117687056e210,       0 },
  { 220, 2.06166058395762107540226850068e211,       0 },
  { 221, 2.44818734065447368884590088393e212,       0 },
  { 222, 4.5768864963859187873930360715e213,        0 },
  { 223, 5.4594577696594763261263589712e214,        0 },
  { 224, 1.02522257519044580837604008002e216,       0 },
  { 225, 1.22837799817338217337843076851e217,       0 },
  { 226, 2.31700301993040752692985058084e218,       0 },
  { 227, 2.78841805585357753356903784452e219,       0 },
  { 228, 5.2827668854413291614000593243e220,        0 },
  { 229, 6.3854773479046925518730966640e221,        0 },
  { 230, 1.21503638365150570712201364459e223,       0 },
  { 231, 1.47504526736598397948268532937e224,       0 },
  { 232, 2.81888441007149324052307165546e225,       0 },
  { 233, 3.4368554729627426721946568174e226,        0 },
  { 234, 6.5961895195672941828239876738e227,        0 },
  { 235, 8.0766103614624452796574435210e228,        0 },
  { 236, 1.55670072661788142714646109101e230,       0 },
  { 237, 1.91415665566659953127881411447e231,       0 },
  { 238, 3.7049477293505577966085773966e232,        0 },
  { 239, 4.5748344070431728797563657336e233,        0 },
  { 240, 8.8918745504413387118605857518e234,        0 },
  { 241, 1.10253509209740466402128414180e236,       0 },
  { 242, 2.15183364120680396827026175195e237,       0 },
  { 243, 2.67916027379669333357172046456e238,       0 },
  { 244, 5.2504740845446016825794386748e239,        0 },
  { 245, 6.5639426708018986672507151382e240,        0 },
  { 246, 1.29161662479797201391454191399e242,       0 },
  { 247, 1.62129383968806897081092663913e243,       0 },
  { 248, 3.2032092294989705945080639467e244,        0 },
  { 249, 4.0370216608232917373192073314e245,        0 },
  { 250, 8.0080230737474264862701598667e246,        0 },
  { 251, 1.01329243686664622606712104019e248,       0 },
  { 252, 2.01802181458435147454008028642e249,       0 },
  { 253, 2.56362986527261495194981623168e250,       0 },
  { 254, 5.1257754090442527453318039275e251,        0 },
  { 255, 6.5372561564451681274720313908e252,        0 },
  { 256, 1.31219850471532870280494180544e254,       0 },
  { 257, 1.68007483220640820876031206743e255,       0 },
  { 258, 3.3854721421655480532367498580e256,        0 },
  { 259, 4.3513938154145972606892082546e257,        0 },
  { 260, 8.8022275696304249384155496309e258,        0 },
  { 261, 1.13571378582320988503988335446e260,       0 },
  { 262, 2.30618362324317133386487400329e261,       0 },
  { 263, 2.98692725671504199765489322224e262,       0 },
  { 264, 6.0883247653619723214032673687e263,        0 },
  { 265, 7.9153572302948612937854670389e264,        0 },
  { 266, 1.61949438758628463749326912007e266,       0 },
  { 267, 2.11340038048872796544071969939e267,       0 },
  { 268, 4.3402449587312428284819612418e268,        0 },
  { 269, 5.6850470235146782270355359914e269,        0 },
  { 270, 1.17186613885743556369012953528e271,       0 },
  { 271, 1.54064774337247779952663025366e272,       0 },
  { 272, 3.1874758976922247332371523360e273,        0 },
  { 273, 4.2059683394068643927077005925e274,        0 },
  { 274, 8.7336839596766957690697974006e275,        0 },
  { 275, 1.15664129333688770799461766294e277,       0 },
  { 276, 2.41049677287076803226326408256e278,       0 },
  { 277, 3.2038963825431789511450909263e279,        0 },
  { 278, 6.7011810285807351296918741495e280,        0 },
  { 279, 8.9388709072954692736948036845e281,        0 },
  { 280, 1.87633068800260583631372476186e283,       0 },
  { 281, 2.51182272495002686590823983534e284,       0 },
  { 282, 5.2912525401673484584047038284e285,        0 },
  { 283, 7.1084583116085760305203187340e286,        0 },
  { 284, 1.50271572140752696218693588728e288,       0 },
  { 285, 2.02591061880844416869829083919e289,       0 },
  { 286, 4.2977669632255271118546366376e290,        0 },
  { 287, 5.8143634759802347641640947085e291,        0 },
  { 288, 1.23775688540895180821413535163e293,       0 },
  { 289, 1.68035104455828784684342337075e294,       0 },
  { 290, 3.5894949676859602438209925197e295,        0 },
  { 291, 4.8898215396646176343143620089e296,        0 },
  { 292, 1.04813253056430039119572981576e298,       0 },
  { 293, 1.43271771112173296685410806860e299,       0 },
  { 294, 3.08150963985904315011544565835e300,       0 },
  { 295, 4.2265172478091122522196188024e301,        0 },
  { 296, 9.1212685339827677243417191487e302,        0 },
  { 297, 1.25527562259930633890922678431e304,       0 },
  /*
  { 298, 2.71813802312686478185383230631e305,       0 },
  { 299, 3.7532741115719259533385880851e306,        0 },
  { 300, 8.1544140693805943455614969189e307,  }
  */
};


/* Chebyshev coefficients for Gamma*(3/4(t+1)+1/2), -1<t<1
 */
static double gstar_a_data[30] = {
  2.16786447866463034423060819465,
 -0.05533249018745584258035832802,
  0.01800392431460719960888319748,
 -0.00580919269468937714480019814,
  0.00186523689488400339978881560,
 -0.00059746524113955531852595159,
  0.00019125169907783353925426722,
 -0.00006124996546944685735909697,
  0.00001963889633130842586440945,
 -6.3067741254637180272515795142e-06,
  2.0288698405861392526872789863e-06,
 -6.5384896660838465981983750582e-07,
  2.1108698058908865476480734911e-07,
 -6.8260714912274941677892994580e-08,
  2.2108560875880560555583978510e-08,
 -7.1710331930255456643627187187e-09,
  2.3290892983985406754602564745e-09,
 -7.5740371598505586754890405359e-10,
  2.4658267222594334398525312084e-10,
 -8.0362243171659883803428749516e-11,
  2.6215616826341594653521346229e-11,
 -8.5596155025948750540420068109e-12,
  2.7970831499487963614315315444e-12,
 -9.1471771211886202805502562414e-13,
  2.9934720198063397094916415927e-13,
 -9.8026575909753445931073620469e-14,
  3.2116773667767153777571410671e-14,
 -1.0518035333878147029650507254e-14,
  3.4144405720185253938994854173e-15,
 -1.0115153943081187052322643819e-15
};
static cheb_series gstar_a_cs = {
  gstar_a_data,
  29,
  -1, 1,
  17
};


/* Chebyshev coefficients for
 * x^2(Gamma*(x) - 1 - 1/(12x)), x = 4(t+1)+2, -1 < t < 1
 */
static double gstar_b_data[] = {
  0.0057502277273114339831606096782,
  0.0004496689534965685038254147807,
 -0.0001672763153188717308905047405,
  0.0000615137014913154794776670946,
 -0.0000223726551711525016380862195,
  8.0507405356647954540694800545e-06,
 -2.8671077107583395569766746448e-06,
  1.0106727053742747568362254106e-06,
 -3.5265558477595061262310873482e-07,
  1.2179216046419401193247254591e-07,
 -4.1619640180795366971160162267e-08,
  1.4066283500795206892487241294e-08,
 -4.6982570380537099016106141654e-09,
  1.5491248664620612686423108936e-09,
 -5.0340936319394885789686867772e-10,
  1.6084448673736032249959475006e-10,
 -5.0349733196835456497619787559e-11,
  1.5357154939762136997591808461e-11,
 -4.5233809655775649997667176224e-12,
  1.2664429179254447281068538964e-12,
 -3.2648287937449326771785041692e-13,
  7.1528272726086133795579071407e-14,
 -9.4831735252566034505739531258e-15,
 -2.3124001991413207293120906691e-15,
  2.8406613277170391482590129474e-15,
 -1.7245370321618816421281770927e-15,
  8.6507923128671112154695006592e-16,
 -3.9506563665427555895391869919e-16,
  1.6779342132074761078792361165e-16,
 -6.0483153034414765129837716260e-17
};
static cheb_series gstar_b_cs = {
  gstar_b_data,
  29,
  -1, 1,
  18
};


/* coefficients for gamma=7, kmax=8  Lanczos method */
static double lanczos_7_c[9] = {
  0.99999999999980993227684700473478,
  676.520368121885098567009190444019,
 -1259.13921672240287047156078755283,
  771.3234287776530788486528258894,
 -176.61502916214059906584551354,
  12.507343278686904814458936853,
 -0.13857109526572011689554707,
  9.984369578019570859563e-6,
  1.50563273514931155834e-7
};

/* complex version of Lanczos method; this is not safe for export
 * since it becomes bad in the left half-plane
 */
static
int
lngamma_lanczos_complex(double zr, double zi, gsl_sf_result * yr, gsl_sf_result * yi)
{
  int k;
  gsl_sf_result log1_r,    log1_i;
  gsl_sf_result logAg_r,   logAg_i;
  double Ag_r, Ag_i;
  double yi_tmp_val, yi_tmp_err;

  zr -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag_r = lanczos_7_c[0];
  Ag_i = 0.0;
  for(k=1; k<=8; k++) {
    double R = zr + k;
    double I = zi;
    double a = lanczos_7_c[k] / (R*R + I*I);
    Ag_r +=  a * R;
    Ag_i -=  a * I;
  }

  gsl_sf_complex_log_e(zr + 7.5, zi, &log1_r,  &log1_i);
  gsl_sf_complex_log_e(Ag_r, Ag_i,   &logAg_r, &logAg_i);

  /* (z+0.5)*log(z+7.5) - (z+7.5) + LogRootTwoPi_ + log(Ag(z)) */
  yr->val = (zr+0.5)*log1_r.val - zi*log1_i.val - (zr+7.5) + LogRootTwoPi_ + logAg_r.val;
  yi->val = zi*log1_r.val + (zr+0.5)*log1_i.val - zi + logAg_i.val;
  yr->err = 4.0 * GSL_DBL_EPSILON * fabs(yr->val);
  yi->err = 4.0 * GSL_DBL_EPSILON * fabs(yi->val);
  yi_tmp_val = yi->val;
  yi_tmp_err = yi->err;
  gsl_sf_angle_restrict_symm_err_e(yi_tmp_val, yi);
  yi->err += yi_tmp_err;
  return GSL_SUCCESS;
}


/* Lanczos method for real x > 0;
 * gamma=7, truncated at 1/(z+8)
 * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
 */
static
int
lngamma_lanczos(double x, gsl_sf_result * result)
{
  int k;
  double Ag;
  double term1, term2;

  x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

  Ag = lanczos_7_c[0];
  for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

  /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
  term1 = (x+0.5)*log((x+7.5)/M_E);
  term2 = LogRootTwoPi_ + log(Ag);
  result->val  = term1 + (term2 - 7.0);
  result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(term1) + fabs(term2) + 7.0);
  result->err += GSL_DBL_EPSILON * fabs(result->val);

  return GSL_SUCCESS;
}

/* x = eps near zero
 * gives double-precision for |eps| < 0.02
 */
static
int
lngamma_sgn_0(double eps, gsl_sf_result * lng, double * sgn)
{
  /* calculate series for g(eps) = Gamma(eps) eps - 1/(1+eps) - eps/2 */
  const double c1  = -0.07721566490153286061;
  const double c2  = -0.01094400467202744461;
  const double c3  =  0.09252092391911371098;
  const double c4  = -0.01827191316559981266;
  const double c5  =  0.01800493109685479790;
  const double c6  = -0.00685088537872380685;
  const double c7  =  0.00399823955756846603;
  const double c8  = -0.00189430621687107802;
  const double c9  =  0.00097473237804513221;
  const double c10 = -0.00048434392722255893;
  const double g6  = c6+eps*(c7+eps*(c8 + eps*(c9 + eps*c10)));
  const double g   = eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*g6)))));

  /* calculate Gamma(eps) eps, a positive quantity */
  const double gee = g + 1.0/(1.0+eps) + 0.5*eps;

  lng->val = log(gee/fabs(eps));
  lng->err = 4.0 * GSL_DBL_EPSILON * fabs(lng->val);
  *sgn = GSL_SIGN(eps);

  return GSL_SUCCESS;
}


/* x near a negative integer
 * Calculates sign as well as log(|gamma(x)|).
 * x = -N + eps
 * assumes N >= 1
 */
static
int
lngamma_sgn_sing(int N, double eps, gsl_sf_result * lng, double * sgn)
{
  if(eps == 0.0) {
    lng->val = 0.0;
    lng->err = 0.0;
    *sgn = 0.0;
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(N == 1) {
    /* calculate series for
     * g = eps gamma(-1+eps) + 1 + eps/2 (1+3eps)/(1-eps^2)
     * double-precision for |eps| < 0.02
     */
    const double c0 =  0.07721566490153286061;
    const double c1 =  0.08815966957356030521;
    const double c2 = -0.00436125434555340577;
    const double c3 =  0.01391065882004640689;
    const double c4 = -0.00409427227680839100;
    const double c5 =  0.00275661310191541584;
    const double c6 = -0.00124162645565305019;
    const double c7 =  0.00065267976121802783;
    const double c8 = -0.00032205261682710437;
    const double c9 =  0.00016229131039545456;
    const double g5 = c5 + eps*(c6 + eps*(c7 + eps*(c8 + eps*c9)));
    const double g  = eps*(c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*g5)))));

    /* calculate eps gamma(-1+eps), a negative quantity */
    const double gam_e = g - 1.0 - 0.5*eps*(1.0+3.0*eps)/(1.0 - eps*eps);

    lng->val = log(fabs(gam_e)/fabs(eps));
    lng->err = 2.0 * GSL_DBL_EPSILON * fabs(lng->val);
    *sgn = ( eps > 0.0 ? -1.0 : 1.0 );
    return GSL_SUCCESS;
  }
  else {
    double g;

    /* series for sin(Pi(N+1-eps))/(Pi eps) modulo the sign
     * double-precision for |eps| < 0.02
     */
    const double cs1 = -1.6449340668482264365;
    const double cs2 =  0.8117424252833536436;
    const double cs3 = -0.1907518241220842137;
    const double cs4 =  0.0261478478176548005;
    const double cs5 = -0.0023460810354558236;
    const double e2  = eps*eps;
    const double sin_ser = 1.0 + e2*(cs1+e2*(cs2+e2*(cs3+e2*(cs4+e2*cs5))));

    /* calculate series for ln(gamma(1+N-eps))
     * double-precision for |eps| < 0.02
     */
    double aeps = fabs(eps);
    double c1, c2, c3, c4, c5, c6, c7;
    double lng_ser;
    gsl_sf_result c0;
    gsl_sf_result psi_0;
    gsl_sf_result psi_1;
    gsl_sf_result psi_2;
    gsl_sf_result psi_3;
    gsl_sf_result psi_4;
    gsl_sf_result psi_5;
    gsl_sf_result psi_6;
    psi_2.val = 0.0;
    psi_3.val = 0.0;
    psi_4.val = 0.0;
    psi_5.val = 0.0;
    psi_6.val = 0.0;
    gsl_sf_lnfact_e(N, &c0);
    gsl_sf_psi_int_e(N+1, &psi_0);
    gsl_sf_psi_1_int_e(N+1, &psi_1);
    if(aeps > 0.00001) gsl_sf_psi_n_e(2, N+1.0, &psi_2);
    if(aeps > 0.0002)  gsl_sf_psi_n_e(3, N+1.0, &psi_3);
    if(aeps > 0.001)   gsl_sf_psi_n_e(4, N+1.0, &psi_4);
    if(aeps > 0.005)   gsl_sf_psi_n_e(5, N+1.0, &psi_5);
    if(aeps > 0.01)    gsl_sf_psi_n_e(6, N+1.0, &psi_6);
    c1 = psi_0.val;
    c2 = psi_1.val/2.0;
    c3 = psi_2.val/6.0;
    c4 = psi_3.val/24.0;
    c5 = psi_4.val/120.0;
    c6 = psi_5.val/720.0;
    c7 = psi_6.val/5040.0;
    lng_ser = c0.val-eps*(c1-eps*(c2-eps*(c3-eps*(c4-eps*(c5-eps*(c6-eps*c7))))));

    /* calculate
     * g = ln(|eps gamma(-N+eps)|)
     *   = -ln(gamma(1+N-eps)) + ln(|eps Pi/sin(Pi(N+1+eps))|)
     */
    g = -lng_ser - log(sin_ser);

    lng->val = g - log(fabs(eps));
    lng->err = c0.err + 2.0 * GSL_DBL_EPSILON * (fabs(g) + fabs(lng->val));

    *sgn = ( GSL_IS_ODD(N) ? -1.0 : 1.0 ) * ( eps > 0.0 ? 1.0 : -1.0 );

    return GSL_SUCCESS;
  }
}


/* This gets bad near the negative half axis. However, this
 * region can be avoided by use of the reflection formula, as usual.
 * Only the first two terms of the series are kept.
 */
#if 0
static
int
lngamma_complex_stirling(const double zr, const double zi, double * lg_r, double * arg)
{
  double re_zinv,  im_zinv;
  double re_zinv2, im_zinv2;
  double re_zinv3, im_zinv3;
  double re_zhlnz, im_zhlnz;
  double r, lnr, theta;
  gsl_sf_complex_log_e(zr, zi, &lnr, &theta);  /* z = r e^{i theta} */
  r = exp(lnr);
  re_zinv =  (zr/r)/r;
  im_zinv = -(zi/r)/r;
  re_zinv2 = re_zinv*re_zinv - im_zinv*im_zinv;
  re_zinv2 = 2.0*re_zinv*im_zinv;
  re_zinv3 = re_zinv2*re_zinv - im_zinv2*im_zinv;
  re_zinv3 = re_zinv2*im_zinv + im_zinv2*re_zinv;
  re_zhlnz = (zr - 0.5)*lnr - zi*theta;
  im_zhlnz = zi*lnr + zr*theta;
  *lg_r = re_zhlnz - zr + 0.5*(M_LN2+M_LNPI) + re_zinv/12.0 - re_zinv3/360.0;
  *arg  = im_zhlnz - zi + 1.0/12.0*im_zinv - im_zinv3/360.0;
  return GSL_SUCCESS;
}
#endif /* 0 */


inline
static
int
lngamma_1_pade(const double eps, gsl_sf_result * result)
{
  /* Use (2,2) Pade for Log[Gamma[1+eps]]/eps
   * plus a correction series.
   */
  const double n1 = -1.0017419282349508699871138440;
  const double n2 =  1.7364839209922879823280541733;
  const double d1 =  1.2433006018858751556055436011;
  const double d2 =  5.0456274100274010152489597514;
  const double num = (eps + n1) * (eps + n2);
  const double den = (eps + d1) * (eps + d2);
  const double pade = 2.0816265188662692474880210318 * num / den;
  const double c0 =  0.004785324257581753;
  const double c1 = -0.01192457083645441;
  const double c2 =  0.01931961413960498;
  const double c3 = -0.02594027398725020;
  const double c4 =  0.03141928755021455;
  const double eps5 = eps*eps*eps*eps*eps;
  const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
  result->val = eps * (pade + corr);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

inline
static
int
lngamma_2_pade(const double eps, gsl_sf_result * result)
{
  /* Use (2,2) Pade for Log[Gamma[2+eps]]/eps
   * plus a correction series.
   */
  const double n1 = 1.000895834786669227164446568;
  const double n2 = 4.209376735287755081642901277;
  const double d1 = 2.618851904903217274682578255;
  const double d2 = 10.85766559900983515322922936;
  const double num = (eps + n1) * (eps + n2);
  const double den = (eps + d1) * (eps + d2);
  const double pade = 2.85337998765781918463568869 * num/den;
  const double c0 =  0.0001139406357036744;
  const double c1 = -0.0001365435269792533;
  const double c2 =  0.0001067287169183665;
  const double c3 = -0.0000693271800931282;
  const double c4 =  0.0000407220927867950;
  const double eps5 = eps*eps*eps*eps*eps;
  const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
  result->val = eps * (pade + corr);
  result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}


/* series for gammastar(x)
 * double-precision for x > 10.0
 */
static
int
gammastar_ser(const double x, gsl_sf_result * result)
{
  /* Use the Stirling series for the correction to Log(Gamma(x)),
   * which is better behaved and easier to compute than the
   * regular Stirling series for Gamma(x).
   */
  const double y = 1.0/(x*x);
  const double c0 =  1.0/12.0;
  const double c1 = -1.0/360.0;
  const double c2 =  1.0/1260.0;
  const double c3 = -1.0/1680.0;
  const double c4 =  1.0/1188.0;
  const double c5 = -691.0/360360.0;
  const double c6 =  1.0/156.0;
  const double c7 = -3617.0/122400.0;
  const double ser = c0 + y*(c1 + y*(c2 + y*(c3 + y*(c4 + y*(c5 + y*(c6 + y*c7))))));
  result->val = exp(ser/x);
  result->err = 2.0 * GSL_DBL_EPSILON * result->val * GSL_MAX_DBL(1.0, ser/x);
  return GSL_SUCCESS;
}


/* Chebyshev expansion for log(gamma(x)/gamma(8))
 * 5 < x < 10
 * -1 < t < 1
 */
static double gamma_5_10_data[24] = {
 -1.5285594096661578881275075214,
  4.8259152300595906319768555035,
  0.2277712320977614992970601978,
 -0.0138867665685617873604917300,
  0.0012704876495201082588139723,
 -0.0001393841240254993658962470,
  0.0000169709242992322702260663,
 -2.2108528820210580075775889168e-06,
  3.0196602854202309805163918716e-07,
 -4.2705675000079118380587357358e-08,
  6.2026423818051402794663551945e-09,
 -9.1993973208880910416311405656e-10,
  1.3875551258028145778301211638e-10,
 -2.1218861491906788718519522978e-11,
  3.2821736040381439555133562600e-12,
 -5.1260001009953791220611135264e-13,
  8.0713532554874636696982146610e-14,
 -1.2798522376569209083811628061e-14,
  2.0417711600852502310258808643e-15,
 -3.2745239502992355776882614137e-16,
  5.2759418422036579482120897453e-17,
 -8.5354147151695233960425725513e-18,
  1.3858639703888078291599886143e-18,
 -2.2574398807738626571560124396e-19
};
static const cheb_series gamma_5_10_cs = {
  gamma_5_10_data,
  23,
  -1, 1,
  11
};


/* gamma(x) for x >= 1/2
 * assumes x >= 1/2
 */
static
int
gamma_xgthalf(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.5) {
    result->val = 1.77245385090551602729817;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  } else if (x <= (GSL_SF_FACT_NMAX + 1.0) && x == floor(x)) {
    int n = (int) floor (x);
    result->val = fact_table[n - 1].f;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(fabs(x - 1.0) < 0.01) {
    /* Use series for Gamma[1+eps] - 1/(1+eps).
     */
    const double eps = x - 1.0;
    const double c1 =  0.4227843350984671394;
    const double c2 = -0.01094400467202744461;
    const double c3 =  0.09252092391911371098;
    const double c4 = -0.018271913165599812664;
    const double c5 =  0.018004931096854797895;
    const double c6 = -0.006850885378723806846;
    const double c7 =  0.003998239557568466030;
    result->val = 1.0/x + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*c7))))));
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(fabs(x - 2.0) < 0.01) {
    /* Use series for Gamma[1 + eps].
     */
    const double eps = x - 2.0;
    const double c1 =  0.4227843350984671394;
    const double c2 =  0.4118403304264396948;
    const double c3 =  0.08157691924708626638;
    const double c4 =  0.07424901075351389832;
    const double c5 = -0.00026698206874501476832;
    const double c6 =  0.011154045718130991049;
    const double c7 = -0.002852645821155340816;
    const double c8 =  0.0021039333406973880085;
    result->val = 1.0 + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(x < 5.0) {
    /* Exponentiating the logarithm is fine, as
     * long as the exponential is not so large
     * that it greatly amplifies the error.
     */
    gsl_sf_result lg;
    lngamma_lanczos(x, &lg);
    result->val = exp(lg.val);
    result->err = result->val * (lg.err + 2.0 * GSL_DBL_EPSILON);
    return GSL_SUCCESS;
  }
  else if(x < 10.0) {
    /* This is a sticky area. The logarithm
     * is too large and the gammastar series
     * is not good.
     */
    const double gamma_8 = 5040.0;
    const double t = (2.0*x - 15.0)/5.0;
    gsl_sf_result c;
    cheb_eval_e(&gamma_5_10_cs, t, &c);
    result->val  = exp(c.val) * gamma_8;
    result->err  = result->val * c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else if(x < GSL_SF_GAMMA_XMAX) {
    /* We do not want to exponentiate the logarithm
     * if x is large because of the inevitable
     * inflation of the error. So we carefully
     * use pow() and exp() with exact quantities.
     */
    double p = pow(x, 0.5*x);
    double e = exp(-x);
    double q = (p * e) * p;
    double pre = M_SQRT2 * M_SQRTPI * q/sqrt(x);
    gsl_sf_result gstar;
    int stat_gs = gammastar_ser(x, &gstar);
    result->val = pre * gstar.val;
    result->err = (x + 2.5) * GSL_DBL_EPSILON * result->val;
    return stat_gs;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_lngamma_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(fabs(x - 1.0) < 0.01) {
    /* Note that we must amplify the errors
     * from the Pade evaluations because of
     * the way we must pass the argument, i.e.
     * writing (1-x) is a loss of precision
     * when x is near 1.
     */
    int stat = lngamma_1_pade(x - 1.0, result);
    result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
    return stat;
  }
  else if(fabs(x - 2.0) < 0.01) {
    int stat = lngamma_2_pade(x - 2.0, result);
    result->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
    return stat;
  }
  else if(x >= 0.5) {
    return lngamma_lanczos(x, result);
  }
  else if(x == 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(fabs(x) < 0.02) {
    double sgn;
    return lngamma_sgn_0(x, result, &sgn);
  }
  else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
    /* Try to extract a fractional
     * part from x.
     */
    double z  = 1.0 - x;
    double s  = sin(M_PI*z);
    double as = fabs(s);
    if(s == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(as < M_PI*0.015) {
      /* x is near a negative integer, -N */
      if(x < INT_MIN + 2.0) {
        result->val = 0.0;
        result->err = 0.0;
        GSL_ERROR ("error", GSL_EROUND);
      }
      else {
        int N = -(int)(x - 0.5);
        double eps = x + N;
        double sgn;
        return lngamma_sgn_sing(N, eps, result, &sgn);
      }
    }
    else {
      gsl_sf_result lg_z;
      lngamma_lanczos(z, &lg_z);
      result->val = M_LNPI - (log(as) + lg_z.val);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg_z.err;
      return GSL_SUCCESS;
    }
  }
  else {
    /* |x| was too large to extract any fractional part */
    result->val = 0.0;
    result->err = 0.0;
    GSL_ERROR ("error", GSL_EROUND);
  }
}


int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double * sgn)
{
  if(fabs(x - 1.0) < 0.01) {
    int stat = lngamma_1_pade(x - 1.0, result_lg);
    result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 1.0));
    *sgn = 1.0;
    return stat;
  }
  else if(fabs(x - 2.0) < 0.01) {
   int stat = lngamma_2_pade(x - 2.0, result_lg);
    result_lg->err *= 1.0/(GSL_DBL_EPSILON + fabs(x - 2.0));
    *sgn = 1.0;
    return stat;
  }
  else if(x >= 0.5) {
    *sgn = 1.0;
    return lngamma_lanczos(x, result_lg);
  }
  else if(x == 0.0) {
    *sgn = 0.0;
    DOMAIN_ERROR(result_lg);
  }
  else if(fabs(x) < 0.02) {
    return lngamma_sgn_0(x, result_lg, sgn);
  }
  else if(x > -0.5/(GSL_DBL_EPSILON*M_PI)) {
   /* Try to extract a fractional
     * part from x.
     */
    double z = 1.0 - x;
    double s = sin(M_PI*x);
    double as = fabs(s);
    if(s == 0.0) {
      *sgn = 0.0;
      DOMAIN_ERROR(result_lg);
    }
    else if(as < M_PI*0.015) {
      /* x is near a negative integer, -N */
      if(x < INT_MIN + 2.0) {
        result_lg->val = 0.0;
        result_lg->err = 0.0;
        *sgn = 0.0;
        GSL_ERROR ("error", GSL_EROUND);
      }
      else {
        int N = -(int)(x - 0.5);
        double eps = x + N;
        return lngamma_sgn_sing(N, eps, result_lg, sgn);
      }
    }
    else {
      gsl_sf_result lg_z;
      lngamma_lanczos(z, &lg_z);
      *sgn = (s > 0.0 ? 1.0 : -1.0);
      result_lg->val = M_LNPI - (log(as) + lg_z.val);
      result_lg->err = 2.0 * GSL_DBL_EPSILON * fabs(result_lg->val) + lg_z.err;
      return GSL_SUCCESS;
    }
  }
  else {
    /* |x| was too large to extract any fractional part */
    result_lg->val = 0.0;
    result_lg->err = 0.0;
    *sgn = 0.0;
    GSL_ERROR ("x too large to extract fraction part", GSL_EROUND);
  }
}


int
gsl_sf_gamma_e(const double x, gsl_sf_result * result)
{
  if(x < 0.5) {
    int rint_x = (int)floor(x+0.5);
    double f_x = x - rint_x;
    double sgn_gamma = ( GSL_IS_EVEN(rint_x) ? 1.0 : -1.0 );
    double sin_term = sgn_gamma * sin(M_PI * f_x) / M_PI;

    if(sin_term == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(x > -169.0) {
      gsl_sf_result g;
      gamma_xgthalf(1.0-x, &g);
      if(fabs(sin_term) * g.val * GSL_DBL_MIN < 1.0) {
        result->val  = 1.0/(sin_term * g.val);
        result->err  = fabs(g.err/g.val) * fabs(result->val);
        result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
      else {
        UNDERFLOW_ERROR(result);
      }
    }
    else {
      /* It is hard to control it here.
       * We can only exponentiate the
       * logarithm and eat the loss of
       * precision.
       */
      gsl_sf_result lng;
      double sgn;
      int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
      int stat_e   = gsl_sf_exp_mult_err_e(lng.val, lng.err, sgn, 0.0, result);
      return GSL_ERROR_SELECT_2(stat_e, stat_lng);
    }
  }
  else {
    return gamma_xgthalf(x, result);
  }
}


int
gsl_sf_gammastar_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x < 0.5) {
    gsl_sf_result lg;
    const int stat_lg = gsl_sf_lngamma_e(x, &lg);
    const double lx = log(x);
    const double c  = 0.5*(M_LN2+M_LNPI);
    const double lnr_val = lg.val - (x-0.5)*lx + x - c;
    const double lnr_err = lg.err + 2.0 * GSL_DBL_EPSILON *((x+0.5)*fabs(lx) + c);
    const int stat_e  = gsl_sf_exp_err_e(lnr_val, lnr_err, result);
    return GSL_ERROR_SELECT_2(stat_lg, stat_e);
  }
  else if(x < 2.0) {
    const double t = 4.0/3.0*(x-0.5) - 1.0;
    return cheb_eval_e(&gstar_a_cs, t, result);
  }
  else if(x < 10.0) {
    const double t = 0.25*(x-2.0) - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&gstar_b_cs, t, &c);
    result->val  = c.val/(x*x) + 1.0 + 1.0/(12.0*x);
    result->err  = c.err/(x*x);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 1.0/GSL_ROOT4_DBL_EPSILON) {
    return gammastar_ser(x, result);
  }
  else if(x < 1.0/GSL_DBL_EPSILON) {
    /* Use Stirling formula for Gamma(x).
     */
    const double xi = 1.0/x;
    result->val = 1.0 + xi/12.0*(1.0 + xi/24.0*(1.0 - xi*(139.0/180.0 + 571.0/8640.0*xi)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 1.0;
    result->err = 1.0/x;
    return GSL_SUCCESS;
  }
}


int
gsl_sf_gammainv_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if (x <= 0.0 && x == floor(x)) { /* negative integer */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  } else if(x < 0.5) {
    gsl_sf_result lng;
    double sgn;
    int stat_lng = gsl_sf_lngamma_sgn_e(x, &lng, &sgn);
    if(stat_lng == GSL_EDOM) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(stat_lng != GSL_SUCCESS) {
      result->val = 0.0;
      result->err = 0.0;
      return stat_lng;
    }
    else {
      return gsl_sf_exp_mult_err_e(-lng.val, lng.err, sgn, 0.0, result);
    }
  }
  else {
    gsl_sf_result g;
    int stat_g = gamma_xgthalf(x, &g);
    if(stat_g == GSL_EOVRFLW) {
      UNDERFLOW_ERROR(result);
    }
    else {
      result->val  = 1.0/g.val;
      result->err  = fabs(g.err/g.val) * fabs(result->val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      CHECK_UNDERFLOW(result);
      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg)
{
  if(zr <= 0.5) {
    /* Transform to right half plane using reflection;
     * in fact we do a little better by stopping at 1/2.
     */
    double x = 1.0-zr;
    double y = -zi;
    gsl_sf_result a, b;
    gsl_sf_result lnsin_r, lnsin_i;

    int stat_l = lngamma_lanczos_complex(x, y, &a, &b);
    int stat_s = gsl_sf_complex_logsin_e(M_PI*zr, M_PI*zi, &lnsin_r, &lnsin_i);

    if(stat_s == GSL_SUCCESS) {
      int stat_r;
      lnr->val = M_LNPI - lnsin_r.val - a.val;
      lnr->err = lnsin_r.err + a.err + 2.0 * GSL_DBL_EPSILON * fabs(lnr->val);
      arg->val = -lnsin_i.val - b.val;
      arg->err = lnsin_i.err + b.err + 2.0 * GSL_DBL_EPSILON * fabs(arg->val);
      stat_r = gsl_sf_angle_restrict_symm_e(&(arg->val));
      return GSL_ERROR_SELECT_2(stat_r, stat_l);
    }
    else {
      DOMAIN_ERROR_2(lnr,arg);
    }
  }
  else {
    /* otherwise plain vanilla Lanczos */
    return lngamma_lanczos_complex(zr, zi, lnr, arg);
  }
}


int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < 0.0 || n < 0) {
    DOMAIN_ERROR(result);
  }
  else if(n == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n == 1) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    const double log2pi = M_LNPI + M_LN2;
    const double ln_test = n*(log(x)+1.0) + 1.0 - (n+0.5)*log(n+1.0) + 0.5*log2pi;

    if(ln_test < GSL_LOG_DBL_MIN+1.0) {
      UNDERFLOW_ERROR(result);
    }
    else if(ln_test > GSL_LOG_DBL_MAX-1.0) {
      OVERFLOW_ERROR(result);
    }
    else {
      double product = 1.0;
      int k;
      for(k=1; k<=n; k++) {
        product *= (x/k);
      }
      result->val = product;
      result->err = n * GSL_DBL_EPSILON * product;
      CHECK_UNDERFLOW(result);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 18) {
    result->val = fact_table[n].f;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n <= GSL_SF_FACT_NMAX){
    result->val = fact_table[n].f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 26) {
    result->val = doub_fact_table[n].f;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(n <= GSL_SF_DOUBLEFACT_NMAX){
    result->val = doub_fact_table[n].f;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    OVERFLOW_ERROR(result);
  }
}


int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= GSL_SF_FACT_NMAX){
    result->val = log(fact_table[n].f);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_lngamma_e(n+1.0, result);
    return GSL_SUCCESS;
  }
}


int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= GSL_SF_DOUBLEFACT_NMAX){
    result->val = log(doub_fact_table[n].f);
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(GSL_IS_ODD(n)) {
    gsl_sf_result lg;
    gsl_sf_lngamma_e(0.5*(n+2.0), &lg);
    result->val = 0.5*(n+1.0) * M_LN2 - 0.5*M_LNPI + lg.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result lg;
    gsl_sf_lngamma_e(0.5*n+1.0, &lg);
    result->val = 0.5*n*M_LN2 + lg.val;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val) + lg.err;
    return GSL_SUCCESS;
  }
}


int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result nf;
    gsl_sf_result mf;
    gsl_sf_result nmmf;
    if(m*2 > n) m = n-m;
    gsl_sf_lnfact_e(n, &nf);
    gsl_sf_lnfact_e(m, &mf);
    gsl_sf_lnfact_e(n-m, &nmmf);
    result->val  = nf.val - mf.val - nmmf.val;
    result->err  = nf.err + mf.err + nmmf.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result)
{
  if(m > n) {
    DOMAIN_ERROR(result);
  }
  else if(m == n || m == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if (n <= GSL_SF_FACT_NMAX) {
    result->val = (fact_table[n].f / fact_table[m].f) / fact_table[n-m].f;
    result->err = 6.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  } else {
    if(m*2 < n) m = n-m;

    if (n - m < 64)  /* compute product for a manageable number of terms */
      {
        double prod = 1.0;
        unsigned int k;

        for(k=n; k>=m+1; k--) {
          double tk = (double)k / (double)(k-m);
          if(tk > GSL_DBL_MAX/prod) {
            OVERFLOW_ERROR(result);
          }
          prod *= tk;
        }
        result->val = prod;
        result->err = 2.0 * GSL_DBL_EPSILON * prod * fabs((double) (n-m));
        return GSL_SUCCESS;
      }
    else
      {
        gsl_sf_result lc;
        const int stat_lc = gsl_sf_lnchoose_e (n, m, &lc);
        const int stat_e  = gsl_sf_exp_err_e(lc.val, lc.err, result);
        return GSL_ERROR_SELECT_2(stat_lc, stat_e);
      }
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_fact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_fact_e(n, &result));
}

double gsl_sf_lnfact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_lnfact_e(n, &result));
}

double gsl_sf_doublefact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_doublefact_e(n, &result));
}

double gsl_sf_lndoublefact(const unsigned int n)
{
  EVAL_RESULT(gsl_sf_lndoublefact_e(n, &result));
}

double gsl_sf_lngamma(const double x)
{
  EVAL_RESULT(gsl_sf_lngamma_e(x, &result));
}

double gsl_sf_gamma(const double x)
{
  EVAL_RESULT(gsl_sf_gamma_e(x, &result));
}

double gsl_sf_gammastar(const double x)
{
  EVAL_RESULT(gsl_sf_gammastar_e(x, &result));
}

double gsl_sf_gammainv(const double x)
{
  EVAL_RESULT(gsl_sf_gammainv_e(x, &result));
}

double gsl_sf_taylorcoeff(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_taylorcoeff_e(n, x, &result));
}

double gsl_sf_choose(unsigned int n, unsigned int m)
{
  EVAL_RESULT(gsl_sf_choose_e(n, m, &result));
}

double gsl_sf_lnchoose(unsigned int n, unsigned int m)
{
  EVAL_RESULT(gsl_sf_lnchoose_e(n, m, &result));
}
/* end inlined source: specfunc/gamma.c */
/* begin inlined source: specfunc/elementary.c */
/* specfunc/elementary.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_elementary.h */
/* specfunc/gsl_sf_elementary.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

/* Miscellaneous elementary functions and operations.
 */
#ifndef __GSL_SF_ELEMENTARY_H__
#define __GSL_SF_ELEMENTARY_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Multiplication.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_multiply(const double x, const double y);


/* Multiplication of quantities with associated errors.
 */
int gsl_sf_multiply_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_ELEMENTARY_H__ */

/* end inlined header: gsl/gsl_sf_elementary.h */
/* inlined: error.h */
/* inlined: check.h */

int
gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result)
{
  const double ax = fabs(x);
  const double ay = fabs(y);

  if(x == 0.0 || y == 0.0) {
    /* It is necessary to eliminate this immediately.
     */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if((ax <= 1.0 && ay >= 1.0) || (ay <= 1.0 && ax >= 1.0)) {
    /* Straddling 1.0 is always safe.
     */
    result->val = x*y;
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double f = 1.0 - 2.0 * GSL_DBL_EPSILON;
    const double min = GSL_MIN_DBL(fabs(x), fabs(y));
    const double max = GSL_MAX_DBL(fabs(x), fabs(y));
    if(max < 0.9 * GSL_SQRT_DBL_MAX || min < (f * DBL_MAX)/max) {
      result->val = GSL_COERCE_DBL(x*y);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      CHECK_UNDERFLOW(result);
      return GSL_SUCCESS;
    }
    else {
      OVERFLOW_ERROR(result);
    }
  }
}


int
gsl_sf_multiply_err_e(const double x, const double dx,
                         const double y, const double dy,
                         gsl_sf_result * result)
{
  int status = gsl_sf_multiply_e(x, y, result);
  result->err += fabs(dx*y) + fabs(dy*x);
  return status;
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_multiply(const double x, const double y)
{
  EVAL_RESULT(gsl_sf_multiply_e(x, y, &result));
}

/* end inlined source: specfunc/elementary.c */
/* begin inlined source: specfunc/zeta.c */
/* specfunc/zeta.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_elementary.h */
/* specfunc/gsl_sf_elementary.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

/* Miscellaneous elementary functions and operations.
 */
#ifndef __GSL_SF_ELEMENTARY_H__
#define __GSL_SF_ELEMENTARY_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Multiplication.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_multiply_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_multiply(const double x, const double y);


/* Multiplication of quantities with associated errors.
 */
int gsl_sf_multiply_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_ELEMENTARY_H__ */

/* end inlined header: gsl/gsl_sf_elementary.h */
/* begin inlined header: gsl/gsl_sf_exp.h */
/* specfunc/gsl_sf_exp.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_EXP_H__
#define __GSL_SF_EXP_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
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

/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e(const double x, gsl_sf_result * result);
double gsl_sf_exp(const double x);


/* Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_exp_mult(const double x, const double y);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
double gsl_sf_expm1(const double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
double gsl_sf_exprel(const double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
double gsl_sf_exprel_2(const double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_exprel_n(const int n, const double x);

int gsl_sf_exprel_n_CF_e(const double n, const double x, gsl_sf_result * result);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);

/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);

__END_DECLS

#endif /* __GSL_SF_EXP_H__ */

/* end inlined header: gsl/gsl_sf_exp.h */
/* begin inlined header: gsl/gsl_sf_gamma.h */
/* specfunc/gsl_sf_gamma.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_GAMMA_H__
#define __GSL_SF_GAMMA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
double gsl_sf_lngamma(const double x);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
double gsl_sf_gamma(const double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
double gsl_sf_gammastar(const double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
double gsl_sf_gammainv(const double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_taylorcoeff(const int n, const double x);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_fact(const unsigned int n);


/* n!! = n(n-2)(n-4) ...
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_doublefact(const unsigned int n);


/* log(n!)
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lnfact(const unsigned int n);


/* log(n!!)
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lndoublefact(const unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochhammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_lnpoch(const double a, const double x);


/* Logarithm of Pochhammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_e(const double a, const double x, gsl_sf_result * result, double * sgn);


/* Pochhammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_poch(const double a, const double x);


/* Relative Pochhammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_pochrel(const double a, const double x);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a >= 0, x >= 0
 *   Q(a,0) := 1
 *   Q(0,x) := 0, x != 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_Q(const double a, const double x);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_P(const double a, const double x);


/* Non-normalized Incomplete Gamma Function
 *
 * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * x >= 0.0
 *   Gamma(a, 0) := Gamma(a)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc(const double a, const double x);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_lnbeta(const double a, const double b);

int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_beta(const double a, const double b);


/* Normalized Incomplete Beta Function
 * B_x(a,b)/B(a,b)
 *
 * a > 0, b > 0, 0 <= x <= 1
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_beta_inc(const double a, const double b, const double x);


/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297

__END_DECLS

#endif /* __GSL_SF_GAMMA_H__ */

/* end inlined header: gsl/gsl_sf_gamma.h */
/* begin inlined header: gsl/gsl_sf_pow_int.h */
/* specfunc/gsl_sf_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_POW_INT_H__
#define __GSL_SF_POW_INT_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Calculate x^n.
 * Does not check for overflow/underflow.
 */
int     gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result);
double  gsl_sf_pow_int(const double x, const int n);


__END_DECLS

#endif /* __GSL_SF_POW_INT_H__ */

/* end inlined header: gsl/gsl_sf_pow_int.h */
/* begin inlined header: gsl/gsl_sf_zeta.h */
/* specfunc/gsl_sf_zeta.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_ZETA_H__
#define __GSL_SF_ZETA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Riemann Zeta Function
 * zeta(n) = Sum[ k^(-n), {k,1,Infinity} ]
 *
 * n=integer, n != 1
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result);
double gsl_sf_zeta_int(const int n);


/* Riemann Zeta Function
 * zeta(x) = Sum[ k^(-s), {k,1,Infinity} ], s != 1.0
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_e(const double s, gsl_sf_result * result);
double gsl_sf_zeta(const double s);


/* Riemann Zeta Function minus 1
 *   useful for evaluating the fractional part
 *   of Riemann zeta for large argument
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zetam1_e(const double s, gsl_sf_result * result);
double gsl_sf_zetam1(const double s);


/* Riemann Zeta Function minus 1 for integer arg
 *   useful for evaluating the fractional part
 *   of Riemann zeta for large argument
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zetam1_int_e(const int s, gsl_sf_result * result);
double gsl_sf_zetam1_int(const int s);


/* Hurwitz Zeta Function
 * zeta(s,q) = Sum[ (k+q)^(-s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result);
double gsl_sf_hzeta(const double s, const double q);


/* Eta Function
 * eta(n) = (1-2^(1-n)) zeta(n)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_int_e(int n, gsl_sf_result * result);
double gsl_sf_eta_int(const int n);


/* Eta Function
 * eta(s) = (1-2^(1-s)) zeta(s)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_e(const double s, gsl_sf_result * result);
double gsl_sf_eta(const double s);


__END_DECLS

#endif /* __GSL_SF_ZETA_H__ */

/* end inlined header: gsl/gsl_sf_zeta.h */
/* inlined: error.h */

/* inlined: chebyshev.h */
/* inlined once already: specfunc/cheb_eval.c */

#define LogTwoPi_  1.8378770664093454835606594728111235279723


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* chebyshev fit for (s(t)-1)Zeta[s(t)]
 * s(t)= (t+1)/2
 * -1 <= t <= 1
 */
static double zeta_xlt1_data[14] = {
  1.48018677156931561235192914649,
  0.25012062539889426471999938167,
  0.00991137502135360774243761467,
 -0.00012084759656676410329833091,
 -4.7585866367662556504652535281e-06,
  2.2229946694466391855561441361e-07,
 -2.2237496498030257121309056582e-09,
 -1.0173226513229028319420799028e-10,
  4.3756643450424558284466248449e-12,
 -6.2229632593100551465504090814e-14,
 -6.6116201003272207115277520305e-16,
  4.9477279533373912324518463830e-17,
 -1.0429819093456189719660003522e-18,
  6.9925216166580021051464412040e-21,
};
static cheb_series zeta_xlt1_cs = {
  zeta_xlt1_data,
  13,
  -1, 1,
  8
};

/* chebyshev fit for (s(t)-1)Zeta[s(t)]
 * s(t)= (19t+21)/2
 * -1 <= t <= 1
 */
static double zeta_xgt1_data[30] = {
  19.3918515726724119415911269006,
   9.1525329692510756181581271500,
   0.2427897658867379985365270155,
  -0.1339000688262027338316641329,
   0.0577827064065028595578410202,
  -0.0187625983754002298566409700,
   0.0039403014258320354840823803,
  -0.0000581508273158127963598882,
  -0.0003756148907214820704594549,
   0.0001892530548109214349092999,
  -0.0000549032199695513496115090,
   8.7086484008939038610413331863e-6,
   6.4609477924811889068410083425e-7,
  -9.6749773915059089205835337136e-7,
   3.6585400766767257736982342461e-7,
  -8.4592516427275164351876072573e-8,
   9.9956786144497936572288988883e-9,
   1.4260036420951118112457144842e-9,
  -1.1761968823382879195380320948e-9,
   3.7114575899785204664648987295e-10,
  -7.4756855194210961661210215325e-11,
   7.8536934209183700456512982968e-12,
   9.9827182259685539619810406271e-13,
  -7.5276687030192221587850302453e-13,
   2.1955026393964279988917878654e-13,
  -4.1934859852834647427576319246e-14,
   4.6341149635933550715779074274e-15,
   2.3742488509048340106830309402e-16,
  -2.7276516388124786119323824391e-16,
   7.8473570134636044722154797225e-17
};
static cheb_series zeta_xgt1_cs = {
  zeta_xgt1_data,
  29,
  -1, 1,
  17
};


/* chebyshev fit for Ln[Zeta[s(t)] - 1 - 2^(-s(t))]
 * s(t)= 10 + 5t
 * -1 <= t <= 1; 5 <= s <= 15
 */
static double zetam1_inter_data[24] = {
  -21.7509435653088483422022339374,
  -5.63036877698121782876372020472,
   0.0528041358684229425504861579635,
  -0.0156381809179670789342700883562,
   0.00408218474372355881195080781927,
  -0.0010264867349474874045036628282,
   0.000260469880409886900143834962387,
  -0.0000676175847209968878098566819447,
   0.0000179284472587833525426660171124,
  -4.83238651318556188834107605116e-6,
   1.31913788964999288471371329447e-6,
  -3.63760500656329972578222188542e-7,
   1.01146847513194744989748396574e-7,
  -2.83215225141806501619105289509e-8,
   7.97733710252021423361012829496e-9,
  -2.25850168553956886676250696891e-9,
   6.42269392950164306086395744145e-10,
  -1.83363861846127284505060843614e-10,
   5.25309763895283179960368072104e-11,
  -1.50958687042589821074710575446e-11,
   4.34997545516049244697776942981e-12,
  -1.25597782748190416118082322061e-12,
   3.61280740072222650030134104162e-13,
  -9.66437239205745207188920348801e-14
};
static cheb_series zetam1_inter_cs = {
  zetam1_inter_data,
  22,
  -1, 1,
  12
};



/* assumes s >= 0 and s != 1.0 */
inline
static int
riemann_zeta_sgt0(double s, gsl_sf_result * result)
{
  if(s < 1.0) {
    gsl_sf_result c;
    cheb_eval_e(&zeta_xlt1_cs, 2.0*s - 1.0, &c);
    result->val = c.val / (s - 1.0);
    result->err = c.err / fabs(s-1.0) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(s <= 20.0) {
    double x = (2.0*s - 21.0)/19.0;
    gsl_sf_result c;
    cheb_eval_e(&zeta_xgt1_cs, x, &c);
    result->val = c.val / (s - 1.0);
    result->err = c.err / (s - 1.0) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double f2 = 1.0 - pow(2.0,-s);
    double f3 = 1.0 - pow(3.0,-s);
    double f5 = 1.0 - pow(5.0,-s);
    double f7 = 1.0 - pow(7.0,-s);
    result->val = 1.0/(f2*f3*f5*f7);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


inline
static int
riemann_zeta1ms_slt0(double s, gsl_sf_result * result)
{
  if(s > -19.0) {
    double x = (-19 - 2.0*s)/19.0;
    gsl_sf_result c;
    cheb_eval_e(&zeta_xgt1_cs, x, &c);
    result->val = c.val / (-s);
    result->err = c.err / (-s) + GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double f2 = 1.0 - pow(2.0,-(1.0-s));
    double f3 = 1.0 - pow(3.0,-(1.0-s));
    double f5 = 1.0 - pow(5.0,-(1.0-s));
    double f7 = 1.0 - pow(7.0,-(1.0-s));
    result->val = 1.0/(f2*f3*f5*f7);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


/* works for 5 < s < 15*/
static int
riemann_zeta_minus_1_intermediate_s(double s, gsl_sf_result * result)
{
  double t = (s - 10.0)/5.0;
  gsl_sf_result c;
  cheb_eval_e(&zetam1_inter_cs, t, &c);
  result->val = exp(c.val) + pow(2.0, -s);
  result->err = (c.err + 2.0*GSL_DBL_EPSILON)*result->val;
  return GSL_SUCCESS;
}


/* assumes s is large and positive
 * write: zeta(s) - 1 = zeta(s) * (1 - 1/zeta(s))
 * and expand a few terms of the product formula to evaluate 1 - 1/zeta(s)
 *
 * works well for s > 15
 */
static int
riemann_zeta_minus1_large_s(double s, gsl_sf_result * result)
{
  double a = pow( 2.0,-s);
  double b = pow( 3.0,-s);
  double c = pow( 5.0,-s);
  double d = pow( 7.0,-s);
  double e = pow(11.0,-s);
  double f = pow(13.0,-s);
  double t1 = a + b + c + d + e + f;
  double t2 = a*(b+c+d+e+f) + b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f;
  /*
  double t3 = a*(b*(c+d+e+f) + c*(d+e+f) + d*(e+f) + e*f) +
              b*(c*(d+e+f) + d*(e+f) + e*f) +
              c*(d*(e+f) + e*f) +
              d*e*f;
  double t4 = a*(b*(c*(d + e + f) + d*(e + f) + e*f) + c*(d*(e+f) + e*f) + d*e*f) +
              b*(c*(d*(e+f) + e*f) + d*e*f) +
              c*d*e*f;
  double t5 = b*c*d*e*f + a*c*d*e*f+ a*b*d*e*f+ a*b*c*e*f+ a*b*c*d*f+ a*b*c*d*e;
  double t6 = a*b*c*d*e*f;
  */
  double numt = t1 - t2 /* + t3 - t4 + t5 - t6 */;
  double zeta = 1.0/((1.0-a)*(1.0-b)*(1.0-c)*(1.0-d)*(1.0-e)*(1.0-f));
  result->val = numt*zeta;
  result->err = (15.0/s + 1.0) * 6.0*GSL_DBL_EPSILON*result->val;
  return GSL_SUCCESS;
}


#if 0
/* zeta(n) */
#define ZETA_POS_TABLE_NMAX   100
static double zeta_pos_int_table_OLD[ZETA_POS_TABLE_NMAX+1] = {
 -0.50000000000000000000000000000,       /* zeta(0) */
  0.0 /* FIXME: DirectedInfinity() */,   /* zeta(1) */
  1.64493406684822643647241516665,       /* ...     */
  1.20205690315959428539973816151,
  1.08232323371113819151600369654,
  1.03692775514336992633136548646,
  1.01734306198444913971451792979,
  1.00834927738192282683979754985,
  1.00407735619794433937868523851,
  1.00200839282608221441785276923,
  1.00099457512781808533714595890,
  1.00049418860411946455870228253,
  1.00024608655330804829863799805,
  1.00012271334757848914675183653,
  1.00006124813505870482925854511,
  1.00003058823630702049355172851,
  1.00001528225940865187173257149,
  1.00000763719763789976227360029,
  1.00000381729326499983985646164,
  1.00000190821271655393892565696,
  1.00000095396203387279611315204,
  1.00000047693298678780646311672,
  1.00000023845050272773299000365,
  1.00000011921992596531107306779,
  1.00000005960818905125947961244,
  1.00000002980350351465228018606,
  1.00000001490155482836504123466,
  1.00000000745071178983542949198,
  1.00000000372533402478845705482,
  1.00000000186265972351304900640,
  1.00000000093132743241966818287,
  1.00000000046566290650337840730,
  1.00000000023283118336765054920,
  1.00000000011641550172700519776,
  1.00000000005820772087902700889,
  1.00000000002910385044497099687,
  1.00000000001455192189104198424,
  1.00000000000727595983505748101,
  1.00000000000363797954737865119,
  1.00000000000181898965030706595,
  1.00000000000090949478402638893,
  1.00000000000045474737830421540,
  1.00000000000022737368458246525,
  1.00000000000011368684076802278,
  1.00000000000005684341987627586,
  1.00000000000002842170976889302,
  1.00000000000001421085482803161,
  1.00000000000000710542739521085,
  1.00000000000000355271369133711,
  1.00000000000000177635684357912,
  1.00000000000000088817842109308,
  1.00000000000000044408921031438,
  1.00000000000000022204460507980,
  1.00000000000000011102230251411,
  1.00000000000000005551115124845,
  1.00000000000000002775557562136,
  1.00000000000000001387778780973,
  1.00000000000000000693889390454,
  1.00000000000000000346944695217,
  1.00000000000000000173472347605,
  1.00000000000000000086736173801,
  1.00000000000000000043368086900,
  1.00000000000000000021684043450,
  1.00000000000000000010842021725,
  1.00000000000000000005421010862,
  1.00000000000000000002710505431,
  1.00000000000000000001355252716,
  1.00000000000000000000677626358,
  1.00000000000000000000338813179,
  1.00000000000000000000169406589,
  1.00000000000000000000084703295,
  1.00000000000000000000042351647,
  1.00000000000000000000021175824,
  1.00000000000000000000010587912,
  1.00000000000000000000005293956,
  1.00000000000000000000002646978,
  1.00000000000000000000001323489,
  1.00000000000000000000000661744,
  1.00000000000000000000000330872,
  1.00000000000000000000000165436,
  1.00000000000000000000000082718,
  1.00000000000000000000000041359,
  1.00000000000000000000000020680,
  1.00000000000000000000000010340,
  1.00000000000000000000000005170,
  1.00000000000000000000000002585,
  1.00000000000000000000000001292,
  1.00000000000000000000000000646,
  1.00000000000000000000000000323,
  1.00000000000000000000000000162,
  1.00000000000000000000000000081,
  1.00000000000000000000000000040,
  1.00000000000000000000000000020,
  1.00000000000000000000000000010,
  1.00000000000000000000000000005,
  1.00000000000000000000000000003,
  1.00000000000000000000000000001,
  1.00000000000000000000000000001,
  1.00000000000000000000000000000,
  1.00000000000000000000000000000,
  1.00000000000000000000000000000
};
#endif /* 0 */


/* zeta(n) - 1 */
#define ZETA_POS_TABLE_NMAX   100
static double zetam1_pos_int_table[ZETA_POS_TABLE_NMAX+1] = {
 -1.5,                               /* zeta(0) */
  0.0,       /* FIXME: Infinity */   /* zeta(1) - 1 */
  0.644934066848226436472415166646,  /* zeta(2) - 1 */
  0.202056903159594285399738161511,
  0.082323233711138191516003696541,
  0.036927755143369926331365486457,
  0.017343061984449139714517929790,
  0.008349277381922826839797549849,
  0.004077356197944339378685238508,
  0.002008392826082214417852769232,
  0.000994575127818085337145958900,
  0.000494188604119464558702282526,
  0.000246086553308048298637998047,
  0.000122713347578489146751836526,
  0.000061248135058704829258545105,
  0.000030588236307020493551728510,
  0.000015282259408651871732571487,
  7.6371976378997622736002935630e-6,
  3.8172932649998398564616446219e-6,
  1.9082127165539389256569577951e-6,
  9.5396203387279611315203868344e-7,
  4.7693298678780646311671960437e-7,
  2.3845050272773299000364818675e-7,
  1.1921992596531107306778871888e-7,
  5.9608189051259479612440207935e-8,
  2.9803503514652280186063705069e-8,
  1.4901554828365041234658506630e-8,
  7.4507117898354294919810041706e-9,
  3.7253340247884570548192040184e-9,
  1.8626597235130490064039099454e-9,
  9.3132743241966818287176473502e-10,
  4.6566290650337840729892332512e-10,
  2.3283118336765054920014559759e-10,
  1.1641550172700519775929738354e-10,
  5.8207720879027008892436859891e-11,
  2.9103850444970996869294252278e-11,
  1.4551921891041984235929632245e-11,
  7.2759598350574810145208690123e-12,
  3.6379795473786511902372363558e-12,
  1.8189896503070659475848321007e-12,
  9.0949478402638892825331183869e-13,
  4.5474737830421540267991120294e-13,
  2.2737368458246525152268215779e-13,
  1.1368684076802278493491048380e-13,
  5.6843419876275856092771829675e-14,
  2.8421709768893018554550737049e-14,
  1.4210854828031606769834307141e-14,
  7.1054273952108527128773544799e-15,
  3.5527136913371136732984695340e-15,
  1.7763568435791203274733490144e-15,
  8.8817842109308159030960913863e-16,
  4.4408921031438133641977709402e-16,
  2.2204460507980419839993200942e-16,
  1.1102230251410661337205445699e-16,
  5.5511151248454812437237365905e-17,
  2.7755575621361241725816324538e-17,
  1.3877787809725232762839094906e-17,
  6.9388939045441536974460853262e-18,
  3.4694469521659226247442714961e-18,
  1.7347234760475765720489729699e-18,
  8.6736173801199337283420550673e-19,
  4.3368086900206504874970235659e-19,
  2.1684043449972197850139101683e-19,
  1.0842021724942414063012711165e-19,
  5.4210108624566454109187004043e-20,
  2.7105054312234688319546213119e-20,
  1.3552527156101164581485233996e-20,
  6.7762635780451890979952987415e-21,
  3.3881317890207968180857031004e-21,
  1.6940658945097991654064927471e-21,
  8.4703294725469983482469926091e-22,
  4.2351647362728333478622704833e-22,
  2.1175823681361947318442094398e-22,
  1.0587911840680233852265001539e-22,
  5.2939559203398703238139123029e-23,
  2.6469779601698529611341166842e-23,
  1.3234889800848990803094510250e-23,
  6.6174449004244040673552453323e-24,
  3.3087224502121715889469563843e-24,
  1.6543612251060756462299236771e-24,
  8.2718061255303444036711056167e-25,
  4.1359030627651609260093824555e-25,
  2.0679515313825767043959679193e-25,
  1.0339757656912870993284095591e-25,
  5.1698788284564313204101332166e-26,
  2.5849394142282142681277617708e-26,
  1.2924697071141066700381126118e-26,
  6.4623485355705318034380021611e-27,
  3.2311742677852653861348141180e-27,
  1.6155871338926325212060114057e-27,
  8.0779356694631620331587381863e-28,
  4.0389678347315808256222628129e-28,
  2.0194839173657903491587626465e-28,
  1.0097419586828951533619250700e-28,
  5.0487097934144756960847711725e-29,
  2.5243548967072378244674341938e-29,
  1.2621774483536189043753999660e-29,
  6.3108872417680944956826093943e-30,
  3.1554436208840472391098412184e-30,
  1.5777218104420236166444327830e-30,
  7.8886090522101180735205378276e-31
};


#define ZETA_NEG_TABLE_NMAX  99
#define ZETA_NEG_TABLE_SIZE  50
static double zeta_neg_int_table[ZETA_NEG_TABLE_SIZE] = {
 -0.083333333333333333333333333333,     /* zeta(-1) */
  0.008333333333333333333333333333,     /* zeta(-3) */
 -0.003968253968253968253968253968,     /* ...      */
  0.004166666666666666666666666667,
 -0.007575757575757575757575757576,
  0.021092796092796092796092796093,
 -0.083333333333333333333333333333,
  0.44325980392156862745098039216,
 -3.05395433027011974380395433027,
  26.4562121212121212121212121212,
 -281.460144927536231884057971014,
  3607.5105463980463980463980464,
 -54827.583333333333333333333333,
  974936.82385057471264367816092,
 -2.0052695796688078946143462272e+07,
  4.7238486772162990196078431373e+08,
 -1.2635724795916666666666666667e+10,
  3.8087931125245368811553022079e+11,
 -1.2850850499305083333333333333e+13,
  4.8241448354850170371581670362e+14,
 -2.0040310656516252738108421663e+16,
  9.1677436031953307756992753623e+17,
 -4.5979888343656503490437943262e+19,
  2.5180471921451095697089023320e+21,
 -1.5001733492153928733711440151e+23,
  9.6899578874635940656497942895e+24,
 -6.7645882379292820990945242302e+26,
  5.0890659468662289689766332916e+28,
 -4.1147288792557978697665486068e+30,
  3.5666582095375556109684574609e+32,
 -3.3066089876577576725680214670e+34,
  3.2715634236478716264211227016e+36,
 -3.4473782558278053878256455080e+38,
  3.8614279832705258893092720200e+40,
 -4.5892974432454332168863989006e+42,
  5.7775386342770431824884825688e+44,
 -7.6919858759507135167410075972e+46,
  1.0813635449971654696354033351e+49,
 -1.6029364522008965406067102346e+51,
  2.5019479041560462843656661499e+53,
 -4.1067052335810212479752045004e+55,
  7.0798774408494580617452972433e+57,
 -1.2804546887939508790190849756e+60,
  2.4267340392333524078020892067e+62,
 -4.8143218874045769355129570066e+64,
  9.9875574175727530680652777408e+66,
 -2.1645634868435185631335136160e+69,
  4.8962327039620553206849224516e+71,    /* ...        */
 -1.1549023923963519663954271692e+74,    /* zeta(-97)  */
  2.8382249570693706959264156336e+76     /* zeta(-99)  */
};


/* coefficients for Maclaurin summation in hzeta()
 * B_{2j}/(2j)!
 */
static double hzeta_c[15] = {
  1.00000000000000000000000000000,
  0.083333333333333333333333333333,
 -0.00138888888888888888888888888889,
  0.000033068783068783068783068783069,
 -8.2671957671957671957671957672e-07,
  2.0876756987868098979210090321e-08,
 -5.2841901386874931848476822022e-10,
  1.3382536530684678832826980975e-11,
 -3.3896802963225828668301953912e-13,
  8.5860620562778445641359054504e-15,
 -2.1748686985580618730415164239e-16,
  5.5090028283602295152026526089e-18,
 -1.3954464685812523340707686264e-19,
  3.5347070396294674716932299778e-21,
 -8.9535174270375468504026113181e-23
};

#define ETA_POS_TABLE_NMAX  100
static double eta_pos_int_table[ETA_POS_TABLE_NMAX+1] = {
0.50000000000000000000000000000,  /* eta(0) */
M_LN2,                            /* eta(1) */
0.82246703342411321823620758332,  /* ...    */
0.90154267736969571404980362113,
0.94703282949724591757650323447,
0.97211977044690930593565514355,
0.98555109129743510409843924448,
0.99259381992283028267042571313,
0.99623300185264789922728926008,
0.99809429754160533076778303185,
0.99903950759827156563922184570,
0.99951714349806075414409417483,
0.99975768514385819085317967871,
0.99987854276326511549217499282,
0.99993917034597971817095419226,
0.99996955121309923808263293263,
0.99998476421490610644168277496,
0.99999237829204101197693787224,
0.99999618786961011347968922641,
0.99999809350817167510685649297,
0.99999904661158152211505084256,
0.99999952325821554281631666433,
0.99999976161323082254789720494,
0.99999988080131843950322382485,
0.99999994039889239462836140314,
0.99999997019885696283441513311,
0.99999998509923199656878766181,
0.99999999254955048496351585274,
0.99999999627475340010872752767,
0.99999999813736941811218674656,
0.99999999906868228145397862728,
0.99999999953434033145421751469,
0.99999999976716989595149082282,
0.99999999988358485804603047265,
0.99999999994179239904531592388,
0.99999999997089618952980952258,
0.99999999998544809143388476396,
0.99999999999272404460658475006,
0.99999999999636202193316875550,
0.99999999999818101084320873555,
0.99999999999909050538047887809,
0.99999999999954525267653087357,
0.99999999999977262633369589773,
0.99999999999988631316532476488,
0.99999999999994315658215465336,
0.99999999999997157829090808339,
0.99999999999998578914539762720,
0.99999999999999289457268000875,
0.99999999999999644728633373609,
0.99999999999999822364316477861,
0.99999999999999911182158169283,
0.99999999999999955591079061426,
0.99999999999999977795539522974,
0.99999999999999988897769758908,
0.99999999999999994448884878594,
0.99999999999999997224442439010,
0.99999999999999998612221219410,
0.99999999999999999306110609673,
0.99999999999999999653055304826,
0.99999999999999999826527652409,
0.99999999999999999913263826204,
0.99999999999999999956631913101,
0.99999999999999999978315956551,
0.99999999999999999989157978275,
0.99999999999999999994578989138,
0.99999999999999999997289494569,
0.99999999999999999998644747284,
0.99999999999999999999322373642,
0.99999999999999999999661186821,
0.99999999999999999999830593411,
0.99999999999999999999915296705,
0.99999999999999999999957648353,
0.99999999999999999999978824176,
0.99999999999999999999989412088,
0.99999999999999999999994706044,
0.99999999999999999999997353022,
0.99999999999999999999998676511,
0.99999999999999999999999338256,
0.99999999999999999999999669128,
0.99999999999999999999999834564,
0.99999999999999999999999917282,
0.99999999999999999999999958641,
0.99999999999999999999999979320,
0.99999999999999999999999989660,
0.99999999999999999999999994830,
0.99999999999999999999999997415,
0.99999999999999999999999998708,
0.99999999999999999999999999354,
0.99999999999999999999999999677,
0.99999999999999999999999999838,
0.99999999999999999999999999919,
0.99999999999999999999999999960,
0.99999999999999999999999999980,
0.99999999999999999999999999990,
0.99999999999999999999999999995,
0.99999999999999999999999999997,
0.99999999999999999999999999999,
0.99999999999999999999999999999,
1.00000000000000000000000000000,
1.00000000000000000000000000000,
1.00000000000000000000000000000,
};


#define ETA_NEG_TABLE_NMAX  99
#define ETA_NEG_TABLE_SIZE  50
static double eta_neg_int_table[ETA_NEG_TABLE_SIZE] = {
 0.25000000000000000000000000000,   /* eta(-1) */
-0.12500000000000000000000000000,   /* eta(-3) */
 0.25000000000000000000000000000,   /* ...      */
-1.06250000000000000000000000000,
 7.75000000000000000000000000000,
-86.3750000000000000000000000000,
 1365.25000000000000000000000000,
-29049.0312500000000000000000000,
 800572.750000000000000000000000,
-2.7741322625000000000000000000e+7,
 1.1805291302500000000000000000e+9,
-6.0523980051687500000000000000e+10,
 3.6794167785377500000000000000e+12,
-2.6170760990658387500000000000e+14,
 2.1531418140800295250000000000e+16,
-2.0288775575173015930156250000e+18,
 2.1708009902623770590275000000e+20,
-2.6173826968455814932120125000e+22,
 3.5324148876863877826668602500e+24,
-5.3042033406864906641493838981e+26,
 8.8138218364311576767253114668e+28,
-1.6128065107490778547354654864e+31,
 3.2355470001722734208527794569e+33,
-7.0876727476537493198506645215e+35,
 1.6890450341293965779175629389e+38,
-4.3639690731216831157655651358e+40,
 1.2185998827061261322605065672e+43,
-3.6670584803153006180101262324e+45,
 1.1859898526302099104271449748e+48,
-4.1120769493584015047981746438e+50,
 1.5249042436787620309090168687e+53,
-6.0349693196941307074572991901e+55,
 2.5437161764210695823197691519e+58,
-1.1396923802632287851130360170e+61,
 5.4180861064753979196802726455e+63,
-2.7283654799994373847287197104e+66,
 1.4529750514918543238511171663e+69,
-8.1705519371067450079777183386e+71,
 4.8445781606678367790247757259e+74,
-3.0246694206649519336179448018e+77,
 1.9858807961690493054169047970e+80,
-1.3694474620720086994386818232e+83,
 9.9070382984295807826303785989e+85,
-7.5103780796592645925968460677e+88,
 5.9598418264260880840077992227e+91,
-4.9455988887500020399263196307e+94,
 4.2873596927020241277675775935e+97,
-3.8791952037716162900707994047e+100,
 3.6600317773156342245401829308e+103,
-3.5978775704117283875784869570e+106    /* eta(-99)  */
};


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/


int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s <= 1.0 || q <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    const double max_bits = 54.0;
    const double ln_term0 = -s * log(q);

    if(ln_term0 < GSL_LOG_DBL_MIN + 1.0) {
      UNDERFLOW_ERROR(result);
    }
    else if(ln_term0 > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR (result);
    }
    else if((s > max_bits && q < 1.0) || (s > 0.5*max_bits && q < 0.25)) {
      result->val = pow(q, -s);
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(s > 0.5*max_bits && q < 1.0) {
      const double p1 = pow(q, -s);
      const double p2 = pow(q/(1.0+q), s);
      const double p3 = pow(q/(2.0+q), s);
      result->val = p1 * (1.0 + p2 + p3);
      result->err = GSL_DBL_EPSILON * (0.5*s + 2.0) * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* Euler-Maclaurin summation formula
       * [Moshier, p. 400, with several typo corrections]
       */
      const int jmax = 12;
      const int kmax = 10;
      int j, k;
      const double pmax  = pow(kmax + q, -s);
      double scp = s;
      double pcp = pmax / (kmax + q);
      double ans = pmax*((kmax+q)/(s-1.0) + 0.5);

      for(k=0; k<kmax; k++) {
        ans += pow(k + q, -s);
      }

      for(j=0; j<=jmax; j++) {
        double delta = hzeta_c[j+1] * scp * pcp;
        ans += delta;
        if(fabs(delta/ans) < 0.5*GSL_DBL_EPSILON) break;
        scp *= (s+2*j+1)*(s+2*j+2);
        pcp /= (kmax + q)*(kmax + q);
      }

      result->val = ans;
      result->err = 2.0 * (jmax + 1.0) * GSL_DBL_EPSILON * fabs(ans);
      return GSL_SUCCESS;
    }
  }
}


int gsl_sf_zeta_e(const double s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s == 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(s >= 0.0) {
    return riemann_zeta_sgt0(s, result);
  }
  else {
    /* reflection formula, [Abramowitz+Stegun, 23.2.5] */

    gsl_sf_result zeta_one_minus_s;
    const int stat_zoms = riemann_zeta1ms_slt0(s, &zeta_one_minus_s);
    const double sin_term = (fmod(s,2.0) == 0.0) ? 0.0 : sin(0.5*M_PI*fmod(s,4.0))/M_PI;

    if(sin_term == 0.0) {
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(s > -170) {
      /* We have to be careful about losing digits
       * in calculating pow(2 Pi, s). The gamma
       * function is fine because we were careful
       * with that implementation.
       * We keep an array of (2 Pi)^(10 n).
       */
      const double twopi_pow[18] = { 1.0,
                                     9.589560061550901348e+007,
                                     9.195966217409212684e+015,
                                     8.818527036583869903e+023,
                                     8.456579467173150313e+031,
                                     8.109487671573504384e+039,
                                     7.776641909496069036e+047,
                                     7.457457466828644277e+055,
                                     7.151373628461452286e+063,
                                     6.857852693272229709e+071,
                                     6.576379029540265771e+079,
                                     6.306458169130020789e+087,
                                     6.047615938853066678e+095,
                                     5.799397627482402614e+103,
                                     5.561367186955830005e+111,
                                     5.333106466365131227e+119,
                                     5.114214477385391780e+127,
                                     4.904306689854036836e+135
                                    };
      const int n = floor((-s)/10.0);
      const double fs = s + 10.0*n;
      const double p = pow(2.0*M_PI, fs) / twopi_pow[n];

      gsl_sf_result g;
      const int stat_g = gsl_sf_gamma_e(1.0-s, &g);
      result->val  = p * g.val * sin_term * zeta_one_minus_s.val;
      result->err  = fabs(p * g.val * sin_term) * zeta_one_minus_s.err;
      result->err += fabs(p * sin_term * zeta_one_minus_s.val) * g.err;
      result->err += GSL_DBL_EPSILON * (fabs(s)+2.0) * fabs(result->val);
      return GSL_ERROR_SELECT_2(stat_g, stat_zoms);
    }
    else {
      /* The actual zeta function may or may not
       * overflow here. But we have no easy way
       * to calculate it when the prefactor(s)
       * overflow. Trying to use log's and exp
       * is no good because we loose a couple
       * digits to the exp error amplification.
       * When we gather a little more patience,
       * we can implement something here. Until
       * then just give up.
       */
      OVERFLOW_ERROR(result);
    }
  }
}


int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n < 0) {
    if(!GSL_IS_ODD(n)) {
      result->val = 0.0; /* exactly zero at even negative integers */
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ZETA_NEG_TABLE_NMAX) {
      result->val = zeta_neg_int_table[-(n+1)/2];
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      return gsl_sf_zeta_e((double)n, result);
    }
  }
  else if(n == 1){
    DOMAIN_ERROR(result);
  }
  else if(n <= ZETA_POS_TABLE_NMAX){
    result->val = 1.0 + zetam1_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
}


int gsl_sf_zetam1_e(const double s, gsl_sf_result * result)
{
  if(s <= 5.0)
  {
    int stat = gsl_sf_zeta_e(s, result);
    result->val = result->val - 1.0;
    return stat;
  }
  else if(s < 15.0)
  {
    return riemann_zeta_minus_1_intermediate_s(s, result);
  }
  else
  {
    return riemann_zeta_minus1_large_s(s, result);
  }
}


int gsl_sf_zetam1_int_e(const int n, gsl_sf_result * result)
{
  if(n < 0) {
    if(!GSL_IS_ODD(n)) {
      result->val = -1.0; /* at even negative integers zetam1 == -1 since zeta is exactly zero */
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ZETA_NEG_TABLE_NMAX) {
      result->val = zeta_neg_int_table[-(n+1)/2] - 1.0;
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      /* could use gsl_sf_zetam1_e here but subtracting 1 makes no difference
         for such large values, so go straight to the result */
      return gsl_sf_zeta_e((double)n, result);
    }
  }
  else if(n == 1){
    DOMAIN_ERROR(result);
  }
  else if(n <= ZETA_POS_TABLE_NMAX){
    result->val = zetam1_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return gsl_sf_zetam1_e(n, result);
  }
}


int gsl_sf_eta_int_e(int n, gsl_sf_result * result)
{
  if(n > ETA_POS_TABLE_NMAX) {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(n >= 0) {
    result->val = eta_pos_int_table[n];
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* n < 0 */

    if(!GSL_IS_ODD(n)) {
      /* exactly zero at even negative integers */
      result->val = 0.0;
      result->err = 0.0;
      return GSL_SUCCESS;
    }
    else if(n > -ETA_NEG_TABLE_NMAX) {
      result->val = eta_neg_int_table[-(n+1)/2];
      result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else {
      gsl_sf_result z;
      gsl_sf_result p;
      int stat_z = gsl_sf_zeta_int_e(n, &z);
      int stat_p = gsl_sf_exp_e((1.0-n)*M_LN2, &p);
      int stat_m = gsl_sf_multiply_e(-p.val, z.val, result);
      result->err  = fabs(p.err * (M_LN2*(1.0-n)) * z.val) + z.err * fabs(p.val);
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_ERROR_SELECT_3(stat_m, stat_p, stat_z);
    }
  }
}


int gsl_sf_eta_e(const double s, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(s > 100.0) {
    result->val = 1.0;
    result->err = GSL_DBL_EPSILON;
    return GSL_SUCCESS;
  }
  else if(fabs(s-1.0) < 10.0*GSL_ROOT5_DBL_EPSILON) {
    double del = s-1.0;
    double c0  = M_LN2;
    double c1  = M_LN2 * (M_EULER - 0.5*M_LN2);
    double c2  = -0.0326862962794492996;
    double c3  =  0.0015689917054155150;
    double c4  =  0.00074987242112047532;
    result->val = c0 + del * (c1 + del * (c2 + del * (c3 + del * c4)));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result z;
    gsl_sf_result p;
    int stat_z = gsl_sf_zeta_e(s, &z);
    int stat_p = gsl_sf_exp_e((1.0-s)*M_LN2, &p);
    int stat_m = gsl_sf_multiply_e(1.0-p.val, z.val, result);
    result->err  = fabs(p.err * (M_LN2*(1.0-s)) * z.val) + z.err * fabs(p.val);
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_3(stat_m, stat_p, stat_z);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_zeta(const double s)
{
  EVAL_RESULT(gsl_sf_zeta_e(s, &result));
}

double gsl_sf_hzeta(const double s, const double a)
{
  EVAL_RESULT(gsl_sf_hzeta_e(s, a, &result));
}

double gsl_sf_zeta_int(const int s)
{
  EVAL_RESULT(gsl_sf_zeta_int_e(s, &result));
}

double gsl_sf_zetam1(const double s)
{
  EVAL_RESULT(gsl_sf_zetam1_e(s, &result));
}

double gsl_sf_zetam1_int(const int s)
{
  EVAL_RESULT(gsl_sf_zetam1_int_e(s, &result));
}

double gsl_sf_eta_int(const int s)
{
  EVAL_RESULT(gsl_sf_eta_int_e(s, &result));
}

double gsl_sf_eta(const double s)
{
  EVAL_RESULT(gsl_sf_eta_e(s, &result));
}
/* end inlined source: specfunc/zeta.c */
/* begin inlined source: specfunc/psi.c */
/* specfunc/psi.c
 *
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2005, 2006 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author: G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_exp.h */
/* specfunc/gsl_sf_exp.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_EXP_H__
#define __GSL_SF_EXP_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
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

/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e(const double x, gsl_sf_result * result);
double gsl_sf_exp(const double x);


/* Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_exp_mult(const double x, const double y);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
double gsl_sf_expm1(const double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
double gsl_sf_exprel(const double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
double gsl_sf_exprel_2(const double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_exprel_n(const int n, const double x);

int gsl_sf_exprel_n_CF_e(const double n, const double x, gsl_sf_result * result);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);

/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);

__END_DECLS

#endif /* __GSL_SF_EXP_H__ */

/* end inlined header: gsl/gsl_sf_exp.h */
/* begin inlined header: gsl/gsl_sf_gamma.h */
/* specfunc/gsl_sf_gamma.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_GAMMA_H__
#define __GSL_SF_GAMMA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
double gsl_sf_lngamma(const double x);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
double gsl_sf_gamma(const double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
double gsl_sf_gammastar(const double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
double gsl_sf_gammainv(const double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_taylorcoeff(const int n, const double x);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_fact(const unsigned int n);


/* n!! = n(n-2)(n-4) ...
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_doublefact(const unsigned int n);


/* log(n!)
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lnfact(const unsigned int n);


/* log(n!!)
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lndoublefact(const unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochhammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_lnpoch(const double a, const double x);


/* Logarithm of Pochhammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_e(const double a, const double x, gsl_sf_result * result, double * sgn);


/* Pochhammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_poch(const double a, const double x);


/* Relative Pochhammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_pochrel(const double a, const double x);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a >= 0, x >= 0
 *   Q(a,0) := 1
 *   Q(0,x) := 0, x != 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_Q(const double a, const double x);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_P(const double a, const double x);


/* Non-normalized Incomplete Gamma Function
 *
 * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * x >= 0.0
 *   Gamma(a, 0) := Gamma(a)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc(const double a, const double x);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_lnbeta(const double a, const double b);

int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_beta(const double a, const double b);


/* Normalized Incomplete Beta Function
 * B_x(a,b)/B(a,b)
 *
 * a > 0, b > 0, 0 <= x <= 1
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_beta_inc(const double a, const double b, const double x);


/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297

__END_DECLS

#endif /* __GSL_SF_GAMMA_H__ */

/* end inlined header: gsl/gsl_sf_gamma.h */
/* begin inlined header: gsl/gsl_sf_zeta.h */
/* specfunc/gsl_sf_zeta.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_ZETA_H__
#define __GSL_SF_ZETA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Riemann Zeta Function
 * zeta(n) = Sum[ k^(-n), {k,1,Infinity} ]
 *
 * n=integer, n != 1
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_int_e(const int n, gsl_sf_result * result);
double gsl_sf_zeta_int(const int n);


/* Riemann Zeta Function
 * zeta(x) = Sum[ k^(-s), {k,1,Infinity} ], s != 1.0
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zeta_e(const double s, gsl_sf_result * result);
double gsl_sf_zeta(const double s);


/* Riemann Zeta Function minus 1
 *   useful for evaluating the fractional part
 *   of Riemann zeta for large argument
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zetam1_e(const double s, gsl_sf_result * result);
double gsl_sf_zetam1(const double s);


/* Riemann Zeta Function minus 1 for integer arg
 *   useful for evaluating the fractional part
 *   of Riemann zeta for large argument
 *
 * s != 1.0
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_zetam1_int_e(const int s, gsl_sf_result * result);
double gsl_sf_zetam1_int(const int s);


/* Hurwitz Zeta Function
 * zeta(s,q) = Sum[ (k+q)^(-s), {k,0,Infinity} ]
 *
 * s > 1.0, q > 0.0
 * exceptions: GSL_EDOM, GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_hzeta_e(const double s, const double q, gsl_sf_result * result);
double gsl_sf_hzeta(const double s, const double q);


/* Eta Function
 * eta(n) = (1-2^(1-n)) zeta(n)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_int_e(int n, gsl_sf_result * result);
double gsl_sf_eta_int(const int n);


/* Eta Function
 * eta(s) = (1-2^(1-s)) zeta(s)
 *
 * exceptions: GSL_EUNDRFLW, GSL_EOVRFLW
 */
int gsl_sf_eta_e(const double s, gsl_sf_result * result);
double gsl_sf_eta(const double s);


__END_DECLS

#endif /* __GSL_SF_ZETA_H__ */

/* end inlined header: gsl/gsl_sf_zeta.h */
/* begin inlined header: gsl/gsl_sf_psi.h */
/* specfunc/gsl_sf_psi.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_PSI_H__
#define __GSL_SF_PSI_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Poly-Gamma Functions
 *
 * psi(m,x) := (d/dx)^m psi(0,x) = (d/dx)^{m+1} log(gamma(x))
 */


/* Di-Gamma Function  psi(n) = psi(0,n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_int(const int n);


/* Di-Gamma Function psi(x) = psi(0, x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi(const double x);


/* Di-Gamma Function Re[psi(1 + I y)]
 *
 * exceptions: none
 */
int     gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result);
double  gsl_sf_psi_1piy(const double y);


/* Di-Gamma Function psi(z) for general complex argument z = x + iy
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  );


/* Tri-Gamma Function psi^(1)(n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_1_int(const int n);


/* Tri-Gamma Function psi^(1)(x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_1_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi_1(const double x);


/* Poly-Gamma Function psi^(n)(x)
 *
 * n >= 0, x > 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result);
double  gsl_sf_psi_n(const int n, const double x);


__END_DECLS

#endif /* __GSL_SF_PSI_H__ */

/* end inlined header: gsl/gsl_sf_psi.h */
/* begin inlined header: gsl/gsl_complex_math.h */
/* complex/gsl_complex_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Jorma Olavi T�htinen, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_MATH_H__
#define __GSL_COMPLEX_MATH_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_complex.h */
/* complex/gsl_complex.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 * Copyright (C) 2020, 2021 Patrick Alken
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_COMPLEX_H__
#define __GSL_COMPLEX_H__

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


/* two consecutive built-in types as a complex number */
typedef double *       gsl_complex_packed ;
typedef float *        gsl_complex_packed_float  ;
typedef long double *  gsl_complex_packed_long_double ;

typedef const double *       gsl_const_complex_packed ;
typedef const float *        gsl_const_complex_packed_float  ;
typedef const long double *  gsl_const_complex_packed_long_double ;


/* 2N consecutive built-in types as N complex numbers */
typedef double *       gsl_complex_packed_array ;
typedef float *        gsl_complex_packed_array_float  ;
typedef long double *  gsl_complex_packed_array_long_double ;

typedef const double *       gsl_const_complex_packed_array ;
typedef const float *        gsl_const_complex_packed_array_float  ;
typedef const long double *  gsl_const_complex_packed_array_long_double ;


/* Yes... this seems weird. Trust us. The point is just that
   sometimes you want to make it obvious that something is
   an output value. The fact that it lacks a 'const' may not
   be enough of a clue for people in some contexts.
 */
typedef double *       gsl_complex_packed_ptr ;
typedef float *        gsl_complex_packed_float_ptr  ;
typedef long double *  gsl_complex_packed_long_double_ptr ;

typedef const double *       gsl_const_complex_packed_ptr ;
typedef const float *        gsl_const_complex_packed_float_ptr  ;
typedef const long double *  gsl_const_complex_packed_long_double_ptr ;

/*
 * If <complex.h> is included, use the C99 complex type.  Otherwise
 * define a type bit-compatible with C99 complex. The GSL_REAL and GSL_IMAG
 * macros require C11 functionality also (_Generic)
 */

/* older gcc compilers claim to be C11 compliant but do not support _Generic */
#if defined(__GNUC__) && (__GNUC__ < 7)
#  define GSL_COMPLEX_LEGACY 1
#endif

#if !defined(GSL_COMPLEX_LEGACY) &&          \
     defined(_Complex_I) &&                  \
     defined(complex) &&                     \
     defined(I) &&                           \
     defined(__STDC__) && (__STDC__ == 1) && \
     defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L) /* C11 */

#  define GSL_COMPLEX_DEFINE(R, C) typedef R _Complex C ;

#  define GSL_COMPLEX_P(zp)        (&(zp))
#  define GSL_COMPLEX_EQ(z1,z2)    ((z1) == (z2))
#  define GSL_SET_COMPLEX(zp,x,y)  (*(zp) = (x) + I * (y))

#  define GSL_REAL(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[0])

#  define GSL_IMAG(z)              (_Generic((z),                             \
                                     complex float       : ((float *) &(z)),  \
                                     complex double      : ((double *) &(z)), \
                                     complex long double : ((long double *) &(z)))[1])

#  define GSL_COMPLEX_P_REAL(zp)   GSL_REAL(*(zp))
#  define GSL_COMPLEX_P_IMAG(zp)   GSL_IMAG(*(zp))
#  define GSL_SET_REAL(zp,x)       do { GSL_REAL(*(zp)) = (x); } while(0)
#  define GSL_SET_IMAG(zp,x)       do { GSL_IMAG(*(zp)) = (x); } while(0)

#else /* legacy complex definitions */

/*
 * According to the C17 standard, 6.2.5 paragraph 13:
 *
 * "Each complex type has the same representation and alignment requirements
 * as an array type containing exactly two elements of the corresponding real
 * type; the first element is equal to the real part, and the second element to
 * the imaginary part, of the complex number."
 */

/*#  define GSL_COMPLEX_DEFINE(R, C) typedef R C[2]*/
#  define GSL_COMPLEX_DEFINE(R, C) typedef struct { R dat[2]; } C ;

#  define GSL_REAL(z)              ((z).dat[0])
#  define GSL_IMAG(z)              ((z).dat[1])
#  define GSL_COMPLEX_P(zp)        ((zp)->dat)
#  define GSL_COMPLEX_P_REAL(zp)   ((zp)->dat[0])
#  define GSL_COMPLEX_P_IMAG(zp)   ((zp)->dat[1])
#  define GSL_COMPLEX_EQ(z1,z2)    (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

#  define GSL_SET_COMPLEX(zp,x,y)  do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#  define GSL_SET_REAL(zp,x)       do {(zp)->dat[0]=(x);} while(0)
#  define GSL_SET_IMAG(zp,y)       do {(zp)->dat[1]=(y);} while(0)

#endif

GSL_COMPLEX_DEFINE(double, gsl_complex)
GSL_COMPLEX_DEFINE(long double, gsl_complex_long_double)
GSL_COMPLEX_DEFINE(float, gsl_complex_float)

#define GSL_SET_COMPLEX_PACKED(zp,n,x,y) do {*((zp)+2*(n))=(x); *((zp)+(2*(n)+1))=(y);} while(0)

__END_DECLS

#endif /* __GSL_COMPLEX_H__ */

/* end inlined header: gsl/gsl_complex.h */
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS }
#else
#define __BEGIN_DECLS           /* empty */
#define __END_DECLS             /* empty */
#endif

__BEGIN_DECLS

/* Complex numbers */

gsl_complex gsl_complex_polar (double r, double theta); /* r= r e^(i theta) */

INLINE_DECL gsl_complex gsl_complex_rect (double x, double y);  /* r= real+i*imag */

#ifdef HAVE_INLINE
INLINE_FUN gsl_complex
gsl_complex_rect (double x, double y)
{                               /* return z = x + i y */
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
#endif

#define GSL_COMPLEX_ONE (gsl_complex_rect(1.0,0.0))
#define GSL_COMPLEX_ZERO (gsl_complex_rect(0.0,0.0))
#define GSL_COMPLEX_NEGONE (gsl_complex_rect(-1.0,0.0))

/* Properties of complex numbers */

double gsl_complex_arg (gsl_complex z); /* return arg(z), -pi< arg(z) <=+pi */
double gsl_complex_abs (gsl_complex z);   /* return |z|   */
double gsl_complex_abs2 (gsl_complex z);  /* return |z|^2 */
double gsl_complex_logabs (gsl_complex z); /* return log|z| */

/* Complex arithmetic operators */

gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b);  /* r=a+b */
gsl_complex gsl_complex_sub (gsl_complex a, gsl_complex b);  /* r=a-b */
gsl_complex gsl_complex_mul (gsl_complex a, gsl_complex b);  /* r=a*b */
gsl_complex gsl_complex_div (gsl_complex a, gsl_complex b);  /* r=a/b */

gsl_complex gsl_complex_add_real (gsl_complex a, double x);  /* r=a+x */
gsl_complex gsl_complex_sub_real (gsl_complex a, double x);  /* r=a-x */
gsl_complex gsl_complex_mul_real (gsl_complex a, double x);  /* r=a*x */
gsl_complex gsl_complex_div_real (gsl_complex a, double x);  /* r=a/x */

gsl_complex gsl_complex_add_imag (gsl_complex a, double y);  /* r=a+iy */
gsl_complex gsl_complex_sub_imag (gsl_complex a, double y);  /* r=a-iy */
gsl_complex gsl_complex_mul_imag (gsl_complex a, double y);  /* r=a*iy */
gsl_complex gsl_complex_div_imag (gsl_complex a, double y);  /* r=a/iy */

gsl_complex gsl_complex_conjugate (gsl_complex z);  /* r=conj(z) */
gsl_complex gsl_complex_inverse (gsl_complex a);    /* r=1/a */
gsl_complex gsl_complex_negative (gsl_complex a);    /* r=-a */

/* Elementary Complex Functions */

gsl_complex gsl_complex_sqrt (gsl_complex z);  /* r=sqrt(z) */
gsl_complex gsl_complex_sqrt_real (double x);  /* r=sqrt(x) (x<0 ok) */

gsl_complex gsl_complex_pow (gsl_complex a, gsl_complex b);  /* r=a^b */
gsl_complex gsl_complex_pow_real (gsl_complex a, double b);  /* r=a^b */

gsl_complex gsl_complex_exp (gsl_complex a);    /* r=exp(a) */
gsl_complex gsl_complex_log (gsl_complex a);    /* r=log(a) (base e) */
gsl_complex gsl_complex_log10 (gsl_complex a);  /* r=log10(a) (base 10) */
gsl_complex gsl_complex_log_b (gsl_complex a, gsl_complex b);   /* r=log_b(a) (base=b) */

/* Complex Trigonometric Functions */

gsl_complex gsl_complex_sin (gsl_complex a);  /* r=sin(a) */
gsl_complex gsl_complex_cos (gsl_complex a);  /* r=cos(a) */
gsl_complex gsl_complex_sec (gsl_complex a);  /* r=sec(a) */
gsl_complex gsl_complex_csc (gsl_complex a);  /* r=csc(a) */
gsl_complex gsl_complex_tan (gsl_complex a);  /* r=tan(a) */
gsl_complex gsl_complex_cot (gsl_complex a);  /* r=cot(a) */

/* Inverse Complex Trigonometric Functions */

gsl_complex gsl_complex_arcsin (gsl_complex a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arcsin_real (double a);  /* r=arcsin(a) */
gsl_complex gsl_complex_arccos (gsl_complex a);  /* r=arccos(a) */
gsl_complex gsl_complex_arccos_real (double a);  /* r=arccos(a) */
gsl_complex gsl_complex_arcsec (gsl_complex a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arcsec_real (double a);  /* r=arcsec(a) */
gsl_complex gsl_complex_arccsc (gsl_complex a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arccsc_real (double a);  /* r=arccsc(a) */
gsl_complex gsl_complex_arctan (gsl_complex a);  /* r=arctan(a) */
gsl_complex gsl_complex_arccot (gsl_complex a);  /* r=arccot(a) */

/* Complex Hyperbolic Functions */

gsl_complex gsl_complex_sinh (gsl_complex a);  /* r=sinh(a) */
gsl_complex gsl_complex_cosh (gsl_complex a);  /* r=coshh(a) */
gsl_complex gsl_complex_sech (gsl_complex a);  /* r=sech(a) */
gsl_complex gsl_complex_csch (gsl_complex a);  /* r=csch(a) */
gsl_complex gsl_complex_tanh (gsl_complex a);  /* r=tanh(a) */
gsl_complex gsl_complex_coth (gsl_complex a);  /* r=coth(a) */

/* Inverse Complex Hyperbolic Functions */

gsl_complex gsl_complex_arcsinh (gsl_complex a);  /* r=arcsinh(a) */
gsl_complex gsl_complex_arccosh (gsl_complex a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arccosh_real (double a);  /* r=arccosh(a) */
gsl_complex gsl_complex_arcsech (gsl_complex a);  /* r=arcsech(a) */
gsl_complex gsl_complex_arccsch (gsl_complex a);  /* r=arccsch(a) */
gsl_complex gsl_complex_arctanh (gsl_complex a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arctanh_real (double a);  /* r=arctanh(a) */
gsl_complex gsl_complex_arccoth (gsl_complex a);  /* r=arccoth(a) */

__END_DECLS

#endif /* __GSL_COMPLEX_MATH_H__ */

/* end inlined header: gsl/gsl_complex_math.h */
#include <stdio.h>

/* inlined: error.h */

/* inlined: chebyshev.h */
/* inlined once already: specfunc/cheb_eval.c */

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* Chebyshev fit for f(y) = Re(Psi(1+Iy)) + M_EULER - y^2/(1+y^2) - y^2/(2(4+y^2))
 * 1 < y < 10
 *   ==>
 * y(x) = (9x + 11)/2,  -1 < x < 1
 * x(y) = (2y - 11)/9
 *
 * g(x) := f(y(x))
 */
static double r1py_data[] = {
   1.59888328244976954803168395603,
   0.67905625353213463845115658455,
  -0.068485802980122530009506482524,
  -0.005788184183095866792008831182,
   0.008511258167108615980419855648,
  -0.004042656134699693434334556409,
   0.001352328406159402601778462956,
  -0.000311646563930660566674525382,
   0.000018507563785249135437219139,
   0.000028348705427529850296492146,
  -0.000019487536014574535567541960,
   8.0709788710834469408621587335e-06,
  -2.2983564321340518037060346561e-06,
   3.0506629599604749843855962658e-07,
   1.3042238632418364610774284846e-07,
  -1.2308657181048950589464690208e-07,
   5.7710855710682427240667414345e-08,
  -1.8275559342450963966092636354e-08,
   3.1020471300626589420759518930e-09,
   6.8989327480593812470039430640e-10,
  -8.7182290258923059852334818997e-10,
   4.4069147710243611798213548777e-10,
  -1.4727311099198535963467200277e-10,
   2.7589682523262644748825844248e-11,
   4.1871826756975856411554363568e-12,
  -6.5673460487260087541400767340e-12,
   3.4487900886723214020103638000e-12,
  -1.1807251417448690607973794078e-12,
   2.3798314343969589258709315574e-13,
   2.1663630410818831824259465821e-15
};
static cheb_series r1py_cs = {
  r1py_data,
  29,
  -1,1,
  18
};


/* Chebyshev fits from SLATEC code for psi(x)

 Series for PSI        on the interval  0.         to  1.00000D+00
                                       with weighted error   2.03E-17
                                        log weighted error  16.69
                              significant figures required  16.39
                                   decimal places required  17.37

 Series for APSI       on the interval  0.         to  2.50000D-01
                                       with weighted error   5.54E-17
                                        log weighted error  16.26
                              significant figures required  14.42
                                   decimal places required  16.86

*/

static double psics_data[23] = {
  -.038057080835217922,
   .491415393029387130,
  -.056815747821244730,
   .008357821225914313,
  -.001333232857994342,
   .000220313287069308,
  -.000037040238178456,
   .000006283793654854,
  -.000001071263908506,
   .000000183128394654,
  -.000000031353509361,
   .000000005372808776,
  -.000000000921168141,
   .000000000157981265,
  -.000000000027098646,
   .000000000004648722,
  -.000000000000797527,
   .000000000000136827,
  -.000000000000023475,
   .000000000000004027,
  -.000000000000000691,
   .000000000000000118,
  -.000000000000000020
};
static cheb_series psi_cs = {
  psics_data,
  22,
  -1, 1,
  17
};

static double apsics_data[16] = {
  -.0204749044678185,
  -.0101801271534859,
   .0000559718725387,
  -.0000012917176570,
   .0000000572858606,
  -.0000000038213539,
   .0000000003397434,
  -.0000000000374838,
   .0000000000048990,
  -.0000000000007344,
   .0000000000001233,
  -.0000000000000228,
   .0000000000000045,
  -.0000000000000009,
   .0000000000000002,
  -.0000000000000000
};
static cheb_series apsi_cs = {
  apsics_data,
  15,
  -1, 1,
  9
};

#define PSI_TABLE_NMAX 100
static double psi_table[PSI_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(0) */
 -M_EULER,                          /* psi(1) */
  0.42278433509846713939348790992,  /* ...    */
  0.92278433509846713939348790992,
  1.25611766843180047272682124325,
  1.50611766843180047272682124325,
  1.70611766843180047272682124325,
  1.87278433509846713939348790992,
  2.01564147795560999653634505277,
  2.14064147795560999653634505277,
  2.25175258906672110764745616389,
  2.35175258906672110764745616389,
  2.44266167997581201673836525479,
  2.52599501330914535007169858813,
  2.60291809023222227314862166505,
  2.67434666166079370172005023648,
  2.74101332832746036838671690315,
  2.80351332832746036838671690315,
  2.86233685773922507426906984432,
  2.91789241329478062982462539988,
  2.97052399224214905087725697883,
  3.02052399224214905087725697883,
  3.06814303986119666992487602645,
  3.11359758531574212447033057190,
  3.15707584618530734186163491973,
  3.1987425128519740085283015864,
  3.2387425128519740085283015864,
  3.2772040513135124700667631249,
  3.3142410883505495071038001619,
  3.3499553740648352213895144476,
  3.3844381326855248765619282407,
  3.4177714660188582098952615740,
  3.4500295305349872421533260902,
  3.4812795305349872421533260902,
  3.5115825608380175451836291205,
  3.5409943255438998981248055911,
  3.5695657541153284695533770196,
  3.5973435318931062473311547974,
  3.6243705589201332743581818244,
  3.6506863483938174848844976139,
  3.6763273740348431259101386396,
  3.7013273740348431259101386396,
  3.7257176179372821503003825420,
  3.7495271417468059598241920658,
  3.7727829557002943319172153216,
  3.7955102284275670591899425943,
  3.8177324506497892814121648166,
  3.8394715810845718901078169905,
  3.8607481768292527411716467777,
  3.8815815101625860745049801110,
  3.9019896734278921969539597029,
  3.9219896734278921969539597029,
  3.9415975165651470989147440166,
  3.9608282857959163296839747858,
  3.9796962103242182164764276160,
  3.9982147288427367349949461345,
  4.0163965470245549168131279527,
  4.0342536898816977739559850956,
  4.0517975495308205809735289552,
  4.0690389288411654085597358518,
  4.0859880813835382899156680552,
  4.1026547480502049565823347218,
  4.1190481906731557762544658694,
  4.1351772229312202923834981274,
  4.1510502388042361653993711433,
  4.1666752388042361653993711433,
  4.1820598541888515500147557587,
  4.1972113693403667015299072739,
  4.2121367424746950597388624977,
  4.2268426248276362362094507330,
  4.2413353784508246420065521823,
  4.2556210927365389277208378966,
  4.2697055997787924488475984600,
  4.2835944886676813377364873489,
  4.2972931188046676391063503626,
  4.3108066323181811526198638761,
  4.3241399656515144859531972094,
  4.3372978603883565912163551041,
  4.3502848733753695782293421171,
  4.3631053861958823987421626300,
  4.3757636140439836645649474401,
  4.3882636140439836645649474401,
  4.4006092930563293435772931191,
  4.4128044150075488557724150703,
  4.4248526077786331931218126607,
  4.4367573696833950978837174226,
  4.4485220755657480390601880108,
  4.4601499825424922251066996387,
  4.4716442354160554434975042364,
  4.4830078717796918071338678728,
  4.4942438268358715824147667492,
  4.5053549379469826935258778603,
  4.5163439489359936825368668713,
  4.5272135141533849868846929582,
  4.5379662023254279976373811303,
  4.5486045001977684231692960239,
  4.5591308159872421073798223397,
  4.5695474826539087740464890064,
  4.5798567610044242379640147796,
  4.5900608426370772991885045755,
  4.6001618527380874001986055856
};


#define PSI_1_TABLE_NMAX 100
static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
  0.0,  /* Infinity */              /* psi(1,0) */
  M_PI*M_PI/6.0,                    /* psi(1,1) */
  0.644934066848226436472415,       /* ...      */
  0.394934066848226436472415,
  0.2838229557371153253613041,
  0.2213229557371153253613041,
  0.1813229557371153253613041,
  0.1535451779593375475835263,
  0.1331370146940314251345467,
  0.1175120146940314251345467,
  0.1051663356816857461222010,
  0.0951663356816857461222010,
  0.0869018728717683907503002,
  0.0799574284273239463058557,
  0.0740402686640103368384001,
  0.0689382278476838062261552,
  0.0644937834032393617817108,
  0.0605875334032393617817108,
  0.0571273257907826143768665,
  0.0540409060376961946237801,
  0.0512708229352031198315363,
  0.0487708229352031198315363,
  0.0465032492390579951149830,
  0.0444371335365786562720078,
  0.0425467743683366902984728,
  0.0408106632572255791873617,
  0.0392106632572255791873617,
  0.0377313733163971768204978,
  0.0363596312039143235969038,
  0.0350841209998326909438426,
  0.0338950603577399442137594,
  0.0327839492466288331026483,
  0.0317433665203020901265817,
  0.03076680402030209012658168,
  0.02984853037475571730748159,
  0.02898347847164153045627052,
  0.02816715194102928555831133,
  0.02739554700275768062003973,
  0.02666508681283803124093089,
  0.02597256603721476254286995,
  0.02531510384129102815759710,
  0.02469010384129102815759710,
  0.02409521984367056414807896,
  0.02352832641963428296894063,
  0.02298749353699501850166102,
  0.02247096461137518379091722,
  0.02197713745088135663042339,
  0.02150454765882086513703965,
  0.02105185413233829383780923,
  0.02061782635456051606003145,
  0.02020133322669712580597065,
  0.01980133322669712580597065,
  0.01941686571420193164987683,
  0.01904704322899483105816086,
  0.01869104465298913508094477,
  0.01834810912486842177504628,
  0.01801753061247172756017024,
  0.01769865306145131939690494,
  0.01739086605006319997554452,
  0.01709360088954001329302371,
  0.01680632711763538818529605,
  0.01652854933985761040751827,
  0.01625980437882562975715546,
  0.01599965869724394401313881,
  0.01574770606433893015574400,
  0.01550356543933893015574400,
  0.01526687904880638577704578,
  0.01503731063741979257227076,
  0.01481454387422086185273411,
  0.01459828089844231513993134,
  0.01438824099085987447620523,
  0.01418415935820681325171544,
  0.01398578601958352422176106,
  0.01379288478501562298719316,
  0.01360523231738567365335942,
  0.01342261726990576130858221,
  0.01324483949212798353080444,
  0.01307170929822216635628920,
  0.01290304679189732236910755,
  0.01273868124291638877278934,
  0.01257845051066194236996928,
  0.01242220051066194236996928,
  0.01226978472038606978956995,
  0.01212106372098095378719041,
  0.01197590477193174490346273,
  0.01183418141592267460867815,
  0.01169577311142440471248438,
  0.01156056489076458859566448,
  0.01142844704164317229232189,
  0.01129931481023821361463594,
  0.01117306812421372175754719,
  0.01104961133409026496742374,
  0.01092885297157366069257770,
  0.01081070552355853781923177,
  0.01069508522063334415522437,
  0.01058191183901270133041676,
  0.01047110851491297833872701,
  0.01036260157046853389428257,
  0.01025632035036012704977199,  /* ...        */
  0.01015219706839427948625679,  /* psi(1,99)  */
  0.01005016666333357139524567   /* psi(1,100) */
};


/* digamma for x both positive and negative; we do both
 * cases here because of the way we use even/odd parts
 * of the function
 */
static int
psi_x(const double x, gsl_sf_result * result)
{
  const double y = fabs(x);

  if(x == 0.0 || x == -1.0 || x == -2.0) {
    DOMAIN_ERROR(result);
  }
  else if(y >= 2.0) {
    const double t = 8.0/(y*y)-1.0;
    gsl_sf_result result_c;
    cheb_eval_e(&apsi_cs, t, &result_c);
    if(x < 0.0) {
      const double s = sin(M_PI*x);
      const double c = cos(M_PI*x);
      if(fabs(s) < 2.0*GSL_SQRT_DBL_MIN) {
        DOMAIN_ERROR(result);
      }
      else {
        result->val  = log(y) - 0.5/x + result_c.val - M_PI * c/s;
        result->err  = M_PI*fabs(x)*GSL_DBL_EPSILON/(s*s);
        result->err += result_c.err;
        result->err += GSL_DBL_EPSILON * fabs(result->val);
        return GSL_SUCCESS;
      }
    }
    else {
      result->val  = log(y) - 0.5/x + result_c.val;
      result->err  = result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
  }
  else { /* -2 < x < 2 */
    gsl_sf_result result_c;

    if(x < -1.0) { /* x = -2 + v */
      const double v  = x + 2.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/(x+1.0);
      const double t3 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);

      result->val  = -(t1 + t2 + t3) + result_c.val;
      result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)) + fabs(x/(t3*t3)));
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < 0.0) { /* x = -1 + v */
      const double v  = x + 1.0;
      const double t1 = 1.0/x;
      const double t2 = 1.0/v;
      cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);

      result->val  = -(t1 + t2) + result_c.val;
      result->err  = GSL_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)));
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x < 1.0) { /* x = v */
      const double t1 = 1.0/x;
      cheb_eval_e(&psi_cs, 2.0*x-1.0, &result_c);

      result->val  = -t1 + result_c.val;
      result->err  = GSL_DBL_EPSILON * t1;
      result->err += result_c.err;
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else { /* x = 1 + v */
      const double v = x - 1.0;
      return cheb_eval_e(&psi_cs, 2.0*v-1.0, result);
    }
  }
}


/* psi(z) for large |z| in the right half-plane; [Abramowitz + Stegun, 6.3.18] */
static
gsl_complex
psi_complex_asymp(gsl_complex z)
{
  /* coefficients in the asymptotic expansion for large z;
   * let w = z^(-2) and write the expression in the form
   *
   *   ln(z) - 1/(2z) - 1/12 w (1 + c1 w + c2 w + c3 w + ... )
   */
  static const double c1 = -0.1;
  static const double c2 =  1.0/21.0;
  static const double c3 = -0.05;

  gsl_complex zi = gsl_complex_inverse(z);
  gsl_complex w  = gsl_complex_mul(zi, zi);
  gsl_complex cs;

  /* Horner method evaluation of term in parentheses */
  gsl_complex sum;
  sum = gsl_complex_mul_real(w, c3/c2);
  sum = gsl_complex_add_real(sum, 1.0);
  sum = gsl_complex_mul_real(sum, c2/c1);
  sum = gsl_complex_mul(sum, w);
  sum = gsl_complex_add_real(sum, 1.0);
  sum = gsl_complex_mul_real(sum, c1);
  sum = gsl_complex_mul(sum, w);
  sum = gsl_complex_add_real(sum, 1.0);

  /* correction added to log(z) */
  cs = gsl_complex_mul(sum, w);
  cs = gsl_complex_mul_real(cs, -1.0/12.0);
  cs = gsl_complex_add(cs, gsl_complex_mul_real(zi, -0.5));

  return gsl_complex_add(gsl_complex_log(z), cs);
}



/* psi(z) for complex z in the right half-plane */
static int
psi_complex_rhp(
  gsl_complex z,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  )
{
  int n_recurse = 0;
  int i;
  gsl_complex a;

  if(GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0)
  {
    result_re->val = 0.0;
    result_im->val = 0.0;
    result_re->err = 0.0;
    result_im->err = 0.0;
    return GSL_EDOM;
  }

  /* compute the number of recurrences to apply */
  if(GSL_REAL(z) < 20.0 && fabs(GSL_IMAG(z)) < 20.0)
  {
    const double sp = sqrt(20.0 + GSL_IMAG(z));
    const double sn = sqrt(20.0 - GSL_IMAG(z));
    const double rhs = sp*sn - GSL_REAL(z);
    if(rhs > 0.0) n_recurse = ceil(rhs);
  }

  /* compute asymptotic at the large value z + n_recurse */
  a = psi_complex_asymp(gsl_complex_add_real(z, n_recurse));

  result_re->err = 2.0 * GSL_DBL_EPSILON * fabs(GSL_REAL(a));
  result_im->err = 2.0 * GSL_DBL_EPSILON * fabs(GSL_IMAG(a));

  /* descend recursively, if necessary */
  for(i = n_recurse; i >= 1; --i)
  {
    gsl_complex zn = gsl_complex_add_real(z, i - 1.0);
    gsl_complex zn_inverse = gsl_complex_inverse(zn);
    a = gsl_complex_sub(a, zn_inverse);

    /* accumulate the error, to catch cancellations */
    result_re->err += 2.0 * GSL_DBL_EPSILON * fabs(GSL_REAL(zn_inverse));
    result_im->err += 2.0 * GSL_DBL_EPSILON * fabs(GSL_IMAG(zn_inverse));
  }

  result_re->val = GSL_REAL(a);
  result_im->val = GSL_IMAG(a);

  result_re->err += 2.0 * GSL_DBL_EPSILON * fabs(result_re->val);
  result_im->err += 2.0 * GSL_DBL_EPSILON * fabs(result_im->val);

  return GSL_SUCCESS;
}



/* generic polygamma; assumes n >= 0 and x > 0
 */
static int
psi_n_xg0(const int n, const double x, gsl_sf_result * result)
{
  if(n == 0) {
    return gsl_sf_psi_e(x, result);
  }
  else {
    /* Abramowitz + Stegun 6.4.10 */
    gsl_sf_result ln_nf;
    gsl_sf_result hzeta;
    int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
    int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
    int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                           hzeta.val, hzeta.err,
                                           result);
    if(GSL_IS_EVEN(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
  }
}



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_psi_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n <= 0) {
    DOMAIN_ERROR(result);
  }
  else if(n <= PSI_TABLE_NMAX) {
    result->val = psi_table[n];
    result->err = GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.3.18 */
    const double c2 = -1.0/12.0;
    const double c3 =  1.0/120.0;
    const double c4 = -1.0/252.0;
    const double c5 =  1.0/240.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2*c5)));
    result->val  = log(n) - 0.5/n + ser;
    result->err  = GSL_DBL_EPSILON * (fabs(log(n)) + fabs(0.5/n) + fabs(ser));
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  return psi_x(x, result);
}


int
gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result)
{
  const double ay = fabs(y);

  /* CHECK_POINTER(result) */

  if(ay > 1000.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    const double yi2 = 1.0/(ay*ay);
    const double lny = log(ay);
    const double sum = yi2 * (1.0/12.0 + 1.0/120.0 * yi2 + 1.0/252.0 * yi2*yi2);
    result->val = lny + sum;
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lny) + fabs(sum));
    return GSL_SUCCESS;
  }
  else if(ay > 10.0) {
    /* [Abramowitz+Stegun, 6.3.19] */
    const double yi2 = 1.0/(ay*ay);
    const double lny = log(ay);
    const double sum = yi2 * (1.0/12.0 +
                         yi2 * (1.0/120.0 +
                           yi2 * (1.0/252.0 +
                             yi2 * (1.0/240.0 +
                               yi2 * (1.0/132.0 + 691.0/32760.0 * yi2)))));
    result->val = lny + sum;
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lny) + fabs(sum));
    return GSL_SUCCESS;
  }
  else if(ay > 1.0){
    const double y2 = ay*ay;
    const double x  = (2.0*ay - 11.0)/9.0;
    const double v  = y2*(1.0/(1.0+y2) + 0.5/(4.0+y2));
    gsl_sf_result result_c;
    cheb_eval_e(&r1py_cs, x, &result_c);
    result->val  = result_c.val - M_EULER + v;
    result->err  = result_c.err;
    result->err += 2.0 * GSL_DBL_EPSILON * (fabs(v) + M_EULER + fabs(result_c.val));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    result->err *= 5.0; /* FIXME: losing a digit somewhere... maybe at x=... ? */
    return GSL_SUCCESS;
  }
  else {
    /* [Abramowitz+Stegun, 6.3.17]
     *
     * -M_EULER + y^2 Sum[1/n 1/(n^2 + y^2), {n,1,M}]
     *   +     Sum[1/n^3, {n,M+1,Infinity}]
     *   - y^2 Sum[1/n^5, {n,M+1,Infinity}]
     *   + y^4 Sum[1/n^7, {n,M+1,Infinity}]
     *   - y^6 Sum[1/n^9, {n,M+1,Infinity}]
     *   + O(y^8)
     *
     * We take M=50 for at least 15 digit precision.
     */
    const int M = 50;
    const double y2 = y*y;
    const double c0 = 0.00019603999466879846570;
    const double c2 = 3.8426659205114376860e-08;
    const double c4 = 1.0041592839497643554e-11;
    const double c6 = 2.9516743763500191289e-15;
    const double p  = c0 + y2 *(-c2 + y2*(c4 - y2*c6));
    double sum = 0.0;
    double v;

    int n;
    for(n=1; n<=M; n++) {
      sum += 1.0/(n * (n*n + y*y));
    }

    v = y2 * (sum + p);
    result->val  = -M_EULER + v;
    result->err  = GSL_DBL_EPSILON * (M_EULER + fabs(v));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  if(n <= 0) {
    DOMAIN_ERROR(result);
  }
  else if(n <= PSI_1_TABLE_NMAX) {
    result->val = psi_1_table[n];
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
  else {
    /* Abramowitz+Stegun 6.4.12
     * double-precision for n > 100
     */
    const double c0 = -1.0/30.0;
    const double c1 =  1.0/42.0;
    const double c2 = -1.0/30.0;
    const double ni2 = (1.0/n)*(1.0/n);
    const double ser =  ni2*ni2 * (c0 + ni2*(c1 + c2*ni2));
    result->val = (1.0 + 0.5/n + 1.0/(6.0*n*n) + ser) / n;
    result->err = GSL_DBL_EPSILON * result->val;
    return GSL_SUCCESS;
  }
}


int gsl_sf_psi_1_e(const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x == 0.0 || x == -1.0 || x == -2.0) {
    DOMAIN_ERROR(result);
  }
  else if(x > 0.0)
  {
    return psi_n_xg0(1, x, result);
  }
  else if(x > -5.0)
  {
    /* Abramowitz + Stegun 6.4.6 */
    int M = -floor(x);
    double fx = x + M;
    double sum = 0.0;
    int m;

    if(fx == 0.0)
      DOMAIN_ERROR(result);

    for(m = 0; m < M; ++m)
      sum += 1.0/((x+m)*(x+m));

    {
      int stat_psi = psi_n_xg0(1, fx, result);
      result->val += sum;
      result->err += M * GSL_DBL_EPSILON * sum;
      return stat_psi;
    }
  }
  else
  {
    /* Abramowitz + Stegun 6.4.7 */
    const double sin_px = sin(M_PI * x);
    const double d = M_PI*M_PI/(sin_px*sin_px);
    gsl_sf_result r;
    int stat_psi = psi_n_xg0(1, 1.0-x, &r);
    result->val = d - r.val;
    result->err = r.err + 2.0*GSL_DBL_EPSILON*d;
    return stat_psi;
  }
}


int gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(n == 0)
  {
    return gsl_sf_psi_e(x, result);
  }
  else if(n == 1)
  {
    return gsl_sf_psi_1_e(x, result);
  }
  else if(n < 0 || x <= 0.0) {
    DOMAIN_ERROR(result);
  }
  else {
    gsl_sf_result ln_nf;
    gsl_sf_result hzeta;
    int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
    int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
    int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
                                           hzeta.val, hzeta.err,
                                           result);
    if(GSL_IS_EVEN(n)) result->val = -result->val;
    return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
  }
}


int
gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  )
{
  if(x >= 0.0)
  {
    gsl_complex z = gsl_complex_rect(x, y);
    return psi_complex_rhp(z, result_re, result_im);
  }
  else
  {
    /* reflection formula [Abramowitz+Stegun, 6.3.7] */
    gsl_complex z = gsl_complex_rect(x, y);
    gsl_complex omz = gsl_complex_rect(1.0 - x, -y);
    gsl_complex zpi = gsl_complex_mul_real(z, M_PI);
    gsl_complex cotzpi = gsl_complex_cot(zpi);
    int ret_val = psi_complex_rhp(omz, result_re, result_im);

    if(GSL_IS_REAL(GSL_REAL(cotzpi)) && GSL_IS_REAL(GSL_IMAG(cotzpi)))
    {
      result_re->val -= M_PI * GSL_REAL(cotzpi);
      result_im->val -= M_PI * GSL_IMAG(cotzpi);
      return ret_val;
    }
    else
    {
      GSL_ERROR("singularity", GSL_EDOM);
    }
  }
}



/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_psi_int(const int n)
{
  EVAL_RESULT(gsl_sf_psi_int_e(n, &result));
}

double gsl_sf_psi(const double x)
{
  EVAL_RESULT(gsl_sf_psi_e(x, &result));
}

double gsl_sf_psi_1piy(const double x)
{
  EVAL_RESULT(gsl_sf_psi_1piy_e(x, &result));
}

double gsl_sf_psi_1_int(const int n)
{
  EVAL_RESULT(gsl_sf_psi_1_int_e(n, &result));
}

double gsl_sf_psi_1(const double x)
{
  EVAL_RESULT(gsl_sf_psi_1_e(x, &result));
}

double gsl_sf_psi_n(const int n, const double x)
{
  EVAL_RESULT(gsl_sf_psi_n_e(n, x, &result));
}
/* end inlined source: specfunc/psi.c */
/* begin inlined source: specfunc/airy.c */
/* specfunc/airy.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_trig.h */
/* specfunc/gsl_sf_trig.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_TRIG_H__
#define __GSL_SF_TRIG_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Sin(x) with GSL semantics. This is actually important
 * because we want to control the error estimate, and trying
 * to guess the error for the standard library implementation
 * every time it is used would be a little goofy.
 */
int gsl_sf_sin_e(double x, gsl_sf_result * result);
double gsl_sf_sin(const double x);


/* Cos(x) with GSL semantics.
 */
int gsl_sf_cos_e(double x, gsl_sf_result * result);
double gsl_sf_cos(const double x);


/* Hypot(x,y) with GSL semantics.
 */
int gsl_sf_hypot_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_hypot(const double x, const double y);


/* Sin(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_sin_e(const double zr, const double zi, gsl_sf_result * szr, gsl_sf_result * szi);


/* Cos(z) for complex z
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_complex_cos_e(const double zr, const double zi, gsl_sf_result * czr, gsl_sf_result * czi);


/* Log(Sin(z)) for complex z
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_complex_logsin_e(const double zr, const double zi, gsl_sf_result * lszr, gsl_sf_result * lszi);


/* Sinc(x) = sin(pi x) / (pi x)
 *
 * exceptions: none
 */
int gsl_sf_sinc_e(double x, gsl_sf_result * result);
double gsl_sf_sinc(const double x);


/* Log(Sinh(x)), x > 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnsinh_e(const double x, gsl_sf_result * result);
double gsl_sf_lnsinh(const double x);


/* Log(Cosh(x))
 *
 * exceptions: none
 */
int gsl_sf_lncosh_e(const double x, gsl_sf_result * result);
double gsl_sf_lncosh(const double x);


/* Convert polar to rectlinear coordinates.
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_polar_to_rect(const double r, const double theta, gsl_sf_result * x, gsl_sf_result * y);

/* Convert rectilinear to polar coordinates.
 * return argument in range [-pi, pi]
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_rect_to_polar(const double x, const double y, gsl_sf_result * r, gsl_sf_result * theta);

/* Sin(x) for quantity with an associated error.
 */
int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result);


/* Cos(x) for quantity with an associated error.
 */
int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result);


/* Force an angle to lie in the range (-pi,pi].
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_symm_e(double * theta);
double gsl_sf_angle_restrict_symm(const double theta);


/* Force an angle to lie in the range [0, 2pi)
 *
 * exceptions: GSL_ELOSS
 */
int gsl_sf_angle_restrict_pos_e(double * theta);
double gsl_sf_angle_restrict_pos(const double theta);


int gsl_sf_angle_restrict_symm_err_e(const double theta, gsl_sf_result * result);

int gsl_sf_angle_restrict_pos_err_e(const double theta, gsl_sf_result * result);


__END_DECLS

#endif /* __GSL_SF_TRIG_H__ */

/* end inlined header: gsl/gsl_sf_trig.h */
/* begin inlined header: gsl/gsl_sf_airy.h */
/* specfunc/gsl_sf_airy.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_AIRY_H__
#define __GSL_SF_AIRY_H__
/* begin inlined header: gsl/gsl_mode.h */
/* gsl_mode.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_MODE_H__
#define __GSL_MODE_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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


/* Some functions can take a mode argument. This
 * is a rough method to do things like control
 * the precision of the algorithm. This mainly
 * occurs in special functions, but we figured
 * it was ok to have a general facility.
 *
 * The mode type is 32-bit field. Most of
 * the fields are currently unused. Users
 * '|' various predefined constants to get
 * a desired mode.
 */
typedef unsigned int gsl_mode_t;


/* Here are the predefined constants.
 * Note that the precision constants
 * are special because they are used
 * to index arrays, so do not change
 * them. The precision information is
 * in the low order 3 bits of gsl_mode_t
 * (the third bit is currently unused).
 */

/* Note that "0" is double precision,
 * so that you get that by default if
 * you forget a flag.
 */
#define GSL_PREC_DOUBLE  0
#define GSL_PREC_SINGLE  1
#define GSL_PREC_APPROX  2

#ifdef HAVE_INLINE
INLINE_FUN unsigned int GSL_MODE_PREC(gsl_mode_t mt);

INLINE_FUN unsigned int
GSL_MODE_PREC(gsl_mode_t mt)
{ return  (mt & (unsigned int)7); }
#else  /* HAVE_INLINE */
#define GSL_MODE_PREC(mt) ((mt) & (unsigned int)7)
#endif /* HAVE_INLINE */


/* Here are some predefined generic modes.
 */
#define GSL_MODE_DEFAULT  0


__END_DECLS

#endif /* __GSL_MODE_H__ */

/* end inlined header: gsl/gsl_mode.h */
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Airy function Ai(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_e(const double x, const gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai(const double x, gsl_mode_t mode);


/* Airy function Bi(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi(const double x, gsl_mode_t mode);


/* scaled Ai(x):
 *                     Ai(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai_scaled(const double x, gsl_mode_t mode);


/* scaled Bi(x):
 *                     Bi(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi_scaled(const double x, gsl_mode_t mode);


/* derivative Ai'(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai_deriv(const double x, gsl_mode_t mode);


/* derivative Bi'(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi_deriv(const double x, gsl_mode_t mode);


/* scaled derivative Ai'(x):
 *                     Ai'(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai_deriv_scaled(const double x, gsl_mode_t mode);


/* scaled derivative:
 *                     Bi'(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi_deriv_scaled(const double x, gsl_mode_t mode);


/* Zeros of Ai(x)
 */
int gsl_sf_airy_zero_Ai_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Ai(unsigned int s);


/* Zeros of Bi(x)
 */
int gsl_sf_airy_zero_Bi_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Bi(unsigned int s);


/* Zeros of Ai'(x)
 */
int gsl_sf_airy_zero_Ai_deriv_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Ai_deriv(unsigned int s);


/* Zeros of Bi'(x)
 */
int gsl_sf_airy_zero_Bi_deriv_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Bi_deriv(unsigned int s);


__END_DECLS

#endif /* __GSL_SF_AIRY_H__ */

/* end inlined header: gsl/gsl_sf_airy.h */
/* inlined: error.h */
/* inlined: check.h */

/* inlined: chebyshev.h */
/* begin inlined source: specfunc/cheb_eval_mode.c */
#ifndef MARLEY_BUILTIN_GSL_CHEB_EVAL_MODE_C
#define MARLEY_BUILTIN_GSL_CHEB_EVAL_MODE_C

static inline int
cheb_eval_mode_e(const cheb_series * cs,
                 const double x,
                 gsl_mode_t mode,
                 gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  int eval_order;

  if(GSL_MODE_PREC(mode) == GSL_PREC_DOUBLE)
    eval_order = cs->order;
  else
    eval_order = cs->order_sp;

  for(j = eval_order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  result->val = y*d - dd + 0.5 * cs->c[0];
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->c[eval_order]);
  return GSL_SUCCESS;
}

#endif /* MARLEY_BUILTIN_GSL_CHEB_EVAL_MODE_C */
/* end inlined source: specfunc/cheb_eval_mode.c */

/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/


/* chebyshev expansions for Airy modulus and phase
   based on SLATEC r9aimp()

 Series for AM21       on the interval -1.25000D-01 to  0.
                                        with weighted error   2.89E-17
                                         log weighted error  16.54
                               significant figures required  14.15
                                    decimal places required  17.34

 Series for ATH1       on the interval -1.25000D-01 to  0.
                                        with weighted error   2.53E-17
                                         log weighted error  16.60
                               significant figures required  15.15
                                    decimal places required  17.38

 Series for AM22       on the interval -1.00000D+00 to -1.25000D-01
                                        with weighted error   2.99E-17
                                         log weighted error  16.52
                               significant figures required  14.57
                                    decimal places required  17.28

 Series for ATH2       on the interval -1.00000D+00 to -1.25000D-01
                                        with weighted error   2.57E-17
                                         log weighted error  16.59
                               significant figures required  15.07
                                    decimal places required  17.34
*/

static double am21_data[37] = {
  0.0065809191761485,
  0.0023675984685722,
  0.0001324741670371,
  0.0000157600904043,
  0.0000027529702663,
  0.0000006102679017,
  0.0000001595088468,
  0.0000000471033947,
  0.0000000152933871,
  0.0000000053590722,
  0.0000000020000910,
  0.0000000007872292,
  0.0000000003243103,
  0.0000000001390106,
  0.0000000000617011,
  0.0000000000282491,
  0.0000000000132979,
  0.0000000000064188,
  0.0000000000031697,
  0.0000000000015981,
  0.0000000000008213,
  0.0000000000004296,
  0.0000000000002284,
  0.0000000000001232,
  0.0000000000000675,
  0.0000000000000374,
  0.0000000000000210,
  0.0000000000000119,
  0.0000000000000068,
  0.0000000000000039,
  0.0000000000000023,
  0.0000000000000013,
  0.0000000000000008,
  0.0000000000000005,
  0.0000000000000003,
  0.0000000000000001,
  0.0000000000000001
};
static cheb_series am21_cs = {
  am21_data,
  36,
  -1, 1,
  20
};

static double ath1_data[36] = {
  -0.07125837815669365,
  -0.00590471979831451,
  -0.00012114544069499,
  -0.00000988608542270,
  -0.00000138084097352,
  -0.00000026142640172,
  -0.00000006050432589,
  -0.00000001618436223,
  -0.00000000483464911,
  -0.00000000157655272,
  -0.00000000055231518,
  -0.00000000020545441,
  -0.00000000008043412,
  -0.00000000003291252,
  -0.00000000001399875,
  -0.00000000000616151,
  -0.00000000000279614,
  -0.00000000000130428,
  -0.00000000000062373,
  -0.00000000000030512,
  -0.00000000000015239,
  -0.00000000000007758,
  -0.00000000000004020,
  -0.00000000000002117,
  -0.00000000000001132,
  -0.00000000000000614,
  -0.00000000000000337,
  -0.00000000000000188,
  -0.00000000000000105,
  -0.00000000000000060,
  -0.00000000000000034,
  -0.00000000000000020,
  -0.00000000000000011,
  -0.00000000000000007,
  -0.00000000000000004,
  -0.00000000000000002
};
static cheb_series ath1_cs = {
  ath1_data,
  35,
  -1, 1,
  15
};

static double am22_data[33] = {
 -0.01562844480625341,
  0.00778336445239681,
  0.00086705777047718,
  0.00015696627315611,
  0.00003563962571432,
  0.00000924598335425,
  0.00000262110161850,
  0.00000079188221651,
  0.00000025104152792,
  0.00000008265223206,
  0.00000002805711662,
  0.00000000976821090,
  0.00000000347407923,
  0.00000000125828132,
  0.00000000046298826,
  0.00000000017272825,
  0.00000000006523192,
  0.00000000002490471,
  0.00000000000960156,
  0.00000000000373448,
  0.00000000000146417,
  0.00000000000057826,
  0.00000000000022991,
  0.00000000000009197,
  0.00000000000003700,
  0.00000000000001496,
  0.00000000000000608,
  0.00000000000000248,
  0.00000000000000101,
  0.00000000000000041,
  0.00000000000000017,
  0.00000000000000007,
  0.00000000000000002
};
static cheb_series am22_cs = {
  am22_data,
  32,
  -1, 1,
  15
};

static double ath2_data[32] = {
   0.00440527345871877,
  -0.03042919452318455,
  -0.00138565328377179,
  -0.00018044439089549,
  -0.00003380847108327,
  -0.00000767818353522,
  -0.00000196783944371,
  -0.00000054837271158,
  -0.00000016254615505,
  -0.00000005053049981,
  -0.00000001631580701,
  -0.00000000543420411,
  -0.00000000185739855,
  -0.00000000064895120,
  -0.00000000023105948,
  -0.00000000008363282,
  -0.00000000003071196,
  -0.00000000001142367,
  -0.00000000000429811,
  -0.00000000000163389,
  -0.00000000000062693,
  -0.00000000000024260,
  -0.00000000000009461,
  -0.00000000000003716,
  -0.00000000000001469,
  -0.00000000000000584,
  -0.00000000000000233,
  -0.00000000000000093,
  -0.00000000000000037,
  -0.00000000000000015,
  -0.00000000000000006,
  -0.00000000000000002
};
static cheb_series ath2_cs = {
  ath2_data,
  31,
  -1, 1,
  16
};


/* Airy modulus and phase for x < -1 */
static
int
airy_mod_phase(const double x, gsl_mode_t mode, gsl_sf_result * mod, gsl_sf_result * phase)
{
  gsl_sf_result result_m;
  gsl_sf_result result_p;
  double m, p;
  double sqx;

  if(x < -2.0) {
    double z = 16.0/(x*x*x) + 1.0;
    cheb_eval_mode_e(&am21_cs, z, mode, &result_m);
    cheb_eval_mode_e(&ath1_cs, z, mode, &result_p);
  }
  else if(x <= -1.0) {
    double z = (16.0/(x*x*x) + 9.0)/7.0;
    cheb_eval_mode_e(&am22_cs, z, mode, &result_m);
    cheb_eval_mode_e(&ath2_cs, z, mode, &result_p);
  }
  else {
    mod->val = 0.0;
    mod->err = 0.0;
    phase->val = 0.0;
    phase->err = 0.0;
    GSL_ERROR ("x is greater than 1.0", GSL_EDOM);
  }

  m =  0.3125 + result_m.val;
  p = -0.625  + result_p.val;

  sqx = sqrt(-x);

  mod->val   = sqrt(m/sqx);
  mod->err  = fabs(mod->val) * (GSL_DBL_EPSILON + fabs(result_m.err/result_m.val));
  phase->val = M_PI_4 - x*sqx * p;
  phase->err = fabs(phase->val) * (GSL_DBL_EPSILON + fabs(result_p.err/result_p.val));

  return GSL_SUCCESS;
}



/* Chebyshev series for Ai(x) with x in [-1,1]
   based on SLATEC ai(x)

 series for aif        on the interval -1.00000d+00 to  1.00000d+00
                                   with weighted error   1.09e-19
                                    log weighted error  18.96
                          significant figures required  17.76
                               decimal places required  19.44

 series for aig        on the interval -1.00000d+00 to  1.00000d+00
                                   with weighted error   1.51e-17
                                    log weighted error  16.82
                          significant figures required  15.19
                               decimal places required  17.27
 */
static double ai_data_f[9] = {
  -0.03797135849666999750,
   0.05919188853726363857,
   0.00098629280577279975,
   0.00000684884381907656,
   0.00000002594202596219,
   0.00000000006176612774,
   0.00000000000010092454,
   0.00000000000000012014,
   0.00000000000000000010
};
static cheb_series aif_cs = {
  ai_data_f,
  8,
  -1, 1,
  8
};

static double ai_data_g[8] = {
   0.01815236558116127,
   0.02157256316601076,
   0.00025678356987483,
   0.00000142652141197,
   0.00000000457211492,
   0.00000000000952517,
   0.00000000000001392,
   0.00000000000000001
};
static cheb_series aig_cs = {
  ai_data_g,
  7,
  -1, 1,
  7
};


/* Chebvyshev series for Bi(x) with x in [-1,1]
   based on SLATEC bi(x)

 series for bif        on the interval -1.00000d+00 to  1.00000d+00
                                        with weighted error   1.88e-19
                                         log weighted error  18.72
                               significant figures required  17.74
                                    decimal places required  19.20

 series for big        on the interval -1.00000d+00 to  1.00000d+00
                                        with weighted error   2.61e-17
                                         log weighted error  16.58
                               significant figures required  15.17
                                    decimal places required  17.03
 */
static double data_bif[9] = {
  -0.01673021647198664948,
   0.10252335834249445610,
   0.00170830925073815165,
   0.00001186254546774468,
   0.00000004493290701779,
   0.00000000010698207143,
   0.00000000000017480643,
   0.00000000000000020810,
   0.00000000000000000018
};
static cheb_series bif_cs = {
  data_bif,
  8,
  -1, 1,
  8
};

static double data_big[8] = {
   0.02246622324857452,
   0.03736477545301955,
   0.00044476218957212,
   0.00000247080756363,
   0.00000000791913533,
   0.00000000001649807,
   0.00000000000002411,
   0.00000000000000002
};
static cheb_series big_cs = {
  data_big,
  7,
  -1, 1,
  7
};


/* Chebyshev series for Bi(x) with x in [1,8]
   based on SLATEC bi(x)
 */
static double data_bif2[10] = {
  0.0998457269381604100,
  0.4786249778630055380,
  0.0251552119604330118,
  0.0005820693885232645,
  0.0000074997659644377,
  0.0000000613460287034,
  0.0000000003462753885,
  0.0000000000014288910,
  0.0000000000000044962,
  0.0000000000000000111
};
static cheb_series bif2_cs = {
  data_bif2,
  9,
  -1, 1,
  9
};

static double data_big2[10] = {
  0.033305662145514340,
  0.161309215123197068,
  0.0063190073096134286,
  0.0001187904568162517,
  0.0000013045345886200,
  0.0000000093741259955,
  0.0000000000474580188,
  0.0000000000001783107,
  0.0000000000000005167,
  0.0000000000000000011
};
static cheb_series big2_cs = {
  data_big2,
  9,
  -1, 1,
  9
};


/* chebyshev for Ai(x) asymptotic factor
   based on SLATEC aie()

 Series for AIP        on the interval  0.          to  1.00000D+00
                   with weighted error   5.10E-17
                    log weighted error  16.29
          significant figures required  14.41
               decimal places required  17.06

 [GJ] Sun Apr 19 18:14:31 EDT 1998
 There was something wrong with these coefficients. I was getting
 errors after 3 or 4 digits. So I recomputed this table. Now I get
 double precision agreement with Mathematica. But it does not seem
 possible that the small differences here would account for the
 original discrepancy. There must have been something wrong with my
 original usage...
*/
static double data_aip[36] = {
 -0.0187519297793867540198,
 -0.0091443848250055004725,
  0.0009010457337825074652,
 -0.0001394184127221491507,
  0.0000273815815785209370,
 -0.0000062750421119959424,
  0.0000016064844184831521,
 -0.0000004476392158510354,
  0.0000001334635874651668,
 -0.0000000420735334263215,
  0.0000000139021990246364,
 -0.0000000047831848068048,
  0.0000000017047897907465,
 -0.0000000006268389576018,
  0.0000000002369824276612,
 -0.0000000000918641139267,
  0.0000000000364278543037,
 -0.0000000000147475551725,
  0.0000000000060851006556,
 -0.0000000000025552772234,
  0.0000000000010906187250,
 -0.0000000000004725870319,
  0.0000000000002076969064,
 -0.0000000000000924976214,
  0.0000000000000417096723,
 -0.0000000000000190299093,
  0.0000000000000087790676,
 -0.0000000000000040927557,
  0.0000000000000019271068,
 -0.0000000000000009160199,
  0.0000000000000004393567,
 -0.0000000000000002125503,
  0.0000000000000001036735,
 -0.0000000000000000509642,
  0.0000000000000000252377,
 -0.0000000000000000125793
/*
  -.0187519297793868
  -.0091443848250055,
   .0009010457337825,
  -.0001394184127221,
   .0000273815815785,
  -.0000062750421119,
   .0000016064844184,
  -.0000004476392158,
   .0000001334635874,
  -.0000000420735334,
   .0000000139021990,
  -.0000000047831848,
   .0000000017047897,
  -.0000000006268389,
   .0000000002369824,
  -.0000000000918641,
   .0000000000364278,
  -.0000000000147475,
   .0000000000060851,
  -.0000000000025552,
   .0000000000010906,
  -.0000000000004725,
   .0000000000002076,
  -.0000000000000924,
   .0000000000000417,
  -.0000000000000190,
   .0000000000000087,
  -.0000000000000040,
   .0000000000000019,
  -.0000000000000009,
   .0000000000000004,
  -.0000000000000002,
   .0000000000000001,
  -.0000000000000000
*/
};
static cheb_series aip_cs = {
  data_aip,
  35,
  -1, 1,
  17
};


/* chebyshev for Bi(x) asymptotic factor
   based on SLATEC bie()

 Series for BIP        on the interval  1.25000D-01 to  3.53553D-01
                   with weighted error   1.91E-17
                    log weighted error  16.72
          significant figures required  15.35
               decimal places required  17.41

 Series for BIP2       on the interval  0.          to  1.25000D-01
                   with weighted error   1.05E-18
                    log weighted error  17.98
          significant figures required  16.74
               decimal places required  18.71
*/
static double data_bip[24] = {
  -0.08322047477943447,
   0.01146118927371174,
   0.00042896440718911,
  -0.00014906639379950,
  -0.00001307659726787,
   0.00000632759839610,
  -0.00000042226696982,
  -0.00000019147186298,
   0.00000006453106284,
  -0.00000000784485467,
  -0.00000000096077216,
   0.00000000070004713,
  -0.00000000017731789,
   0.00000000002272089,
   0.00000000000165404,
  -0.00000000000185171,
   0.00000000000059576,
  -0.00000000000012194,
   0.00000000000001334,
   0.00000000000000172,
  -0.00000000000000145,
   0.00000000000000049,
  -0.00000000000000011,
   0.00000000000000001
};
static cheb_series bip_cs = {
  data_bip,
  23,
  -1, 1,
  14
};

static double data_bip2[29] = {
  -0.113596737585988679,
   0.0041381473947881595,
   0.0001353470622119332,
   0.0000104273166530153,
   0.0000013474954767849,
   0.0000001696537405438,
  -0.0000000100965008656,
  -0.0000000167291194937,
  -0.0000000045815364485,
   0.0000000003736681366,
   0.0000000005766930320,
   0.0000000000621812650,
  -0.0000000000632941202,
  -0.0000000000149150479,
   0.0000000000078896213,
   0.0000000000024960513,
  -0.0000000000012130075,
  -0.0000000000003740493,
   0.0000000000002237727,
   0.0000000000000474902,
  -0.0000000000000452616,
  -0.0000000000000030172,
   0.0000000000000091058,
  -0.0000000000000009814,
  -0.0000000000000016429,
   0.0000000000000005533,
   0.0000000000000002175,
  -0.0000000000000001737,
  -0.0000000000000000010
};
static cheb_series bip2_cs = {
  data_bip2,
  28,
  -1, 1,
  10
};


/* assumes x >= 1.0 */
inline static int
airy_aie(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  double sqx = sqrt(x);
  double z   = 2.0/(x*sqx) - 1.0;
  double y   = sqrt(sqx);
  gsl_sf_result result_c;
  cheb_eval_mode_e(&aip_cs, z, mode, &result_c);
  result->val = (0.28125 + result_c.val)/y;
  result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}

/* assumes x >= 2.0 */
static int airy_bie(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  const double ATR =  8.7506905708484345;
  const double BTR = -2.0938363213560543;

  if(x < 4.0) {
    double sqx = sqrt(x);
    double z   = ATR/(x*sqx) + BTR;
    double y   = sqrt(sqx);
    gsl_sf_result result_c;
    cheb_eval_mode_e(&bip_cs, z, mode, &result_c);
    result->val = (0.625 + result_c.val)/y;
    result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    double sqx = sqrt(x);
    double z   = 16.0/(x*sqx) - 1.0;
    double y   = sqrt(sqx);
    gsl_sf_result result_c;
    cheb_eval_mode_e(&bip2_cs, z, mode, &result_c);
    result->val = (0.625 + result_c.val)/y;
    result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_airy_Ai_e(const double x, const gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result cos_result;
    int stat_mp  = airy_mod_phase(x, mode, &mod, &theta);
    int stat_cos = gsl_sf_cos_err_e(theta.val, theta.err, &cos_result);
    result->val  = mod.val * cos_result.val;
    result->err  = fabs(mod.val * cos_result.err) + fabs(cos_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_cos);
  }
  else if(x <= 1.0) {
    const double z = x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs, z, mode, &result_c1);
    result->val  = 0.375 + (result_c0.val - x*(0.25 + result_c1.val));
    result->err  = result_c0.err + fabs(x*result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double x32 = x * sqrt(x);
    double s   = exp(-2.0*x32/3.0);
    gsl_sf_result result_aie;
    int stat_aie = airy_aie(x, mode, &result_aie);
    result->val  = result_aie.val * s;
    result->err  = result_aie.err * s + result->val * x32 * GSL_DBL_EPSILON;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return stat_aie;
  }
}


int
gsl_sf_airy_Ai_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result cos_result;
    int stat_mp  = airy_mod_phase(x, mode, &mod, &theta);
    int stat_cos = gsl_sf_cos_err_e(theta.val, theta.err, &cos_result);
    result->val  = mod.val * cos_result.val;
    result->err  = fabs(mod.val * cos_result.err) + fabs(cos_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_cos);
  }
  else if(x <= 1.0) {
    const double z = x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs, z, mode, &result_c1);
    result->val  = 0.375 + (result_c0.val - x*(0.25 + result_c1.val));
    result->err  = result_c0.err + fabs(x*result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    if(x > 0.0) {
      const double scale = exp(2.0/3.0 * sqrt(z));
      result->val *= scale;
      result->err *= scale;
    }

    return GSL_SUCCESS;
  }
  else {
    return airy_aie(x, mode, result);
  }
}


int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result sin_result;
    int stat_mp  = airy_mod_phase(x, mode, &mod, &theta);
    int stat_sin = gsl_sf_sin_err_e(theta.val, theta.err, &sin_result);
    result->val  = mod.val * sin_result.val;
    result->err  = fabs(mod.val * sin_result.err) + fabs(sin_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_sin);
  }
  else if(x < 1.0) {
    const double z = x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big_cs, z, mode, &result_c1);
    result->val  = 0.625 + result_c0.val + x*(0.4375 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    const double z = (2.0*x*x*x - 9.0)/7.0;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c1);
    result->val  = 1.125 + result_c0.val + x*(0.625 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double y = 2.0*x*sqrt(x)/3.0;
    const double s = exp(y);

    if(y > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR(result);
    }
    else {
      gsl_sf_result result_bie;
      int stat_bie = airy_bie(x, mode, &result_bie);
      result->val  = result_bie.val * s;
      result->err  = result_bie.err * s + fabs(1.5*y * (GSL_DBL_EPSILON * result->val));
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return stat_bie;
    }
  }
}


int
gsl_sf_airy_Bi_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result sin_result;
    int stat_mp  = airy_mod_phase(x, mode, &mod, &theta);
    int stat_sin = gsl_sf_sin_err_e(theta.val, theta.err, &sin_result);
    result->val  = mod.val * sin_result.val;
    result->err  = fabs(mod.val * sin_result.err) + fabs(sin_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_sin);
  }
  else if(x < 1.0) {
    const double z = x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big_cs, z, mode, &result_c1);
    result->val  = 0.625 + result_c0.val + x*(0.4375 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    if(x > 0.0) {
      const double scale = exp(-2.0/3.0 * sqrt(z));
      result->val *= scale;
      result->err *= scale;
    }
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    const double x3 = x*x*x;
    const double z  = (2.0*x3 - 9.0)/7.0;
    const double s  = exp(-2.0/3.0 * sqrt(x3));
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c1);
    result->val  = s * (1.125 + result_c0.val + x*(0.625 + result_c1.val));
    result->err  = s * (result_c0.err + fabs(x * result_c1.err));
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    return airy_bie(x, mode, result);
  }
}


/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* inlined: eval.h */

double gsl_sf_airy_Ai(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_e(x, mode, &result));
}

double gsl_sf_airy_Ai_scaled(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_scaled_e(x, mode, &result));
}

double gsl_sf_airy_Bi(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_e(x, mode, &result));
}

double gsl_sf_airy_Bi_scaled(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_scaled_e(x, mode, &result));
}


/* end inlined source: specfunc/airy.c */
/* begin inlined source: specfunc/coulomb.c */
/* specfunc/coulomb.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

/* Evaluation of Coulomb wave functions F_L(eta, x), G_L(eta, x),
 * and their derivatives. A combination of Steed's method, asymptotic
 * results, and power series.
 *
 * Steed's method:
 *  [Barnett, CPC 21, 297 (1981)]
 * Power series and other methods:
 *  [Biedenharn et al., PR 97, 542 (1954)]
 *  [Bardin et al., CPC 3, 73 (1972)]
 *  [Abad+Sesma, CPC 71, 110 (1992)]
 */
/* begin inlined header: config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
#define MARLEY_BUILTIN_GSL_CONFIG_WRAPPER_H
/* begin inlined header: gsl_config.h */
#ifndef MARLEY_BUILTIN_GSL_CONFIG_H
#define MARLEY_BUILTIN_GSL_CONFIG_H

#define HAVE_INLINE 1
#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_DECL_ISINF 1
#define GSL_COERCE_DBL(x) (x)

#define GSL_DISABLE_DEPRECATED 0
#define PACKAGE_VERSION "2.8"
#define VERSION "2.8"

#endif

/* end inlined header: gsl_config.h */
#endif

/* end inlined header: config.h */
/* begin inlined header: gsl/gsl_math.h */
/* gsl_math.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MATH_H__
#define __GSL_MATH_H__
#include <math.h>
/* begin inlined header: gsl/gsl_sys.h */
/* sys/gsl_sys.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_SYS_H__
#define __GSL_SYS_H__

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

double gsl_log1p (const double x);
double gsl_expm1 (const double x);
double gsl_hypot (const double x, const double y);
double gsl_hypot3 (const double x, const double y, const double z);
double gsl_acosh (const double x);
double gsl_asinh (const double x);
double gsl_atanh (const double x);

int gsl_isnan (const double x);
int gsl_isinf (const double x);
int gsl_finite (const double x);

double gsl_nan (void);
double gsl_posinf (void);
double gsl_neginf (void);
double gsl_fdiv (const double x, const double y);

double gsl_coerce_double (const double x);
float gsl_coerce_float (const float x);
long double gsl_coerce_long_double (const long double x);

double gsl_ldexp(const double x, const int e);
double gsl_frexp(const double x, int * e);

int gsl_fcmp (const double x1, const double x2, const double epsilon);

__END_DECLS

#endif /* __GSL_SYS_H__ */

/* end inlined header: gsl/gsl_sys.h */
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
/* begin inlined header: gsl/gsl_machine.h */
/* Author:  B. Gough and G. Jungman */
#ifndef __GSL_MACHINE_H__
#define __GSL_MACHINE_H__

#include <limits.h>
#include <float.h>

/* magic constants; mostly for the benefit of the implementation */

/* -*-MACHINE CONSTANTS-*-
 *
 * PLATFORM: Whiz-O-Matic 9000
 * FP_PLATFORM: IEEE-Virtual
 * HOSTNAME: nnn.lanl.gov
 * DATE: Fri Nov 20 17:53:26 MST 1998
 */
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define GSL_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define GSL_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define GSL_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_SQRT_DBL_MIN   1.4916681462400413e-154
#define GSL_ROOT3_DBL_MIN  2.8126442852362996e-103
#define GSL_ROOT4_DBL_MIN  1.2213386697554620e-77
#define GSL_ROOT5_DBL_MIN  2.9476022969691763e-62
#define GSL_ROOT6_DBL_MIN  5.3034368905798218e-52
#define GSL_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define GSL_DBL_MAX        1.7976931348623157e+308
#define GSL_SQRT_DBL_MAX   1.3407807929942596e+154
#define GSL_ROOT3_DBL_MAX  5.6438030941222897e+102
#define GSL_ROOT4_DBL_MAX  1.1579208923731620e+77
#define GSL_ROOT5_DBL_MAX  4.4765466227572707e+61
#define GSL_ROOT6_DBL_MAX  2.3756689782295612e+51
#define GSL_LOG_DBL_MAX    7.0978271289338397e+02

#define GSL_FLT_EPSILON        1.1920928955078125e-07
#define GSL_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define GSL_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define GSL_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define GSL_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define GSL_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define GSL_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define GSL_FLT_MIN        1.1754943508222875e-38
#define GSL_SQRT_FLT_MIN   1.0842021724855044e-19
#define GSL_ROOT3_FLT_MIN  2.2737367544323241e-13
#define GSL_ROOT4_FLT_MIN  3.2927225399135965e-10
#define GSL_ROOT5_FLT_MIN  2.5944428542140822e-08
#define GSL_ROOT6_FLT_MIN  4.7683715820312542e-07
#define GSL_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define GSL_FLT_MAX        3.4028234663852886e+38
#define GSL_SQRT_FLT_MAX   1.8446743523953730e+19
#define GSL_ROOT3_FLT_MAX  6.9814635196223242e+12
#define GSL_ROOT4_FLT_MAX  4.2949672319999986e+09
#define GSL_ROOT5_FLT_MAX  5.0859007855960041e+07
#define GSL_ROOT6_FLT_MAX  2.6422459233807749e+06
#define GSL_LOG_FLT_MAX    8.8722839052068352e+01

#define GSL_SFLT_EPSILON        4.8828125000000000e-04
#define GSL_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define GSL_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define GSL_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define GSL_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define GSL_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define GSL_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)

/* !MACHINE CONSTANTS! */


/* a little internal backwards compatibility */
#define GSL_MACH_EPS  GSL_DBL_EPSILON



/* Here are the constants related to or derived from
 * machine constants. These are not to be confused with
 * the constants that define various precision levels
 * for the precision/error system.
 *
 * This information is determined at configure time
 * and is platform dependent. Edit at your own risk.
 *
 * PLATFORM: WHIZ-O-MATIC
 * CONFIG-DATE: Thu Nov 19 19:27:18 MST 1998
 * CONFIG-HOST: nnn.lanl.gov
 */

/* machine precision constants */
/* #define GSL_MACH_EPS         1.0e-15 */
#define GSL_SQRT_MACH_EPS       3.2e-08
#define GSL_ROOT3_MACH_EPS      1.0e-05
#define GSL_ROOT4_MACH_EPS      0.000178
#define GSL_ROOT5_MACH_EPS      0.00100
#define GSL_ROOT6_MACH_EPS      0.00316
#define GSL_LOG_MACH_EPS       (-34.54)


#endif /* __GSL_MACHINE_H__ */

/* end inlined header: gsl/gsl_machine.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
/* begin inlined header: gsl/gsl_nan.h */
/* gsl_nan.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_NAN_H__
#define __GSL_NAN_H__

#ifdef INFINITY
# define GSL_POSINF INFINITY
# define GSL_NEGINF (-INFINITY)
#elif defined(HUGE_VAL)
# define GSL_POSINF HUGE_VAL
# define GSL_NEGINF (-HUGE_VAL)
#else
# define GSL_POSINF (gsl_posinf())
# define GSL_NEGINF (gsl_neginf())
#endif

#ifdef NAN
# define GSL_NAN NAN
#elif defined(INFINITY)
# define GSL_NAN (INFINITY/INFINITY)
#else
# define GSL_NAN (gsl_nan())
#endif

#define GSL_POSZERO (+0.0)
#define GSL_NEGZERO (-0.0)

#endif /* __GSL_NAN_H__ */

/* end inlined header: gsl/gsl_nan.h */
/* begin inlined header: gsl/gsl_pow_int.h */
/* gsl_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

INLINE_DECL double gsl_pow_2(const double x);
INLINE_DECL double gsl_pow_3(const double x);
INLINE_DECL double gsl_pow_4(const double x);
INLINE_DECL double gsl_pow_5(const double x);
INLINE_DECL double gsl_pow_6(const double x);
INLINE_DECL double gsl_pow_7(const double x);
INLINE_DECL double gsl_pow_8(const double x);
INLINE_DECL double gsl_pow_9(const double x);

#ifdef HAVE_INLINE
INLINE_FUN double gsl_pow_2(const double x) { return x*x;   }
INLINE_FUN double gsl_pow_3(const double x) { return x*x*x; }
INLINE_FUN double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
INLINE_FUN double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
INLINE_FUN double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
INLINE_FUN double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
INLINE_FUN double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
INLINE_FUN double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#endif

double gsl_pow_int(double x, int n);
double gsl_pow_uint(double x, unsigned int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_pow_int.h */
/* begin inlined header: gsl/gsl_minmax.h */
/* gsl_minmax.h
 *
 * Copyright (C) 2008 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_MINMAX_H__
#define __GSL_MINMAX_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))
#define GSL_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double gsl_max (double a, double b);
double gsl_min (double a, double b);

/* inline-friendly strongly typed versions */
#ifdef HAVE_INLINE

INLINE_FUN int GSL_MAX_INT (int a, int b);
INLINE_FUN int GSL_MIN_INT (int a, int b);
INLINE_FUN double GSL_MAX_DBL (double a, double b);
INLINE_FUN double GSL_MIN_DBL (double a, double b);
INLINE_FUN long double GSL_MAX_LDBL (long double a, long double b);
INLINE_FUN long double GSL_MIN_LDBL (long double a, long double b);

INLINE_FUN int
GSL_MAX_INT (int a, int b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN int
GSL_MIN_INT (int a, int b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN double
GSL_MAX_DBL (double a, double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN double
GSL_MIN_DBL (double a, double b)
{
  return GSL_MIN (a, b);
}

INLINE_FUN long double
GSL_MAX_LDBL (long double a, long double b)
{
  return GSL_MAX (a, b);
}

INLINE_FUN long double
GSL_MIN_LDBL (long double a, long double b)
{
  return GSL_MIN (a, b);
}
#else
#define GSL_MAX_INT(a,b)   GSL_MAX(a,b)
#define GSL_MIN_INT(a,b)   GSL_MIN(a,b)
#define GSL_MAX_DBL(a,b)   GSL_MAX(a,b)
#define GSL_MIN_DBL(a,b)   GSL_MIN(a,b)
#define GSL_MAX_LDBL(a,b)  GSL_MAX(a,b)
#define GSL_MIN_LDBL(a,b)  GSL_MIN(a,b)
#endif /* HAVE_INLINE */

__END_DECLS

#endif /* __GSL_POW_INT_H__ */

/* end inlined header: gsl/gsl_minmax.h */
#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

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

/* other needlessly compulsive abstractions */

#define GSL_IS_ODD(n)  ((n) & 1)
#define GSL_IS_EVEN(n) (!(GSL_IS_ODD(n)))
#define GSL_SIGN(x)    ((x) >= 0.0 ? 1 : -1)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define GSL_IS_REAL(x) (gsl_finite(x))

/* Definition of an arbitrary function with parameters */

struct gsl_function_struct
{
  double (* function) (double x, void * params);
  void * params;
};

typedef struct gsl_function_struct gsl_function ;

#define GSL_FN_EVAL(F,x) (*((F)->function))(x,(F)->params)

/* Definition of an arbitrary function returning two values, r1, r2 */

struct gsl_function_fdf_struct
{
  double (* f) (double x, void * params);
  double (* df) (double x, void * params);
  void (* fdf) (double x, void * params, double * f, double * df);
  void * params;
};

typedef struct gsl_function_fdf_struct gsl_function_fdf ;

#define GSL_FN_FDF_EVAL_F(FDF,x) (*((FDF)->f))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_DF(FDF,x) (*((FDF)->df))(x,(FDF)->params)
#define GSL_FN_FDF_EVAL_F_DF(FDF,x,y,dy) (*((FDF)->fdf))(x,(FDF)->params,(y),(dy))


/* Definition of an arbitrary vector-valued function with parameters */

struct gsl_function_vec_struct
{
  int (* function) (double x, double y[], void * params);
  void * params;
};

typedef struct gsl_function_vec_struct gsl_function_vec ;

#define GSL_FN_VEC_EVAL(F,x,y) (*((F)->function))(x,y,(F)->params)

__END_DECLS

#endif /* __GSL_MATH_H__ */

/* end inlined header: gsl/gsl_math.h */
/* begin inlined header: gsl/gsl_errno.h */
/* err/gsl_errno.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_ERRNO_H__
#define __GSL_ERRNO_H__

#include <stdio.h>
#include <errno.h>
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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
  GSL_FAILURE  = -1,
  GSL_CONTINUE = -2,  /* iteration has not converged */
  GSL_EDOM     = 1,   /* input domain error, e.g sqrt(-1) */
  GSL_ERANGE   = 2,   /* output range error, e.g. exp(1e100) */
  GSL_EFAULT   = 3,   /* invalid pointer */
  GSL_EINVAL   = 4,   /* invalid argument supplied by user */
  GSL_EFAILED  = 5,   /* generic failure */
  GSL_EFACTOR  = 6,   /* factorization failed */
  GSL_ESANITY  = 7,   /* sanity check failed - shouldn't happen */
  GSL_ENOMEM   = 8,   /* malloc failed */
  GSL_EBADFUNC = 9,   /* problem with user-supplied function */
  GSL_ERUNAWAY = 10,  /* iterative process is out of control */
  GSL_EMAXITER = 11,  /* exceeded max number of iterations */
  GSL_EZERODIV = 12,  /* tried to divide by zero */
  GSL_EBADTOL  = 13,  /* user specified an invalid tolerance */
  GSL_ETOL     = 14,  /* failed to reach the specified tolerance */
  GSL_EUNDRFLW = 15,  /* underflow */
  GSL_EOVRFLW  = 16,  /* overflow  */
  GSL_ELOSS    = 17,  /* loss of accuracy */
  GSL_EROUND   = 18,  /* failed because of roundoff error */
  GSL_EBADLEN  = 19,  /* matrix, vector lengths are not conformant */
  GSL_ENOTSQR  = 20,  /* matrix not square */
  GSL_ESING    = 21,  /* apparent singularity detected */
  GSL_EDIVERGE = 22,  /* integral or series is divergent */
  GSL_EUNSUP   = 23,  /* requested feature is not supported by the hardware */
  GSL_EUNIMPL  = 24,  /* requested feature not (yet) implemented */
  GSL_ECACHE   = 25,  /* cache limit exceeded */
  GSL_ETABLE   = 26,  /* table limit exceeded */
  GSL_ENOPROG  = 27,  /* iteration is not making progress towards solution */
  GSL_ENOPROGJ = 28,  /* jacobian evaluations are not improving the solution */
  GSL_ETOLF    = 29,  /* cannot reach the specified tolerance in F */
  GSL_ETOLX    = 30,  /* cannot reach the specified tolerance in X */
  GSL_ETOLG    = 31,  /* cannot reach the specified tolerance in gradient */
  GSL_EOF      = 32   /* end of file */
} ;

void gsl_error (const char * reason, const char * file, int line,
                int gsl_errno);

void gsl_stream_printf (const char *label, const char *file,
                        int line, const char *reason);

typedef void gsl_error_handler_t (const char * reason, const char * file,
                                  int line, int gsl_errno);

gsl_error_handler_t *
gsl_set_error_handler (gsl_error_handler_t * new_handler);

/* GSL_ERROR: call the error handler, and return the error code */

#define GSL_ERROR(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return gsl_errno ; \
       } while (0)

/* GSL_ERROR_VAL: call the error handler, and return the given value */

#define GSL_ERROR_VAL(reason, gsl_errno, value) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return value ; \
       } while (0)

/* GSL_ERROR_VOID: call the error handler, and then return
   (for void functions which still need to generate an error) */

#define GSL_ERROR_VOID(reason, gsl_errno) \
       do { \
       gsl_error (reason, __FILE__, __LINE__, gsl_errno) ; \
       return ; \
       } while (0)

/* GSL_ERROR_NULL suitable for out-of-memory conditions */

#define GSL_ERROR_NULL(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, 0)

/* Sometimes you have several status results returned from
 * function calls and you want to combine them in some sensible
 * way. You cannot produce a "total" status condition, but you can
 * pick one from a set of conditions based on an implied hierarchy.
 *
 * In other words:
 *    you have: status_a, status_b, ...
 *    you want: status = (status_a if it is bad, or status_b if it is bad,...)
 *
 * In this example you consider status_a to be more important and
 * it is checked first, followed by the others in the order specified.
 *
 * Here are some dumb macros to do this.
 */
#define GSL_ERROR_SELECT_2(a,b)       ((a) != GSL_SUCCESS ? (a) : ((b) != GSL_SUCCESS ? (b) : GSL_SUCCESS))
#define GSL_ERROR_SELECT_3(a,b,c)     ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_2(b,c))
#define GSL_ERROR_SELECT_4(a,b,c,d)   ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_3(b,c,d))
#define GSL_ERROR_SELECT_5(a,b,c,d,e) ((a) != GSL_SUCCESS ? (a) : GSL_ERROR_SELECT_4(b,c,d,e))

#define GSL_STATUS_UPDATE(sp, s) do { if ((s) != GSL_SUCCESS) *(sp) = (s);} while(0)

__END_DECLS

#endif /* __GSL_ERRNO_H__ */

/* end inlined header: gsl/gsl_errno.h */
/* begin inlined header: gsl/gsl_sf_exp.h */
/* specfunc/gsl_sf_exp.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_EXP_H__
#define __GSL_SF_EXP_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
/* begin inlined header: gsl/gsl_precision.h */
/* gsl_precision.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_PRECISION_H__
#define __GSL_PRECISION_H__
/* begin inlined header: gsl/gsl_types.h */
/* gsl_types.h
 *
 * Copyright (C) 2001, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_TYPES_H__
#define __GSL_TYPES_H__

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef GSL_DLL
#    ifdef DLL_EXPORT
#      define GSL_VAR extern __declspec(dllexport)
#    else
#      define GSL_VAR extern __declspec(dllimport)
#    endif
#  else
#    define GSL_VAR extern
#  endif
#else
#  define GSL_VAR extern
#endif

#endif

#endif /* __GSL_TYPES_H__ */

/* end inlined header: gsl/gsl_types.h */
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


/* A type for the precision indicator.
 * This is mainly for pedagogy.
 */
typedef  unsigned int  gsl_prec_t;


/* The number of precision types.
 * Remember that precision-mode
 * can index an array.
 */
#define _GSL_PREC_T_NUM 3


/* Arrays containing derived
 * precision constants for the
 * different precision levels.
 */
GSL_VAR const double gsl_prec_eps[];
GSL_VAR const double gsl_prec_sqrt_eps[];
GSL_VAR const double gsl_prec_root3_eps[];
GSL_VAR const double gsl_prec_root4_eps[];
GSL_VAR const double gsl_prec_root5_eps[];
GSL_VAR const double gsl_prec_root6_eps[];


__END_DECLS

#endif /* __GSL_PRECISION_H__ */

/* end inlined header: gsl/gsl_precision.h */
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

/* Provide an exp() function with GSL semantics,
 * i.e. with proper error checking, etc.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e(const double x, gsl_sf_result * result);
double gsl_sf_exp(const double x);


/* Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_e10_e(const double x, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e(const double x, const double y, gsl_sf_result * result);
double gsl_sf_exp_mult(const double x, const double y);


/* Exponentiate and multiply by a given factor:  y * Exp(x)
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_e10_e(const double x, const double y, gsl_sf_result_e10 * result);


/* exp(x)-1
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_expm1_e(const double x, gsl_sf_result * result);
double gsl_sf_expm1(const double x);


/* (exp(x)-1)/x = 1 + x/2 + x^2/(2*3) + x^3/(2*3*4) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_e(const double x, gsl_sf_result * result);
double gsl_sf_exprel(const double x);


/* 2(exp(x)-1-x)/x^2 = 1 + x/3 + x^2/(3*4) + x^3/(3*4*5) + ...
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_exprel_2_e(double x, gsl_sf_result * result);
double gsl_sf_exprel_2(const double x);


/* Similarly for the N-th generalization of
 * the above. The so-called N-relative exponential
 *
 * exprel_N(x) = N!/x^N (exp(x) - Sum[x^k/k!, {k,0,N-1}])
 *             = 1 + x/(N+1) + x^2/((N+1)(N+2)) + ...
 *             = 1F1(1,1+N,x)
 */
int gsl_sf_exprel_n_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_exprel_n(const int n, const double x);

int gsl_sf_exprel_n_CF_e(const double n, const double x, gsl_sf_result * result);


/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e(const double x, const double dx, gsl_sf_result * result);

/* Exponentiate a quantity with an associated error.
 */
int gsl_sf_exp_err_e10_e(const double x, const double dx, gsl_sf_result_e10 * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e(const double x, const double dx, const double y, const double dy, gsl_sf_result * result);


/* Exponentiate and multiply by a given factor:  y * Exp(x),
 * for quantities with associated errors.
 *
 * exceptions: GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_exp_mult_err_e10_e(const double x, const double dx, const double y, const double dy, gsl_sf_result_e10 * result);

__END_DECLS

#endif /* __GSL_SF_EXP_H__ */

/* end inlined header: gsl/gsl_sf_exp.h */
/* begin inlined header: gsl/gsl_sf_psi.h */
/* specfunc/gsl_sf_psi.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_PSI_H__
#define __GSL_SF_PSI_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Poly-Gamma Functions
 *
 * psi(m,x) := (d/dx)^m psi(0,x) = (d/dx)^{m+1} log(gamma(x))
 */


/* Di-Gamma Function  psi(n) = psi(0,n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_int(const int n);


/* Di-Gamma Function psi(x) = psi(0, x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi(const double x);


/* Di-Gamma Function Re[psi(1 + I y)]
 *
 * exceptions: none
 */
int     gsl_sf_psi_1piy_e(const double y, gsl_sf_result * result);
double  gsl_sf_psi_1piy(const double y);


/* Di-Gamma Function psi(z) for general complex argument z = x + iy
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_complex_psi_e(
  const double x,
  const double y,
  gsl_sf_result * result_re,
  gsl_sf_result * result_im
  );


/* Tri-Gamma Function psi^(1)(n)
 *
 * n > 0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_1_int_e(const int n, gsl_sf_result * result);
double  gsl_sf_psi_1_int(const int n);


/* Tri-Gamma Function psi^(1)(x)
 *
 * x != 0.0, -1.0, -2.0, ...
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int     gsl_sf_psi_1_e(const double x, gsl_sf_result * result);
double  gsl_sf_psi_1(const double x);


/* Poly-Gamma Function psi^(n)(x)
 *
 * n >= 0, x > 0.0
 * exceptions: GSL_EDOM
 */
int     gsl_sf_psi_n_e(const int n, const double x, gsl_sf_result * result);
double  gsl_sf_psi_n(const int n, const double x);


__END_DECLS

#endif /* __GSL_SF_PSI_H__ */

/* end inlined header: gsl/gsl_sf_psi.h */
/* begin inlined header: gsl/gsl_sf_airy.h */
/* specfunc/gsl_sf_airy.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_AIRY_H__
#define __GSL_SF_AIRY_H__
/* begin inlined header: gsl/gsl_mode.h */
/* gsl_mode.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_MODE_H__
#define __GSL_MODE_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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


/* Some functions can take a mode argument. This
 * is a rough method to do things like control
 * the precision of the algorithm. This mainly
 * occurs in special functions, but we figured
 * it was ok to have a general facility.
 *
 * The mode type is 32-bit field. Most of
 * the fields are currently unused. Users
 * '|' various predefined constants to get
 * a desired mode.
 */
typedef unsigned int gsl_mode_t;


/* Here are the predefined constants.
 * Note that the precision constants
 * are special because they are used
 * to index arrays, so do not change
 * them. The precision information is
 * in the low order 3 bits of gsl_mode_t
 * (the third bit is currently unused).
 */

/* Note that "0" is double precision,
 * so that you get that by default if
 * you forget a flag.
 */
#define GSL_PREC_DOUBLE  0
#define GSL_PREC_SINGLE  1
#define GSL_PREC_APPROX  2

#ifdef HAVE_INLINE
INLINE_FUN unsigned int GSL_MODE_PREC(gsl_mode_t mt);

INLINE_FUN unsigned int
GSL_MODE_PREC(gsl_mode_t mt)
{ return  (mt & (unsigned int)7); }
#else  /* HAVE_INLINE */
#define GSL_MODE_PREC(mt) ((mt) & (unsigned int)7)
#endif /* HAVE_INLINE */


/* Here are some predefined generic modes.
 */
#define GSL_MODE_DEFAULT  0


__END_DECLS

#endif /* __GSL_MODE_H__ */

/* end inlined header: gsl/gsl_mode.h */
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Airy function Ai(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_e(const double x, const gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai(const double x, gsl_mode_t mode);


/* Airy function Bi(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi(const double x, gsl_mode_t mode);


/* scaled Ai(x):
 *                     Ai(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai_scaled(const double x, gsl_mode_t mode);


/* scaled Bi(x):
 *                     Bi(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi_scaled(const double x, gsl_mode_t mode);


/* derivative Ai'(x)
 *
 * exceptions: GSL_EUNDRFLW
 */
int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai_deriv(const double x, gsl_mode_t mode);


/* derivative Bi'(x)
 *
 * exceptions: GSL_EOVRFLW
 */
int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi_deriv(const double x, gsl_mode_t mode);


/* scaled derivative Ai'(x):
 *                     Ai'(x)   x < 0
 *   exp(+2/3 x^{3/2}) Ai'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Ai_deriv_scaled(const double x, gsl_mode_t mode);


/* scaled derivative:
 *                     Bi'(x)   x < 0
 *   exp(-2/3 x^{3/2}) Bi'(x)   x > 0
 *
 * exceptions: none
 */
int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result);
double gsl_sf_airy_Bi_deriv_scaled(const double x, gsl_mode_t mode);


/* Zeros of Ai(x)
 */
int gsl_sf_airy_zero_Ai_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Ai(unsigned int s);


/* Zeros of Bi(x)
 */
int gsl_sf_airy_zero_Bi_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Bi(unsigned int s);


/* Zeros of Ai'(x)
 */
int gsl_sf_airy_zero_Ai_deriv_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Ai_deriv(unsigned int s);


/* Zeros of Bi'(x)
 */
int gsl_sf_airy_zero_Bi_deriv_e(unsigned int s, gsl_sf_result * result);
double gsl_sf_airy_zero_Bi_deriv(unsigned int s);


__END_DECLS

#endif /* __GSL_SF_AIRY_H__ */

/* end inlined header: gsl/gsl_sf_airy.h */
/* begin inlined header: gsl/gsl_sf_pow_int.h */
/* specfunc/gsl_sf_pow_int.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_POW_INT_H__
#define __GSL_SF_POW_INT_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Calculate x^n.
 * Does not check for overflow/underflow.
 */
int     gsl_sf_pow_int_e(double x, int n, gsl_sf_result * result);
double  gsl_sf_pow_int(const double x, const int n);


__END_DECLS

#endif /* __GSL_SF_POW_INT_H__ */

/* end inlined header: gsl/gsl_sf_pow_int.h */
/* begin inlined header: gsl/gsl_sf_gamma.h */
/* specfunc/gsl_sf_gamma.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_GAMMA_H__
#define __GSL_SF_GAMMA_H__
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method.
 * Returns the real part of Log[Gamma[x]] when x < 0,
 * i.e. Log[|Gamma[x]|].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_e(double x, gsl_sf_result * result);
double gsl_sf_lngamma(const double x);


/* Log[Gamma(x)], x not a negative integer
 * Uses real Lanczos method. Determines
 * the sign of Gamma[x] as well as Log[|Gamma[x]|] for x < 0.
 * So Gamma[x] = sgn * Exp[result_lg].
 *
 * exceptions: GSL_EDOM, GSL_EROUND
 */
int gsl_sf_lngamma_sgn_e(double x, gsl_sf_result * result_lg, double *sgn);


/* Gamma(x), x not a negative integer
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EROUND
 */
int gsl_sf_gamma_e(const double x, gsl_sf_result * result);
double gsl_sf_gamma(const double x);


/* Regulated Gamma Function, x > 0
 * Gamma^*(x) = Gamma(x)/(Sqrt[2Pi] x^(x-1/2) exp(-x))
 *            = (1 + 1/(12x) + ...),  x->Inf
 * A useful suggestion of Temme.
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gammastar_e(const double x, gsl_sf_result * result);
double gsl_sf_gammastar(const double x);


/* 1/Gamma(x)
 * Uses real Lanczos method.
 *
 * exceptions: GSL_EUNDRFLW, GSL_EROUND
 */
int gsl_sf_gammainv_e(const double x, gsl_sf_result * result);
double gsl_sf_gammainv(const double x);


/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
 */
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);


/* x^n / n!
 *
 * x >= 0.0, n >= 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_taylorcoeff_e(const int n, const double x, gsl_sf_result * result);
double gsl_sf_taylorcoeff(const int n, const double x);


/* n!
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_fact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_fact(const unsigned int n);


/* n!! = n(n-2)(n-4) ...
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_doublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_doublefact(const unsigned int n);


/* log(n!)
 * Faster than ln(Gamma(n+1)) for n < 170; defers for larger n.
 *
 * exceptions: none
 */
int gsl_sf_lnfact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lnfact(const unsigned int n);


/* log(n!!)
 *
 * exceptions: none
 */
int gsl_sf_lndoublefact_e(const unsigned int n, gsl_sf_result * result);
double gsl_sf_lndoublefact(const unsigned int n);


/* log(n choose m)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnchoose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_lnchoose(unsigned int n, unsigned int m);


/* n choose m
 *
 * exceptions: GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_choose_e(unsigned int n, unsigned int m, gsl_sf_result * result);
double gsl_sf_choose(unsigned int n, unsigned int m);


/* Logarithm of Pochhammer (Apell) symbol
 *   log( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a > 0, a+x > 0
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_lnpoch(const double a, const double x);


/* Logarithm of Pochhammer (Apell) symbol, with sign information.
 *   result = log( |(a)_x| )
 *   sgn    = sgn( (a)_x )
 *   where (a)_x := Gamma[a + x]/Gamma[a]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_lnpoch_sgn_e(const double a, const double x, gsl_sf_result * result, double * sgn);


/* Pochhammer (Apell) symbol
 *   (a)_x := Gamma[a + x]/Gamma[x]
 *
 * a != neg integer, a+x != neg integer
 *
 * exceptions:  GSL_EDOM, GSL_EOVRFLW
 */
int gsl_sf_poch_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_poch(const double a, const double x);


/* Relative Pochhammer (Apell) symbol
 *   ((a,x) - 1)/x
 *   where (a,x) = (a)_x := Gamma[a + x]/Gamma[a]
 *
 * exceptions:  GSL_EDOM
 */
int gsl_sf_pochrel_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_pochrel(const double a, const double x);


/* Normalized Incomplete Gamma Function
 *
 * Q(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * a >= 0, x >= 0
 *   Q(a,0) := 1
 *   Q(0,x) := 0, x != 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_Q(const double a, const double x);


/* Complementary Normalized Incomplete Gamma Function
 *
 * P(a,x) = 1/Gamma(a) Integral[ t^(a-1) e^(-t), {t,0,x} ]
 *
 * a > 0, x >= 0
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc_P(const double a, const double x);


/* Non-normalized Incomplete Gamma Function
 *
 * Gamma(a,x) := Integral[ t^(a-1) e^(-t), {t,x,Infinity} ]
 *
 * x >= 0.0
 *   Gamma(a, 0) := Gamma(a)
 *
 * exceptions: GSL_EDOM
 */
int gsl_sf_gamma_inc_e(const double a, const double x, gsl_sf_result * result);
double gsl_sf_gamma_inc(const double a, const double x);


/* Logarithm of Beta Function
 * Log[B(a,b)]
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM
 */
int gsl_sf_lnbeta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_lnbeta(const double a, const double b);

int gsl_sf_lnbeta_sgn_e(const double x, const double y, gsl_sf_result * result, double * sgn);


/* Beta Function
 * B(a,b)
 *
 * a > 0, b > 0
 * exceptions: GSL_EDOM, GSL_EOVRFLW, GSL_EUNDRFLW
 */
int gsl_sf_beta_e(const double a, const double b, gsl_sf_result * result);
double gsl_sf_beta(const double a, const double b);


/* Normalized Incomplete Beta Function
 * B_x(a,b)/B(a,b)
 *
 * a > 0, b > 0, 0 <= x <= 1
 * exceptions: GSL_EDOM, GSL_EUNDRFLW
 */
int gsl_sf_beta_inc_e(const double a, const double b, const double x, gsl_sf_result * result);
double gsl_sf_beta_inc(const double a, const double b, const double x);


/* The maximum x such that gamma(x) is not
 * considered an overflow.
 */
#define GSL_SF_GAMMA_XMAX  171.0

/* The maximum n such that gsl_sf_fact(n) does not give an overflow. */
#define GSL_SF_FACT_NMAX 170

/* The maximum n such that gsl_sf_doublefact(n) does not give an overflow. */
#define GSL_SF_DOUBLEFACT_NMAX 297

__END_DECLS

#endif /* __GSL_SF_GAMMA_H__ */

/* end inlined header: gsl/gsl_sf_gamma.h */
/* begin inlined header: gsl/gsl_sf_coulomb.h */
/* specfunc/gsl_sf_coulomb.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_COULOMB_H__
#define __GSL_SF_COULOMB_H__
/* begin inlined header: gsl/gsl_mode.h */
/* gsl_mode.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Gough and G. Jungman */

#ifndef __GSL_MODE_H__
#define __GSL_MODE_H__
/* begin inlined header: gsl/gsl_inline.h */
/* gsl_inline.h
 *
 * Copyright (C) 2008, 2009 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INLINE_H__
#define __GSL_INLINE_H__

/* In recent versiions of GCC, the inline keyword has two different
   forms: GNU and C99.

   In GNU mode we can use 'extern inline' to make inline functions
   work like macros.  The function is only inlined--it is never output
   as a definition in an object file.

   In the new C99 mode 'extern inline' has a different meaning--it
   causes the definition of the function to be output in each object
   file where it is used.  This will result in multiple-definition
   errors on linking.  The 'inline' keyword on its own (without
   extern) has the same behavior as the original GNU 'extern inline'.

   The C99 style is the default with -std=c99 in GCC 4.3.

   This header file allows either form of inline to be used by
   redefining the macros INLINE_DECL and INLINE_FUN.  These are used
   in the public header files as

        INLINE_DECL double gsl_foo (double x);
        #ifdef HAVE_INLINE
        INLINE_FUN double gsl_foo (double x) { return x+1.0; } ;
        #endif

*/

#ifdef HAVE_INLINE
#  if defined(__GNUC_STDC_INLINE__) || defined(GSL_C99_INLINE) || defined(HAVE_C99_INLINE)
#    define INLINE_DECL inline  /* use C99 inline */
#    define INLINE_FUN inline
#  else
#    define INLINE_DECL         /* use GNU extern inline */
#    define INLINE_FUN extern inline
#  endif
#else
#  define INLINE_DECL /* */
#endif

/* Range checking conditions in headers do not require any run-time
   tests of the global variable gsl_check_range.  They are enabled or
   disabled in user code at compile time with GSL_RANGE_CHECK macro.
   See also build.h. */
#define GSL_RANGE_COND(x) (x)

#endif /* __GSL_INLINE_H__ */

/* end inlined header: gsl/gsl_inline.h */
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


/* Some functions can take a mode argument. This
 * is a rough method to do things like control
 * the precision of the algorithm. This mainly
 * occurs in special functions, but we figured
 * it was ok to have a general facility.
 *
 * The mode type is 32-bit field. Most of
 * the fields are currently unused. Users
 * '|' various predefined constants to get
 * a desired mode.
 */
typedef unsigned int gsl_mode_t;


/* Here are the predefined constants.
 * Note that the precision constants
 * are special because they are used
 * to index arrays, so do not change
 * them. The precision information is
 * in the low order 3 bits of gsl_mode_t
 * (the third bit is currently unused).
 */

/* Note that "0" is double precision,
 * so that you get that by default if
 * you forget a flag.
 */
#define GSL_PREC_DOUBLE  0
#define GSL_PREC_SINGLE  1
#define GSL_PREC_APPROX  2

#ifdef HAVE_INLINE
INLINE_FUN unsigned int GSL_MODE_PREC(gsl_mode_t mt);

INLINE_FUN unsigned int
GSL_MODE_PREC(gsl_mode_t mt)
{ return  (mt & (unsigned int)7); }
#else  /* HAVE_INLINE */
#define GSL_MODE_PREC(mt) ((mt) & (unsigned int)7)
#endif /* HAVE_INLINE */


/* Here are some predefined generic modes.
 */
#define GSL_MODE_DEFAULT  0


__END_DECLS

#endif /* __GSL_MODE_H__ */

/* end inlined header: gsl/gsl_mode.h */
/* begin inlined header: gsl/gsl_sf_result.h */
/* specfunc/gsl_sf_result.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#ifndef __GSL_SF_RESULT_H__
#define __GSL_SF_RESULT_H__

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

struct gsl_sf_result_struct {
  double val;
  double err;
};
typedef struct gsl_sf_result_struct gsl_sf_result;

#define GSL_SF_RESULT_SET(r,v,e) do { (r)->val=(v); (r)->err=(e); } while(0)


struct gsl_sf_result_e10_struct {
  double val;
  double err;
  int    e10;
};
typedef struct gsl_sf_result_e10_struct gsl_sf_result_e10;


int gsl_sf_result_smash_e(const gsl_sf_result_e10 * re, gsl_sf_result * r);


__END_DECLS

#endif /* __GSL_SF_RESULT_H__ */

/* end inlined header: gsl/gsl_sf_result.h */
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


/* Normalized hydrogenic bound states, radial dependence. */

/* R_1 := 2Z sqrt(Z) exp(-Z r)
 */
int gsl_sf_hydrogenicR_1_e(const double Z, const double r, gsl_sf_result * result);
double gsl_sf_hydrogenicR_1(const double Z, const double r);

/* R_n := norm exp(-Z r/n) (2Z/n)^l Laguerre[n-l-1, 2l+1, 2Z/n r]
 *
 * normalization such that psi(n,l,r) = R_n Y_{lm}
 */
int gsl_sf_hydrogenicR_e(const int n, const int l, const double Z, const double r, gsl_sf_result * result);
double gsl_sf_hydrogenicR(const int n, const int l, const double Z, const double r);


/* Coulomb wave functions F_{lam_F}(eta,x), G_{lam_G}(eta,x)
 * and their derivatives; lam_G := lam_F - k_lam_G
 *
 * lam_F, lam_G > -0.5
 * x > 0.0
 *
 * Conventions of Abramowitz+Stegun.
 *
 * Because there can be a large dynamic range of values,
 * overflows are handled gracefully. If an overflow occurs,
 * GSL_EOVRFLW is signalled and exponent(s) are returned
 * through exp_F, exp_G. These are such that
 *
 *   F_L(eta,x)  =  fc[k_L] * exp(exp_F)
 *   G_L(eta,x)  =  gc[k_L] * exp(exp_G)
 *   F_L'(eta,x) = fcp[k_L] * exp(exp_F)
 *   G_L'(eta,x) = gcp[k_L] * exp(exp_G)
 */
int
gsl_sf_coulomb_wave_FG_e(const double eta, const double x,
                            const double lam_F,
                            const int  k_lam_G,
                            gsl_sf_result * F, gsl_sf_result * Fp,
                            gsl_sf_result * G, gsl_sf_result * Gp,
                            double * exp_F, double * exp_G);


/* F_L(eta,x), G_L(eta,x) as arrays */
int gsl_sf_coulomb_wave_FG_array(double lam_min, int kmax,
                                double eta, double x,
                                double * fc_array, double * gc_array,
                                double * F_exponent,
                                double * G_exponent
                                );


__END_DECLS

#endif /* __GSL_SF_COULOMB_H__ */

/* end inlined header: gsl/gsl_sf_coulomb.h */
/* inlined: error.h */

/* the full definition of C_L(eta) for any valid L and eta
 * [Abramowitz and Stegun 14.1.7]
 * This depends on the complex gamma function. For large
 * arguments the phase of the complex gamma function is not
 * very accurately determined. However the modulus is, and that
 * is all that we need to calculate C_L.
 *
 * This is not valid for L <= -3/2  or  L = -1.
 */
static
int
CLeta(double L, double eta, gsl_sf_result * result)
{
  gsl_sf_result ln1; /* log of numerator Gamma function */
  gsl_sf_result ln2; /* log of denominator Gamma function */
  double sgn = 1.0;
  double arg_val, arg_err;

  if(fabs(eta/(L+1.0)) < GSL_DBL_EPSILON) {
    gsl_sf_lngamma_e(L+1.0, &ln1);
  }
  else {
    gsl_sf_result p1;                 /* phase of numerator Gamma -- not used */
    gsl_sf_lngamma_complex_e(L+1.0, eta, &ln1, &p1); /* should be ok */
  }

  gsl_sf_lngamma_e(2.0*(L+1.0), &ln2);
  if(L < -1.0) sgn = -sgn;

  arg_val  = L*M_LN2 - 0.5*eta*M_PI + ln1.val - ln2.val;
  arg_err  = ln1.err + ln2.err;
  arg_err += GSL_DBL_EPSILON * (fabs(L*M_LN2) + fabs(0.5*eta*M_PI));
  return gsl_sf_exp_err_e(arg_val, arg_err, result);
}


/* Evaluate the series for Phi_L(eta,x) and Phi_L*(eta,x)
 * [Abramowitz+Stegun 14.1.5]
 * [Abramowitz+Stegun 14.1.13]
 *
 * The sequence of coefficients A_k^L is
 * manifestly well-controlled for L >= -1/2
 * and eta < 10.
 *
 * This makes sense since this is the region
 * away from threshold, and you expect
 * the evaluation to become easier as you
 * get farther from threshold.
 *
 * Empirically, this is quite well-behaved for
 *   L >= -1/2
 *   eta < 10
 *   x   < 10
 */
#if 0
static
int
coulomb_Phi_series(const double lam, const double eta, const double x,
                   double * result, double * result_star)
{
  int kmin =   5;
  int kmax = 200;
  int k;
  double Akm2 = 1.0;
  double Akm1 = eta/(lam+1.0);
  double Ak;

  double xpow = x;
  double sum  = Akm2 + Akm1*x;
  double sump = (lam+1.0)*Akm2 + (lam+2.0)*Akm1*x;
  double prev_abs_del   = fabs(Akm1*x);
  double prev_abs_del_p = (lam+2.0) * prev_abs_del;

  for(k=2; k<kmax; k++) {
    double del;
    double del_p;
    double abs_del;
    double abs_del_p;

    Ak = (2.0*eta*Akm1 - Akm2)/(k*(2.0*lam + 1.0 + k));

    xpow *= x;
    del   = Ak*xpow;
    del_p = (k+lam+1.0)*del;
    sum  += del;
    sump += del_p;

    abs_del   = fabs(del);
    abs_del_p = fabs(del_p);

    if(          abs_del/(fabs(sum)+abs_del)          < GSL_DBL_EPSILON
       &&   prev_abs_del/(fabs(sum)+prev_abs_del)     < GSL_DBL_EPSILON
       &&      abs_del_p/(fabs(sump)+abs_del_p)       < GSL_DBL_EPSILON
       && prev_abs_del_p/(fabs(sump)+prev_abs_del_p)  < GSL_DBL_EPSILON
       && k > kmin
       ) break;

    /* We need to keep track of the previous delta because when
     * eta is near zero the odd terms of the sum are very small
     * and this could lead to premature termination.
     */
    prev_abs_del   = abs_del;
    prev_abs_del_p = abs_del_p;

    Akm2 = Akm1;
    Akm1 = Ak;
  }

  *result      = sum;
  *result_star = sump;

  if(k==kmax) {
    GSL_ERROR ("error", GSL_EMAXITER);
  }
  else {
    return GSL_SUCCESS;
  }
}
#endif /* 0 */


/* Determine the connection phase, phi_lambda.
 * See coulomb_FG_series() below. We have
 * to be careful about sin(phi)->0. Note that
 * there is an underflow condition for large
 * positive eta in any case.
 */
static
int
coulomb_connection(const double lam, const double eta,
                   double * cos_phi, double * sin_phi)
{
  if(eta > -GSL_LOG_DBL_MIN/2.0*M_PI-1.0) {
    *cos_phi = 1.0;
    *sin_phi = 0.0;
    GSL_ERROR ("error", GSL_EUNDRFLW);
  }
  else if(eta > -GSL_LOG_DBL_EPSILON/(4.0*M_PI)) {
    const double eps = 2.0 * exp(-2.0*M_PI*eta);
    const double tpl = tan(M_PI * lam);
    const double dth = eps * tpl / (tpl*tpl + 1.0);
    *cos_phi = -1.0 + 0.5 * dth*dth;
    *sin_phi = -dth;
    return GSL_SUCCESS;
  }
  else {
    double X   = tanh(M_PI * eta) / tan(M_PI * lam);
    double phi = -atan(X) - (lam + 0.5) * M_PI;
    *cos_phi = cos(phi);
    *sin_phi = sin(phi);
    return GSL_SUCCESS;
  }
}


/* Evaluate the Frobenius series for F_lam(eta,x) and G_lam(eta,x).
 * Homegrown algebra. Evaluates the series for F_{lam} and
 * F_{-lam-1}, then uses
 *    G_{lam} = (F_{lam} cos(phi) - F_{-lam-1}) / sin(phi)
 * where
 *    phi = Arg[Gamma[1+lam+I eta]] - Arg[Gamma[-lam + I eta]] - (lam+1/2)Pi
 *        = Arg[Sin[Pi(-lam+I eta)] - (lam+1/2)Pi
 *        = atan2(-cos(lam Pi)sinh(eta Pi), -sin(lam Pi)cosh(eta Pi)) - (lam+1/2)Pi
 *
 *        = -atan(X) - (lam+1/2) Pi,  X = tanh(eta Pi)/tan(lam Pi)
 *
 * Not appropriate for lam <= -1/2, lam = 0, or lam >= 1/2.
 */
static
int
coulomb_FG_series(const double lam, const double eta, const double x,
                  gsl_sf_result * F, gsl_sf_result * G)
{
  const int max_iter = 800;
  gsl_sf_result ClamA;
  gsl_sf_result ClamB;
  int stat_A = CLeta(lam, eta, &ClamA);
  int stat_B = CLeta(-lam-1.0, eta, &ClamB);
  const double tlp1 = 2.0*lam + 1.0;
  const double pow_x = pow(x, lam);
  double cos_phi_lam;
  double sin_phi_lam;

  double uA_mm2 = 1.0;                  /* uA sum is for F_{lam} */
  double uA_mm1 = x*eta/(lam+1.0);
  double uA_m;
  double uB_mm2 = 1.0;                  /* uB sum is for F_{-lam-1} */
  double uB_mm1 = -x*eta/lam;
  double uB_m;
  double A_sum = uA_mm2 + uA_mm1;
  double B_sum = uB_mm2 + uB_mm1;
  double A_abs_del_prev = fabs(A_sum);
  double B_abs_del_prev = fabs(B_sum);
  gsl_sf_result FA, FB;
  int m = 2;

  int stat_conn = coulomb_connection(lam, eta, &cos_phi_lam, &sin_phi_lam);

  if(stat_conn == GSL_EUNDRFLW) {
    F->val = 0.0;  /* FIXME: should this be set to Inf too like G? */
    F->err = 0.0;
    OVERFLOW_ERROR(G);
  }

  while(m < max_iter) {
    double abs_dA;
    double abs_dB;
    uA_m = x*(2.0*eta*uA_mm1 - x*uA_mm2)/(m*(m+tlp1));
    uB_m = x*(2.0*eta*uB_mm1 - x*uB_mm2)/(m*(m-tlp1));
    A_sum += uA_m;
    B_sum += uB_m;
    abs_dA = fabs(uA_m);
    abs_dB = fabs(uB_m);
    if(m > 15) {
      /* Don't bother checking until we have gone out a little ways;
       * a minor optimization. Also make sure to check both the
       * current and the previous increment because the odd and even
       * terms of the sum can have very different behaviour, depending
       * on the value of eta.
       */
      double max_abs_dA = GSL_MAX(abs_dA, A_abs_del_prev);
      double max_abs_dB = GSL_MAX(abs_dB, B_abs_del_prev);
      double abs_A = fabs(A_sum);
      double abs_B = fabs(B_sum);
      if(   max_abs_dA/(max_abs_dA + abs_A) < 4.0*GSL_DBL_EPSILON
         && max_abs_dB/(max_abs_dB + abs_B) < 4.0*GSL_DBL_EPSILON
         ) break;
    }
    A_abs_del_prev = abs_dA;
    B_abs_del_prev = abs_dB;
    uA_mm2 = uA_mm1;
    uA_mm1 = uA_m;
    uB_mm2 = uB_mm1;
    uB_mm1 = uB_m;
    m++;
  }

  FA.val = A_sum * ClamA.val * pow_x * x;
  FA.err = fabs(A_sum) * ClamA.err * pow_x * x + 2.0*GSL_DBL_EPSILON*fabs(FA.val);
  FB.val = B_sum * ClamB.val / pow_x;
  FB.err = fabs(B_sum) * ClamB.err / pow_x + 2.0*GSL_DBL_EPSILON*fabs(FB.val);

  F->val = FA.val;
  F->err = FA.err;

  G->val = (FA.val * cos_phi_lam - FB.val)/sin_phi_lam;
  G->err = (FA.err * fabs(cos_phi_lam) + FB.err)/fabs(sin_phi_lam);

  if(m >= max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_ERROR_SELECT_2(stat_A, stat_B);
}


/* Evaluate the Frobenius series for F_0(eta,x) and G_0(eta,x).
 * See [Bardin et al., CPC 3, 73 (1972), (14)-(17)];
 * note the misprint in (17): nu_0=1 is correct, not nu_0=0.
 */
static
int
coulomb_FG0_series(const double eta, const double x,
                   gsl_sf_result * F, gsl_sf_result * G)
{
  const int max_iter = 800;
  const double x2  = x*x;
  const double tex = 2.0*eta*x;
  gsl_sf_result C0;
  int stat_CL = CLeta(0.0, eta, &C0);
  gsl_sf_result r1pie;
  int psi_stat = gsl_sf_psi_1piy_e(eta, &r1pie);
  double u_mm2 = 0.0;  /* u_0 */
  double u_mm1 = x;    /* u_1 */
  double u_m;
  double v_mm2 = 1.0;                               /* nu_0 */
  double v_mm1 = tex*(2.0*M_EULER-1.0+r1pie.val);   /* nu_1 */
  double v_m;
  double u_sum = u_mm2 + u_mm1;
  double v_sum = v_mm2 + v_mm1;
  double u_abs_del_prev = fabs(u_sum);
  double v_abs_del_prev = fabs(v_sum);
  int m = 2;
  double u_sum_err = 2.0 * GSL_DBL_EPSILON * fabs(u_sum);
  double v_sum_err = 2.0 * GSL_DBL_EPSILON * fabs(v_sum);
  double ln2x = log(2.0*x);

  while(m < max_iter) {
    double abs_du;
    double abs_dv;
    double m_mm1 = m*(m-1.0);
    u_m = (tex*u_mm1 - x2*u_mm2)/m_mm1;
    v_m = (tex*v_mm1 - x2*v_mm2 - 2.0*eta*(2*m-1)*u_m)/m_mm1;
    u_sum += u_m;
    v_sum += v_m;
    abs_du = fabs(u_m);
    abs_dv = fabs(v_m);
    u_sum_err += 2.0 * GSL_DBL_EPSILON * abs_du;
    v_sum_err += 2.0 * GSL_DBL_EPSILON * abs_dv;
    if(m > 15) {
      /* Don't bother checking until we have gone out a little ways;
       * a minor optimization. Also make sure to check both the
       * current and the previous increment because the odd and even
       * terms of the sum can have very different behaviour, depending
       * on the value of eta.
       */
      double max_abs_du = GSL_MAX(abs_du, u_abs_del_prev);
      double max_abs_dv = GSL_MAX(abs_dv, v_abs_del_prev);
      double abs_u = fabs(u_sum);
      double abs_v = fabs(v_sum);
      if(   max_abs_du/(max_abs_du + abs_u) < 40.0*GSL_DBL_EPSILON
         && max_abs_dv/(max_abs_dv + abs_v) < 40.0*GSL_DBL_EPSILON
         ) break;
    }
    u_abs_del_prev = abs_du;
    v_abs_del_prev = abs_dv;
    u_mm2 = u_mm1;
    u_mm1 = u_m;
    v_mm2 = v_mm1;
    v_mm1 = v_m;
    m++;
  }

  F->val  = C0.val * u_sum;
  F->err  = C0.err * fabs(u_sum);
  F->err += fabs(C0.val) * u_sum_err;
  F->err += 2.0 * GSL_DBL_EPSILON * fabs(F->val);

  G->val  = (v_sum + 2.0*eta*u_sum * ln2x) / C0.val;
  G->err  = (fabs(v_sum) + fabs(2.0*eta*u_sum * ln2x)) / fabs(C0.val) * fabs(C0.err/C0.val);
  G->err += (v_sum_err + fabs(2.0*eta*u_sum_err*ln2x)) / fabs(C0.val);
  G->err += 2.0 * GSL_DBL_EPSILON * fabs(G->val);

  if(m == max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return GSL_ERROR_SELECT_2(psi_stat, stat_CL);
}


/* Evaluate the Frobenius series for F_{-1/2}(eta,x) and G_{-1/2}(eta,x).
 * Homegrown algebra.
 */
static
int
coulomb_FGmhalf_series(const double eta, const double x,
                       gsl_sf_result * F, gsl_sf_result * G)
{
  const int max_iter = 800;
  const double rx  = sqrt(x);
  const double x2  = x*x;
  const double tex = 2.0*eta*x;
  gsl_sf_result Cmhalf;
  int stat_CL = CLeta(-0.5, eta, &Cmhalf);
  double u_mm2 = 1.0;                      /* u_0 */
  double u_mm1 = tex * u_mm2;              /* u_1 */
  double u_m;
  double v_mm2, v_mm1, v_m;
  double f_sum, g_sum;
  double tmp1;
  gsl_sf_result rpsi_1pe;
  gsl_sf_result rpsi_1p2e;
  int m = 2;

  gsl_sf_psi_1piy_e(eta,     &rpsi_1pe);
  gsl_sf_psi_1piy_e(2.0*eta, &rpsi_1p2e);

  v_mm2 = 2.0*M_EULER - M_LN2 - rpsi_1pe.val + 2.0*rpsi_1p2e.val;
  v_mm1 = tex*(v_mm2 - 2.0*u_mm2);

  f_sum = u_mm2 + u_mm1;
  g_sum = v_mm2 + v_mm1;

  while(m < max_iter) {
    double m2 = m*m;
    u_m = (tex*u_mm1 - x2*u_mm2)/m2;
    v_m = (tex*v_mm1 - x2*v_mm2 - 2.0*m*u_m)/m2;
    f_sum += u_m;
    g_sum += v_m;
    if(   f_sum != 0.0
       && g_sum != 0.0
       && (fabs(u_m/f_sum) + fabs(v_m/g_sum) < 10.0*GSL_DBL_EPSILON)) break;
    u_mm2 = u_mm1;
    u_mm1 = u_m;
    v_mm2 = v_mm1;
    v_mm1 = v_m;
    m++;
  }

  F->val = Cmhalf.val * rx * f_sum;
  F->err = Cmhalf.err * fabs(rx * f_sum) + 2.0*GSL_DBL_EPSILON*fabs(F->val);

  tmp1 = f_sum*log(x);
  G->val = -rx*(tmp1 + g_sum)/Cmhalf.val;
  G->err = fabs(rx)*(fabs(tmp1) + fabs(g_sum))/fabs(Cmhalf.val) * fabs(Cmhalf.err/Cmhalf.val);

  if(m == max_iter)
    GSL_ERROR ("error", GSL_EMAXITER);
  else
    return stat_CL;
}


/* Evolve the backwards recurrence for F,F'.
 *
 *    F_{lam-1}  = (S_lam F_lam + F_lam') / R_lam
 *    F_{lam-1}' = (S_lam F_{lam-1} - R_lam F_lam)
 * where
 *    R_lam = sqrt(1 + (eta/lam)^2)
 *    S_lam = lam/x + eta/lam
 *
 */
static
int
coulomb_F_recur(double lam_min, int kmax,
                double eta, double x,
                double F_lam_max, double Fp_lam_max,
                double * F_lam_min, double * Fp_lam_min
                )
{
  double x_inv = 1.0/x;
  double fcl = F_lam_max;
  double fpl = Fp_lam_max;
  double lam_max = lam_min + kmax;
  double lam = lam_max;
  int k;

  for(k=kmax-1; k>=0; k--) {
    double el = eta/lam;
    double rl = hypot(1.0, el);
    double sl = el  + lam*x_inv;
    double fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fpl    =  fc_lm1*sl - fcl*rl;
    fcl    =  fc_lm1;
    lam -= 1.0;
  }

  *F_lam_min  = fcl;
  *Fp_lam_min = fpl;
  return GSL_SUCCESS;
}


/* Evolve the forward recurrence for G,G'.
 *
 *   G_{lam+1}  = (S_lam G_lam - G_lam')/R_lam
 *   G_{lam+1}' = R_{lam+1} G_lam - S_lam G_{lam+1}
 *
 * where S_lam and R_lam are as above in the F recursion.
 */
static
int
coulomb_G_recur(const double lam_min, const int kmax,
                const double eta, const double x,
                const double G_lam_min, const double Gp_lam_min,
                double * G_lam_max, double * Gp_lam_max
                )
{
  double x_inv = 1.0/x;
  double gcl = G_lam_min;
  double gpl = Gp_lam_min;
  double lam = lam_min + 1.0;
  int k;

  for(k=1; k<=kmax; k++) {
    double el = eta/lam;
    double rl = hypot(1.0, el);
    double sl = el + lam*x_inv;
    double gcl1 = (sl*gcl - gpl)/rl;
    gpl   = rl*gcl - sl*gcl1;
    gcl   = gcl1;
    lam += 1.0;
  }

  *G_lam_max  = gcl;
  *Gp_lam_max = gpl;
  return GSL_SUCCESS;
}


/* Evaluate the first continued fraction, giving
 * the ratio F'/F at the upper lambda value.
 * We also determine the sign of F at that point,
 * since it is the sign of the last denominator
 * in the continued fraction.
 */
static
int
coulomb_CF1(double lambda,
            double eta, double x,
            double * fcl_sign,
            double * result,
            int * count
            )
{
  const double CF1_small = 1.e-30;
  const double CF1_abort = 1.0e+05;
  const double CF1_acc   = 2.0*GSL_DBL_EPSILON;
  const double x_inv     = 1.0/x;
  const double px        = lambda + 1.0 + CF1_abort;

  double pk = lambda + 1.0;
  double F  = eta/pk + pk*x_inv;
  double D, C;
  double df;

  *fcl_sign = 1.0;
  *count = 0;

  if(fabs(F) < CF1_small) F = CF1_small;
  D = 0.0;
  C = F;

  do {
    double pk1 = pk + 1.0;
    double ek  = eta / pk;
    double rk2 = 1.0 + ek*ek;
    double tk  = (pk + pk1)*(x_inv + ek/pk1);
    D   =  tk - rk2 * D;
    C   =  tk - rk2 / C;
    if(fabs(C) < CF1_small) C = CF1_small;
    if(fabs(D) < CF1_small) D = CF1_small;
    D = 1.0/D;
    df = D * C;
    F  = F * df;
    if(D < 0.0) {
      /* sign of result depends on sign of denominator */
      *fcl_sign = - *fcl_sign;
    }
    pk = pk1;
    if( pk > px ) {
      *result = F;
      GSL_ERROR ("error", GSL_ERUNAWAY);
    }
    ++(*count);
  }
  while(fabs(df-1.0) > CF1_acc);

  *result = F;
  return GSL_SUCCESS;
}


#if 0
static
int
old_coulomb_CF1(const double lambda,
                double eta, double x,
                double * fcl_sign,
                double * result
                )
{
  const double CF1_abort = 1.e5;
  const double CF1_acc   = 10.0*GSL_DBL_EPSILON;
  const double x_inv     = 1.0/x;
  const double px        = lambda + 1.0 + CF1_abort;

  double pk = lambda + 1.0;

  double D;
  double df;

  double F;
  double p;
  double pk1;
  double ek;

  double fcl = 1.0;

  double tk;

  while(1) {
    ek = eta/pk;
    F = (ek + pk*x_inv)*fcl + (fcl - 1.0)*x_inv;
    pk1 = pk + 1.0;
    if(fabs(eta*x + pk*pk1) > CF1_acc) break;
    fcl = (1.0 + ek*ek)/(1.0 + eta*eta/(pk1*pk1));
    pk = 2.0 + pk;
  }

  D  = 1.0/((pk + pk1)*(x_inv + ek/pk1));
  df = -fcl*(1.0 + ek*ek)*D;

  if(fcl != 1.0) fcl = -1.0;
  if(D    < 0.0) fcl = -fcl;

  F = F + df;

  p = 1.0;
  do {
    pk = pk1;
    pk1 = pk + 1.0;
    ek  = eta / pk;
    tk  = (pk + pk1)*(x_inv + ek/pk1);
    D   =  tk - D*(1.0+ek*ek);
    if(fabs(D) < sqrt(CF1_acc)) {
      p += 1.0;
      if(p > 2.0) {
        printf("HELP............\n");
      }
    }
    D = 1.0/D;
    if(D < 0.0) {
      /* sign of result depends on sign of denominator */
      fcl = -fcl;
    }
    df = df*(D*tk - 1.0);
    F  = F + df;
    if( pk > px ) {
      GSL_ERROR ("error", GSL_ERUNAWAY);
    }
  }
  while(fabs(df) > fabs(F)*CF1_acc);

  *fcl_sign = fcl;
  *result = F;
  return GSL_SUCCESS;
}
#endif /* 0 */


/* Evaluate the second continued fraction to
 * obtain the ratio
 *    (G' + i F')/(G + i F) := P + i Q
 * at the specified lambda value.
 */
static
int
coulomb_CF2(const double lambda, const double eta, const double x,
            double * result_P, double * result_Q, int * count
            )
{
  int status = GSL_SUCCESS;

  const double CF2_acc   = 4.0*GSL_DBL_EPSILON;
  const double CF2_abort = 2.0e+05;

  const double wi    = 2.0*eta;
  const double x_inv = 1.0/x;
  const double e2mm1 = eta*eta + lambda*(lambda + 1.0);

  double ar = -e2mm1;
  double ai =  eta;

  double br =  2.0*(x - eta);
  double bi =  2.0;

  double dr =  br/(br*br + bi*bi);
  double di = -bi/(br*br + bi*bi);

  double dp = -x_inv*(ar*di + ai*dr);
  double dq =  x_inv*(ar*dr - ai*di);

  double A, B, C, D;

  double pk =  0.0;
  double P  =  0.0;
  double Q  =  1.0 - eta*x_inv;

  *count = 0;

  do {
    P += dp;
    Q += dq;
    pk += 2.0;
    ar += pk;
    ai += wi;
    bi += 2.0;
    D  = ar*dr - ai*di + br;
    di = ai*dr + ar*di + bi;
    C  = 1.0/(D*D + di*di);
    dr =  C*D;
    di = -C*di;
    A  = br*dr - bi*di - 1.;
    B  = bi*dr + br*di;
    C  = dp*A  - dq*B;
    dq = dp*B  + dq*A;
    dp = C;
    if(pk > CF2_abort) {
      status = GSL_ERUNAWAY;
      break;
    }
    ++(*count);
  }
  while(fabs(dp)+fabs(dq) > (fabs(P)+fabs(Q))*CF2_acc);

  if(Q < CF2_abort*GSL_DBL_EPSILON*fabs(P)) {
    status = GSL_ELOSS;
  }

  *result_P = P;
  *result_Q = Q;
  return status;
}


/* WKB evaluation of F, G. Assumes  0 < x < turning point.
 * Overflows are trapped, GSL_EOVRFLW is signalled,
 * and an exponent is returned such that:
 *
 *   result_F = fjwkb * exp(-exponent)
 *   result_G = gjwkb * exp( exponent)
 *
 * See [Biedenharn et al. Phys. Rev. 97, 542-554 (1955), Section IV]
 *
 * Unfortunately, this is not very accurate in general. The
 * test cases typically have 3-4 digits of precision. One could
 * argue that this is ok for general use because, for instance,
 * F is exponentially small in this region and so the absolute
 * accuracy is still roughly acceptable. But it would be better
 * to have a systematic method for improving the precision. See
 * the Abad+Sesma method discussion below.
 */
static
int
coulomb_jwkb(const double lam, const double eta, const double x,
             gsl_sf_result * fjwkb, gsl_sf_result * gjwkb,
             double * exponent)
{
  const double llp1      = lam*(lam+1.0) + 6.0/35.0;
  const double llp1_eff  = GSL_MAX(llp1, 0.0);
  const double rho_ghalf = sqrt(x*(2.0*eta - x) + llp1_eff);
  const double sinh_arg  = sqrt(llp1_eff/(eta*eta+llp1_eff)) * rho_ghalf / x;
  const double sinh_inv  = log(sinh_arg + hypot(1.0,sinh_arg));

  const double phi = fabs(rho_ghalf - eta*atan2(rho_ghalf,x-eta) - sqrt(llp1_eff) * sinh_inv);

  const double zeta_half = pow(3.0*phi/2.0, 1.0/3.0);
  const double prefactor = sqrt(M_PI*phi*x/(6.0 * rho_ghalf));

  double F = prefactor * 3.0/zeta_half;
  double G = prefactor * 3.0/zeta_half; /* Note the sqrt(3) from Bi normalization */
  double F_exp;
  double G_exp;

  const double airy_scale_exp = phi;
  gsl_sf_result ai;
  gsl_sf_result bi;
  gsl_sf_airy_Ai_scaled_e(zeta_half*zeta_half, GSL_MODE_DEFAULT, &ai);
  gsl_sf_airy_Bi_scaled_e(zeta_half*zeta_half, GSL_MODE_DEFAULT, &bi);
  F *= ai.val;
  G *= bi.val;
  F_exp = log(F) - airy_scale_exp;
  G_exp = log(G) + airy_scale_exp;

  if(G_exp >= GSL_LOG_DBL_MAX) {
    fjwkb->val = F;
    gjwkb->val = G;
    fjwkb->err = 1.0e-3 * fabs(F); /* FIXME: real error here ... could be smaller */
    gjwkb->err = 1.0e-3 * fabs(G);
    *exponent = airy_scale_exp;
    GSL_ERROR ("error", GSL_EOVRFLW);
  }
  else {
    fjwkb->val = exp(F_exp);
    gjwkb->val = exp(G_exp);
    fjwkb->err = 1.0e-3 * fabs(fjwkb->val);
    gjwkb->err = 1.0e-3 * fabs(gjwkb->val);
    *exponent = 0.0;
    return GSL_SUCCESS;
  }
}


/* Asymptotic evaluation of F and G below the minimal turning point.
 *
 * This is meant to be a drop-in replacement for coulomb_jwkb().
 * It uses the expressions in [Abad+Sesma]. This requires some
 * work because I am not sure where it is valid. They mumble
 * something about |x| < |lam|^(-1/2) or 8|eta x| > lam when |x| < 1.
 * This seems true, but I thought the result was based on a uniform
 * expansion and could be controlled by simply using more terms.
 */
#if 0
static
int
coulomb_AS_xlt2eta(const double lam, const double eta, const double x,
                   gsl_sf_result * f_AS, gsl_sf_result * g_AS,
                   double * exponent)
{
  /* no time to do this now... */
}
#endif /* 0 */



/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_coulomb_wave_FG_e(const double eta, const double x,
                            const double lam_F,
                            const int  k_lam_G,      /* lam_G = lam_F - k_lam_G */
                            gsl_sf_result * F, gsl_sf_result * Fp,
                            gsl_sf_result * G, gsl_sf_result * Gp,
                            double * exp_F, double * exp_G)
{
  const double lam_G = lam_F - k_lam_G;

  if(x < 0.0 || lam_F <= -0.5 || lam_G <= -0.5) {
    GSL_SF_RESULT_SET(F,  0.0, 0.0);
    GSL_SF_RESULT_SET(Fp, 0.0, 0.0);
    GSL_SF_RESULT_SET(G,  0.0, 0.0);
    GSL_SF_RESULT_SET(Gp, 0.0, 0.0);
    *exp_F = 0.0;
    *exp_G = 0.0;
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(x == 0.0) {
    gsl_sf_result C0;
    CLeta(0.0, eta, &C0);
    GSL_SF_RESULT_SET(F,  0.0, 0.0);
    GSL_SF_RESULT_SET(Fp, 0.0, 0.0);
    GSL_SF_RESULT_SET(G,  0.0, 0.0); /* FIXME: should be Inf */
    GSL_SF_RESULT_SET(Gp, 0.0, 0.0); /* FIXME: should be Inf */
    *exp_F = 0.0;
    *exp_G = 0.0;
    if(lam_F == 0.0){
      GSL_SF_RESULT_SET(Fp, C0.val, C0.err);
    }
    if(lam_G == 0.0) {
      GSL_SF_RESULT_SET(Gp, 1.0/C0.val, fabs(C0.err/C0.val)/fabs(C0.val));
    }
    GSL_ERROR ("domain error", GSL_EDOM);
    /* After all, since we are asking for G, this is a domain error... */
  }
  else if(x < 1.2 && 2.0*M_PI*eta < 0.9*(-GSL_LOG_DBL_MIN) && fabs(eta*x) < 10.0) {
    /* Reduce to a small lambda value and use the series
     * representations for F and G. We cannot allow eta to
     * be large and positive because the connection formula
     * for G_lam is badly behaved due to an underflow in sin(phi_lam)
     * [see coulomb_FG_series() and coulomb_connection() above].
     * Note that large negative eta is ok however.
     */
    const double SMALL = GSL_SQRT_DBL_EPSILON;
    const int N    = (int)(lam_F + 0.5);
    const int span = GSL_MAX(k_lam_G, N);
    const double lam_min = lam_F - N;    /* -1/2 <= lam_min < 1/2 */
    double F_lam_F, Fp_lam_F;
    double G_lam_G = 0.0, Gp_lam_G = 0.0;
    double F_lam_F_err, Fp_lam_F_err;
    double Fp_over_F_lam_F;
    double F_sign_lam_F;
    double F_lam_min_unnorm, Fp_lam_min_unnorm;
    double Fp_over_F_lam_min;
    gsl_sf_result F_lam_min;
    gsl_sf_result G_lam_min, Gp_lam_min;
    double F_scale;
    double Gerr_frac;
    double F_scale_frac_err;
    double F_unnorm_frac_err;

    /* Determine F'/F at lam_F. */
    int CF1_count;
    int stat_CF1 = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count);

    int stat_ser;
    int stat_Fr;
    int stat_Gr;

    /* Recurse down with unnormalized F,F' values. */
    F_lam_F  = SMALL;
    Fp_lam_F = Fp_over_F_lam_F * F_lam_F;
    if(span != 0) {
      stat_Fr = coulomb_F_recur(lam_min, span, eta, x,
                                F_lam_F, Fp_lam_F,
                                &F_lam_min_unnorm, &Fp_lam_min_unnorm
                                );
    }
    else {
      F_lam_min_unnorm  =  F_lam_F;
      Fp_lam_min_unnorm = Fp_lam_F;
      stat_Fr = GSL_SUCCESS;
    }

    /* Determine F and G at lam_min. */
    if(lam_min == -0.5) {
      stat_ser = coulomb_FGmhalf_series(eta, x, &F_lam_min, &G_lam_min);
    }
    else if(lam_min == 0.0) {
      stat_ser = coulomb_FG0_series(eta, x, &F_lam_min, &G_lam_min);
    }
    else if(lam_min == 0.5) {
      /* This cannot happen. */
      F->val  = F_lam_F;
      F->err  = 2.0 * GSL_DBL_EPSILON * fabs(F->val);
      Fp->val = Fp_lam_F;
      Fp->err = 2.0 * GSL_DBL_EPSILON * fabs(Fp->val);
      G->val  = G_lam_G;
      G->err  = 2.0 * GSL_DBL_EPSILON * fabs(G->val);
      Gp->val = Gp_lam_G;
      Gp->err = 2.0 * GSL_DBL_EPSILON * fabs(Gp->val);
      *exp_F = 0.0;
      *exp_G = 0.0;
      GSL_ERROR ("error", GSL_ESANITY);
    }
    else {
      stat_ser = coulomb_FG_series(lam_min, eta, x, &F_lam_min, &G_lam_min);
    }

    /* Determine remaining quantities. */
    Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm;
    Gp_lam_min.val  = Fp_over_F_lam_min*G_lam_min.val - 1.0/F_lam_min.val;
    Gp_lam_min.err  = fabs(Fp_over_F_lam_min)*G_lam_min.err;
    Gp_lam_min.err += fabs(1.0/F_lam_min.val) * fabs(F_lam_min.err/F_lam_min.val);
    F_scale     = F_lam_min.val / F_lam_min_unnorm;

    /* Apply scale to the original F,F' values. */
    F_scale_frac_err  = fabs(F_lam_min.err/F_lam_min.val);
    F_unnorm_frac_err = 2.0*GSL_DBL_EPSILON*(CF1_count+span+1);
    F_lam_F     *= F_scale;
    F_lam_F_err  = fabs(F_lam_F) * (F_unnorm_frac_err + F_scale_frac_err);
    Fp_lam_F    *= F_scale;
    Fp_lam_F_err = fabs(Fp_lam_F) * (F_unnorm_frac_err + F_scale_frac_err);

    /* Recurse up to get the required G,G' values. */
    stat_Gr = coulomb_G_recur(lam_min, GSL_MAX(N-k_lam_G,0), eta, x,
                              G_lam_min.val, Gp_lam_min.val,
                              &G_lam_G, &Gp_lam_G
                              );

    F->val  = F_lam_F;
    F->err  = F_lam_F_err;
    F->err += 2.0 * GSL_DBL_EPSILON * fabs(F_lam_F);

    Fp->val  = Fp_lam_F;
    Fp->err  = Fp_lam_F_err;
    Fp->err += 2.0 * GSL_DBL_EPSILON * fabs(Fp_lam_F);

    Gerr_frac = fabs(G_lam_min.err/G_lam_min.val) + fabs(Gp_lam_min.err/Gp_lam_min.val);

    G->val  = G_lam_G;
    G->err  = Gerr_frac * fabs(G_lam_G);
    G->err += 2.0 * (CF1_count+1) * GSL_DBL_EPSILON * fabs(G->val);

    Gp->val  = Gp_lam_G;
    Gp->err  = Gerr_frac * fabs(Gp->val);
    Gp->err += 2.0 * (CF1_count+1) * GSL_DBL_EPSILON * fabs(Gp->val);

    *exp_F = 0.0;
    *exp_G = 0.0;

    return GSL_ERROR_SELECT_4(stat_ser, stat_CF1, stat_Fr, stat_Gr);
  }
  else if(x < 2.0*eta) {
    /* Use WKB approximation to obtain F and G at the two
     * lambda values, and use the Wronskian and the
     * continued fractions for F'/F to obtain F' and G'.
     */
    gsl_sf_result F_lam_F, G_lam_F;
    gsl_sf_result F_lam_G, G_lam_G;
    double exp_lam_F, exp_lam_G;
    int stat_lam_F;
    int stat_lam_G;
    int stat_CF1_lam_F;
    [[maybe_unused]] int stat_CF1_lam_G;
    int CF1_count;
    double Fp_over_F_lam_F;
    double Fp_over_F_lam_G;
    double F_sign_lam_F;
    double F_sign_lam_G;

    stat_lam_F = coulomb_jwkb(lam_F, eta, x, &F_lam_F, &G_lam_F, &exp_lam_F);
    if(k_lam_G == 0) {
      stat_lam_G = stat_lam_F;
      F_lam_G = F_lam_F;
      G_lam_G = G_lam_F;
      exp_lam_G = exp_lam_F;
    }
    else {
      stat_lam_G = coulomb_jwkb(lam_G, eta, x, &F_lam_G, &G_lam_G, &exp_lam_G);
    }

    stat_CF1_lam_F = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count);
    if(k_lam_G == 0) {
      stat_CF1_lam_G  = stat_CF1_lam_F;
      F_sign_lam_G    = F_sign_lam_F;
      Fp_over_F_lam_G = Fp_over_F_lam_F;
    }
    else {
      stat_CF1_lam_G = coulomb_CF1(lam_G, eta, x, &F_sign_lam_G, &Fp_over_F_lam_G, &CF1_count);
    }

    F->val = F_lam_F.val;
    F->err = F_lam_F.err;

    G->val = G_lam_G.val;
    G->err = G_lam_G.err;

    Fp->val  = Fp_over_F_lam_F * F_lam_F.val;
    Fp->err  = fabs(Fp_over_F_lam_F) * F_lam_F.err;
    Fp->err += 2.0*GSL_DBL_EPSILON*fabs(Fp->val);

    Gp->val  = Fp_over_F_lam_G * G_lam_G.val - 1.0/F_lam_G.val;
    Gp->err  = fabs(Fp_over_F_lam_G) * G_lam_G.err;
    Gp->err += fabs(1.0/F_lam_G.val) * fabs(F_lam_G.err/F_lam_G.val);

    *exp_F = exp_lam_F;
    *exp_G = exp_lam_G;

    if(stat_lam_F == GSL_EOVRFLW || stat_lam_G == GSL_EOVRFLW) {
      GSL_ERROR ("overflow", GSL_EOVRFLW);
    }
    else {
      return GSL_ERROR_SELECT_2(stat_lam_F, stat_lam_G);
    }
  }
  else {
    /* x > 2 eta, so we know that we can find a lambda value such
     * that x is above the turning point. We do this, evaluate
     * using Steed's method at that oscillatory point, then
     * use recursion on F and G to obtain the required values.
     *
     * lam_0   = a value of lambda such that x is below the turning point
     * lam_min = minimum of lam_0 and the requested lam_G, since
     *           we must go at least as low as lam_G
     */
    const double SMALL = GSL_SQRT_DBL_EPSILON;
    const double C = sqrt(1.0 + 4.0*x*(x-2.0*eta));
    const int N = ceil(lam_F - C + 0.5);
    const double lam_0   = lam_F - GSL_MAX(N, 0);
    const double lam_min = GSL_MIN(lam_0, lam_G);
    double F_lam_F, Fp_lam_F;
    double G_lam_G, Gp_lam_G;
    double F_lam_min_unnorm, Fp_lam_min_unnorm;
    double F_lam_min;
    [[maybe_unused]] double Fp_lam_min;
    double G_lam_min, Gp_lam_min;
    double Fp_over_F_lam_F;
    double Fp_over_F_lam_min;
    double F_sign_lam_F, F_sign_lam_min;
    double P_lam_min, Q_lam_min;
    double alpha;
    double gamma;
    double F_scale;

    int CF1_count;
    int CF2_count;
    int stat_CF1 = coulomb_CF1(lam_F, eta, x, &F_sign_lam_F, &Fp_over_F_lam_F, &CF1_count);
    int stat_CF2;
    int stat_Fr;
    int stat_Gr;

    int F_recur_count;
    int G_recur_count;

    double err_amplify;

    F_lam_F  = F_sign_lam_F * SMALL;  /* unnormalized */
    Fp_lam_F = Fp_over_F_lam_F * F_lam_F;

    /* Backward recurrence to get F,Fp at lam_min */
    F_recur_count = GSL_MAX(k_lam_G, N);
    stat_Fr = coulomb_F_recur(lam_min, F_recur_count, eta, x,
                              F_lam_F, Fp_lam_F,
                              &F_lam_min_unnorm, &Fp_lam_min_unnorm
                              );
    Fp_over_F_lam_min = Fp_lam_min_unnorm / F_lam_min_unnorm;

    /* Steed evaluation to complete evaluation of F,Fp,G,Gp at lam_min */
    stat_CF2 = coulomb_CF2(lam_min, eta, x, &P_lam_min, &Q_lam_min, &CF2_count);
    alpha = Fp_over_F_lam_min - P_lam_min;
    gamma = alpha/Q_lam_min;

    F_sign_lam_min = GSL_SIGN(F_lam_min_unnorm) ;

    F_lam_min  = F_sign_lam_min / sqrt(alpha*alpha/Q_lam_min + Q_lam_min);
    Fp_lam_min = Fp_over_F_lam_min * F_lam_min;
    G_lam_min  = gamma * F_lam_min;
    Gp_lam_min = (P_lam_min * gamma - Q_lam_min) * F_lam_min;

    /* Apply scale to values of F,Fp at lam_F (the top). */
    F_scale = F_lam_min / F_lam_min_unnorm;
    F_lam_F  *= F_scale;
    Fp_lam_F *= F_scale;

    /* Forward recurrence to get G,Gp at lam_G (the top). */
    G_recur_count = GSL_MAX(N-k_lam_G,0);
    stat_Gr = coulomb_G_recur(lam_min, G_recur_count, eta, x,
                              G_lam_min, Gp_lam_min,
                              &G_lam_G, &Gp_lam_G
                              );

    err_amplify = CF1_count + CF2_count + F_recur_count + G_recur_count + 1;

    F->val  = F_lam_F;
    F->err  = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(F->val);

    Fp->val = Fp_lam_F;
    Fp->err = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(Fp->val);

    G->val  = G_lam_G;
    G->err  = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(G->val);

    Gp->val = Gp_lam_G;
    Gp->err = 8.0*err_amplify*GSL_DBL_EPSILON * fabs(Gp->val);

    *exp_F = 0.0;
    *exp_G = 0.0;

    return GSL_ERROR_SELECT_4(stat_CF1, stat_CF2, stat_Fr, stat_Gr);
  }
}


int
gsl_sf_coulomb_wave_FG_array(double lam_min, int kmax,
                                  double eta, double x,
                                  double * fc_array, double * gc_array,
                                  double * F_exp, double * G_exp)
{
  const double x_inv = 1.0/x;
  const double lam_max = lam_min + kmax;
  gsl_sf_result F, Fp;
  gsl_sf_result G, Gp;

  int stat_FG = gsl_sf_coulomb_wave_FG_e(eta, x, lam_max, kmax,
                                            &F, &Fp, &G, &Gp, F_exp, G_exp);

  double fcl  = F.val;
  double fpl = Fp.val;
  double lam = lam_max;
  int k;

  double gcl, gpl;

  fc_array[kmax] = F.val;

  for(k=kmax-1; k>=0; k--) {
    double el = eta/lam;
    double rl = hypot(1.0, el);
    double sl = el  + lam*x_inv;
    double fc_lm1;
    fc_lm1 = (fcl*sl + fpl)/rl;
    fc_array[k] = fc_lm1;
    fpl         =  fc_lm1*sl - fcl*rl;
    fcl         =  fc_lm1;
    lam -= 1.0;
  }

  gcl = G.val;
  gpl = Gp.val;
  lam = lam_min + 1.0;

  gc_array[0] = G.val;

  for(k=1; k<=kmax; k++) {
    double el = eta/lam;
    double rl = hypot(1.0, el);
    double sl = el + lam*x_inv;
    double gcl1 = (sl*gcl - gpl)/rl;
    gc_array[k] = gcl1;
    gpl         = rl*gcl - sl*gcl1;
    gcl         = gcl1;
    lam += 1.0;
  }

  return stat_FG;
}


/* end inlined source: specfunc/coulomb.c */
                                 /* gsl_sf_coulomb_wave_FG_array (future)  */


/* ======================================================================== */
/* cdf/chi-squared support: erfc, gamma_inc, cdf/gamma, cdf/chisq           */
/* ======================================================================== */

/* begin inlined source: specfunc/erfc.c (gsl_sf_erfc_e only) */
/* specfunc/erfc.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  J. Theiler (modifications by G. Jungman) */


/* (stripped: #include directives — already handled above) */
/* (only gsl_sf_erfc_e is included; log_erfc, erf, erf_Z, erf_Q, hazard
 * are omitted as they are not needed by MARLEY) */

#define LogRootPi_  0.57236494292470008706


static double erfc8_sum(double x)
{
  /* estimates erfc(x) valid for 8 < x < 100 */
  /* This is based on index 5725 in Hart et al */

  static double P[] = {
      2.97886562639399288862,
      7.409740605964741794425,
      6.1602098531096305440906,
      5.019049726784267463450058,
      1.275366644729965952479585264,
      0.5641895835477550741253201704
  };
  static double Q[] = {
      3.3690752069827527677,
      9.608965327192787870698,
      17.08144074746600431571095,
      12.0489519278551290360340491,
      9.396034016235054150430579648,
      2.260528520767326969591866945,
      1.0
  };
  double num=0.0, den=0.0;
  int i;

  num = P[5];
  for (i=4; i>=0; --i) {
      num = x*num + P[i];
  }
  den = Q[6];
  for (i=5; i>=0; --i) {
      den = x*den + Q[i];
  }

  return num/den;
}

inline
static double erfc8(double x)
{
  double e;
  e = erfc8_sum(x);
  e *= exp(-x*x);
  return e;
}

#if 0
/* Abramowitz+Stegun, 7.2.14 */
static double erfcasympsum(double x)
{
  int i;
  double e = 1.;
  double coef = 1.;
  for (i=1; i<5; ++i) {
    /* coef *= -(2*i-1)/(2*x*x); ??? [GJ] */
    coef *= -(2*i+1)/(i*(4*x*x*x*x));
    e += coef;
    /*
    if (fabs(coef) < 1.0e-15) break;
    if (fabs(coef) > 1.0e10) break;

    [GJ]: These tests are not useful. This function is only
    used below. Took them out; they gum up the pipeline.
    */
  }
  return e;
}
#endif /* 0 */


/* Chebyshev fit for erfc((t+1)/2), -1 < t < 1
 */
static double erfc_xlt1_data[20] = {
  1.06073416421769980345174155056,
 -0.42582445804381043569204735291,
  0.04955262679620434040357683080,
  0.00449293488768382749558001242,
 -0.00129194104658496953494224761,
 -0.00001836389292149396270416979,
  0.00002211114704099526291538556,
 -5.23337485234257134673693179020e-7,
 -2.78184788833537885382530989578e-7,
  1.41158092748813114560316684249e-8,
  2.72571296330561699984539141865e-9,
 -2.06343904872070629406401492476e-10,
 -2.14273991996785367924201401812e-11,
  2.22990255539358204580285098119e-12,
  1.36250074650698280575807934155e-13,
 -1.95144010922293091898995913038e-14,
 -6.85627169231704599442806370690e-16,
  1.44506492869699938239521607493e-16,
  2.45935306460536488037576200030e-18,
 -9.29599561220523396007359328540e-19
};
static cheb_series erfc_xlt1_cs = {
  erfc_xlt1_data,
  19,
  -1, 1,
  12
};

/* Chebyshev fit for erfc(x) exp(x^2), 1 < x < 5, x = 2t + 3, -1 < t < 1
 */
static double erfc_x15_data[25] = {
  0.44045832024338111077637466616,
 -0.143958836762168335790826895326,
  0.044786499817939267247056666937,
 -0.013343124200271211203618353102,
  0.003824682739750469767692372556,
 -0.001058699227195126547306482530,
  0.000283859419210073742736310108,
 -0.000073906170662206760483959432,
  0.000018725312521489179015872934,
 -4.62530981164919445131297264430e-6,
  1.11558657244432857487884006422e-6,
 -2.63098662650834130067808832725e-7,
  6.07462122724551777372119408710e-8,
 -1.37460865539865444777251011793e-8,
  3.05157051905475145520096717210e-9,
 -6.65174789720310713757307724790e-10,
  1.42483346273207784489792999706e-10,
 -3.00141127395323902092018744545e-11,
  6.22171792645348091472914001250e-12,
 -1.26994639225668496876152836555e-12,
  2.55385883033257575402681845385e-13,
 -5.06258237507038698392265499770e-14,
  9.89705409478327321641264227110e-15,
 -1.90685978789192181051961024995e-15,
  3.50826648032737849245113757340e-16
};
static cheb_series erfc_x15_cs = {
  erfc_x15_data,
  24,
  -1, 1,
  16
};

/* Chebyshev fit for erfc(x) x exp(x^2), 5 < x < 10, x = (5t + 15)/2, -1 < t < 1
 */
static double erfc_x510_data[20] = {
  1.11684990123545698684297865808,
  0.003736240359381998520654927536,
 -0.000916623948045470238763619870,
  0.000199094325044940833965078819,
 -0.000040276384918650072591781859,
  7.76515264697061049477127605790e-6,
 -1.44464794206689070402099225301e-6,
  2.61311930343463958393485241947e-7,
 -4.61833026634844152345304095560e-8,
  8.00253111512943601598732144340e-9,
 -1.36291114862793031395712122089e-9,
  2.28570483090160869607683087722e-10,
 -3.78022521563251805044056974560e-11,
  6.17253683874528285729910462130e-12,
 -9.96019290955316888445830597430e-13,
  1.58953143706980770269506726000e-13,
 -2.51045971047162509999527428316e-14,
  3.92607828989125810013581287560e-15,
 -6.07970619384160374392535453420e-16,
  9.12600607264794717315507477670e-17
};
static cheb_series erfc_x510_cs = {
  erfc_x510_data,
  19,
  -1, 1,
  12
};

#if 0
inline
static double
erfc_asymptotic(double x)
{
  return exp(-x*x)/x * erfcasympsum(x) / M_SQRTPI;
}
inline
static double
log_erfc_asymptotic(double x)
{
  return log(erfcasympsum(x)/x) - x*x - LogRootPi_;
}
#endif /* 0 */


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_sf_erfc_e(double x, gsl_sf_result * result)
{
  const double ax = fabs(x);
  double e_val, e_err;

  /* CHECK_POINTER(result) */

  if(ax <= 1.0) {
    double t = 2.0*ax - 1.0;
    gsl_sf_result c;
    cheb_eval_e(&erfc_xlt1_cs, t, &c);
    e_val = c.val;
    e_err = c.err;
  }
  else if(ax <= 5.0) {
    double ex2 = exp(-x*x);
    double t = 0.5*(ax-3.0);
    gsl_sf_result c;
    cheb_eval_e(&erfc_x15_cs, t, &c);
    e_val = ex2 * c.val;
    e_err = ex2 * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON);
  }
  else if(ax < 10.0) {
    double exterm = exp(-x*x) / ax;
    double t = (2.0*ax - 15.0)/5.0;
    gsl_sf_result c;
    cheb_eval_e(&erfc_x510_cs, t, &c);
    e_val = exterm * c.val;
    e_err = exterm * (c.err + 2.0*fabs(x)*GSL_DBL_EPSILON + GSL_DBL_EPSILON);
  }
  else {
    e_val = erfc8(ax);
    e_err = (x*x + 1.0) * GSL_DBL_EPSILON * fabs(e_val);
  }

  if(x < 0.0) {
    result->val  = 2.0 - e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    result->val  = e_val;
    result->err  = e_err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}

/* end inlined source: specfunc/erfc.c */

/* begin inlined source: specfunc/gamma_inc.c (gsl_sf_gamma_inc_P, gsl_sf_gamma_inc_Q only) */
/* specfunc/gamma_inc.c
 *
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */


/* (stripped: #include directives — already handled above) */
/* (gsl_sf_gamma_inc_e omitted: requires gsl_sf_expint_E1_e not in this subset;
 * that function is not called by MARLEY) */

/* The dominant part,
 * D(a,x) := x^a e^(-x) / Gamma(a+1)
 */
static
int
gamma_inc_D(const double a, const double x, gsl_sf_result * result)
{
  if(a < 10.0) {
    double lnr;
    gsl_sf_result lg;
    gsl_sf_lngamma_e(a+1.0, &lg);
    lnr = a * log(x) - x - lg.val;
    result->val = exp(lnr);
    result->err = 2.0 * GSL_DBL_EPSILON * (fabs(lnr) + 1.0) * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    gsl_sf_result gstar;
    gsl_sf_result ln_term;
    double term1;
    if (x < 0.5*a) {
      double u = x/a;
      double ln_u = log(u);
      ln_term.val = ln_u - u + 1.0;
      ln_term.err = (fabs(ln_u) + fabs(u) + 1.0) * GSL_DBL_EPSILON;
    } else {
      double mu = (x-a)/a;
      gsl_sf_log_1plusx_mx_e(mu, &ln_term);  /* log(1+mu) - mu */
      /* Propagate cancellation error from x-a, since the absolute
         error of mu=x-a is DBL_EPSILON */
      ln_term.err += GSL_DBL_EPSILON * fabs(mu);
    };
    gsl_sf_gammastar_e(a, &gstar);
    term1 = exp(a*ln_term.val)/sqrt(2.0*M_PI*a);
    result->val  = term1/gstar.val;
    result->err  = 2.0 * GSL_DBL_EPSILON * (fabs(a*ln_term.val) + 1.0) * fabs(result->val);
    /* Include propagated error from log term */
    result->err += fabs(a) * ln_term.err * fabs(result->val);
    result->err += gstar.err/fabs(gstar.val) * fabs(result->val);
    return GSL_SUCCESS;
  }

}


/* P series representation.
 */
static
int
gamma_inc_P_series(const double a, const double x, gsl_sf_result * result)
{
  const int nmax = 10000;

  gsl_sf_result D;
  int stat_D = gamma_inc_D(a, x, &D);

  /* Approximating the terms of the series using Stirling's
     approximation gives t_n = (x/a)^n * exp(-n(n+1)/(2a)), so the
     convergence condition is n^2 / (2a) + (1-(x/a) + (1/2a)) n >>
     -log(GSL_DBL_EPS) if we want t_n < O(1e-16) t_0. The condition
     below detects cases where the minimum value of n is > 5000 */

  if (x > 0.995 * a && a > 1e5) { /* Difficult case: try continued fraction */
    gsl_sf_result cf_res;
    int status =  gsl_sf_exprel_n_CF_e(a, x, &cf_res);
    result->val = D.val * cf_res.val;
    result->err = fabs(D.val * cf_res.err) + fabs(D.err * cf_res.val);
    return status;
  }

  /* Series would require excessive number of terms */

  if (x > (a + nmax)) {
    GSL_ERROR ("gamma_inc_P_series x>>a exceeds range", GSL_EMAXITER);
  }

  /* Normal case: sum the series */

  {
    double sum  = 1.0;
    double term = 1.0;
    double remainder;
    int n;

    /* Handle lower part of the series where t_n is increasing, |x| > a+n */

    int nlow = (x > a) ? (x - a): 0;

    for(n=1; n < nlow; n++) {
      term *= x/(a+n);
      sum  += term;
    }

    /* Handle upper part of the series where t_n is decreasing, |x| < a+n */

    for (/* n = previous n */ ; n<nmax; n++)  {
      term *= x/(a+n);
      sum  += term;
      if(fabs(term/sum) < GSL_DBL_EPSILON) break;
    }

    /*  Estimate remainder of series ~ t_(n+1)/(1-x/(a+n+1)) */
    {
      double tnp1 = (x/(a+n)) * term;
      remainder =  tnp1 / (1.0 - x/(a + n + 1.0));
    }

    result->val  = D.val * sum;
    result->err  = D.err * fabs(sum) + fabs(D.val * remainder);
    result->err += (1.0 + n) * GSL_DBL_EPSILON * fabs(result->val);

    if(n == nmax && fabs(remainder/sum) > GSL_SQRT_DBL_EPSILON)
      GSL_ERROR ("gamma_inc_P_series failed to converge", GSL_EMAXITER);
    else
      return stat_D;
  }
}


/* Q large x asymptotic
 */
static
int
gamma_inc_Q_large_x(const double a, const double x, gsl_sf_result * result)
{
  const int nmax = 5000;

  gsl_sf_result D;
  const int stat_D = gamma_inc_D(a, x, &D);

  double sum  = 1.0;
  double term = 1.0;
  double last = 1.0;
  int n;
  for(n=1; n<nmax; n++) {
    term *= (a-n)/x;
    if(fabs(term/last) > 1.0) break;
    if(fabs(term/sum)  < GSL_DBL_EPSILON) break;
    sum  += term;
    last  = term;
  }

  result->val  = D.val * (a/x) * sum;
  result->err  = D.err * fabs((a/x) * sum);
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error in large x asymptotic", GSL_EMAXITER);
  else
    return stat_D;
}


/* Uniform asymptotic for x near a, a and x large.
 * See [Temme, p. 285]
 */
static
int
gamma_inc_Q_asymp_unif(const double a, const double x, gsl_sf_result * result)
{
  const double rta = sqrt(a);
  const double eps = (x-a)/a;

  gsl_sf_result ln_term;
  const int stat_ln = gsl_sf_log_1plusx_mx_e(eps, &ln_term);  /* log(1+eps) - eps */
  const double eta  = GSL_SIGN(eps) * sqrt(-2.0*ln_term.val);

  gsl_sf_result erfc;

  double R;
  double c0, c1;

  /* This used to say erfc(eta*M_SQRT2*rta), which is wrong.
   * The sqrt(2) is in the denominator. Oops.
   * Fixed: [GJ] Mon Nov 15 13:25:32 MST 2004
   */
  gsl_sf_erfc_e(eta*rta/M_SQRT2, &erfc);

  if(fabs(eps) < GSL_ROOT5_DBL_EPSILON) {
    c0 = -1.0/3.0 + eps*(1.0/12.0 - eps*(23.0/540.0 - eps*(353.0/12960.0 - eps*589.0/30240.0)));
    c1 = -1.0/540.0 - eps/288.0;
  }
  else {
    const double rt_term = sqrt(-2.0 * ln_term.val/(eps*eps));
    const double lam = x/a;
    c0 = (1.0 - 1.0/rt_term)/eps;
    c1 = -(eta*eta*eta * (lam*lam + 10.0*lam + 1.0) - 12.0 * eps*eps*eps) / (12.0 * eta*eta*eta*eps*eps*eps);
  }

  R = exp(-0.5*a*eta*eta)/(M_SQRT2*M_SQRTPI*rta) * (c0 + c1/a);

  result->val  = 0.5 * erfc.val + R;
  result->err  = GSL_DBL_EPSILON * fabs(R * 0.5 * a*eta*eta) + 0.5 * erfc.err;
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);

  return stat_ln;
}


/* Continued fraction which occurs in evaluation
 * of Q(a,x) or Gamma(a,x).
 *
 *              1   (1-a)/x  1/x  (2-a)/x   2/x  (3-a)/x
 *   F(a,x) =  ---- ------- ----- -------- ----- -------- ...
 *             1 +   1 +     1 +   1 +      1 +   1 +
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no).
 *
 * Split out from gamma_inc_Q_CF() by GJ [Tue Apr  1 13:16:41 MST 2003].
 * See gamma_inc_Q_CF() below.
 *
 */
static int
gamma_inc_F_CF(const double a, const double x, gsl_sf_result * result)
{
  const int    nmax  =  5000;
  const double small =  gsl_pow_3 (GSL_DBL_EPSILON);

  double hn = 1.0;           /* convergent */
  double Cn = 1.0 / small;
  double Dn = 1.0;
  int n;

  /* n == 1 has a_1, b_1, b_0 independent of a,x,
     so that has been done by hand                */
  for ( n = 2 ; n < nmax ; n++ )
  {
    double an;
    double delta;

    if(GSL_IS_ODD(n))
      an = 0.5*(n-1)/x;
    else
      an = (0.5*n-a)/x;

    Dn = 1.0 + an * Dn;
    if ( fabs(Dn) < small )
      Dn = small;
    Cn = 1.0 + an/Cn;
    if ( fabs(Cn) < small )
      Cn = small;
    Dn = 1.0 / Dn;
    delta = Cn * Dn;
    hn *= delta;
    if(fabs(delta-1.0) < GSL_DBL_EPSILON) break;
  }

  result->val = hn;
  result->err = 2.0*GSL_DBL_EPSILON * fabs(hn);
  result->err += GSL_DBL_EPSILON * (2.0 + 0.5*n) * fabs(result->val);

  if(n == nmax)
    GSL_ERROR ("error in CF for F(a,x)", GSL_EMAXITER);
  else
    return GSL_SUCCESS;
}


/* Continued fraction for Q.
 *
 * Q(a,x) = D(a,x) a/x F(a,x)
 *
 * Hans E. Plesser, 2002-01-22 (hans dot plesser at itf dot nlh dot no):
 *
 * Since the Gautschi equivalent series method for CF evaluation may lead
 * to singularities, I have replaced it with the modified Lentz algorithm
 * given in
 *
 * I J Thompson and A R Barnett
 * Coulomb and Bessel Functions of Complex Arguments and Order
 * J Computational Physics 64:490-509 (1986)
 *
 * In consequence, gamma_inc_Q_CF_protected() is now obsolete and has been
 * removed.
 *
 * Identification of terms between the above equation for F(a, x) and
 * the first equation in the appendix of Thompson&Barnett is as follows:
 *
 *    b_0 = 0, b_n = 1 for all n > 0
 *
 *    a_1 = 1
 *    a_n = (n/2-a)/x    for n even
 *    a_n = (n-1)/(2x)   for n odd
 *
 */
static
int
gamma_inc_Q_CF(const double a, const double x, gsl_sf_result * result)
{
  gsl_sf_result D;
  gsl_sf_result F;
  const int stat_D = gamma_inc_D(a, x, &D);
  const int stat_F = gamma_inc_F_CF(a, x, &F);

  result->val  = D.val * (a/x) * F.val;
  result->err  = D.err * fabs((a/x) * F.val) + fabs(D.val * a/x * F.err);

  return GSL_ERROR_SELECT_2(stat_F, stat_D);
}


/* Useful for small a and x. Handles the subtraction analytically.
 */
static
int
gamma_inc_Q_series(const double a, const double x, gsl_sf_result * result)
{
  double term1;  /* 1 - x^a/Gamma(a+1) */
  double sum;    /* 1 + (a+1)/(a+2)(-x)/2! + (a+1)/(a+3)(-x)^2/3! + ... */
  int stat_sum;
  double term2;  /* a temporary variable used at the end */

  {
    /* Evaluate series for 1 - x^a/Gamma(a+1), small a
     */
    const double pg21 = -2.404113806319188570799476;  /* PolyGamma[2,1] */
    const double lnx  = log(x);
    const double el   = M_EULER+lnx;
    const double c1 = -el;
    const double c2 = M_PI*M_PI/12.0 - 0.5*el*el;
    const double c3 = el*(M_PI*M_PI/12.0 - el*el/6.0) + pg21/6.0;
    const double c4 = -0.04166666666666666667
                       * (-1.758243446661483480 + lnx)
                       * (-0.764428657272716373 + lnx)
                       * ( 0.723980571623507657 + lnx)
                       * ( 4.107554191916823640 + lnx);
    const double c5 = -0.0083333333333333333
                       * (-2.06563396085715900 + lnx)
                       * (-1.28459889470864700 + lnx)
                       * (-0.27583535756454143 + lnx)
                       * ( 1.33677371336239618 + lnx)
                       * ( 5.17537282427561550 + lnx);
    const double c6 = -0.0013888888888888889
                       * (-2.30814336454783200 + lnx)
                       * (-1.65846557706987300 + lnx)
                       * (-0.88768082560020400 + lnx)
                       * ( 0.17043847751371778 + lnx)
                       * ( 1.92135970115863890 + lnx)
                       * ( 6.22578557795474900 + lnx);
    const double c7 = -0.00019841269841269841
                       * (-2.5078657901291800 + lnx)
                       * (-1.9478900888958200 + lnx)
                       * (-1.3194837322612730 + lnx)
                       * (-0.5281322700249279 + lnx)
                       * ( 0.5913834939078759 + lnx)
                       * ( 2.4876819633378140 + lnx)
                       * ( 7.2648160783762400 + lnx);
    const double c8 = -0.00002480158730158730
                       * (-2.677341544966400 + lnx)
                       * (-2.182810448271700 + lnx)
                       * (-1.649350342277400 + lnx)
                       * (-1.014099048290790 + lnx)
                       * (-0.191366955370652 + lnx)
                       * ( 0.995403817918724 + lnx)
                       * ( 3.041323283529310 + lnx)
                       * ( 8.295966556941250 + lnx);
    const double c9 = -2.75573192239859e-6
                       * (-2.8243487670469080 + lnx)
                       * (-2.3798494322701120 + lnx)
                       * (-1.9143674728689960 + lnx)
                       * (-1.3814529102920370 + lnx)
                       * (-0.7294312810261694 + lnx)
                       * ( 0.1299079285269565 + lnx)
                       * ( 1.3873333251885240 + lnx)
                       * ( 3.5857258865210760 + lnx)
                       * ( 9.3214237073814600 + lnx);
    const double c10 = -2.75573192239859e-7
                       * (-2.9540329644556910 + lnx)
                       * (-2.5491366926991850 + lnx)
                       * (-2.1348279229279880 + lnx)
                       * (-1.6741881076349450 + lnx)
                       * (-1.1325949616098420 + lnx)
                       * (-0.4590034650618494 + lnx)
                       * ( 0.4399352987435699 + lnx)
                       * ( 1.7702236517651670 + lnx)
                       * ( 4.1231539047474080 + lnx)
                       * ( 10.342627908148680 + lnx);

    term1 = a*(c1+a*(c2+a*(c3+a*(c4+a*(c5+a*(c6+a*(c7+a*(c8+a*(c9+a*c10)))))))));
  }

  {
    /* Evaluate the sum.
     */
    const int nmax = 5000;
    double t = 1.0;
    int n;
    sum = 1.0;

    for(n=1; n<nmax; n++) {
      t *= -x/(n+1.0);
      sum += (a+1.0)/(a+n+1.0)*t;
      if(fabs(t/sum) < GSL_DBL_EPSILON) break;
    }

    if(n == nmax)
      stat_sum = GSL_EMAXITER;
    else
      stat_sum = GSL_SUCCESS;
  }

  term2 = (1.0 - term1) * a/(a+1.0) * x * sum;
  result->val  = term1 + term2;
  result->err  = GSL_DBL_EPSILON * (fabs(term1) + 2.0*fabs(term2));
  result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
  return stat_sum;
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_gamma_inc_Q_e(const double a, const double x, gsl_sf_result * result)
{
  if(a < 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(a == 0.0)
  {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x <= 0.5*a) {
    /* If the series is quick, do that. It is
     * robust and simple.
     */
    gsl_sf_result P;
    int stat_P = gamma_inc_P_series(a, x, &P);
    result->val  = 1.0 - P.val;
    result->err  = P.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_P;
  }
  else if(a >= 1.0e+06 && (x-a)*(x-a) < a) {
    /* Then try the difficult asymptotic regime.
     * This is the only way to do this region.
     */
    return gamma_inc_Q_asymp_unif(a, x, result);
  }
  else if(a < 0.2 && x < 5.0) {
    /* Cancellations at small a must be handled
     * analytically; x should not be too big
     * either since the series terms grow
     * with x and log(x).
     */
    return gamma_inc_Q_series(a, x, result);
  }
  else if(a <= x) {
    if(x <= 1.0e+06) {
      /* Continued fraction is excellent for x >~ a.
       * We do not let x be too large when x > a since
       * it is somewhat pointless to try this there;
       * the function is rapidly decreasing for
       * x large and x > a, and it will just
       * underflow in that region anyway. We
       * catch that case in the standard
       * large-x method.
       */
      return gamma_inc_Q_CF(a, x, result);
    }
    else {
      return gamma_inc_Q_large_x(a, x, result);
    }
  }
  else {
    if(x > a - sqrt(a)) {
      /* Continued fraction again. The convergence
       * is a little slower here, but that is fine.
       * We have to trade that off against the slow
       * convergence of the series, which is the
       * only other option.
       */
      return gamma_inc_Q_CF(a, x, result);
    }
    else {
      gsl_sf_result P;
      int stat_P = gamma_inc_P_series(a, x, &P);
      result->val  = 1.0 - P.val;
      result->err  = P.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_P;
    }
  }
}


int
gsl_sf_gamma_inc_P_e(const double a, const double x, gsl_sf_result * result)
{
  if(a <= 0.0 || x < 0.0) {
    DOMAIN_ERROR(result);
  }
  else if(x == 0.0) {
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x < 20.0 || x < 0.5*a) {
    /* Do the easy series cases. Robust and quick.
     */
    return gamma_inc_P_series(a, x, result);
  }
  else if(a > 1.0e+06 && (x-a)*(x-a) < a) {
    /* Crossover region. Note that Q and P are
     * roughly the same order of magnitude here,
     * so the subtraction is stable.
     */
    gsl_sf_result Q;
    int stat_Q = gamma_inc_Q_asymp_unif(a, x, &Q);
    result->val  = 1.0 - Q.val;
    result->err  = Q.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_Q;
  }
  else if(a <= x) {
    /* Q <~ P in this area, so the
     * subtractions are stable.
     */
    gsl_sf_result Q;
    int stat_Q;
    if(a > 0.2*x) {
      stat_Q = gamma_inc_Q_CF(a, x, &Q);
    }
    else {
      stat_Q = gamma_inc_Q_large_x(a, x, &Q);
    }
    result->val  = 1.0 - Q.val;
    result->err  = Q.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_Q;
  }
  else {
    if((x-a)*(x-a) < a) {
      /* This condition is meant to insure
       * that Q is not very close to 1,
       * so the subtraction is stable.
       */
      gsl_sf_result Q;
      int stat_Q = gamma_inc_Q_CF(a, x, &Q);
      result->val  = 1.0 - Q.val;
      result->err  = Q.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return stat_Q;
    }
    else {
      return gamma_inc_P_series(a, x, result);
    }
  }
}

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

/* (eval.h macros already defined above) */

double gsl_sf_gamma_inc_P(const double a, const double x)
{
  EVAL_RESULT(gsl_sf_gamma_inc_P_e(a, x, &result));
}

double gsl_sf_gamma_inc_Q(const double a, const double x)
{
  EVAL_RESULT(gsl_sf_gamma_inc_Q_e(a, x, &result));
}

/* end inlined source: specfunc/gamma_inc.c */

/* begin inlined source: cdf/gamma.c (gsl_cdf_gamma_Q only) */
/* cdf/cdf_gamma.c
 *
 * Copyright (C) 2003 Jason Stover.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author: J. Stover */

/* (stripped: #include directives — already handled above) */
/* (gsl_cdf_gamma_P omitted: not needed by MARLEY) */

double
gsl_cdf_gamma_Q (const double x, const double a, const double b)
{
  double Q;
  double y;

  /* Upstream GSL can return boundary values in some invalid-parameter and/or
   * x <= 0 cases. The reduced MARLEY fallback subset intentionally tightens
   * domain handling for unsupported, non-physical parameter domains so that
   * unsupported usage fails loudly via the GSL C error path.
   */
  if (a <= 0.0 || b <= 0.0)
    {
      GSL_ERROR_VAL ("unsupported parameter domain in built-in gsl_cdf_gamma_Q"
        " (require a > 0 and b > 0)", GSL_EDOM, GSL_NAN);
    }

  y = x / b;

  if (x <= 0.0)
    {
      return 1.0;
    }

  if (y < a)
    {
      Q = 1 - gsl_sf_gamma_inc_P (a, y);
    }
  else
    {
      Q = gsl_sf_gamma_inc_Q (a, y);
    }

  return Q;
}
/* end inlined source: cdf/gamma.c */

/* begin inlined source: cdf/chisq.c (gsl_cdf_chisq_Q only) */
/* cdf/cdf_chisq.c
 *
 * Copyright (C) 2003, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* (stripped: #include directives — already handled above) */
/* (gsl_cdf_chisq_P omitted: not needed by MARLEY) */

double
gsl_cdf_chisq_Q (const double x, const double nu)
{
  /* Keep parameter-domain behavior aligned with the built-in gamma-Q guard:
   * reject unsupported/non-physical nu <= 0 through the GSL error path.
   */
  if (nu <= 0.0)
    {
      GSL_ERROR_VAL ("unsupported parameter domain in built-in gsl_cdf_chisq_Q"
        " (require nu > 0)", GSL_EDOM, GSL_NAN);
    }

  return gsl_cdf_gamma_Q (x, nu / 2, 2.0);
}
/* end inlined source: cdf/chisq.c */

#endif // MARLEY_FOUND_GSL
