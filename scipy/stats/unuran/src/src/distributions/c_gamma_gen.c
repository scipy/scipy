/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_gamma_gen.c                                                *
 *                                                                           *
 *   Special generators for Gamma distribution                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <methods/cstd.h>
#include <methods/cstd_struct.h>
#include <methods/x_gen_source.h>
#include <distr/distr_source.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
#include "unur_distributions.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int gamma_gll_init( struct unur_gen *gen );
inline static int gamma_gs_init( struct unur_gen *gen );
inline static int gamma_gd_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define MAX_gen_params  8      /* maximal number of parameters for generator */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define alpha (DISTR.params[0])   /* shape */
#define beta  (DISTR.params[1])   /* scale */
#define gamma (DISTR.params[2])   /* location */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_gamma_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for gamma distribution                  */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case 1:  /* Acceptance Rejection combined with Acceptance Complement */
    if (gen==NULL) return UNUR_SUCCESS; /* test existence only  */
    if (alpha < 1.) {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_gamma_gs );
      return gamma_gs_init( gen );
    }
    else {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_gamma_gd );
      return gamma_gd_init( gen );

    }

  case 2:  /* Rejection with log-logistic envelopes */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_gamma_gll );
    return gamma_gll_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }

} /* end of _unur_stdgen_gamma_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Gamma Distribution: Rejection from log-logistic envelopes                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               Gamma distribution with parameter a.                        *
 *                                                                           *
 * REFERENCE:  - R.C.H. Cheng (1977): The Generation of Gamma Variables      *
 *               with Non-Integral Shape Parameter,                          *
 *               Appl. Statist. 26(1), 71-75                                 *
 *                                                                           *
 * Implemented by Ralf Kremer                                                *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define aa  GEN->gen_param[0]
#define bb  GEN->gen_param[1]
#define cc  GEN->gen_param[2]
/*---------------------------------------------------------------------------*/

inline static int
gamma_gll_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  aa = (alpha > 1.0) ? sqrt(alpha + alpha - 1.0) : alpha;
  bb = alpha - 1.386294361;
  cc = alpha + aa;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of gamma_gll_init() */

double 
_unur_stdgen_sample_gamma_gll( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;
  double u1,u2,v,r,z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    u1 = uniform();
    u2 = uniform();
    v = log(u1 / (1.0 - u1)) / aa;
    X = alpha * exp(v);
    r = bb + cc * v - X;
    z = u1 * u1 * u2;
    if (r + 2.504077397 >= 4.5 * z) break;
    if (r >= log(z)) break;
  }
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==1) ? X : gamma + beta * X );

} /* end of _unur_stdgen_sample_gamma_gll() */

/*---------------------------------------------------------------------------*/
#undef aa
#undef bb
#undef cc
/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Gamma Distribution: Acceptance Rejection combined with                    *
 *                     Acceptance Complement                                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Gamma distribution with    *
 *               parameter alpha > 0.                                        *
 *               Acceptance Rejection  gs  for  alpha < 1,                   *
 *               Acceptance Complement gd  for  alpha >= 1.                  *
 *                                                                           *
 * REFERENCE:  - J.H. Ahrens, U. Dieter (1974): Computer methods for         *
 *               sampling from gamma, beta, Poisson and binomial             *
 *               distributions, Computing 12, 223-246.                       *
 *             - J.H. Ahrens, U. Dieter (1982): Generating gamma variates    *
 *               by a modified rejection technique,                          *
 *               Communications of the ACM 25, 47-54.                        *
 *                                                                           *
 * Implemented by G. Thallinger and E. Stadlober, August 1990                *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define b   GEN->gen_param[0]
/*---------------------------------------------------------------------------*/

inline static int
gamma_gs_init( struct unur_gen *gen )
     /* CASE alpha < 1: Acceptance rejection algorithm gs */
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  b = 1. + 0.36788794412 * alpha;       /* Step 1 */
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of gamma_gs_init() */

double 
_unur_stdgen_sample_gamma_gs( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X, p;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    p = b * uniform();
    if (p <= 1.) {                   /* Step 2. Case gds <= 1 */
      X = exp(log(p) / alpha);
      if ( log(uniform()) <= -X)
	break;
    }
    else {                           /* Step 3. Case gds > 1 */
      X = - log ((b - p) / alpha);
      if ( log(uniform()) <= ((alpha - 1.) * log(X)))
	break;
    }
  }
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==1) ? X : gamma + beta * X );

} /* end of _unur_stdgen_sample_gamma_gs() */

/*---------------------------------------------------------------------------*/
#undef b
/*---------------------------------------------------------------------------*/
#define ss   GEN->gen_param[0]
#define s    GEN->gen_param[1]
#define d    GEN->gen_param[2]
#define r    GEN->gen_param[3]
#define q0   GEN->gen_param[4]
#define b    GEN->gen_param[5]
#define c    GEN->gen_param[6]
#define si   GEN->gen_param[7]

#define q1   0.0416666664
#define q2   0.0208333723
#define q3   0.0079849875
#define q4   0.0015746717
#define q5  -0.0003349403
#define q6   0.0003340332
#define q7   0.0006053049
#define q8  -0.0004701849
#define q9   0.0001710320
#define a1   0.333333333
#define a2  -0.249999949
#define a3   0.199999867
#define a4  -0.166677482
#define a5   0.142873973
#define a6  -0.124385581
#define a7   0.110368310
#define a8  -0.112750886
#define a9   0.104089866
#define e1   1.000000000
#define e2   0.499999994
#define e3   0.166666848
#define e4   0.041664508
#define e5   0.008345522
#define e6   0.001353826
#define e7   0.000247453
/*---------------------------------------------------------------------------*/
#define NORMAL  gen->gen_aux   /* pointer to normal variate generator        */
/*---------------------------------------------------------------------------*/

inline static int
gamma_gd_init( struct unur_gen *gen )
     /* CASE alpha >= 1: Acceptance complement algorithm gd */
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  /* Step 1. Preparations */
  ss = alpha - 0.5;
  s = sqrt(ss);
  d = 5.656854249 - 12. * s;

  /* Step 4. Set-up for hat case */
  r = 1. / alpha;
  q0 = ((((((((q9 * r + q8) * r + q7) * r + q6) * r + q5) * r + q4) *
	  r + q3) * r + q2) * r + q1) * r;
  if (alpha > 3.686) {
    if (alpha > 13.022) {
      b = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }
    else {
      b = 1.654 + 0.0076 * ss;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    }
  }
  else {
    b = 0.463 + s - 0.178 * ss;
    si = 1.235;
    c = 0.195 / s - 0.079 + 0.016 * s;
  }

  /* make a normal variate generator (use default special generator) */
  if (NORMAL==NULL) {
    struct unur_distr *distr = unur_distr_normal(NULL,0);
    struct unur_par *par = unur_cstd_new( distr );
    NORMAL = (par) ? _unur_init(par) : NULL;
    _unur_check_NULL( NULL, NORMAL, UNUR_ERR_NULL );
    /* need same uniform random number generator as slash generator */
    NORMAL->urng = gen->urng;
    /* copy debugging flags */
    NORMAL->debug = gen->debug;
    /* we do not need the distribution object any more */
    _unur_distr_free( distr );
  }
  /* else we are in the re-init mode 
     --> there is no necessity to make the generator object again */

  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of gamma_gd_init() */

double 
_unur_stdgen_sample_gamma_gd( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double U,X,E;
  double q,sign_U,t,v,w,x;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  do {

    /* Step 2. Normal deviate */
    t = _unur_sample_cont(NORMAL);
    x = s + 0.5 * t;
    X = x * x;
    if (t >= 0.)
      break;         /* Immediate acceptance */

    /* Step 3. Uniform random number */
    U = uniform();
    if (d * U <= t * t * t) 
      break;         /* Squeeze acceptance */
    
    /* Step 5. Calculation of q */
    if (x > 0.) {
      /* Step 6. */
      v = t / (s + s);
      if (fabs(v) > 0.25)
	q = q0 - s * t + 0.25 * t * t + (ss + ss) * log(1. + v);
      else
	q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) * v + a6) *
				    v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
      /* Step 7. Quotient acceptance */
      if (log(1. - U) <= q)
	break;
    }
    
    /* Step 8. Double exponential deviate t */
    while (1) {
      do {
	E = -log( uniform() );
	U = uniform();
	U = U + U - 1.;
	sign_U = (U > 0) ? 1. : -1.;
	t = b + (E * si) * sign_U;
      } while (t <= -0.71874483771719);   /* Step 9. Rejection of t */
      
      /* Step 10. New q(t) */
      v = t / (s + s);
      if (fabs(v) > 0.25)
	q = q0 - s * t + 0.25 * t * t + (ss + ss) * log(1. + v);
      else
	q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) * v + a6) *
				    v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
      
      /* Step 11. */
      if (q <= 0.)
	continue; 
      
      if (q > 0.5)
	w = exp(q) - 1.;
      else
	w = ((((((e7 * q + e6) * q + e5) * q + e4) * q + e3) * q + e2) *
	   q + e1) * q;
      
      /* Step 12. Hat acceptance */
      if ( c * U * sign_U <= w * exp(E - 0.5 * t * t)) {
	x = s + 0.5 * t;
	X = x * x;
	break;
      }
    }

  } while (0);

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==1) ? X : gamma + beta * X );

} /* end of _unur_stdgen_sample_gamma_gd() */

/*---------------------------------------------------------------------------*/
#undef ss
#undef s 
#undef d 
#undef r 
#undef q0
#undef b 
#undef c 
#undef si

#undef q1
#undef q2
#undef q3
#undef q4
#undef q5
#undef q6
#undef q7
#undef q8
#undef q9
#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef a7
#undef a8
#undef a9
#undef e1
#undef e2
#undef e3
#undef e4
#undef e5
#undef e6
#undef e7
/*---------------------------------------------------------------------------*/
#undef NORMAL
/*---------------------------------------------------------------------------*/
