/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_beta_gen.c                                                 *
 *                                                                           *
 *   Special generators for Beta distribution                                *
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
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int beta_bc_init( struct unur_gen *gen );
inline static int beta_bb_init( struct unur_gen *gen );

inline static int beta_b00_init( struct unur_gen *gen );
inline static int beta_b01_init( struct unur_gen *gen );
inline static int beta_b1prs_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
 /* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params 22      /* maximal number of parameters for generator */

#define p     (DISTR.params[0])   /* shape parameter */
#define q     (DISTR.params[1])   /* shape parameter */
#define a     (DISTR.params[2])   /* left boundary */
#define b     (DISTR.params[3])   /* right boundary */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_beta_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for gamma distribution.                 */
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
  case 1:  /* Rejection with log-logistic envelopes */
    if (gen==NULL) return UNUR_SUCCESS; /* test existence only  */
    if (p>1. && q>1.) {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_bb );
      return beta_bb_init( gen );
    }
    else {
      _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_bc );
      return beta_bc_init( gen );
    }

  case 2:  /* Stratified Rejection/Patchwork Rejection */
    if (gen==NULL) return UNUR_SUCCESS; /* test existence only  */ 
    if (p>1.)
      if (q>1.) {    /* p > 1 && q > 1 */
	_unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b1prs );
	return beta_b1prs_init( gen );
      }
      else {         /* p > 1 && q <= 1 */
	_unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b01 );
	return beta_b01_init( gen );
      }
    else
      if (q>1.) {    /* p <= 1 && q > 1 */
	_unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b01 );
	return beta_b01_init( gen );
      }
      else {         /* p <= 1 && q <= 1 */
	_unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_beta_b00 );
	return beta_b00_init( gen );
      }

  default: /* no such generator */
    return UNUR_FAILURE;

  }

} /* end of _unur_stdgen_beta_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Beta Distribution: Acceptance/Rejection from log-logistic hats for the    *
 *                    beta prime distribution                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Beta distribution with     *
 *               parameters p,q (p > 0, q > 0).                              *
 *               It combines algorithms bb (p > 1, q > 1) and                *
 *               bc (p <= 1 or q <= 1).                                      *
 *                                                                           *
 * REFERENCE : - R.C.H. Cheng (1978): Generating beta variates with          *
 *               nonintegral shape parameters,                               *
 *               Communications of the ACM 21, 317-322.                      *
 *                                                                           *
 * ROUTINES:   - bb ... Beta generator for p > 1, q > 1                      *
 *             - bc ... Beta generator for p <= 1 or q <= 1                  *
 *                                                                           *
 * Implemented by E. Stadlober, R. Kremer, 1990                              *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define am      (GEN->gen_param[0])
#define bm      (GEN->gen_param[1])
#define al      (GEN->gen_param[2])
#define alnam   (GEN->gen_param[3])
#define be      (GEN->gen_param[4])
#define ga      (GEN->gen_param[5])
#define si      (GEN->gen_param[6])
#define rk1     (GEN->gen_param[7])
#define rk2     (GEN->gen_param[8])
/*---------------------------------------------------------------------------*/

inline static int
beta_bc_init( struct unur_gen *gen )
     /* p <= 1. || q <= 1. */ 
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  am = (p > q) ? p : q;
  bm = (p < q) ? p : q;
  al = am + bm;
  alnam = al * log(al/am) - 1.386294361;
  be = 1.0 / bm;
  si = 1.0 + am - bm;
  rk1 = si * (0.013888889 + 0.041666667 * bm) / (am * be - 0.77777778);
  rk2 = 0.25 + (0.5 + 0.25 / si) * bm;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of beta_bc_init() */

double 
_unur_stdgen_sample_beta_bc(  struct unur_gen *gen )
     /* p <= 1. || q <= 1. */ 
{
  /* -X- generator code -X- */
  double X;
  double u1,u2,v,w,y,z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    /* Step 1 */
    u1 = uniform();
    u2 = uniform();

    if (u1 < 0.5) {
      /* Step 2 */
      y = u1 * u2;
      z = u1 * y;

      if ((0.25 * u2 - y + z) >= rk1) 
	continue;  /* goto 1 */

      /* Step 5 */
      v = be * log(u1 / (1.0 - u1));
      if (v > 80.0) {
	if (alnam < log(z))
	  continue;
	X = (_unur_FP_same(am,p)) ? 1.0 : 0.0;
	break;
      }
      else {
	w = am * exp(v);
	if ((al * (log(al / (bm + w)) + v) - 1.386294361) < log(z))
	  continue;  /* goto 1 */

	/* Step 6_a */
	X = (!_unur_FP_same(am,p)) ? bm / (bm + w) : w / (bm + w);
	break;
      }
    }
    
    else {
      /* Step 3 */
      z = u1 * u1 * u2;
      if (z < 0.25) {
	/* Step 5 */
	v = be * log(u1 / (1.0 - u1));
	if (v > 80.0) {
	  X = (_unur_FP_same(am,p)) ? 1.0 : 0.0;
	  break;
	}

	w = am * exp(v);
	X = (!_unur_FP_same(am,p)) ? bm / (bm + w) : w / (bm + w);
	break;
      }
      else {
	if (z >= rk2)
	  continue;
	v = be * log(u1 / (1.0 - u1));
	if ( v > 80.0) {
	  if (alnam < log(z))
	    continue;
	  X = (_unur_FP_same(am,p)) ? 1.0 : 0.0;
	  break;
	}
	w = am * exp(v);
	if ((al * (log(al / (bm + w)) + v) - 1.386294361) < log(z))
	  continue;  /* goto 1 */

	/* Step 6_b */
	X = (!_unur_FP_same(am,p))? bm / (bm + w) : w / (bm + w);
	break;
      }
    }
  }

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==2) ? X : a + (b-a) * X);

} /* end of _unur_stdgen_sample_beta_bc() */

/*---------------------------------------------------------------------------*/

inline static int
beta_bb_init( struct unur_gen *gen )
     /* p > 1. && q > 1 */ 
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  am = (p < q) ? p : q;
  bm = (p > q) ? p : q;
  al = am + bm;
  be = sqrt((al - 2.0)/(2.0 * p * q - al));
  ga = am + 1.0 / be;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of beta_bb_init() */

double 
_unur_stdgen_sample_beta_bb(  struct unur_gen *gen )
     /* p > 1. && q > 1 */ 
{
  /* -X- generator code -X- */
  double X;
  double u1,u2,v,w,z,r,s,t;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    /* Step 1 */
    u1 = uniform();
    u2 = uniform();
    v = be * log(u1 / (1.0 - u1));
    w = am * exp(v);
    z = u1 * u1 * u2;
    r = ga * v - 1.386294361;
    s = am + r - w;

    /* Step 2 */
    if (s + 2.609437912 < 5.0 * z) {
      /* Step 3 */
      t = log(z);
      if (s < t)
	/* Step 4 */
	if (r + al * log(al/(bm + w)) < t) 
	  continue;
    }

    /* Step 5 */
    X = (_unur_FP_same(am,p)) ? w / (bm + w) : bm / (bm + w);
    break;
  }
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==2) ? X : a + (b-a) * X);

} /* end of _unur_stdgen_sample_beta_bb() */

/*---------------------------------------------------------------------------*/
#undef am
#undef bm
#undef al
#undef alnam
#undef be
#undef ga
#undef si
#undef rk1
#undef rk2
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Beta Distribution: Stratified Rejection/Patchwork Rejection               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Beta distribution with     *
 *               parameters p,q (p > 0, q > 0).                              *
 *                                                                           *
 * REFERENCES: - H. Sakasegawa (1983): Stratified rejection and squeeze      *
 *               method for generating beta random numbers,                  *
 *               Ann. Inst. Statist. Math. 35 B, 291-302.                    *
 *             - H. Zechner, E. Stadlober (1993): Generating beta variates   *
 *               via patchwork rejection, Computing 50, 1-18.                *
 *                                                                           *
 * ROUTINES:   - b00   ... Beta generator for p < 1, q < 1                   *
 *             - b01   ... Beta generator for p < 1 < q or p > 1 > q         *
 *             - b1prs ... Beta generator for p > 1, q > 1                   *
 *                                                                           *
 * Implemented by H. Zechner and F. Niederl, July 1994                       *
 *****************************************************************************
 *                                                                           *
 *   For parameters p < 1 , q < 1  and  p < 1 < q   or  q < 1 < p the        *
 *   stratified rejection methods b00 and b01 of Sakasegawa are used. Both   *
 *   procedures employ suitable two-part power functions from which samples  *
 *   can be obtained by inversion.                                           *
 *   If  p > 1, q > 1 (unimodal case) the patchwork rejection method b1prs   *
 *   of Zechner/Stadlober is utilized: The area below the density function   *
 *   f(x) in its body is rearranged by certain point reflections. Within a   *
 *   large center interval variates are sampled efficiently by rejection     *
 *   from uniform hats. Rectangular immediate acceptance regions speed up    *
 *   the generation. The remaining tails are covered by exponential          *
 *   functions.                                                              *
 *   If (p-1)(q-1) = 0 sampling is done by inversion if either p or q are    *
 *   not equal to one. If  p = q = 1  a uniform random variate is delivered. *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define p_      (GEN->gen_param[0])
#define q_      (GEN->gen_param[1])
#define c       (GEN->gen_param[2])
#define t       (GEN->gen_param[3])
#define fp      (GEN->gen_param[4])
#define fq      (GEN->gen_param[5])
#define p1      (GEN->gen_param[6])
#define p2      (GEN->gen_param[7])
/*---------------------------------------------------------------------------*/

inline static int
beta_b00_init( struct unur_gen *gen )
     /* p < 1. && q < 1 */ 
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  p_ = p - 1.;
  q_ = q - 1.;
  c = (q * q_) / (p * p_);                            /* q(1-q) / p(1-p) */
  t = _unur_FP_same(c,1.) ? 0.5 : (1. - sqrt(c))/(1. - c);   /* t = t_opt       */
  fp = exp(p_ * log(t));
  fq = exp(q_ * log(1. - t));                        /* f(t) = fa * fb  */
  
  p1 = t/p;                                           /* 0 < X < t       */
  p2 = (1. - t)/q + p1;                              /* t < X < 1       */
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of beta_b00_init() */

double 
_unur_stdgen_sample_beta_b00(  struct unur_gen *gen )
     /* p < 1. && q < 1 */
{
  /* -X- generator code -X- */
  double U, V, X, Z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    U = uniform() * p2;
    if (U <= p1) {  /*  X < t                                                */
      Z = exp(log(U/p1) / p);
      X = t * Z;
      /* squeeze accept:   L(x) = 1 + (1 - q)x                               */
      V = uniform() * fq;
      if (V <= 1. - q_ * X)
	break;
      /* squeeze reject:   U(x) = 1 + ((1 - t)^(q-1) - 1)/t * x              */
      if (V <= 1. + (fq - 1.) * Z) {
	/* quotient accept:  quot(x) = (1 - x)^(q-1) / fq                    */
	if (log(V) <= q_ * log(1. - X))
	  break;
      }
    }
    else {          /*  X > t  */
      Z = exp( log( (U-p1)/(p2-p1) ) / q);
      X = 1. - (1. - t)*Z;
      /* squeeze accept:   L(x) = 1 + (1 - p)(1 - x)                         */
      V = uniform() * fp;
      if (V <= 1.0 - p_*(1. - X))
	break;
      /* squeeze reject:   U(x) = 1 + (t^(p-1) - 1)/(1 - t) * (1 - x)        */
      if (V <= 1.0 + (fp - 1.) * Z) {
	/* quotient accept:  quot(x) = x^(p-1) / fp                          */
	if (log(V) <= p_ * log(X))  
	  break;
      }
    }
  }
  /* -X- end of generator code -X- */
  
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);

} /* end of _unur_stdgen_sample_beta_b00() */

/*---------------------------------------------------------------------------*/
#undef p_
#undef q_
#undef c
#undef t
#undef fp
#undef fq
#undef p1
#undef p2
/*---------------------------------------------------------------------------*/
#define pint    (GEN->gen_param[0])
#define qint    (GEN->gen_param[1])
#define p_      (GEN->gen_param[2])
#define q_      (GEN->gen_param[3])
#define t       (GEN->gen_param[4])
#define fp      (GEN->gen_param[5])
#define fq      (GEN->gen_param[6])
#define ml      (GEN->gen_param[7])
#define mu      (GEN->gen_param[8])
#define p1      (GEN->gen_param[9])
#define p2      (GEN->gen_param[10])
/*---------------------------------------------------------------------------*/

inline static int
beta_b01_init( struct unur_gen *gen )
     /* p < 1. < q || p > 1. > q */ 
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  /* internal use of p and q */
  if (p>q) {
    /* swap p and q */
    pint = q;
    qint = p;
  }
  else {
    pint = p;
    qint = q;
  }
  
  p_ = pint - 1.;
  q_ = qint - 1.;
  t = p_/(pint - qint);                   /* one step Newton * start value t   */
  fq = exp((q_ - 1.) * log(1. - t));
  fp = pint - (pint + q_) * t;
  t -= (t - (1. - fp) * (1. - t) * fq / qint) / (1. - fp*fq);
  fp = exp(p_ * log(t));
  fq = exp(q_ * log(1. - t));                         /* f(t) = fa * fb  */
  if (q_ <= 1.0) {
    ml = (1. - fq) / t;                               /*   ml = -m1      */
    mu = q_ * t;                                      /*   mu = -m2 * t  */
  }
  else {
    ml = q_;
    mu = 1. - fq;
  }
  p1 = t/pint;                                           /*  0 < X < t      */
  p2 = fq * (1. - t)/qint + p1;                          /*  t < X < 1      */
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;
  
} /* end of beta_b01_init() */

double 
_unur_stdgen_sample_beta_b01(  struct unur_gen *gen )
     /* p < 1. < q || p > 1. > q */ 
{
  /* -X- generator code -X- */
  double U, V, X, Z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    U = uniform() * p2;
    if (U <= p1) {    /*  X < t                                              */
      Z = exp( log(U/p1) / pint);
      X = t * Z;
      /* squeeze accept:   L(x) = 1 + m1*x,  ml = -m1                        */
      V = uniform();
      if (V <= 1. - ml * X)
	break;
      /* squeeze reject:   U(x) = 1 + m2*x,  mu = -m2 * t                    */
      if (V <= 1. - mu * Z)
	/* quotient accept:  quot(x) = (1 - x)^(q-1)                         */
	if (log(V) <= q_ * log(1. - X))
	  break;
    }
    else {             /*  X > t                                             */
      Z = exp( log((U-p1)/(p2-p1)) / qint);
      X = 1. - (1. - t) * Z;
      /* squeeze accept:   L(x) = 1 + (1 - p)(1 - x)                         */
      V = uniform() * fp;
      if (V <= 1. - p_ * (1. - X))
	break;
      /* squeeze reject:   U(x) = 1 + (t^(p-1) - 1)/(1 - t) * (1 - x)        */
      if (V <= 1. + (fp - 1.) * Z)
	/* quotient accept:  quot(x) = (x)^(p-1) / fp                        */
	if (log(V) <= p_ * log(X))
	  break;
    }
  }
  if (p>q)
    /* p and q has been swapped */
    X = 1. - X;
  /* -X- end of generator code -X- */
  
  return ((DISTR.n_params==2) ? X : a + (b-a) * X);

} /* end of _unur_stdgen_sample_beta_b01() */

/*---------------------------------------------------------------------------*/
#undef pint
#undef qint
#undef p_
#undef q_
#undef t 
#undef fp
#undef fq
#undef ml
#undef mu
#undef p1
#undef p2
/*---------------------------------------------------------------------------*/
#define p_      (GEN->gen_param[0])
#define q_      (GEN->gen_param[1])
#define s       (GEN->gen_param[2])
#define m       (GEN->gen_param[3])
#define D       (GEN->gen_param[4])
#define Dl      (GEN->gen_param[5])
#define x1      (GEN->gen_param[6])
#define x2      (GEN->gen_param[7])
#define x4      (GEN->gen_param[8])
#define x5      (GEN->gen_param[9])
#define f1      (GEN->gen_param[10])
#define f2      (GEN->gen_param[11])
#define f4      (GEN->gen_param[12])
#define f5      (GEN->gen_param[13])
#define ll      (GEN->gen_param[14])
#define lr      (GEN->gen_param[15])
#define z2      (GEN->gen_param[16])
#define z4      (GEN->gen_param[17])
#define p1      (GEN->gen_param[18])
#define p2      (GEN->gen_param[19])
#define p3      (GEN->gen_param[20])
#define p4      (GEN->gen_param[21])
/*---------------------------------------------------------------------------*/

inline static int
beta_b1prs_init( struct unur_gen *gen )
     /* p > 1. && q > 1. */ 
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  p_ = p - 1.0;
  q_ = q - 1.0;
  s = p_ + q_;
  m = p_ / s;

  if (p_ > 1. || q_ > 1.)
    D = sqrt(m * (1. - m) / (s - 1.));

  if (p_ <= 1.) {
    x2 = Dl = m * 0.5;
    x1 = z2 = f1 = ll = 0.;
  }
  else {
    x2 = m - D;
    x1 = x2 - D;
    z2 = x2 * (1. - (1. - x2)/(s * D));
    if (x1 <= 0. || (s - 6.) * x2 - p_ + 3. > 0.) {
      x1 = z2;  x2 = (x1 + m) * 0.5;
      Dl = m - x2;
    }
    else {
      Dl = D;
    }
    f1 = exp( p_ * log(x1/m) + q_ * log((1. - x1)/(1. - m)) );
    ll = x1 * (1.0 - x1) / (s * (m - x1));            /* z1 = x1 - ll   */
  }
  f2 = exp( p_ * log(x2/m) + q_ * log((1. - x2)/(1. - m)) );

  if (q_ <= 1.) {
    D = (1. - m) * 0.5;
    x4 = 1. - D;
    x5 = z4 = 1.;
    f5 = lr = 0.;
  }
  else {
    x4 = m + D;
    x5 = x4 + D;
    z4 = x4 * (1. + (1. - x4)/(s * D));
    if (x5 >= 1. || (s - 6.) * x4 - p_ + 3. < 0.) {
      x5 = z4;
      x4 = (m + x5) * 0.5;
      D = x4 - m;
    }
    f5 = exp( p_ * log(x5/m) + q_ * log((1. - x5)/(1. - m)) );
    lr = x5 * (1. - x5) / (s * (x5 - m));            /* z5 = x5 + lr   */
  }
  f4 = exp( p_ * log(x4/m) + q_ * log((1. - x4)/(1. - m)) );

  p1 = f2 * (Dl + Dl);                                /*  x1 < X < m    */
  p2 = f4 * (D  + D) + p1;                            /*  m  < X < x5   */
  p3 = f1 * ll       + p2;                            /*       X < x1   */
  p4 = f5 * lr       + p3;                            /*  x5 < X        */
  /* -X- end of setup code -X- */
  
  return UNUR_SUCCESS;

} /* end of beta_b1prs_init() */

double 
_unur_stdgen_sample_beta_b1prs(  struct unur_gen *gen )
     /* p > 1. && q > 1. */ 
{
  /* -X- generator code -X- */
  double U, V, W, X, Y;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    U = uniform() * p4;
    if (U <= p1) {
      /* immediate accept:  x2 < X < m, - f(x2) < W < 0                      */
      W = U/Dl - f2;
      if (W <= 0.) {
	X = m - U/f2;
	break;
      }
      /* immediate accept:  x1 < X < x2, 0 < W < f(x1)                       */
      if (W <= f1) {
	X = x2 - W/f1 * Dl;
	break;
      }
      /* candidates for acceptance-rejection-test                            */
      U = uniform();
      V = Dl * U;
      X = x2 - V;
      Y = x2 + V;
      /* squeeze accept:    L(x) = f(x2) (x - z2) / (x2 - z2)                */
      if (W * (x2 - z2) <= f2 * (X - z2))
	break;
      V = f2 + f2 - W;
      if (V < 1.) {
	/* squeeze accept:    L(x) = f(x2) + (1 - f(x2))(x - x2)/(m - x2)    */
	if (V <= f2 + (1. - f2) * U) {
	  X = Y;
	  break;
	}
	/* quotient accept:   x2 < Y < m,   W >= 2f2 - f(Y)                  */
	if (V <= exp( p_ * log(Y/m) + q_ * log((1. - Y)/(1. - m)) ) ) {
	  X = Y;
	  break;
	}
      }
    }
    else 
      if (U <= p2) {
	U -= p1;
	/* immediate accept:  m < X < x4, - f(x4) < W < 0                    */
	W = U/D - f4;
	if (W <= 0.) {
	  X = m + U/f4;
	  break;
	}
	/* immediate accept:  x4 < X < x5, 0 < W < f(x5)                     */
	if (W <= f5) {
	  X = x4 + W/f5 * D;
	  break;
	}
	/* candidates for acceptance-rejection-test                          */
	U = uniform();
	V = D * U;
	X = x4 + V;
	Y = x4 - V;
	/* squeeze accept:    L(x) = f(x4) (z4 - x) / (z4 - x4)              */
	if (W * (z4 - x4) <= f4 * (z4 - X))
	  break;
	V = f4 + f4 - W;
	if (V < 1.) {
	  /* squeeze accept:    L(x) = f(x4) + (1 - f(x4))(x4 - x)/(x4 - m)  */
	  if (V <= f4 + (1.0 - f4) * U) {
	    X = Y;
	    break;
	  }
	  /* quotient accept:   m < Y < x4,   W >= 2f4 - f(Y)                */
	  if (V <= exp( p_ * log(Y/m) + q_ * log((1. - Y)/(1. - m)) ) ) {
	    X = Y;
	    break;
	  }
	}
      }
      else
	if (U <= p3) {                                    /*      X < x1     */
	  U = (U - p2)/(p3 - p2);
	  Y = log(U);
	  X = x1 + ll * Y;
	  if (X <= 0.)                                    /*      X > 0!!    */
	    continue; 
	  W = U * uniform();
	  /* squeeze accept:    L(x) = f(x1) (x - z1) / (x1 - z1)            */
	  /*                    z1 = x1 - ll,   W <= 1 + (X - x1)/ll         */
	  if (W <= 1. + Y)
	    break;
	  W *= f1;
	}
	else {                                            /*    x5 < X       */
	  U = (U - p3)/(p4 - p3);
	  Y = log(U);
	  X = x5 - lr * Y;
	  if (X >= 1.)                                    /*      X < 1!!    */
	    continue;
	  W = U * uniform();
	  /* squeeze accept:    L(x) = f(x5) (z5 - x) / (z5 - x5)            */
	  /*                    z5 = x5 + lr,   W <= 1 + (x5 - X)/lr         */
	  if (W <= 1. + Y)
	    break;
	  W *= f5;
	}
    /* density accept:  f(x) = (x/m)^(p_) ((1 - x)/(1 - m))^(q_)             */
    if (log(W) <= p_ * log(X/m) + q_ * log((1. - X)/(1. - m)))
      break;
  }
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==2) ? X : a + (b-a) * X);

} /* end of _unur_stdgen_sample_beta_b1prs() */

/*---------------------------------------------------------------------------*/
#undef p_
#undef q_
#undef s 
#undef m 
#undef D 
#undef Dl
#undef x1
#undef x2
#undef x4
#undef x5
#undef f1
#undef f2
#undef f4
#undef f5
#undef ll
#undef lr
#undef z2
#undef z4
#undef p1
#undef p2
#undef p3
#undef p4
/*---------------------------------------------------------------------------*/
