/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_normal_gen.c                                               *
 *                                                                           *
 *   Special generators for Normal distribution                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010  Wolfgang Hoermann and Josef Leydold            *
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

inline static int normal_bm_init( struct unur_gen *gen );
inline static int normal_pol_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params  1      /* maximal number of parameters for generator */

#define mu    (DISTR.params[0])   /* location */
#define sigma (DISTR.params[1])   /* scale    */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_normal_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for normal distribution                 */
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

  case 1:    /* Box-Muller method */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_bm );
    return normal_bm_init( gen );

  case 2:    /* Polarmethod with rejection */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_pol );
    return normal_pol_init( gen );

  case 3:    /* Kindermann-Ramage method */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_kr );
    return UNUR_SUCCESS;

  case 0:    /* DEFAULT */
  case 4:    /* Acceptance-complement ratio */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_acr );
    return UNUR_SUCCESS;

  case 5:    /* "Naive" ratio-of-uniforms */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_nquo );
    return UNUR_SUCCESS;

  case 6:    /* Ratio-of-uniforms with squeeze */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_quo );
    return UNUR_SUCCESS;

  case 7:    /* Ratio-of-uniforms with quadratic bounding curves */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_leva );
    return UNUR_SUCCESS;

  case 99:   /* infamous sum-of-12-uniforms method. DO NOT USE */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_normal_sum );
    return UNUR_SUCCESS;

  default: /* no such generator */
    return UNUR_FAILURE;
  }

} /* end of _unur_stdgen_normal_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Sinus-Cosinus or Box/Muller Method                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * This method is based on the transformation                                *
 * x = sqrt(-2 ln(u)) cos(2pi*v), y = sqrt(-2 ln(u)) sin(2pi*v)              *
 * which converts two independent (0,1)-Uniforms u and v to two              *
 * independent standard Normal variates x and y.                             *
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - G.E.P. Box, M.E. Muller (1958):                             *
 *               A note on the generation of random normal deviates,         *
 *               Annals Math. Statist. 29, 610-611.                          *
 *                                                                           *
 * Implemented by W. Hoermann, April 1992                                    *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

#define Xstore  GEN->gen_param[0]
#define flag    GEN->flag

int
normal_bm_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  Xstore = 0.;
  flag = 1;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of normal_bm_init() */

double
_unur_stdgen_sample_normal_bm( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;
  double u,v,s;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  do {
    flag = -flag;
    if (flag > 0) { X = Xstore; break; }

    u = uniform();
    v = uniform();
    s = sqrt(-2.0 * log(u));
    Xstore = s * sin(2 * M_PI * v);
    X = s * cos(2 * M_PI * v);
  } while(0);
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_bm() */

#undef Xstore
#undef flag

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Polarmethod with rejection                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the standard                   *
 *               Normal distribution  N(0,1).                                *
 *                                                                           *
 * REFERENCE:  - G. Marsaglia (1962): Improving the Polar Method for         *
 *               Generating a Pair of Random Variables,                      *
 *               Boeing Sci. Res. Lab., Seattle, Washington.                 *
 *                                                                           *
 * Implemented by W. Hoermann, April 1992                                    *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

#define Xstore  GEN->gen_param[0]
#define flag    GEN->flag

int
normal_pol_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  Xstore = 0.;
  flag = 1;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of normal_pol_init() */

double
_unur_stdgen_sample_normal_pol( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;
  double s,x,y,tmp;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  do {
    flag = -flag;
    if (flag > 0) { X = Xstore; break; }

    while(1) {
      x = 2. * uniform() - 1.;
      y = 2. * uniform() - 1.;
      s = x*x + y*y;
      if( s < 1. ) {
	tmp = sqrt( -2. * log(s) / s );
	Xstore = y * tmp;
	X = x * tmp;
	break;
      }
    }
  } while(0);
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_pol() */

#undef Xstore
#undef flag

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: "Naive" Ratio of uniforms                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - A.J. Kindermann, F.J.Monahan (1977): Computing generation   *
 *               of random variables using the ratio of uniform deviates,    *
 *               ACM TOMS 3(3), 257-260                                      *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

double
_unur_stdgen_sample_normal_nquo( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;
  double u,v;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    u = uniform();
    if (_unur_iszero(u)) u = 1.;
    v = (uniform() - 0.5) * 0.857763885 * 2;
    X = v/u;
    if (X*X <= -4. * log(u)) 
      break;
  }
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_nquo() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Ratio of uniforms with squeeze                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - L. Barabesi (1993): Random variate generation               *
 *               by using the ratio-of-uniforms method, p. 133               *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

double
_unur_stdgen_sample_normal_quo( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;
  double r,w;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    r = uniform();
    X = (2.101083837941101 * uniform() - 1.050541918970551) / sqrt(r);
    w = X * X;
    if (4. - 4.186837275258269 * r < w) {
      if (1.5/r - 0.920558458320164 < w) 
	continue;
      if (-3.*log(r) < w )
	continue;
    }
    break;
  }
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_quo() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Ratio of uniforms with quadratic bounding curves     *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - J.L. Leva (1992):                                           *
 *               Algorithm 712; a normal random number generator,            *
 *               ACM TOMS 18(4), 454-455                                     *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

double
_unur_stdgen_sample_normal_leva( struct unur_gen *gen )
{
  /* -X- generator code -X- */
#define S    0.449871
#define T   -0.386595
#define A    0.19600
#define B    0.25472
#define RA   0.27597
#define RB   0.27846

  double X;
  double u,v,x,y,q;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  while (1) {
    u = uniform();
    v = uniform();
    v = 1.7156 * (v - 0.5);
    x = u - S;
    y = fabs(v) - T;
    q = x * x + y * (A * y - B * x);
    X = v/u;
    if( q < RA ) break;
    if( q > RB ) continue;
    if (v*v > -4.*log(u)*u*u) continue;
    break;
  }

#undef S
#undef T  
#undef A  
#undef B  
#undef RA 
#undef RB 
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_leva() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Kindermann-Ramage (patchwork) method                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - Kinderman A.J., Ramage J.G. (1976):                         *
 *               Computer Generation of Normal Random Variables,             *
 *               J. Am. Stat. Assoc. 71(356), 893 - 898.                     *
 *                                                                           *
 *             - Tirler G., Dalgaard P., Hoermann W., and Leydold J. (2004): *
 *               An Error in the {Kinderman-Ramage} Method and How to Fix It,*
 *               Comp. Stat. Data Anal. 47(3), 433 - 440.                    *
 *                                                                           *
 * Remark: This is the fixed version of the algorithm.                       *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

double
_unur_stdgen_sample_normal_kr( struct unur_gen *gen )
{
  /* -X- generator code -X- */
#define XI 2.216035867166471
#define PIhochK 0.3989422804
  
  double U, V, W, X;
  double t, z;
  
  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);
  
  U = uniform();
  
  if (U < 0.884070402298758) {
    V = uniform();
    X = XI * (1.131131635444180 * U + V - 1.);
  }
  
  else if (U >= 0.973310954173898) {
    do {
      V = uniform();
      W = uniform();
      if (_unur_iszero(W)) { t=0.; continue; }
      t = XI * XI/2. - log(W);
    } while ( (V*V*t) > (XI*XI/2.) );
    X = (U < 0.986655477086949) ? pow(2*t,0.5) : -pow(2*t,0.5);
  }
  
  else if (U>=0.958720824790463) {
    do {
      V = uniform();
      W = uniform();
      z = V - W;
      t = XI - 0.630834801921960 * _unur_min(V,W);
    } while (_unur_max(V,W) > 0.755591531667601 &&
	     0.034240503750111 * fabs(z) > (PIhochK * exp(t*t/(-2.)) - 0.180025191068563*(XI-fabs(t))) );
    X = (z<0) ? t : -t;
  }
  
  else if (U>=0.911312780288703) {
    do {
      V = uniform();
      W = uniform();
      z = V - W;
      t = 0.479727404222441 + 1.105473661022070 * _unur_min(V,W);
    } while (_unur_max(V,W) > 0.872834976671790 &&
	     0.049264496373128*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
    X = (z<0) ? t : -t;
  }
  
  else {
    do {
      V = uniform();
      W = uniform(); 
      z = V - W;
      t = 0.479727404222441 - 0.595507138015940 * _unur_min(V,W);
      /* 
       * The following line rejects negative values for t 
       * (the KR algorithm only generates random variates from
       * the half-normal distribution).
       * However, it is missing in the original algorithm as 
       * compiled in the paper.
       */
      if (t<=0.) continue;
    } while (_unur_max(V,W)>0.805777924423817 &&
	     0.053377549506886*fabs(z) > (PIhochK * exp(t*t/(-2)) -0.180025191068563*(XI-fabs(t))) );
    X = (z<0) ? t : -t;
  }
  
#undef XI
#undef PIhochK 
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );
  
} /* end of _unur_stdgen_sample_normal_kr() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: Acceptance-complement ratio                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - W. Hoermann and G. Derflinger (1990):                       *
 *               The ACR Methodfor generating normal random variables,       *
 *               OR Spektrum 12 (1990), 181-185.                             *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

double
_unur_stdgen_sample_normal_acr( struct unur_gen *gen )
{
  /* -X- generator code -X- */
#define c1 1.448242853
#define c2 3.307147487
#define c3 1.46754004
#define d1 1.036467755
#define d2 5.295844968
#define d3 3.631288474
#define hm 0.483941449
#define zm 0.107981933
#define hp 4.132731354
#define zp 18.52161694
#define phln 0.4515827053
#define hm1 0.516058551
#define hp1 3.132731354
#define hzm 0.375959516
#define hzmp 0.591923442
/*zhm 0.967882898*/

#define as 0.8853395638
#define bs 0.2452635696
#define cs 0.2770276848
#define b  0.5029324303
#define x0 0.4571828819
#define ym 0.187308492 
#define s  0.7270572718 
#define t  0.03895759111

  double X;
  double rn,x,y,z;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  do {
    y = uniform();

    if (y>hm1) {
      X = hp*y-hp1; break; }
  
    else if (y<zm) {  
      rn = zp*y-1;
      X = (rn>0) ? (1+rn) : (-1+rn);
      break;
    } 

    else if (y<hm) {  
      rn = uniform();
      rn = rn-1+rn;
      z = (rn>0) ? 2-rn : -2-rn;
      if ((c1-y)*(c3+fabs(z))<c2) {
	X = z; break; }
      else {  
	x = rn*rn;
	if ((y+d1)*(d3+x)<d2) {
	  X = rn; break; }
	else if (hzmp-y<exp(-(z*z+phln)/2)) {
	  X = z; break; }
	else if (y+hzm<exp(-(x+phln)/2)) {
	  X = rn; break; }
      }
    }

    while (1) {
      x = uniform();
      y = ym * uniform();
      z = x0 - s*x - y;
      if (z>0) 
	rn = 2+y/x;
      else {
	x = 1-x;
	y = ym-y;
	rn = -(2+y/x);
      }
      if ((y-as+x)*(cs+x)+bs<0) {
	X = rn; break; }
      else if (y<x+t)
	if (rn*rn<4*(b-log(x))) {
	  X = rn; break; }
    }
  } while(0);

  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_acr() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Normal Distribution: infamous sum-of-12-uniforms method.                  *
 *                                                                           *
 *                    NEVER use this method!!                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               standard Normal distribution  N(0,1).                       *
 *                                                                           *
 * REFERENCE:  - W. Hoermann and G. Derflinger (1990):                       *
 *               The ACR Methodfor generating normal random variables,       *
 *               OR Spektrum 12 (1990), 181-185.                             *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

double 
_unur_stdgen_sample_normal_sum( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  X = ( uniform() + uniform() + uniform() + uniform() + uniform() + uniform() +
	uniform() + uniform() + uniform() + uniform() + uniform() + uniform()
	- 6 );
  /* -X- end of generator code -X- */

  return ((DISTR.n_params==0) ? X : mu + sigma * X );

} /* end of _unur_stdgen_sample_normal_sum() */

/*---------------------------------------------------------------------------*/
#undef mu
#undef sigma
/*---------------------------------------------------------------------------*/

