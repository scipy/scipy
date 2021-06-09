/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_binomial_gen.c                                             *
 *                                                                           *
 *   Special generators for Binomial distribution                            *
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
#include <methods/cstd.h>   /* for the definition of `UNUR_STDGEN_INVERSION' */
#include <methods/dstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int binomial_bruec_init( struct unur_gen *gen );
static int _unur_stdgen_sample_binomial_bruec( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_dstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params  (11)   /* maximal number of parameters for generator */
#define MAX_gen_iparams  (3)   /* maximal number of integer param. for gen.  */

/* parameters */
#define par_n  (DISTR.params[0])
#define par_p  (DISTR.params[1])

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_binomial_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Binomial distribution               */
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
  case 1:  /* Ratio of Uniforms/Inversion */
     if (gen==NULL) return UNUR_SUCCESS; /* test existence only  */
     _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_binomial_bruec );
     return binomial_bruec_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_binomial_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

#define flogfak(k) _unur_SF_ln_factorial(k)

/*****************************************************************************
 *                                                                           *
 *  Binomial Distribution - Ratio of Uniforms/Inversion                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Ratio of Uniforms method combined with Inversion for sampling from       *
 *  Binomial distributions with parameters n (number of trials) and          *
 *  p (probability of success).                                              *
 *  For  min(n*p,n*(1-p)) < 5  the Inversion method is applied:              *
 *  The random numbers are generated via sequential search, starting at the  *
 *  lowest index k=0. The cumulative probabilities are avoided by using the  *
 *  technique of chop-down.                                                  *
 *  For  min(n*p,n*(1-p)) >=5  Ratio of Uniforms is employed:                *
 *  A table mountain hat function h(x) with optimal scale parameter s for    *
 *  fixed location parameter  a = mu+1/2  is used. If the candidate k is     *
 *  near the mode and k>29 (respectively n-k > 29) f(k) is computed          *
 *  recursively starting at the mode m. The algorithm is valid for           *
 *  np > 0, n > =1, p <= 0.5. For p > 0.5, p is replaced by 1-p and          *
 *  k is replaced by n-k.                                                    *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:    - bruec samples a random number from the Binomial            *
 *                distribution with parameters n and p . It is valid for     *
 *                n * min(p,1-p) > 0.                                        *
 *                                                                           *
 * REFERENCE:   - E. Stadlober (1989): Sampling from Poisson, binomial and   *
 *                hypergeometric distributions: ratio of uniforms as a       *
 *                simple and fast alternative,                               *
 *                Bericht 303, Math. Stat. Sektion, Forschungsgesellschaft   *
 *                Joanneum, Graz.                                            *
 *                                                                           *
 * Implemented by R.Kremer 1990, revised by P.Busswald, July 1992            *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#define p       (DISTR.params[1])
#define n       (GEN->gen_iparam[0])

#define b       (GEN->gen_iparam[1])
#define m       (GEN->gen_iparam[2])

#define par     (GEN->gen_param[0])
#define q1      (GEN->gen_param[1])
#define np      (GEN->gen_param[3])
#define a       (GEN->gen_param[4])
#define h       (GEN->gen_param[5])
#define g       (GEN->gen_param[6])
#define r       (GEN->gen_param[7])
#define t       (GEN->gen_param[8])
#define r1      (GEN->gen_param[9])
#define p0      (GEN->gen_param[10])

/*---------------------------------------------------------------------------*/

int
binomial_bruec_init( struct unur_gen *gen )
{
  int bh,k1;
  double c,x; 

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
    GEN->n_gen_iparam = MAX_gen_iparams;
    GEN->gen_iparam = _unur_xmalloc(GEN->n_gen_iparam * sizeof(int));
  }

  /* convert integer parameters that are stored in an array of type 'double' */
  /* into those of type 'int' and store it in working array GEN->gen_iparam. */
  n = (int) par_n;

  /* -X- setup code -X- */

  par = _unur_min(p, 1.-p);
  q1 = 1.0 - par;
  np = n*par;                                /*np=min(n*p,n*(1-p))*/

  if (np < 5) {
    p0 = exp(n*log(q1));                     /* Set-up for Inversion */
    bh = (int)(np + 10.0*sqrt(np*q1));
    b = _unur_min(n,bh);                     /* safety-bound */
  }

  else {                                     /* Set-up for Ratio of Uniforms */
    m = (int)(np + par);                     /* mode */
    a = np + 0.5;                            /* shift parameter */
    c = sqrt(2.0 * a * q1);
    r = par/q1;
    t = (n+1) * r;
    r1 = log(r);
    bh = (int)(a + 7.0*c);
    b = _unur_min(n,bh);                     /* safety-bound */
    g = flogfak(m) + flogfak(n-m);           /* binomial const. */
    k1 = (int)(a-c);
    x = (a-k1-1.0)/(a-k1);
    if((n-k1)*par*x*x > (k1+1)*q1)
      k1++;                                  /* h=2*s */
    h = (a-k1) * exp(.5*((k1-m)*r1+g-flogfak(k1)-flogfak(n-k1))+M_LN2);
  }
	 
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of binomial_bruec_init() */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_binomial_bruec( struct unur_gen *gen )
{
  /* -X- generator code -X- */

  int i,k;
  double u,f,x,lf;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);
  
  if (np<5) {
    /* Inversion/Chop-down */

    double pk;

    k = 0;
    pk = p0;
    u = uniform();
    while (u>pk) {
      ++k;
      if (k>b) {
	u = uniform();
	k = 0;
	pk=p0;
      }
      else {
	u -= pk;
	pk=(double)(((n-k+1)*par*pk)/(k*q1));
      }
    }
    return ((p>0.5) ? n-k:k);
  }

  /* Ratio of Uniforms */
  for (;;) {
    do {
      u = uniform();
      x = a+h*(uniform()-0.5)/u;
    } while (x < 0 || ((k=(int)x) > b));    /* check, if k is valed candidate */

    if ((abs(m-k)<=15) && ((k>29)||(n-k>29)) ) {
      f = 1.0;                              /* compute f(k) recursively */
      if (m<k) {
	for (i=m;i<k;) 
	  f *= t / (double)++i-r;           /* f - f(k) */
	if (u*u <= f) break;                /* u^2<=f   */
      }
      else {
	for (i=k;i<m;)
	  f *= t / (double)++i-r;           /* f - 1/f(k) */
	if (u*u*f <= 1.0)
	  break;                            /* u^2<=f     */
      }
    }
    else {
      lf = (k-m)*r1+g-flogfak(k)-flogfak(n-k);       /* lf - ln(f(k)) */
      if ( u * (4.0 - u) - 3.0 <= lf) break;         /* lower squeeze */
      if (u*(u-lf) <= 1.0)                           /* upper squeeze */
	if (2.0*log(u) <= lf) break;                 /* final acceptance */
    }
  }
  return((p > 0.5) ? n-k : k);

  /* -X- end of generator code -X- */

} /* end of _unur_stdgen_sample_binomial_bruec() */

/*---------------------------------------------------------------------------*/

#undef p
#undef n

#undef m
#undef b

#undef par
#undef q1
#undef np
#undef al
#undef h
#undef g
#undef r
#undef t
#undef r1
#undef p0

#undef par_n
#undef par_p

/*---------------------------------------------------------------------------*/
