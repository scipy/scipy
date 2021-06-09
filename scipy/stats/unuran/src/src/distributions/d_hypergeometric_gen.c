/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_hypergeometric_gen.c                                       *
 *                                                                           *
 *   Special generators for Hypergeometric distribution                      *
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

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_dstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params  (8)   /* maximal number of parameters for generator  */
#define MAX_gen_iparams (9)   /* maximal number of integer parameters for generator */

/* parameters */
#define par_N  (DISTR.params[0])
#define par_M  (DISTR.params[1])
#define par_n  (DISTR.params[2])

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int hypergeometric_hruec_init( struct unur_gen *gen );
static int _unur_stdgen_sample_hypergeometric_hruec( struct unur_gen *gen );
inline static int h_util(int N_, int M_, int n_, int k_);

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_hypergeometric_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Hypergeometric distribution         */
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
  case 1:  /* HRUEC  method */
     _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_hypergeometric_hruec );
     return hypergeometric_hruec_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_hypergeometric_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 *   Hypergeometric Distribution - Ratio of Uniforms/Inversion               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Ratio of Uniforms combined with Inversion for sampling from              *
 *  Hypergeometric distributions with parameters N, M and n. The algorithm   *
 *  is valid for M <= N/2, n <= N/2. Otherwise parameters (at the beginning  *
 *  of the algorithm) and random numbers k are adapted in function h_util(). *
 *  For mode m < 5 Inversion is applied:                                     *
 *  The random numbers are generated via sequential search, starting at the  *
 *  lowest index k=0. The cumulative probabilities are avoided by using the  *
 *  technique of chop-down.                                                  *
 *  For mode  m >=5  Ratio of Uniforms is employed: A table mountain hat     *
 *  function h(x) with optimal scale parameter s for fixed location          *
 *  parameter  a = mu+1/2  is used.                                          *
 *  If the mode m <= 20 and the candidate k is near the mode f(k) is         *
 *  computed recursively starting at the mode  m.                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:    - hruec samples a random number from the Hypergeometric      *
 *                distribution with parameters N (number of red and          *
 *                black balls), M (number of red balls) and n (number of     *
 *                trials) valid for N >= 2, M,n <= N.                        *
 * REFERENCE:   - E. Stadlober (1989): Sampling from Poisson, binomial and   *
 *                hypergeometric distributions: ratio of uniforms as a       *
 *                simple and fast alternative,                               *
 *                Bericht 303, Math. Stat. Sektion,                          *
 *                Forschungsgesellschaft Joanneum, Graz.                     *
 *                                                                           *
 * Implemented by R.Kremer 1990, revised by P.Busswald, July 1992            *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

#define flogfak(k) (_unur_SF_ln_factorial(k))
#define delta(k) (flogfak(k)+flogfak(Mc-k)+flogfak(nc-k)+flogfak(NMn+k))

#define N       (GEN->gen_iparam[0])
#define M       (GEN->gen_iparam[1])
#define n       (GEN->gen_iparam[2])

#define b       (GEN->gen_iparam[3])
#define m       (GEN->gen_iparam[4])
#define NMn     (GEN->gen_iparam[5])
#define Mc      (GEN->gen_iparam[6])
#define nc      (GEN->gen_iparam[7])
#define N_half  (GEN->gen_iparam[8])

#define NMnp    (GEN->gen_param[0])
#define Np      (GEN->gen_param[1])
#define Mp      (GEN->gen_param[2])
#define np      (GEN->gen_param[3])
#define g       (GEN->gen_param[4])
#define a       (GEN->gen_param[5])
#define h       (GEN->gen_param[6])
#define p0      (GEN->gen_param[7])

/*---------------------------------------------------------------------------*/

int
hypergeometric_hruec_init( struct unur_gen *gen )
{
  int k1,bh;
  double x,p,q,c,my;

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
  N = (int) par_N;
  M = (int) par_M;
  n = (int) par_n;

  /* -X- setup code -X- */

  N_half = N/2;                      /* Preparations of the parameters  */
  Mc = (M<=N_half) ? M : N-M;        /* if M<=N/2, M is replaced by N-M */
  nc = (n<=N_half) ? n : N-n;        /* if n<=N/2, n is replaced by n-M */
  
  /* Set-up */
  Np = (double) N;
  Mp = (double) Mc;
  np = (double) nc;

  NMn = N - Mc - nc;
  NMnp = Np - Mp - np;
  p = Mp / Np;
  q = 1.0 - p;
  my = np * p;
  bh = _unur_min(nc,Mc);
  m = (int) ((np+1.0)*(Mp+1.0)/(Np+2.0));       /* mode */

  if (m < 5) {
    /* Set-up for Inversion */
    c = my + 10.0*sqrt(my*q*(1.0-np/Np));
    b = _unur_min(bh,(int)c);                   /* safety-bound */
    p0 = exp(flogfak(N-Mc)+flogfak(N-nc)-flogfak(NMn)-flogfak(N));
  }

  else {
    /* Set-up for Ratio of Uniforms */
    a = my+0.5;
    c = sqrt(2.0*a*q*(1.0-np/Np));
    b = _unur_min(bh,(int)(a+7.0*c));           /* safety-bound */
    g = delta(m);
    k1 = (int)(a-c);
    x = (a-k1-1.0)/(a-k1);
    if((np-k1)*(p-(double)k1/Np)*x*x > (k1+1)*(q-(np-k1-1.0)/Np))
      k1++;
    h = (a-k1)*exp(0.5*(g-delta(k1))+M_LN2);    /* h=2*s */
  }

  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of hypergeometric_hruec_init() */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_hypergeometric_hruec( struct unur_gen *gen )
{
  int k,i;
  double x,u,f,lf;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);

  /* -X- generator code -X- */
  
  if (m<5) {                                     /* Inversion/Chop-down */
    double pk;

    k = 0;
    pk = p0;
    u = uniform();
    while (u>pk) {
      ++k;
      if (k>b) {
	u = uniform();
	k = 0;
	pk = p0;
      }
      else {
	u -= pk;
	pk *= ((Mp-k+1.0)*(np-k+1.0)) / ((double)k*(NMnp+k));
      }
    }

    return (h_util(N,M,n,k));
  }

  for (;;) {                                    /* Ratio of Uniforms */
    do {
      u = uniform();
      x = a + h*(uniform()-0.5) / u;
    } while (x < 0 || ((k=(int)x) > b));        /* check, if k is valid candidate */

    if (m <= 20 || abs(m-k) <= 15) {           /* compute f(k) recursively */
      f = 1.0;
      if (m<k) {
	for (i=m+1;i<=k;i++)
	  f *= ((Mp-i+1.0)*(np-i+1.0)) / ((double)i*(NMnp+i));
	if (u*u <= f) break;                    /* f - f(k), u^2<=f */
      }
      else {
	for (i=k+1;i<=m;i++)
	  f *= ((Mp-i+1.0)*(np-i+1.0)) / ((double)i*(NMnp+i));
	if (u*u*f <= 1.0) break;                /* f - 1/f(k), u^2<=f */
      }
    }

    else {
      lf = g - delta(k);                        /* lf - ln(f(k)) */
      if ( u * (4.0 - u) - 3.0 <= lf) break;    /* lower squeeze */
      if (u*(u-lf) <= 1.0)                      /* upper squeeze */
	if (2.0*log(u) <= lf) break;            /* final acceptance */
    }
  }

  return (h_util(N,M,n,k));

  /* -X- end of generator code -X- */

} /* end of _unur_stdgen_sample_hypergeometric_hruec() */

/*---------------------------------------------------------------------------*/

#undef N
#undef M
#undef n

#undef b    
#undef m   
#undef NMn 
#undef Mc
#undef nc
#undef N_half

#undef NMnp
#undef Np  
#undef Mp  
#undef np  
#undef g   
#undef a   
#undef h   
#undef p0  

#undef delta
#undef flogfak

#undef par_N
#undef par_M
#undef par_n

/*---------------------------------------------------------------------------*/

int
h_util(int N_, int M_, int n_, int k)
     /* Transformation to variate k from H(N,M,n) */
{
 int N_half;

 N_half = N_/2;

 if (n_ <= N_half)
	 {
	  if (M_ <= N_half) return(k);   /* no replacements */
	  else return(n_ - k);           /* M has been replaced by N-M, therefore */
	 }                            /* k has to be replaced by n-k           */
 else
	 {
	  if (M_ <= N_half) return(M_ - k); /* n h.b.r. by N-n, therefore k by M-k */
	  else return(M_ - N_ + n_ + k);       /* M h.b.r. by N-M and n by N-n,       */
	 }                            /* therefore k by M-N+n+k              */
} /* end of h_util() */

/*---------------------------------------------------------------------------*/
