/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_zipf_gen.c                                                 *
 *                                                                           *
 *   Special generators for Zipf (or Zeta) distribution                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
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

#include <limits.h>
#include <unur_source.h>
#include <methods/cstd.h>   /* for the definition of `UNUR_STDGEN_INVERSION' */
#include <methods/dstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int zipf_zet_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_dstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params  2      /* maximal number of parameters for generator */

/* parameters */
#define rho  (DISTR.params[0])    /* shape */
#define tau  (DISTR.params[1])

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_zipf_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Zipf distribution                    */
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
  case 1:  /* Acceptance Rejection */
    _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_zipf_zet );
    return zipf_zet_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_zipf_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Zipf Distribution: Acceptance Rejection                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Zipf distribution with     *    
 *               parameters rho > 0 and tau >= 0.                            *
 *                                                                           *
 * REFERENCE : - J. Dagpunar (1988): Principles of Random Variate Generation,*
 *               Clarendon Press, Oxford.                                    *
 *                                                                           *
 * Implemented by P. Busswald, September 1992                                *
 *****************************************************************************
 *                                                                           *
 *  To sample from the Zipf (or Zeta) distribution with parameters rho and   *
 *  tau it suffices to sample variates x from the distribution with density  *
 *  function                                                                 *
 *     f(x) = B * {[x+0.5]+tau}^(-(1+rho)) ( x > 0.5 )                       *
 *  and then deliver k=[x+0.5].                                              *
 *  1/B = Sum[(j+tau)^-(rho+1)]  (j=1,2,...) converges for rho >= 0.5.       *
 *  It is not necessary to compute B, because variates x are generated by    *
 *  acceptance rejection using density function                              *
 *     g(x) = rho * (c+0.5)^rho * (c+x)^-(rho+1).                            *
 *  Integer overflow is possible, when roh is small (rho <= .5) and tau      *
 *  large. In this case a new sample is generated.                           *
 *  If rho and tau satisfy the inequality                                    *
 *     rho > 0.14 + tau * 1.85e-8 + 0.02 * ln(tau)                           *
 *  the percentage of overflow is less than 1%, so that the result is        *
 *  reliable.                                                                *
 *  If either rho > 100  or  k > 10000 numerical problems in computing the   *
 *  theoretical moments arise, therefore rho<=100 and k<=10000 are           *
 *  recommended.                                                             *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define c   (GEN->gen_param[0])
#define d   (GEN->gen_param[1])
/*---------------------------------------------------------------------------*/

inline static int
zipf_zet_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  if (rho<tau) {
    c = tau - 0.5;
    d = 0.;
  }
  else {
    c = rho - 0.5;
    d = (1. + rho) * log((1. + tau)/(1. + rho));
  }
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of zipf_zet_init() */

int
_unur_stdgen_sample_zipf_zet( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double U, V, E, X;
  int K;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);

  do {
    do {
      U = uniform();
      V = uniform();
      X = (c+0.5) * exp( -log(U)/rho ) - c;
    } while (X <= 0.5 || X >= (double) INT_MAX);
    K = (long int) (X+0.5);
    E = -log(V);
  } while ( E < (1.+rho) * log( (K+tau)/(X+c)) - d );
  /* -X- end of generator code -X- */

  return K;
  
} /* end of _unur_stdgen_sample_zipf_zet() */

/*---------------------------------------------------------------------------*/
#undef c
#undef d
/*---------------------------------------------------------------------------*/
