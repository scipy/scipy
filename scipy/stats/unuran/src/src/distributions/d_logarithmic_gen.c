/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      d_logarithmic_gen.c                                          *
 *                                                                           *
 *   Special generators for Logarithmic distribution                         *
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

#include <unur_source.h>
#include <methods/cstd.h>   /* for the definition of `UNUR_STDGEN_INVERSION' */
#include <methods/dstd_struct.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int logarithmic_lsk_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_dstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

#define MAX_gen_params  2      /* maximal number of parameters for generator */

/* parameters */
#define theta  (DISTR.params[0])    /* shape */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_logarithmic_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Logarithmic distribution            */
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
  case 1:  /* Inversion/Transformation */
    _unur_dstd_set_sampling_routine(gen, _unur_stdgen_sample_logarithmic_lsk );
    return logarithmic_lsk_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_logarithmic_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Logarithmic Distribution: Inversion/Transformation                        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Logarithmic distribution   *
 *               with parameter 0 < theta < 1.                               *
 *                                                                           *
 * REFERENCE:  - A.W. Kemp (1981): Efficient generation of logarithmically   *
 *               distributed pseudo-random variables,                        *
 *               Appl. Statist. 30, 249-253.                                 *
 *                                                                           *
 * Implemented by R.Kremer 1990, revised by P.Busswald, July 1992            *
 *****************************************************************************
 *                                                                           *
 * The algorithm combines Inversion and Transformation.                      *
 * It is based on the following fact: A random variable X from the           *
 * Logarithmic distribution has the property that X for fixed Y = y is       *
 * Geometric distributed with                                                *
 *  (*)   P(X=x|Y=y) = (1-y) * y^(x-1)                                       *
 * where Y has distribution function F(y) = ln(1-y) / ln(1-p).               *
 * So first random numbers y are generated by simple Inversion, then         *
 * k = (int)(1+ln(u) / ln(y)) is a Geometric random number and because of    *
 * (*) a Logarithmic one.                                                    *
 * To speed up the algorithm squeezes are used as well as the fact, that     *
 * many of the random numbers are 1 or 2 (depending on special               *
 * circumstances).                                                           *
 * On an IBM/PC 486 optimal performance is achieved, if for p<0.97 simple    *
 * inversion is used and otherwise the transformation.                       *
 * On an IBM/PC 286 inversion should be restricted to p<0.90.                *
 *                                                                           *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define t   (GEN->gen_param[0])
#define h   (GEN->gen_param[1])

#define theta_limit  0.97
/* theta <  theta_limit --> Inversion
   theta >= theta_limit --> Transformation  */
/*---------------------------------------------------------------------------*/

inline static int
logarithmic_lsk_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_DSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  if (theta < theta_limit) 
    t = -theta / log(1.0 - theta);
  else 
    h=log(1.0 - theta);
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of logarithmic_lsk_init() */


int
_unur_stdgen_sample_logarithmic_lsk( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double U, V, p, q;
  int K;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);
  COOKIE_CHECK(gen,CK_DSTD_GEN,INT_MAX);

  U = uniform();

  if (theta < theta_limit) {
    /* Inversion/Chop-down */
    K = 1;
    p = t;
    while (U > p) {
      U -= p;
      K++;
      p *= theta * (K - 1.)/((double) K);
    }
    return K;
  }

  else {
    /* Transformation  */
    if (U > theta) 
      return 1;

    V = uniform();
    q = 1. - exp(V * h);
    if ( U <= q * q) {
      K = 1 + (int)(log(U)/log(q));
      return K;
    }
    
    return ((U > q) ? 1 : 2);
  }

  /* -X- end of generator code -X- */
  
} /* end of _unur_stdgen_sample_logarithmic_lsk() */

/*---------------------------------------------------------------------------*/
#undef t
#undef h
#undef theta_limit
/*---------------------------------------------------------------------------*/
