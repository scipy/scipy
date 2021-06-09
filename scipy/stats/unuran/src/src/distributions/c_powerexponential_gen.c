/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_powerexponential_gen.c                                     *
 *                                                                           *
 *   Special generators for Power-exponential distribution                   *
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

inline static int powerexponential_epd_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define MAX_gen_params  2      /* maximal number of parameters for generator */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

/*---------------------------------------------------------------------------*/
#define tau    (DISTR.params[0])        /* shape */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_powerexponential_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for power-exponential distribution      */
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
  case 1:  /* Transformed density rejection */
    { /* check parameters of distribution */
      double d_tau = (par) ? par->distr->data.cont.params[0] : tau;
      if (d_tau < 1.) {
	_unur_error(NULL,UNUR_ERR_GEN_CONDITION,"");
	return UNUR_ERR_GEN_CONDITION;
      }
    }

    /* tau >= 1 !!!! */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_powerexponential_epd );
    return powerexponential_epd_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }

} /* end of _unur_stdgen_powerexponential_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Powerexponential Distribution: Transformed density rejection              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Power-exponential          *
 *               distribution with parameter tau >= 1 by using the           *
 *               non-universal rejection method for logconcave densities.    *
 *                                                                           *
 * REFERENCE:  - L. Devroye (1986): Non-Uniform Random Variate Generation,   *
 *               Springer Verlag, New York.                                  *
 *                                                                           *
 * Implemented by K. Lehner, 1990                                            *
 * Revised by F. Niederl, August 1992                                        *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define s    GEN->gen_param[0]
#define sm1  GEN->gen_param[1]
/*---------------------------------------------------------------------------*/

inline static int
powerexponential_epd_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  if (GEN->gen_param == NULL) {
    GEN->n_gen_param = MAX_gen_params;
    GEN->gen_param = _unur_xmalloc(GEN->n_gen_param * sizeof(double));
  }

  /* -X- setup code -X- */
  s = 1. / tau;
  sm1 = 1. - s;
  /* -X- end of setup code -X- */

  return UNUR_SUCCESS;

} /* end of powerexponential_epd_init() */

double 
_unur_stdgen_sample_powerexponential_epd( struct unur_gen *gen )
{
  /* -X- generator code -X- */
  double U,u1,V,X,y;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  do {
    U = 2. * uniform() - 1.;                                  /* U(-1.0/1.0) */
    u1 = fabs(U);                                             /* u1=|u|      */
    V = uniform();                                            /* U(0/1)      */

    if (u1 <= sm1)
      /* Uniform hat-function for x <= (1-1/tau)   */
      X = u1;
    else {                       
      /* Exponential hat-function for x > (1-1/tau) */
      y = tau * (1. - u1);                                         /* U(0/1) */
      X = sm1 - s * log(y);
      V *= y;
    }
  } while (log(V) > -exp(log(X)*tau));               /* Acceptance/Rejection */
  
  /* Random sign */
  if (U > 0.)
    X = -X;

  /* -X- end of generator code -X- */

  return X;

} /* end of _unur_stdgen_sample_powerexponential_epd() */

/*---------------------------------------------------------------------------*/
#undef s
#undef sm1
/*---------------------------------------------------------------------------*/
