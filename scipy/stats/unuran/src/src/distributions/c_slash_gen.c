/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      c_slash_gen.c                                                *
 *                                                                           *
 *   Special generators for Slash distribution                               *
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
#include <methods/x_gen_source.h>
#include <distr/distr_source.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions_source.h"
#include "unur_distributions.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

inline static int slash_slash_init( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define uniform()  _unur_call_urng(gen->urng) /* call for uniform prng       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_slash_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for Slash distribution                  */
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
  case 1:  /* Ratio of normal and uniform random variates */
    _unur_cstd_set_sampling_routine(gen, _unur_stdgen_sample_slash_slash );
    return slash_slash_init( gen );

  default: /* no such generator */
    return UNUR_FAILURE;
  }
  
} /* end of _unur_stdgen_slash_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * Slash Distribution:   Z/U  (Z from N(0,1), U from U(0,1))                 *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the Slash distribution.        *
 *                                                                           *
 * Implemented by R. Kremer, 1990                                            *
 *****************************************************************************
 *    WinRand (c) 1995 Ernst Stadlober, Institut fuer Statistitk, TU Graz    *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#define NORMAL  gen->gen_aux   /* pointer to normal variate generator        */
/*---------------------------------------------------------------------------*/

inline static int
slash_slash_init( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_CSTD_GEN,UNUR_ERR_COOKIE);

  /* -X- setup code -X- */

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

} /* end of slash_slash_init() */

double
_unur_stdgen_sample_slash_slash( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,INFINITY);
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  /* -X- generator code -X- */
  return (_unur_sample_cont(NORMAL) / uniform());
  /* -X- end of generator code -X- */
  
} /* end of _unur_stdgen_sample_slash_slash() */

/*---------------------------------------------------------------------------*/
#undef NORMAL
/*---------------------------------------------------------------------------*/




