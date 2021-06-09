/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_default.c                                                    *
 *                                                                           *
 *   routines to set, change and get the pointers to the                     *
 *   UNURAN default uniform random number generators.                        *
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
#include "urng.h"
#include <uniform/urng_builtin.h>
#include <uniform/urng_fvoid.h>
#include <uniform/urng_randomshift.h>

#if defined(UNURAN_HAS_GSL) && defined(UNUR_URNG_UNURAN)
#  include <uniform/urng_gsl.h>
#  include <uniform/urng_gslqrng.h>
#endif

#if defined(UNURAN_HAS_PRNG) && defined(UNUR_URNG_UNURAN)
#  include <uniform/urng_prng.h>
#endif

#if defined(UNURAN_HAS_RNGSTREAM) && defined(UNUR_URNG_UNURAN)
#  include <uniform/urng_rngstreams.h>
#endif

/*---------------------------------------------------------------------------*/
/* pointer to default uniform random number generator */

static UNUR_URNG *urng_default = NULL;
static UNUR_URNG *urng_aux_default = NULL;

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Main uniform random number generator                                   **/
/**                                                                         **/
/*****************************************************************************/

UNUR_URNG *
unur_get_default_urng( void )
     /*----------------------------------------------------------------------*/
     /* return default uniform random number generator                       */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to default generator                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* default generator already running ? */
  if( urng_default == NULL ) {
    /* have to initialize default generator first */
    urng_default = UNUR_URNG_DEFAULT;

    if( urng_default == NULL ) {
      /* some parameters invalid! */
      _unur_error("URNG",UNUR_ERR_NULL,"Cannot set default URNG. EXIT !!!");
      /* we cannot recover from this error */
      exit(EXIT_FAILURE);
    }
  }

  /* return default generator */
  return (urng_default);
} /* end of unur_get_default_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_set_default_urng( UNUR_URNG *urng_new )
     /*----------------------------------------------------------------------*/
     /* set default uniform random number generator and return old one       */
     /*                                                                      */
     /* parameters: pointer to new default uniform random number generator   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to old  uniform random number generator                    */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng_old = urng_default;

  /* NULL pointer not allowed */
  _unur_check_NULL("URNG", urng_new, urng_default);

  urng_default = urng_new;     /* set urng */

  /* return old default generator */
  return (urng_old);
} /* end of unur_set_default_urng() */


/*****************************************************************************/
/**                                                                         **/
/**  Auxiliary uniform random number generator                              **/
/**                                                                         **/
/*****************************************************************************/

UNUR_URNG *
unur_get_default_urng_aux( void )
     /*----------------------------------------------------------------------*/
     /* return default auxilliary uniform random number generator            */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to default auxilliary uniform generator                    */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* default generator already running ? */
  if( urng_aux_default == NULL ) {
    /* have to initialize default generator first */
    urng_aux_default = UNUR_URNG_AUX_DEFAULT;

    if( urng_aux_default == NULL ) {
      /* some parameters invalid! */
      _unur_error("URNG",UNUR_ERR_NULL,"Cannot set default auxilliary URNG. EXIT !!!");
      /* we cannot recover from this error */
      exit(EXIT_FAILURE);
    }
  }

  /* return default generator */
  return (urng_aux_default);
} /* end of unur_get_default_urng_aux() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_set_default_urng_aux( UNUR_URNG *urng_aux_new )
     /*----------------------------------------------------------------------*/
     /* set default auxilliary uniform RNG and return old one.               */
     /*                                                                      */
     /* parameters: pointer to new default auxilliary uniform RNG            */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to old auxilliary uniform RNG                              */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng_aux_old = urng_aux_default;

  /* NULL pointer not allowed */
  _unur_check_NULL("URNG", urng_aux_new, urng_aux_default);

  urng_aux_default = urng_aux_new;     /* set auxilliary urng */

  /* return old default generator */
  return (urng_aux_old);
} /* end of unur_set_default_urng_aux() */

/*---------------------------------------------------------------------------*/

