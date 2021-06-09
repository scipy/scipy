/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_set.c                                                        *
 *                                                                           *
 *   routines to set, change and get the pointers to the                     *
 *   uniform random number generators used in generator objects.
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

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Main uniform random number generator                                   **/
/**                                                                         **/
/*****************************************************************************/

int
unur_set_urng( struct unur_par *par, UNUR_URNG *urng )
     /*----------------------------------------------------------------------*/
     /* set uniform random number generator                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   urng    ... pointer to uniform random number generator             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, par, UNUR_ERR_NULL );
  _unur_check_NULL("URNG", urng, UNUR_ERR_NULL);

  /* main generator */
  par->urng = urng;

  /* overwrite auxilliary generator */
  if (par->urng_aux) par->urng_aux = urng;

  return UNUR_SUCCESS;
} /* end of unur_set_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_get_urng( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get uniform random number generator                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old uniform RNG                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->urng;
} /* end of unur_get_urng() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_chg_urng( struct unur_gen *gen, UNUR_URNG *urng )
     /*----------------------------------------------------------------------*/
     /* set uniform random number generator                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   urng    ... pointer to uniform random number generator             */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old uniform RNG                                         */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng_old;

  /* check arguments */
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng,NULL);

  urng_old = gen->urng;

  /* set pointer to main URNG */
  gen->urng = urng;

  /* also set pointer in auxiliary generator objects */
  if (gen->gen_aux)
    unur_chg_urng(gen->gen_aux,urng);

  if (gen->gen_aux_list && gen->n_gen_aux_list) {
    int i;
    for (i=0; i<gen->n_gen_aux_list; i++) {
      if (gen->gen_aux_list[i])
	unur_chg_urng(gen->gen_aux_list[i],urng);
    }
  }

  /* overwrite auxilliary URNG */
  if (gen->urng_aux) gen->urng_aux = urng;

  return urng_old;
} /* end of unur_chg_urng() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Auxiliary uniform random number generator                              **/
/**                                                                         **/
/*****************************************************************************/

int
unur_set_urng_aux( struct unur_par *par, UNUR_URNG *urng_aux )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par_aux ... pointer to parameter for building generator object     */
     /*   urng    ... pointer to auxilliary uniform random number generator  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, par, UNUR_ERR_NULL );
  _unur_check_NULL("URNGaux", urng_aux, UNUR_ERR_NULL);

  if (par->urng_aux == NULL)
    /* no auxilliary generator is required */
    return UNUR_ERR_GENERIC;

  par->urng_aux = urng_aux;

  return UNUR_SUCCESS;
} /* end of unur_set_urng_aux() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_get_urng_aux( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get auxilliary uniform random number generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old auxilliary uniform RNG                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,NULL);

  return gen->urng_aux;
} /* end of unur_get_urng_aux() */

/*---------------------------------------------------------------------------*/

UNUR_URNG *
unur_chg_urng_aux( struct unur_gen *gen, UNUR_URNG *urng_aux )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   urng_aux ... pointer to auxilliary uniform random number generator */
     /*                                                                      */
     /* return:                                                              */
     /*   Pointer to old auxilliary uniform RNG                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng_aux_old;

  /* check arguments */
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng_aux,NULL);

  if (gen->urng_aux == NULL) 
    /* no auxilliary generator is required */
    return NULL;

  urng_aux_old = gen->urng_aux;

  /* set pointer to main URNG */
  gen->urng_aux = urng_aux;

  /* also set pointer in auxiliary generator objects */
  if (gen->gen_aux)
    unur_chg_urng_aux(gen->gen_aux,urng_aux);

  if (gen->gen_aux_list && gen->n_gen_aux_list) {
    int i;
    for (i=0; i<gen->n_gen_aux_list; i++) {
      if (gen->gen_aux_list[i])
	unur_chg_urng_aux(gen->gen_aux_list[i],urng_aux);
    }
  }

  return urng_aux_old;
} /* end of unur_chg_urng_aux() */

/*---------------------------------------------------------------------------*/

int
unur_use_urng_aux_default( UNUR_PAR *par )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator to default            */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  if (par->urng_aux == NULL)
    /* no auxilliary generator is required */
    return UNUR_ERR_GENERIC;

  /* set aux URNG */
  par->urng_aux = unur_get_default_urng_aux();

  return UNUR_SUCCESS;

} /* end of unur_use_urng_aux_default() */

/*---------------------------------------------------------------------------*/

int
unur_chgto_urng_aux_default( UNUR_GEN *gen )
     /*----------------------------------------------------------------------*/
     /* set auxilliary uniform random number generator to default            */
     /* (initialize generator if necessary)                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  if (gen->urng_aux == NULL)
    /* no auxilliary generator is required */
    return UNUR_ERR_GENERIC;

  /* set aux URNG */
  gen->urng_aux = unur_get_default_urng_aux();

  return UNUR_SUCCESS;

} /* end of unur_chgto_urng_aux_default() */

/*---------------------------------------------------------------------------*/


