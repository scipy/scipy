/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: urng_unuran.c                                                     *
 *                                                                           *
 *   Unified interface for UNURAN uniform random number generators.          *
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
#ifdef UNUR_URNG_UNURAN
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Create new URNG object                                                 **/
/**                                                                         **/
/*****************************************************************************/

UNUR_URNG *
unur_urng_new( double (*sampleunif)(void *state), void *state )
     /*----------------------------------------------------------------------*/
     /* create a new URNG object for uniform random number generator.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   sampleunif ... pointer to sampling routine                         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to URNG object                                             */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  UNUR_URNG *urng = NULL;         /* pointer to URNG object */

  /* check arguments */
  _unur_check_NULL( "URNG", sampleunif, NULL );

  /* allocate memory for URNG object */
  urng = _unur_xmalloc( sizeof(struct unur_urng) );

  /* copy parameters into object */
  urng->sampleunif = sampleunif;
  urng->state      = state;

  /* initialize optional functions (set to not available) */
  urng->samplearray = NULL;
  urng->sync     = NULL;
  urng->seed     = ULONG_MAX;
  urng->setseed  = NULL;
  urng->delete   = NULL;
  urng->reset    = NULL;
  urng->nextsub  = NULL;
  urng->resetsub = NULL;
  urng->anti     = NULL;

  /* set magic cookie */
  COOKIE_SET(urng,CK_URNG);

  /* return object */
  return urng;
} /* end of unur_urng_new() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_sample_array( UNUR_URNG *urng, 
			   unsigned int (*samplearray)(void *state, double *X, int dim) )
     /*----------------------------------------------------------------------*/
     /* Set function to sample random point                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng        ... pointer to URNG object                             */
     /*   samplearray ... function for sampling random point                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->samplearray = samplearray;
  return UNUR_SUCCESS;

} /* end of unur_urng_set_sample_array() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_sync( UNUR_URNG *urng, void (*sync)(void *state) )
     /*----------------------------------------------------------------------*/
     /* Set function for jumping into defined state ("sync")                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   sync   ... function for syncing generator object                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->sync = sync;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_reset() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_seed( UNUR_URNG *urng, void (*setseed)(void *state, unsigned long seed) )
     /*----------------------------------------------------------------------*/
     /* Set function to seed.                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng    ... pointer to URNG object                                 */
     /*   setseed ... function for seeding generator                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->setseed = setseed;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_seed() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_anti( UNUR_URNG *urng, void (*setanti)(void *state, int anti) )
     /*----------------------------------------------------------------------*/
     /* Set function to switch to antithetic random numbers.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng    ... pointer to URNG object                                 */
     /*   setanti ... function for switching antithetic flag                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->anti = setanti;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_anti() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_reset( UNUR_URNG *urng, void (*reset)(void *state) )
     /*----------------------------------------------------------------------*/
     /* Set function for reseting URNG.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   reset  ... function for reseting generator object                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->reset = reset;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_reset() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_nextsub( UNUR_URNG *urng, void (*nextsub)(void *state) )
     /*----------------------------------------------------------------------*/
     /* Set function for jumping to start of the next substream of URNG.     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng    ... pointer to URNG object                                 */
     /*   nextsub ... function for jumping to next substream                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->nextsub = nextsub;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_nextsub() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_resetsub( UNUR_URNG *urng, void (*resetsub)(void *state) )
     /*----------------------------------------------------------------------*/
     /* Set function for jumping to start of the current substream of URNG.  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng     ... pointer to URNG object                                */
     /*   resetsub ... function for reseting current substream               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->resetsub = resetsub;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_resetsub() */

/*---------------------------------------------------------------------------*/

int
unur_urng_set_delete( UNUR_URNG *urng, void (*delete)(void *state) )
     /*----------------------------------------------------------------------*/
     /* Set function for destroying URNG.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   delete ... function for destroying generator object                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  urng->delete  = delete;
  return UNUR_SUCCESS;
} /* end of unur_urng_set_delete() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Handle a URNG object                                                   **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
unur_urng_sync (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Jump into defined state ("sync")                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  /* check whether we can reset the URNG object */
  if (urng->sync == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"sync");
    return UNUR_ERR_URNG_MISS;
  }

  /* jump to next substream */
  urng->sync (urng->state);

  return UNUR_SUCCESS;
} /* end of unur_urng_sync() */ 

/*---------------------------------------------------------------------------*/

int
unur_urng_seed (UNUR_URNG *urng, unsigned long seed)
     /*----------------------------------------------------------------------*/
     /* (Re-)Seed URNG.                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   seed   ... new seed for generator                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  /* check whether we can set the antithetic flag */
  if (urng->setseed == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"seeding function");
    return UNUR_ERR_URNG_MISS;
  }

  /* set seed */
  urng->setseed (urng->state,seed);

  /* store seed */
  urng->seed = seed;

  return UNUR_SUCCESS;
} /* end of unur_urng_seed() */

/*---------------------------------------------------------------------------*/

int
unur_urng_anti (UNUR_URNG *urng, int anti)
     /*----------------------------------------------------------------------*/
     /* Set antithetic flag for URNG.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*   anti   ... antithetic flag                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  /* check whether we can set the antithetic flag */
  if (urng->anti == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"antithetic flag");
    return UNUR_ERR_URNG_MISS;
  }

  /* set flag */
  urng->anti (urng->state,anti);

  return UNUR_SUCCESS;
} /* end of unur_urng_anti() */ 

/*---------------------------------------------------------------------------*/

int
unur_urng_nextsub (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Jump to start of next substream in URNG.                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  /* check whether we can reset the URNG object */
  if (urng->nextsub == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"next substream");
    return UNUR_ERR_URNG_MISS;
  }

  /* jump to next substream */
  urng->nextsub (urng->state);

  return UNUR_SUCCESS;
} /* end of unur_urng_nextsub() */ 

/*---------------------------------------------------------------------------*/

int
unur_urng_resetsub (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Reset current substream of the URNG object.                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) 
    /* use default generator */
    urng = unur_get_default_urng();

  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);

  /* check whether we can reset the URNG object */
  if (urng->resetsub == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"reset substream");
    return UNUR_ERR_URNG_MISS;
  }

  /* reset substream */
  urng->resetsub (urng->state);

  return UNUR_SUCCESS;
} /* end of unur_urng_resetsub() */ 

/*---------------------------------------------------------------------------*/

void
unur_urng_free (UNUR_URNG *urng)
     /*----------------------------------------------------------------------*/
     /* Destroy URNG object.                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   urng   ... pointer to URNG object                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  if (urng == NULL) return;  /* nothing to do */
  COOKIE_CHECK(urng,CK_URNG,RETURN_VOID);

  if (urng->delete != NULL) urng->delete (urng->state);
  free (urng);
  urng = NULL;

  return;
} /* end of unur_urng_free() */ 

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Handle URNG object in a generator object                               **/
/**                                                                         **/
/*****************************************************************************/

int
unur_gen_sync (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Jump into defined state ("sync") of uniform generator in             */
     /* generator object.                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_sync(gen->urng);
} /* end of unur_gen_sync() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_seed (UNUR_GEN *gen, unsigned long seed)
     /*----------------------------------------------------------------------*/
     /* (Re-) Seed uniform generator in generator object.                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*   seed ... new seed for generator                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_seed(gen->urng, seed);
} /* end of unur_gen_seed() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_anti (UNUR_GEN *gen, int anti)
     /*----------------------------------------------------------------------*/
     /* Set antithetic flag of uniform generator in generator object.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*   anti ... antithetic flag                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_anti(gen->urng, anti);
} /* end of unur_gen_anti() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_reset (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Reset uniform generator in generator object.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_reset(gen->urng);
} /* end of unur_gen_reset() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_nextsub (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Jump uniform generator in generator object to start of next          */
     /* substream.                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_nextsub(gen->urng);
} /* end of unur_gen_nextsub() */ 

/*---------------------------------------------------------------------------*/

int
unur_gen_resetsub (UNUR_GEN *gen)
     /*----------------------------------------------------------------------*/
     /* Reset current substream of uniform generator in generator object.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... generator object                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check argument */
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );

  return unur_urng_resetsub(gen->urng);
} /* end of unur_gen_resetsub() */ 

/*---------------------------------------------------------------------------*/
#endif    /* UNUR_URNG_UNURAN */
/*---------------------------------------------------------------------------*/

