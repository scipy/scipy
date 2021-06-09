/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      gibbs.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    Gibbs sampler using full conditional distributions.          *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF                                                            *
 *      Produce a value x consistent with its density                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function and its derivatives                  *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      mode of the density                                                  *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <distr/condi.h>
#include <distr/cvec.h>
#include <urng/urng.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "ars.h"
#include "tdr.h"

#include "gibbs.h"
#include "gibbs_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define GIBBS_VARMASK_VARIANT     0x000fu    /* indicates variant            */
#define GIBBS_VARIANT_COORD       0x0001u    /* coordinate sampler           */
#define GIBBS_VARIANT_RANDOMDIR   0x0002u    /* random direction sampler     */

#define GIBBS_VARMASK_T           0x00f0u    /* indicates transformation     */
#define GIBBS_VAR_T_SQRT          0x0010u    /* T(x) = -1/sqrt(x)            */
#define GIBBS_VAR_T_LOG           0x0020u    /* T(x) = log(x)                */
#define GIBBS_VAR_T_POW           0x0030u    /* T(x) = -x^c                  */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define GIBBS_DEBUG_CONDI   0x01000000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define GIBBS_SET_C          0x001u    /* set parameter for transformation T */
#define GIBBS_SET_X0         0x002u    /* set starting point                 */
#define GIBBS_SET_THINNING   0x004u    /* set thinning factor                */
#define GIBBS_SET_BURNIN     0x008u    /* set length of burn-in              */

/*---------------------------------------------------------------------------*/

#define GENTYPE "GIBBS"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_init( struct unur_par *par );
static int _unur_gibbs_coord_init( struct unur_gen *gen );
static int _unur_gibbs_randomdir_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_gibbs_coord_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_gibbs_randomdir_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_gibbs_normalgen( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create a normal random variate generator                                  */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_random_unitvector( struct unur_gen *gen, double *direction );
/*---------------------------------------------------------------------------*/
/* generate a random direction vector                                        */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_debug_init_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before init of generator starts.                                    */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_debug_init_condi( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print list of conditional generators.                                     */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_debug_burnin_failed( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after burnin has failed.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_gibbs_debug_init_finished( const struct unur_gen *gen, int success );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_gibbs_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_gibbs_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_gibbs_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

/* generators for conditional distributions */
#define GEN_CONDI     gen->gen_aux_list     

/* an auxiliary generator for standard normal variates */
#define GEN_NORMAL    gen->gen_aux

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CVEC *
_unur_gibbs_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_RANDOMDIR:
    return _unur_gibbs_randomdir_sample_cvec;
  case GIBBS_VARIANT_COORD:
  default:
    return _unur_gibbs_coord_sample_cvec;
  }
} /* end of _unur_gibbs_getSAMPLE() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_gibbs_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (DISTR_IN.logpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"logPDF");
    return NULL;
  }
  if (DISTR_IN.dlogpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dlogPDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_gibbs_par) );
  COOKIE_SET(par,CK_GIBBS_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR->c_T      = 0.;        /* parameter for transformation (-1. <= c < 0.) */

  par->method   = UNUR_METH_GIBBS ;       /* method                          */
  par->variant  = GIBBS_VARIANT_COORD;    /* default variant                 */
  par->set      = 0u;                     /* inidicate default parameters    */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  PAR->thinning = 1;                /* thinning factor                       */
  PAR->burnin   = 0;                /* length of burn-in for chain           */
  PAR->x0       = NULL;             /* starting point of chain, default is 0 */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_gibbs_init;

  return par;

} /* end of unur_gibbs_new() */

/*****************************************************************************/

int
unur_gibbs_set_variant_coordinate( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Coordinate Sampler :                                                 */
     /* Sampling along the coordinate directions (cyclic).                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );

  /* we use a bit in variant */
  par->variant = (par->variant & ~GIBBS_VARMASK_VARIANT) | GIBBS_VARIANT_COORD;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_gibbs_set_variant_coordinate() */

/*---------------------------------------------------------------------------*/

int
unur_gibbs_set_variant_random_direction( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* Random Direction Sampler :                                           */
     /* Sampling along the random directions.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );
            
  /* we use a bit in variant */
  par->variant = (par->variant & ~GIBBS_VARMASK_VARIANT) | GIBBS_VARIANT_RANDOMDIR;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_gibbs_set_variant_coordinate() */

/*---------------------------------------------------------------------------*/

int
unur_gibbs_set_c( struct unur_par *par, double c )
     /*----------------------------------------------------------------------*/
     /* set parameter c for transformation T_c used for sampling from        */
     /* conditional distributions                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   c    ... parameter c                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );

  /* check new parameter for generator */
  if (c > 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c > 0");
    return UNUR_ERR_PAR_SET;
  }
  /** TODO: ... **/
  /*    if (c <= -1.) { */
  /*      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"c <= -1 only if domain is bounded. Use `TABL' method then."); */
  /*      return 0; */
  /*    } */
  /** TODO: ... **/
  if (c < -0.5) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_SET,"c < -0.5 not implemented yet");
    return UNUR_ERR_PAR_SET;
  }
  if (!_unur_iszero(c) && c > -0.5) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"-0.5 < c < 0 not recommended. using c = -0.5 instead.");
    c = -0.5;
  }
    
  /* store date */
  PAR->c_T = c;

  /* changelog */
  par->set |= GIBBS_SET_C;

  return UNUR_SUCCESS;

} /* end of unur_gibbs_set_c() */

/*---------------------------------------------------------------------------*/

int 
unur_gibbs_set_startingpoint( struct unur_par *par, const double *x0)
     /*----------------------------------------------------------------------*/
     /* set starting point for chain                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   x0       ... starting point of chain                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );

  /* store data */
  PAR->x0 = x0;

  /* changelog */
  par->set |= GIBBS_SET_X0;

  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_gibbs_set_startingpoint() */

/*---------------------------------------------------------------------------*/

int
unur_gibbs_set_thinning( struct unur_par *par, int thinning )
     /*----------------------------------------------------------------------*/
     /* set thinning factor for chain                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   thinning ... thinning factor                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );

  /* check new parameter for generator */
  if (thinning < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"thinning < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->thinning = thinning;

  /* changelog */
  par->set |= GIBBS_SET_THINNING;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_gibbs_set_thinning() */

/*---------------------------------------------------------------------------*/

int
unur_gibbs_set_burnin( struct unur_par *par, int burnin )
     /*----------------------------------------------------------------------*/
     /* set length of burn-in for chain                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   burnin ... length of burn-in                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, GIBBS );

  /* check new parameter for generator */
  if (burnin < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"burnin < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->burnin = burnin;

  /* changelog */
  par->set |= GIBBS_SET_BURNIN;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_gibbs_set_burnin() */

/*---------------------------------------------------------------------------*/

const double *
unur_gibbs_get_state( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get current state of the Gibbs chain                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to chain ... on success                                    */
     /*   NULL             ... on error                                      */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, NULL );
  if (gen->method != UNUR_METH_GIBBS) {
    _unur_error(gen->genid, UNUR_ERR_GEN_INVALID,"");
    return NULL;
  }

  return GEN->state;
} /* end of unur_gibbs_get_state() */

/*---------------------------------------------------------------------------*/

int 
unur_gibbs_chg_state( struct unur_gen *gen, const double *state )
     /*----------------------------------------------------------------------*/
     /* chg current state of the Gibbs chain                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   state ... new state of chain                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, GIBBS, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( gen->genid, state, UNUR_ERR_NULL );

  /* copy state */
  memcpy( GEN->state, state, GEN->dim * sizeof(double));

  return UNUR_SUCCESS;
} /* end of unur_gibbs_chg_state() */

/*---------------------------------------------------------------------------*/

int 
unur_gibbs_reset_state( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* reset current state of the Gibbs chain to starting point             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, GIBBS, UNUR_ERR_GEN_INVALID );

  /* copy state */
  memcpy( GEN->state, GEN->x0, GEN->dim * sizeof(double));

  if (gen->variant & GIBBS_VARIANT_COORD)
    GEN->coord = (GEN->dim)-1;

  return UNUR_SUCCESS;
} /* end of unur_gibbs_reset_state() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_gibbs_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_GIBBS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_GIBBS_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_gibbs_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_gibbs_debug_init_start(gen);
#endif

  /* make generators for conditional distributions */
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:  /* --- coordinate direction sampling --- */
    if (_unur_gibbs_coord_init(gen)!=UNUR_SUCCESS) {
      _unur_gibbs_free(gen); return NULL;
    }
    break;
  case GIBBS_VARIANT_RANDOMDIR:  /* --- random direction sampling --- */
    if (_unur_gibbs_randomdir_init(gen)!=UNUR_SUCCESS) {
      _unur_gibbs_free(gen); return NULL;
    }
    break;
  default:
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_gibbs_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_gibbs_debug_init_condi(gen);
#endif

  /* run burn-in */
  if (GEN->burnin > 0 ) {
    int thinning, burnin;
    double *X;

    /* allocate memory for random vector */
    X = _unur_xmalloc( GEN->dim * sizeof(double) );

    /* store thinning factor; we use 1 for burn-in */
    thinning = GEN->thinning;
    GEN->thinning = 1;

    for (burnin = GEN->burnin; burnin>0; --burnin) {
      if ( _unur_sample_vec(gen,X) != UNUR_SUCCESS ) {
#ifdef UNUR_ENABLE_LOGGING
	_unur_gibbs_debug_burnin_failed(gen);
	if (gen->debug) _unur_gibbs_debug_init_finished(gen,FALSE);
#endif
	_unur_gibbs_free(gen); free (X); return NULL;
      }
    }

    /* restore thinning factor */
    GEN->thinning = thinning;
    /* free memory for random vector */
    free (X);
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_gibbs_debug_init_finished(gen,TRUE);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_gibbs_init() */

/*---------------------------------------------------------------------------*/

int
_unur_gibbs_coord_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for coordinate sampler                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par_condi;
  struct unur_gen *gen_condi;
  int i;
  int errorcode = UNUR_SUCCESS;

  /* conditional distribution object */
  GEN->distr_condi = unur_distr_condi_new( gen->distr, GEN->state, NULL, 0);

  for (i=0; i<GEN->dim; i++) {

    if ( (errorcode=unur_distr_condi_set_condition(GEN->distr_condi,GEN->state,NULL,i))
	 != UNUR_SUCCESS )
      break;

    /* make parameter object */
    switch( gen->variant & GIBBS_VARMASK_T ) {
    case GIBBS_VAR_T_LOG:
      /* use more robust method ARS for T = log */
      par_condi = unur_ars_new(GEN->distr_condi);
      unur_ars_set_reinit_percentiles(par_condi,2,NULL);
      break;
      
    case GIBBS_VAR_T_SQRT:
      /* we only have method TDR for T = -1/sqrt */
      par_condi = unur_tdr_new(GEN->distr_condi);
      unur_tdr_set_reinit_percentiles(par_condi,2,NULL);
      unur_tdr_set_c(par_condi,-0.5);
      unur_tdr_set_usedars(par_condi,FALSE);
      unur_tdr_set_variant_gw(par_condi);
      break;
      
    case GIBBS_VAR_T_POW:
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    
    /* we do not need a private copy since otherwise we cannot update the generator object */
    unur_set_use_distr_privatecopy( par_condi, FALSE );
    /* debugging messages from the conditional generator should be rare */
    /* otherwise the size if the LOG file explodes                      */
    unur_set_debug( par_condi, (gen->debug&GIBBS_DEBUG_CONDI)?gen->debug:(gen->debug?1u:0u));
  
    /* we use the same URNG for all auxiliary generators */
    unur_set_urng( par_condi, gen->urng );
  
    /* init generator object for sampling from conditional distributions */
    gen_condi = unur_init(par_condi);
  
    if (gen_condi == NULL) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		  "Cannot create generator for conditional distributions");
#ifdef UNUR_ENABLE_LOGGING
      if (gen->debug) _unur_gibbs_debug_init_finished(gen,FALSE);
#endif
      errorcode = UNUR_ERR_GEN_CONDITION;
      break;
    }

    /* store pointer to generator of conditional distribution */
    GEN_CONDI[i] = gen_condi;

    if (i==0 && DISTR.domainrect==NULL) {
      /* when the domain is not bounded by a rectangle we only  */
      /* need a clone for each dimension (except the first one) */
      for (i=1; i<GEN->dim; i++)
	GEN_CONDI[i] = unur_gen_clone(gen_condi);
      break;
    }
  }

  return errorcode;
} /* end of _unur_gibbs_coord_init() */

/*---------------------------------------------------------------------------*/

int
_unur_gibbs_randomdir_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize new generator for random direction sampler                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par_condi;
  struct unur_gen *gen_condi;

  /* we need an auxiliary generator for normal random variates */
  GEN_NORMAL = _unur_gibbs_normalgen( gen );
  if ( GEN_NORMAL == NULL ) return UNUR_FAILURE;
  
  /* conditional distribution object */
  _unur_gibbs_random_unitvector( gen, GEN->direction );
  GEN->distr_condi = unur_distr_condi_new( gen->distr, GEN->state, GEN->direction, 0);
  
  /* make parameter object */
  switch( gen->variant & GIBBS_VARMASK_T ) {
  case GIBBS_VAR_T_LOG:
    /* use more robust method ARS for T = log */
    par_condi = unur_ars_new(GEN->distr_condi);
    unur_ars_set_reinit_percentiles(par_condi,2,NULL);
    break;

  case GIBBS_VAR_T_SQRT:
    /* we only have method TDR for T = -1/sqrt */
    par_condi = unur_tdr_new(GEN->distr_condi);
    unur_tdr_set_reinit_percentiles(par_condi,2,NULL);
    unur_tdr_set_c(par_condi,-0.5);
    unur_tdr_set_usedars(par_condi,FALSE);
    break;

  case GIBBS_VAR_T_POW:
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  /* we do not need a private copy since otherwise we cannot update the generator object */
  unur_set_use_distr_privatecopy( par_condi, FALSE );
  /* debugging messages from the conditional generator should be rare */
  /* otherwise the size if the LOG file explodes                      */
  unur_set_debug( par_condi, (gen->debug&GIBBS_DEBUG_CONDI)?gen->debug:(gen->debug?1u:0u));
  
  /* we use the same URNG for all auxiliary generators */
  unur_set_urng( par_condi, gen->urng );
  
  /* init generator object for sampling from conditional distributions */
  gen_condi = unur_init(par_condi);
  
  if (gen_condi == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,
		"Cannot create generator for conditional distributions");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_gibbs_debug_init_finished(gen,FALSE);
#endif
    return UNUR_ERR_GEN_CONDITION;
  }

  /* store generator in structure. we only need one such generator */
  *GEN_CONDI = gen_condi;
    
  return UNUR_SUCCESS;
} /* end of _unur_gibbs_randomdir_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_gibbs_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_GIBBS_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_gibbs_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_GIBBS_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim;

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* which transformation for conditional distributions */
  if ( _unur_iszero(PAR->c_T) )
    par->variant = (par->variant & (~GIBBS_VARMASK_T)) | GIBBS_VAR_T_LOG;
  else if (_unur_FP_same(PAR->c_T, -0.5))
    par->variant = (par->variant & (~GIBBS_VARMASK_T)) | GIBBS_VAR_T_SQRT;
  else
    par->variant = (par->variant & (~GIBBS_VARMASK_T)) | GIBBS_VAR_T_POW;

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_gibbs_getSAMPLE(gen);
  gen->destroy = _unur_gibbs_free;
  gen->clone = _unur_gibbs_clone;

  /* variant of sampling method */
  gen->variant = par->variant;        
  
  /* copy parameters into generator object */
  GEN->thinning = PAR->thinning;           /* thinning factor                */
  GEN->burnin = PAR->burnin;               /* length of burnin               */
  GEN->c_T = PAR->c_T;                     /* parameter for transformation   */
  
  /* allocate memory for state */
  GEN->state = _unur_xmalloc( GEN->dim * sizeof(double));
  GEN->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
  if (PAR->x0 == NULL) 
    PAR->x0 = unur_distr_cvec_get_center(gen->distr);
  memcpy( GEN->state, PAR->x0, GEN->dim * sizeof(double));
  memcpy( GEN->x0, PAR->x0, GEN->dim * sizeof(double));

  /* generator(s) for conditional distributions */
  GEN->distr_condi = NULL;
  GEN_CONDI = _unur_xmalloc( GEN->dim * sizeof(struct unur_gen *) );
  gen->n_gen_aux_list = GEN->dim;   /* size of array GEN_CONDI */
  for (i=0; i<GEN->dim; i++) GEN_CONDI[i] = NULL;

  /* allocate memory for random direction */
  GEN->direction = _unur_xmalloc( GEN->dim * sizeof(double));

  /* defaults */
  GEN->coord = (GEN->dim)-1;      /* current coordinate of GIBBS chain.
				     we want to start with coordinate 0. */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_gibbs_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_gibbs_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_gibbs_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
#define CLONE         ((struct unur_gibbs_gen*)clone->datap)
#define CLONE_CONDI   clone->gen_aux_list     

  int i;
  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_GIBBS_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy state */
  CLONE->state = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy( CLONE->state, GEN->state, GEN->dim * sizeof(double));
  CLONE->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
  memcpy( CLONE->x0, GEN->x0, GEN->dim * sizeof(double));
  
  /* copy generators for conditional distributions */
  if (GEN->distr_condi) CLONE->distr_condi = _unur_distr_clone( GEN->distr_condi );

  /* GEN_CONDI is cloned by _unur_generic_clone */
  /* however, these use a pointer to GEN->distr which must be updated, too */
  if (CLONE_CONDI) {
    for (i=0; i<GEN->dim; i++)
      if (CLONE_CONDI[i])
	CLONE_CONDI[i]->distr = CLONE->distr_condi;
  }

  /* allocate memory for random direction */
  CLONE->direction = _unur_xmalloc( GEN->dim * sizeof(double));

  return clone;

#undef CLONE
#undef CLONE_CONDI
} /* end of _unur_gibbs_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_gibbs_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_GIBBS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free vectors */
  if (GEN->state) free (GEN->state);
  if (GEN->x0) free (GEN->x0);
  if (GEN->direction) free (GEN->direction);

  /* free conditional distribution object */
  if (GEN->distr_condi) _unur_distr_free (GEN->distr_condi);

  /* GEN_CONDI is freed by _unur_generic_free */

  _unur_generic_free(gen);

} /* end of _unur_gibbs_free() */

/*****************************************************************************/

int
_unur_gibbs_coord_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  double X;
  int thinning;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_GIBBS_GEN,UNUR_ERR_COOKIE);

  for (thinning = GEN->thinning; thinning > 0; --thinning) {

    /* update coordinate direction */
    GEN->coord = (GEN->coord + 1) % GEN->dim;
    
    /* check state of chain */
    if (!_unur_isfinite(GEN->state[GEN->coord]))
      /* there has been a fatal error during the last sampling from
	 the conditional distribution. This probably has been caused
	 by improper target distribution for which the method
	 does not work (i.e., it is not T-concave). 
      */
      continue;

    /* update conditional distribution */
    unur_distr_condi_set_condition( GEN->distr_condi, GEN->state, NULL, GEN->coord);

    /* reinit generator object */
    if (unur_reinit(GEN_CONDI[GEN->coord]) == UNUR_SUCCESS) {

      /* sample from distribution */
      X = unur_sample_cont(GEN_CONDI[GEN->coord]);

      if (_unur_isfinite(X)) {
	/* update state */
	GEN->state[GEN->coord] = X;
	continue;
      }
    }
    
    /* ERROR: reset to starting point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"reset chain");
    unur_gibbs_reset_state(gen);
    return UNUR_FAILURE;

  }
  
  /* copy current state into given vector */
  memcpy(vec, GEN->state, GEN->dim * sizeof(double)); 

  return UNUR_SUCCESS;

} /* end of _unur_gibbs_coord_sample_cvec() */

/*---------------------------------------------------------------------------*/

int
_unur_gibbs_randomdir_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  int i;
  double X;
  int thinning;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_GIBBS_GEN,UNUR_ERR_COOKIE);

  for (thinning = GEN->thinning; thinning > 0; --thinning) {

    /* check state of chain */
    if (!_unur_isfinite(GEN->state[0]))
      /* there has been a fatal error during the last sampling from
	 the conditional distribution. This probably has been caused
	 by improper target distribution for which the method
	 does not work (i.e., it is not T-concave). 
      */
      break;

    /* new random direction */
    _unur_gibbs_random_unitvector( gen, GEN->direction );
    
    /* update conditional distribution */
    unur_distr_condi_set_condition( GEN->distr_condi, GEN->state, GEN->direction, 0);

    /* reinit generator object */
    if (unur_reinit(*GEN_CONDI) == UNUR_SUCCESS) {

      /* sample from distribution */
      X = unur_sample_cont(*GEN_CONDI);

      if (_unur_isfinite(X)) {
	/* update state */
	for (i=0; i<GEN->dim; i++)
	  GEN->state[i] += X * GEN->direction[i];	  
	continue;
      }
    }

    /* ERROR: reset to starting point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"reset chain");
    unur_gibbs_reset_state(gen);
    return UNUR_FAILURE;

  }

  /* copy current state into given vector */
  memcpy(vec, GEN->state, GEN->dim * sizeof(double)); 

  return UNUR_SUCCESS;

} /* end of _unur_gibbs_randomdir_sample_cvec() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

struct unur_gen *
_unur_gibbs_normalgen( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create a normal random variate generator                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to GIBBS generator object                    */
     /*                                                                      */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen   *normalgen;
  struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
  struct unur_par   *normalpar = unur_arou_new( normaldistr );

  unur_arou_set_usedars( normalpar, TRUE );
  normalgen = unur_init( normalpar );
  _unur_distr_free( normaldistr );
  if (normalgen == NULL) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,
		"Cannot create aux Gaussian generator");
    return NULL;
  }

  /* uniform random number generator and debugging flags */
  normalgen->urng = gen->urng;
  normalgen->debug = gen->debug;

  return normalgen;

} /* end of _unur_gibbs_normalgen() */

/*---------------------------------------------------------------------------*/

void
_unur_gibbs_random_unitvector( struct unur_gen *gen, double *direction )
     /*----------------------------------------------------------------------*/
     /* generate a random direction vector                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*   direction ... random vector (result)                               */
     /*----------------------------------------------------------------------*/
{
  int i;

  do {
    for (i=0; i<GEN->dim; i++) 
      direction[i] = unur_sample_cont(GEN_NORMAL);
    /* normalize direction vector */
    _unur_vector_normalize(GEN->dim, direction);

    /* there is an extremely small change that direction is the null before
       normalizing. In this case non of its coordinates are finite. */
  } while (!_unur_isfinite(direction[0]));

} /* end of _unur_gibbs_random_unitvector() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_gibbs_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = GIBBS (Markov Chain - GIBBS sampler)\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:
    fprintf(LOG,"coordinate sampling (original Gibbs sampler)  [default]\n"); break;
  case GIBBS_VARIANT_RANDOMDIR:
    fprintf(LOG,"random directions\n"); break;
  }
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: transformation T_c(x) for TDR method = ",gen->genid);
  switch( gen->variant & GIBBS_VARMASK_T ) {
  case GIBBS_VAR_T_LOG:
    fprintf(LOG,"log(x)  ... c = 0");                   break;
  case GIBBS_VAR_T_SQRT:
    fprintf(LOG,"-1/sqrt(x)  ... c = -1/2");            break;
  case GIBBS_VAR_T_POW:
    fprintf(LOG,"-x^(%g)  ... c = %g",GEN->c_T,GEN->c_T); break;
  }
  _unur_print_if_default(gen,GIBBS_SET_C);
  fprintf(LOG,"\n%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:
    fprintf(LOG,"%s: sampling routine = _unur_gibbs_coord_sample()\n",gen->genid);
    break;
  case GIBBS_VARIANT_RANDOMDIR:
    fprintf(LOG,"%s: sampling routine = _unur_gibbs_randomdir_sample()\n",gen->genid);
    break;
  }

  fprintf(LOG,"%s: thinning = %d",gen->genid,GEN->thinning);
  _unur_print_if_default(gen,GIBBS_SET_THINNING);
  fprintf(LOG,"\n%s: burn-in = %d",gen->genid,GEN->burnin);
  _unur_print_if_default(gen,GIBBS_SET_BURNIN);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_matrix_print_vector( GEN->dim, GEN->x0, "starting point = ", LOG, gen->genid, "\t   ");

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_gibbs_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_gibbs_debug_init_finished( const struct unur_gen *gen, int success )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   success ... whether init has failed or not                         */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (success) 
    fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  else
    fprintf(LOG,"%s: INIT failed **********************\n",gen->genid);

} /* end of _unur_gibbs_debug_init_finished() */

/*---------------------------------------------------------------------------*/

void
_unur_gibbs_debug_burnin_failed( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info after burnin has failed                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: Burn-in failed --> INIT failed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_gibbs_debug_burnin_failed() */

/*---------------------------------------------------------------------------*/

void
_unur_gibbs_debug_init_condi( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write list of conditional generators into LOG file                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*----------------------------------------------------------------------*/
{
  int i;
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_GIBBS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  switch (gen->variant & GIBBS_VARMASK_VARIANT) {
  case GIBBS_VARIANT_COORD:
    fprintf(LOG,"%s: generators for full conditional distributions = \n",gen->genid);
    fprintf(LOG,"%s:\t",gen->genid);
    for (i=0; i<GEN->dim; i++)
      fprintf(LOG,"[%s] ", GEN_CONDI[i]->genid);
    fprintf(LOG,"\n%s:\n",gen->genid);
    break;
  case GIBBS_VARIANT_RANDOMDIR:
    fprintf(LOG,"%s: generators for full conditional distributions = [%s]\n",gen->genid,
	    GEN_CONDI[0]->genid);
    fprintf(LOG,"%s: generator for random directions = [%s]\n",gen->genid,
	    GEN_NORMAL->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    break;
  }

} /* end of _unur_gibbs_debug_init_condi() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_gibbs_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_distr_cvec_info_domain(gen);
  _unur_string_append(info,"   center    = ");
  _unur_distr_info_vector( gen, unur_distr_cvec_get_center(gen->distr), GEN->dim);
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]");
    else
      _unur_string_append(info,"  [default]");
  }
  _unur_string_append(info,"\n\n");
  
  /*   if (help) { */
  /*   _unur_string_append(info,"\n"); */
  /*   } */

  /* method */
  _unur_string_append(info,"method: GIBBS (GIBBS sampler [MCMC])\n");
  _unur_string_append(info,"   variant = %s\n",
		      ((gen->variant & GIBBS_VARMASK_VARIANT)==GIBBS_VARIANT_COORD)
		      ? "coordinate sampling [default]" : "random direction sampling");

  /* used transformation */
  _unur_string_append(info,"   T_c(x) = ");
  switch( gen->variant & GIBBS_VARMASK_T ) {
  case GIBBS_VAR_T_LOG:
    _unur_string_append(info,"log(x)  ... c = 0\n"); break;
  case GIBBS_VAR_T_SQRT:
    _unur_string_append(info,"-1/sqrt(x)  ... c = -1/2\n"); break;
  case GIBBS_VAR_T_POW:
    _unur_string_append(info,"-x^(%g)  ... c = %g\n",GEN->c_T,GEN->c_T); break;
  }

  _unur_string_append(info,"   thinning = %d\n", GEN->thinning);
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   rejection constant = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");

    switch (gen->variant & GIBBS_VARMASK_VARIANT) {
    case GIBBS_VARIANT_COORD:
      _unur_string_append(info,"   variant_coordinate  [default]\n"); break;
    case GIBBS_VARIANT_RANDOMDIR:
      _unur_string_append(info,"   variant_random_direction\n"); break;
    }

    _unur_string_append(info,"   c = %g  %s\n", GEN->c_T,
 			(gen->set & GIBBS_SET_C) ? "" : "[default]");
    _unur_string_append(info,"   thinning = %d  %s\n", GEN->thinning,
 			(gen->set & GIBBS_SET_THINNING) ? "" : "[default]");
    _unur_string_append(info,"   burnin = %d  %s\n", GEN->burnin,
 			(gen->set & GIBBS_SET_THINNING) ? "" : "[default]");
    
    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_gibbs_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
    */
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
} /* end of _unur_gibbs_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
