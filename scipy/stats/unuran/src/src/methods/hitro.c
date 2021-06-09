/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      hitro.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    HITRO sampler using full conditional distributions.          *
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
#include <distr/cvec.h>
#include <urng/urng.h>
#include <utils/matrix_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#include "unur_methods_source.h"
#include "arou.h"
#include "x_gen.h"
#include "x_gen_source.h"

#include "hitro.h"
#include "hitro_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif


/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* minimum value for multiplier (same number also appears in error message)  */
#define HITRO_MIN_MULTIPLIER  (1.0001)

/* starting value for u and v bounds when not provided by user */
#define HITRO_START_UVMIN  (1.e-3)

/* increase size of adaptive bounding rectangle by this factor */
#define HITRO_DEFAULT_ADAPTIVE_MULTIPLIER  (1.1)

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define HITRO_VARMASK_VARIANT     0x000fu   /* indicates variant             */
#define HITRO_VARIANT_COORD       0x0001u   /* coordinate sampler            */
#define HITRO_VARIANT_RANDOMDIR   0x0002u   /* random direction sampler      */

#define HITRO_VARFLAG_ADAPTLINE   0x0010u   /* use adaptive line sampling    */
#define HITRO_VARFLAG_ADAPTRECT   0x0020u   /* use adaptive bounding rectangle */
#define HITRO_VARFLAG_BOUNDRECT   0x0040u   /* use entire bounding rectangle */
#define HITRO_VARFLAG_BOUNDDOMAIN 0x0080u   /* use bounded domain (if given) */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define HITRO_SET_R          0x0001u   /* set parameter r                    */
#define HITRO_SET_X0         0x0002u   /* set starting point                 */
#define HITRO_SET_THINNING   0x0004u   /* set thinning factor                */
#define HITRO_SET_BURNIN     0x0008u   /* set length of burn-in              */
#define HITRO_SET_U          0x0010u   /* set u-boundaries of bounding rectangles */
#define HITRO_SET_V          0x0020u   /* set upper v-bound of bounding rectangles */

#define HITRO_SET_ADAPTLINE  0x0100u   /* set adaptive line sampling         */
#define HITRO_SET_ADAPTRECT  0x0200u   /* set adaptive bounding rectangle    */
#define HITRO_SET_BOUNDRECT  0x0400u   /* set entire bounding rectangle      */
#define HITRO_SET_ADAPTMULT  0x0800u   /* set multiplier for adaptive rectangles */

/*---------------------------------------------------------------------------*/

#define GENTYPE "HITRO"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitro_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitro_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitro_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_hitro_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_hitro_coord_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_hitro_randomdir_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_hitro_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

static void _unur_hitro_xy_to_vu( const struct unur_gen *gen, const double *x, double y, double *vu );
static void _unur_hitro_vu_to_x( const struct unur_gen *gen, const double *vu, double *x );
/*---------------------------------------------------------------------------*/
/* transforming between coordinates xy in original scale and                 */
/* coordinates vu in Ratio-of-Unforms scale.                                 */
/*---------------------------------------------------------------------------*/

static double _unur_hitro_xv_to_u( const struct unur_gen *gen, double x, double v, int k );
/*---------------------------------------------------------------------------*/
/* transform point where we are given the v-coordinate (in RoU scale)        */
/* and the k-th x-coordinate (in xy scale).                                  */
/*---------------------------------------------------------------------------*/

static int _unur_hitro_vu_is_inside_region( const struct unur_gen *gen, const double *vu );
/*---------------------------------------------------------------------------*/
/* check whether point with vu-coordinates is inside acceptance region       */
/* of RoU method.                                                            */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hitro_normalgen( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create a normal random variate generator                                  */
/*---------------------------------------------------------------------------*/

static void _unur_hitro_random_unitvector( struct unur_gen *gen, double *direction );
/*---------------------------------------------------------------------------*/
/* generate a random direction vector                                        */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_hitro_debug_init_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before init of generator starts.                                    */
/*---------------------------------------------------------------------------*/

static void _unur_hitro_debug_init_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_hitro_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_hitro_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_hitro_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_hitro_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/* an auxiliary generator for standard normal variates */
#define GEN_NORMAL    gen->gen_aux

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CVEC *
_unur_hitro_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant & HITRO_VARMASK_VARIANT) {
  case HITRO_VARIANT_COORD:
    return _unur_hitro_coord_sample_cvec;
  case HITRO_VARIANT_RANDOMDIR:
  default:
    return _unur_hitro_randomdir_sample_cvec;
  }
} /* end of _unur_hitro_getSAMPLE() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_hitro_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_hitro_par) );
  COOKIE_SET(par,CK_HITRO_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_HITRO ;             /* method                    */
  par->variant  = ( HITRO_VARIANT_COORD |       /* default variant           */
		    HITRO_VARFLAG_ADAPTLINE );  /*   see also _unur_hitro_init() */
  par->set      = 0u;                     /* inidicate default parameters    */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  PAR->r        = 1.;               /* parameter r                           */
  PAR->thinning = 1;                /* thinning factor                       */
  PAR->burnin   = 0;                /* length of burn-in for chain           */
  PAR->x0       = NULL;             /* starting point of chain, default is 0 */
  PAR->adaptive_mult = HITRO_DEFAULT_ADAPTIVE_MULTIPLIER; /* multiplier for adaptive rectangles  */

  PAR->vmax     = -1.;        /* v-boundary of bounding rectangle (unknown)  */
  PAR->umin     = NULL;       /* u-boundary of bounding rectangle (unknown)  */
  PAR->umax     = NULL;       /* u-boundary of bounding rectangle (unknown)  */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_hitro_init;

  return par;

} /* end of unur_hitro_new() */

/*****************************************************************************/

int
unur_hitro_set_variant_coordinate( struct unur_par *par )
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
  _unur_check_par_object( par, HITRO );

  /* we use a bit in variant */
  par->variant = (par->variant & ~HITRO_VARMASK_VARIANT) | HITRO_VARIANT_COORD;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitro_set_variant_coordinate() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_variant_random_direction( struct unur_par *par )
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
  _unur_check_par_object( par, HITRO );
            
  /* we use a bit in variant */
  par->variant = (par->variant & ~HITRO_VARMASK_VARIANT) | HITRO_VARIANT_RANDOMDIR;
  
  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitro_set_variant_coordinate() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_use_adaptiveline( struct unur_par *par, int adaptive )
     /*----------------------------------------------------------------------*/
     /* Enable/Disable adaptive line sampling                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   adaptive ... whether adaptive line sampling is enabled             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
            
  /* we use a bit in variant */
  par->variant = (adaptive) 
    ? (par->variant | HITRO_VARFLAG_ADAPTLINE) 
    : (par->variant & (~HITRO_VARFLAG_ADAPTLINE));

  /* changelog */
  par->set |= HITRO_SET_ADAPTLINE;

  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitro_set_use_adaptiveline() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_use_adaptiverectangle( struct unur_par *par, int adaptive )
     /*----------------------------------------------------------------------*/
     /* Enable/Disable adaptive bounding rectangle                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   adaptive ... whether adaptive rectangle is enabled                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
            
  /* we use a bit in variant */
  par->variant = (adaptive) 
    ? (par->variant | HITRO_VARFLAG_ADAPTRECT) 
    : (par->variant & (~HITRO_VARFLAG_ADAPTRECT));

  /* changelog */
  par->set |= HITRO_SET_ADAPTRECT;

  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitro_set_use_adaptiverectangle() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_use_boundingrectangle( struct unur_par *par, int rectangle )
     /*----------------------------------------------------------------------*/
     /* Whether to use entire bounding rectangle or just upper bound         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   rectangle ... whether using entire bounding rectangle is enabled   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
            
  /* we use a bit in variant */
  par->variant = (rectangle) 
    ? (par->variant | HITRO_VARFLAG_BOUNDRECT) 
    : (par->variant & (~HITRO_VARFLAG_BOUNDRECT));

  /* changelog */
  par->set |= HITRO_SET_BOUNDRECT;

  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitro_set_use_adaptiverectangle() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* set parameter r for ratio-of-uniforms                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   r    ... parameter r                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );

  /* check new parameter for generator */
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->r = r;

  /* changelog */
  par->set |= HITRO_SET_R;

  return UNUR_SUCCESS;

} /* end of unur_hitro_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_u( struct unur_par *par, const double *umin, const double *umax )
     /*----------------------------------------------------------------------*/
     /* u-boundaries for bounding rectangle                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   umin ... lower left boundary of rectangle                          */
     /*   umax ... upper right boundary of rectangle                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int d;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );

  /* check new parameter for generator */
  for (d = 0; d < par->distr->dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
    if (! (_unur_isfinite(umax[d]) && _unur_isfinite(umin[d])) ) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"rectangle not bounded");
      return UNUR_ERR_PAR_SET;
    }
  }

  /* set values */
  PAR->umin = umin;
  PAR->umax = umax;

  /* changelog */
  par->set |= HITRO_SET_U;

  /* ok */
  return UNUR_SUCCESS;

} /* end of unur_hitro_set_u() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_v( struct unur_par *par, double vmax )
     /*----------------------------------------------------------------------*/
     /* upper v-bound for bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   vmax ... upper bound of rectangle                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  if (! _unur_isfinite(vmax) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"rectangle not bounded");
    return UNUR_ERR_PAR_SET;
  }

  /* store values */
  PAR->vmax = vmax;

  /* changelog */
  par->set |= HITRO_SET_V;

  /* ok */
  return UNUR_SUCCESS;

} /* end of unur_hitro_set_v() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_adaptive_multiplier( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set multiplier for adaptive rectangles                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   factor   ... multiplier                                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, HITRO );

  /* check new parameter for generator */
  if (factor < HITRO_MIN_MULTIPLIER) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"multiplier too small (<= 1.0001)");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->adaptive_mult = factor;

  /* changelog */
  par->set |= HITRO_SET_ADAPTMULT;

  /* ok */
  return UNUR_SUCCESS;

} /* end of unur_hitro_set_adaptive_multiplier() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_startingpoint( struct unur_par *par, const double *x0)
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
  _unur_check_par_object( par, HITRO );

  /* store data */
  PAR->x0 = x0;

  /* changelog */
  par->set |= HITRO_SET_X0;

  /* ok */
  return UNUR_SUCCESS;
} /* end of unur_hitro_set_startingpoint() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_thinning( struct unur_par *par, int thinning )
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
  _unur_check_par_object( par, HITRO );

  /* check new parameter for generator */
  if (thinning < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"thinning < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->thinning = thinning;

  /* changelog */
  par->set |= HITRO_SET_THINNING;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitro_set_thinning() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_set_burnin( struct unur_par *par, int burnin )
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
  _unur_check_par_object( par, HITRO );

  /* check new parameter for generator */
  if (burnin < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"burnin < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->burnin = burnin;

  /* changelog */
  par->set |= HITRO_SET_BURNIN;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hitro_set_burnin() */

/*---------------------------------------------------------------------------*/

const double *
unur_hitro_get_state( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get current state of the HITRO chain                                 */
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
  if (gen->method != UNUR_METH_HITRO) {
    _unur_error(gen->genid, UNUR_ERR_GEN_INVALID,"");
    return NULL;
  }

  return GEN->state;
} /* end of unur_hitro_get_state() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_chg_state( struct unur_gen *gen, const double *state )
     /*----------------------------------------------------------------------*/
     /* chg current state of the HITRO chain                                 */
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
  _unur_check_gen_object( gen, HITRO, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( gen->genid, state, UNUR_ERR_NULL );

  /* check input */
  if ( ! _unur_hitro_vu_is_inside_region(gen,state) ) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"invalid state");
    return UNUR_ERR_PAR_SET;
  }

  /* copy state */
  memcpy( GEN->state, state, GEN->dim * sizeof(double));

  return UNUR_SUCCESS;
} /* end of unur_hitro_chg_state() */

/*---------------------------------------------------------------------------*/

int
unur_hitro_reset_state( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* reset current state of the HITRO chain to starting point             */
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
  _unur_check_gen_object( gen, HITRO, UNUR_ERR_GEN_INVALID );

  /* copy state */
  memcpy( GEN->state, GEN->x0, GEN->dim * sizeof(double));

  /* set state to start from */
  _unur_hitro_xy_to_vu(gen, GEN->x0, GEN->fx0/2., GEN->state );
  memcpy( GEN->vu, GEN->state, (GEN->dim + 1) * sizeof(double) );
  GEN->vumax[0] = pow(GEN->fx0, 1./(GEN->r * GEN->dim + 1.)) * (1. + DBL_EPSILON); 

  if (gen->variant & HITRO_VARIANT_COORD) GEN->coord = 0;

  return UNUR_SUCCESS;
} /* end of unur_hitro_reset_state() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_hitro_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_HITRO ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HITRO_PAR,NULL);

  /* check variant flags */
  /* ( we assume here that 
     par->variant & (HITRO_VARFLAG_ADAPTRECT | HITRO_VARFLAG_BOUNDRECT) 
     are set to 0 in unur_hitro_new() 
  */
  if ( par->variant & HITRO_VARIANT_COORD ) {
    /* coordinate direction sampling; not all variants are allowed */

    if (_unur_distr_cvec_has_boundeddomain(par->distr))
      /* we must have a bounded rectangular domain ... */
      par->variant |= HITRO_VARFLAG_BOUNDDOMAIN;
    else
      /* ... or a bounding rectangle for acceptance region */
      par->variant |= HITRO_VARFLAG_BOUNDRECT;

    /* we do not override the flag set by user for adaptive rectangle */
    if (!(par->set & HITRO_SET_ADAPTRECT) )
      par->variant |= HITRO_VARFLAG_ADAPTRECT;
  }

  /* create a new empty generator object */
  gen = _unur_hitro_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_hitro_debug_init_start(gen);
#endif

  /* we need a starting point for our chain */
  GEN->fx0 = PDF(GEN->x0);
  if ( (GEN->fx0 / 2.) <= 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"x0 not in support of PDF");
    _unur_hitro_free(gen); return NULL;
  }
  /* set state to start from */
  _unur_hitro_xy_to_vu(gen, GEN->x0, GEN->fx0/2., GEN->state );
  memcpy( GEN->vu, GEN->state, (GEN->dim + 1) * sizeof(double) );
  GEN->vumax[0] = pow(GEN->fx0, 1./(GEN->r * GEN->dim + 1.)) * (1. + DBL_EPSILON); 

  /* routines for sampling and destroying generator */
  if (gen->variant & HITRO_VARIANT_RANDOMDIR ) {
    /* we need an auxiliary generator for normal random variates */
    GEN_NORMAL = _unur_hitro_normalgen( gen );
    if ( GEN_NORMAL == NULL ) {
      _unur_hitro_free(gen); return NULL;
    }
  }

  /* computing required data about bounding rectangle */
  if ( !(gen->variant & HITRO_VARFLAG_ADAPTRECT) )
    if (_unur_hitro_rectangle(gen) != UNUR_SUCCESS ) {
      _unur_hitro_free(gen); return NULL;
    }
  /* else 
     we use an adaptive bounding rectangle --> nothing to do here */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_hitro_debug_init_finished(gen);
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

    for (burnin = GEN->burnin; burnin>0; --burnin)
      _unur_sample_vec(gen,X);

    /* restore thinning factor */
    GEN->thinning = thinning;
    /* free memory for random vector */
    free (X);
  }

  /* creation of generator object successfull */
  gen->status = UNUR_SUCCESS;

  /* o.k. */
  return gen;

} /* end of _unur_hitro_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_hitro_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HITRO_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_hitro_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_HITRO_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim;

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_hitro_getSAMPLE(gen);
  gen->destroy = _unur_hitro_free;
  gen->clone = _unur_hitro_clone;

  /* variant of sampling method */
  gen->variant = par->variant;

  /* copy parameters into generator object */
  GEN->thinning = PAR->thinning; /* thinning factor                          */
  GEN->burnin = PAR->burnin;     /* length of burnin                         */
  GEN->r = PAR->r;               /* parameter r for RoU                      */
  GEN->adaptive_mult = PAR->adaptive_mult; /* multiplier for adaptive rectangles */

  /* get center of the distribution */
  GEN->center = unur_distr_cvec_get_center(gen->distr);

  /* store starting point / get default starting point */
  GEN->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
  if (PAR->x0 == NULL)
    PAR->x0 = unur_distr_cvec_get_center(gen->distr);
  memcpy( GEN->x0, PAR->x0, GEN->dim * sizeof(double));

  /* bounding rectangle */
  GEN->vumin = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
  GEN->vumax = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
  /* bounding rectangle v-coordinate */
  GEN->vumin[0] = 0.;
  GEN->vumax[0] = (PAR->vmax > 0.) ? PAR->vmax : HITRO_START_UVMIN;  
  /* lower left and upper right vertex of bounding rectangle in u-hyperplane */
  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    if (PAR->umin && PAR->umax) {
      memcpy (GEN->vumin+1, PAR->umin, GEN->dim * sizeof(double) );
      memcpy (GEN->vumax+1, PAR->umax, GEN->dim * sizeof(double) );
    }
    else {
      for (i=1; i<GEN->dim+1; i++) GEN->vumin[i] = -HITRO_START_UVMIN;
      for (i=1; i<GEN->dim+1; i++) GEN->vumax[i] =  HITRO_START_UVMIN;
    }
  }
  /* else: we do not need the u-bounds */

  /* initialize remaining pointers */
  GEN->state = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
  GEN->x     = _unur_xmalloc( GEN->dim * sizeof(double) );
  GEN->vu    = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );

  /* allocate memory for random direction */
  GEN->direction = _unur_xmalloc( (1 + GEN->dim) * sizeof(double));

  /* defaults */
  GEN->coord = 0;           /* current coordinate of HITRO chain. */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_hitro_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_hitro_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hitro_clone( const struct unur_gen *gen )
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
#define CLONE         ((struct unur_hitro_gen*)clone->datap)

/*   int i; */
  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HITRO_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* get center of the distribution */
  CLONE->center = unur_distr_cvec_get_center(clone->distr);

  /* copy state */
  if (GEN->state) {
    CLONE->state = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
    memcpy( CLONE->state, GEN->state, (1 + GEN->dim) * sizeof(double) );
  }
  if (GEN->vumin) {
    CLONE->vumin = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
    memcpy( CLONE->vumin, GEN->vumin, (GEN->dim+1) * sizeof(double) );
  }
  if (GEN->vumax) {
    CLONE->vumax = _unur_xmalloc( (GEN->dim+1) * sizeof(double) );
    memcpy( CLONE->vumax, GEN->vumax, (GEN->dim+1) * sizeof(double) );
  }

  /* copy starting points */
  if (GEN->x0) {
    CLONE->x0 = _unur_xmalloc( GEN->dim * sizeof(double));
    memcpy( CLONE->x0, GEN->x0, GEN->dim * sizeof(double));
  }

  /* working array */
  if (GEN->x) {
    CLONE->x = _unur_xmalloc( GEN->dim * sizeof(double) );
    memcpy( CLONE->x, GEN->x, GEN->dim * sizeof(double) );
  }
  if (GEN->vu) {
    CLONE->vu = _unur_xmalloc( (1 + GEN->dim) * sizeof(double) );
    memcpy( CLONE->vu, GEN->vu, (1 + GEN->dim) * sizeof(double) );
  }
  if (GEN->direction) {
    CLONE->direction = _unur_xmalloc( (1 + GEN->dim) * sizeof(double));
    memcpy( CLONE->direction, GEN->direction, (1 + GEN->dim) * sizeof(double));
  }

  return clone;

#undef CLONE
} /* end of _unur_hitro_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_hitro_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_HITRO ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_hitro_debug_free(gen);
#endif

  /* free vectors */
  if (GEN->state) free (GEN->state);
  if (GEN->x0) free (GEN->x0);
  if (GEN->x) free (GEN->x);
  if (GEN->vu) free (GEN->vu);
  if (GEN->direction) free (GEN->direction);
  if (GEN->vumin) free (GEN->vumin);
  if (GEN->vumax) free (GEN->vumax);

  _unur_generic_free(gen);

} /* end of _unur_hitro_free() */

/*****************************************************************************/

int
_unur_hitro_coord_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  int thinning;
  double lmin, lmax, lmid;   /* l.h.s., r.h.s. endpoint and midpoint of line segment */
  double *vuaux;  /* pointer to auxiliary working array */
  int coord;      /* direction */
  double U;       /* uniform random number */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_HITRO_GEN,UNUR_ERR_COOKIE);

  /* pointer to auxiliary working array */
  vuaux = GEN->vu;

  for (thinning = GEN->thinning; thinning > 0; --thinning) {

    /* update coordinate direction */
    coord = GEN->coord = (GEN->coord + 1) % (GEN->dim + 1);

    /* --- l.h.s. and r.h.s. endpoints of line segment --- */
    if (! (gen->variant & HITRO_VARFLAG_BOUNDDOMAIN) || coord == 0) {
      /* no bounded domain or v coordinate */
      lmin = GEN->vumin[coord];
      lmax = GEN->vumax[coord];
    }
    else {
      /* bounded domain; we look at one of the u-coordinates */
      int k = coord-1;
      double *domain = DISTR.domainrect;
      lmin = _unur_hitro_xv_to_u(gen, domain[2*k], vuaux[0], k );
      lmax = _unur_hitro_xv_to_u(gen, domain[2*k+1], vuaux[0], k );
      if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
	/* we also have a bounded domain; look whether this leads to shorter intervals */
	lmin = _unur_max(lmin,GEN->vumin[coord]);
	lmax = _unur_min(lmax,GEN->vumax[coord]);
      }
    }

    /* --- adaptive bounding rectangle --- */
    if ( gen->variant & HITRO_VARFLAG_ADAPTRECT ) {
      lmid = 0.5 * (lmin + lmax);
      /* check whether endpoint is outside of region */
      vuaux[coord] = lmax;
      while ( _unur_hitro_vu_is_inside_region(gen,vuaux) ) {
	lmax = lmid + (lmax-lmid) * GEN->adaptive_mult;
	GEN->vumax[coord] = vuaux[coord] = lmax;
      }
      vuaux[coord] = lmin;
      /* no update of lmin in case of coord == 0 since v-coordinate always > 0 */
      while ( coord!=0 && _unur_hitro_vu_is_inside_region(gen,vuaux) ) {
	lmin = lmid + (lmin-lmid) * GEN->adaptive_mult;
	GEN->vumin[coord] = vuaux[coord] = lmin;
      }
    }

    /* --- new coordinate for next point --- */
    while (1) {
      U = _unur_call_urng(gen->urng);
      vuaux[coord] = U * lmin + (1.-U) * lmax;
      if (_unur_hitro_vu_is_inside_region(gen,vuaux) )
	break;

      /* --- adaptive line sampling --- */
      if ( gen->variant & HITRO_VARFLAG_ADAPTLINE ) {
	if (GEN->state[coord] < vuaux[coord])
	  lmax = vuaux[coord];
	else
	  lmin = vuaux[coord];
      }
    }

    /* store new state */
    GEN->state[coord] = vuaux[coord];
  }
    
  /* transform current state into point in original scale */
  _unur_hitro_vu_to_x( gen, GEN->state, vec );
  
  return UNUR_SUCCESS;

} /* end of _unur_hitro_coord_sample_cvec() */

/*---------------------------------------------------------------------------*/

int
_unur_hitro_randomdir_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define new_point(ll)  { int j; for (j=0;j<dim+1;j++) vuaux[j] = GEN->state[j]+(ll)*GEN->direction[j]; }

  int thinning;
  int i, d, k;
  double lambda, lb[2];  /* line segments, endpoints */
  double *vuaux;       /* pointer to auxiliary working array */
  double U;
  int update;
  int dim = GEN->dim;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_HITRO_GEN,UNUR_ERR_COOKIE);

  /* entire bounding rectangle or just covering plate ? */
  d = (gen->variant & HITRO_VARFLAG_BOUNDRECT) ? dim+1 : 1;

  /* pointer to auxiliary working array */
  vuaux = GEN->vu;

  for (thinning = GEN->thinning; thinning > 0; --thinning) {

    /* new random direction */
    _unur_hitro_random_unitvector( gen, GEN->direction );
    
    /* --- l.h.s. and r.h.s. endpoint of line segment --- */
    lb[1] = INFINITY;  
    lb[0] = -INFINITY;
    for (i=0; i<d; i++) {
      lambda = (GEN->vumin[i] - GEN->state[i]) / GEN->direction[i];
      if (lambda>0 && lambda<lb[1]) lb[1] = lambda;
      if (lambda<0 && lambda>lb[0]) lb[0] = lambda;
      lambda = (GEN->vumax[i] - GEN->state[i]) / GEN->direction[i];
      if (lambda>0 && lambda<lb[1]) lb[1] = lambda;
      if (lambda<0 && lambda>lb[0]) lb[0] = lambda;
    }
    if (! (_unur_isfinite(lb[0]) && _unur_isfinite(lb[1])) ) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"line segment not bounded, try again");
      continue;
    }

    /* --- adaptive bounding rectangle --- */
    if ( gen->variant & HITRO_VARFLAG_ADAPTRECT ) {
      /* check whether endpoint is outside of region */
      for (k=0; k<2; k++) {
	update = FALSE;
	while (1) {
	  new_point(lb[k]);
	  if (! _unur_hitro_vu_is_inside_region(gen,vuaux) )
	    break;
	  update = TRUE;
	  lb[k] *= GEN->adaptive_mult;
	}
	if (update) {
	  /* we have to update vumin and vmax */
	  new_point(lb[k]);
	  for (i=0; i<d; i++) {
	    if (vuaux[i] < GEN->vumin[i] && i!=0) GEN->vumin[i] =  vuaux[i];
	    if (vuaux[i] > GEN->vumax[i])         GEN->vumax[i] =  vuaux[i];
	    /* remark: vmin is always 0 */
	  }			
	}
      }
    }		

    /* --- new coordinate for next point --- */
    while (1) {
      U = _unur_call_urng(gen->urng); 
      lambda = U * lb[0] + (1.-U) * lb[1];
      new_point(lambda);
      if (_unur_hitro_vu_is_inside_region(gen,vuaux) )
	break;

      /* --- adaptive line sampling --- */
      if ( gen->variant & HITRO_VARFLAG_ADAPTLINE ) {
	if (lambda < 0) lb[0] = lambda;
	else            lb[1] = lambda;
      }
    }

    /* store new state */
    memcpy( GEN->state, vuaux, (dim+1)*sizeof(double) );
  }
    
  /* transform current state into point in original scale */
  _unur_hitro_vu_to_x( gen, GEN->state, vec );
  
  return UNUR_SUCCESS;

#undef new_point

} /* end of _unur_hitro_randomdir_sample_cvec() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_hitro_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute bounding rectangle                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  int d; /* index used in dimension loops (0 <= d < dim) */
  struct MROU_RECTANGLE *rr;

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_HITRO_GEN, UNUR_ERR_COOKIE );

  /* Boundary rectangle is already set ? */
  if ((gen->set & HITRO_SET_U) && (gen->set & HITRO_SET_V))
    return UNUR_SUCCESS;

  /* Allocating and filling mrou_rectangle struct */
  rr = _unur_mrou_rectangle_new();

  rr->distr  = gen->distr;
  rr->dim    = GEN->dim;
  rr->umin   = GEN->vumin+1;
  rr->umax   = GEN->vumax+1;
  rr->r      = GEN->r;
  rr->center = GEN->center;
  rr->genid  = gen->genid;
  rr->bounding_rectangle = 
    ( (gen->variant & HITRO_VARFLAG_BOUNDRECT) && !(gen->set & HITRO_SET_U) )
    ? 1 : 0;
  
  /* calculate bounding rectangle */
  if ( _unur_mrou_rectangle_compute(rr) != UNUR_SUCCESS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"Cannot compute bounding rectangle, try adaptive");
    gen->variant &= HITRO_VARFLAG_ADAPTRECT;
    free(rr); return UNUR_ERR_GEN_CONDITION;
  }

  if (!(gen->set & HITRO_SET_V)) {
    /* user has not provided any upper bound for v */
    GEN->vumax[0] = rr->vmax;
  }

  if (rr->bounding_rectangle) {
    /* user has not provided required bounds for u */
    for (d=0; d<GEN->dim; d++) GEN->vumin[d+1] = rr->umin[d];
    for (d=0; d<GEN->dim; d++) GEN->vumax[d+1] = rr->umax[d];
  }

  /* free working space */
  free(rr);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_hitro_rectangle() */

/*---------------------------------------------------------------------------*/

double
_unur_hitro_xv_to_u( const struct unur_gen *gen, double x, double v, int k )
     /*----------------------------------------------------------------------*/
     /* transform point where we are given the v-coordinate (in RoU scale)   */
     /* and the k-th x-coordinate (in xy scale)                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   x     ... x coordinates of point (in original scale)               */
     /*   v     ... v coordinate of point (in RoU scale)                     */
     /*   k     ... index of coordinate                                      */ 
     /*----------------------------------------------------------------------*/
{
  if (_unur_isone(GEN->r)) 
    return (x - GEN->center[k]) * v;
  else
    return (x - GEN->center[k]) * pow(v,GEN->r) ;
} /* end of _unur_hitro_x_to_u() */

/*---------------------------------------------------------------------------*/

void 
_unur_hitro_xy_to_vu( const struct unur_gen *gen, const double *x, double y, double *vu )
     /*----------------------------------------------------------------------*/
     /* transform point with coordinates xy in the original scale into       */
     /* point in RoU-scale with coordinates vu.                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   x    ... x coordinates of point                                    */
     /*   y    ... y coordinate of point                                     */
     /*   vu   ... coordinates of transformed point                          */
     /*----------------------------------------------------------------------*/
{
  int d;
  double v;
  double *u = vu+1;
   
  vu[0] = v = pow(y, 1./(GEN->r * GEN->dim + 1.));

  if (_unur_isone(GEN->r)) 
    for (d=0; d<GEN->dim; d++)  u[d] = (x[d] - GEN->center[d]) * v;
  else
    for (d=0; d<GEN->dim; d++)  u[d] = (x[d] - GEN->center[d]) * pow(v,GEN->r) ;

} /* end of _unur_hitro_xy_to_vu() */

/*---------------------------------------------------------------------------*/

void 
_unur_hitro_vu_to_x( const struct unur_gen *gen, const double *vu, double *x )
     /*----------------------------------------------------------------------*/
     /* transform point with coordinates vu in RoU-scale into original       */
     /* scale with coordinates xy. The y coordinate is not of interest here  */
     /* and thus omitted.                                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   vu   ... coordinates of transformed point                          */
     /*   x    ... x coordinates of point                                    */
     /*----------------------------------------------------------------------*/
{
  int d;
  double v = vu[0];
  const double *u = vu+1;

  if (v<=0.) {
    /* we are outside of the domain --> return 0 */
    for (d=0; d<GEN->dim; d++)  x[d] = 0.;
    return;
  }
    
  if (_unur_isone(GEN->r))
    for (d=0; d<GEN->dim; d++)  x[d] = u[d]/v + GEN->center[d];
  else
    for (d=0; d<GEN->dim; d++)  x[d] = u[d]/pow(v,GEN->r) + GEN->center[d];

} /* end of _unur_hitro_vu_to_x() */

/*---------------------------------------------------------------------------*/

int
_unur_hitro_vu_is_inside_region( const struct unur_gen *gen, const double *vu )
     /*----------------------------------------------------------------------*/
     /* check whether point is inside RoU acceptance region                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   vu   ... point in uv space                                         */
     /*----------------------------------------------------------------------*/
{
  double y;

  /* v-coordinate */
  double v = vu[0];

  /* transform into original scale */
  _unur_hitro_vu_to_x( gen, vu, GEN->x );

  /* PDF at x */
  y = PDF(GEN->x);

  /* we are outside of domain if x is not in support of PDF or */
  /* v is non-positive.                                        */
  if (y <= 0. || v <= 0.) return FALSE;

  /* now check whether point is in region below PDF */
  return ( (v < pow(y,1./(GEN->r * GEN->dim + 1.))) ? TRUE : FALSE );

} /* end of _unur_hitro_vu_is_inside_region() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hitro_normalgen( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create a normal random variate generator                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to HITRO generator object                    */
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

} /* end of _unur_hitro_normalgen() */

/*---------------------------------------------------------------------------*/

void
_unur_hitro_random_unitvector( struct unur_gen *gen, double *direction )
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
    for (i=0; i<GEN->dim+1; i++)
      direction[i] = unur_sample_cont(GEN_NORMAL);
    /* normalize direction vector */
    _unur_vector_normalize(GEN->dim+1, direction);

    /* there is an extremely small change that direction is the null before
       normalizing. In this case non of its coordinates are finite. */
  } while (!_unur_isfinite(direction[0]));

} /* end of _unur_hitro_random_unitvector() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_hitro_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = HITRO (Markov Chain - HITRO sampler)\n",gen->genid);
  fprintf(LOG,"%s: variant = ",gen->genid);
  switch (gen->variant & HITRO_VARMASK_VARIANT) {
  case HITRO_VARIANT_COORD:
    fprintf(LOG,"coordinate sampling [default]\n"); break;
  case HITRO_VARIANT_RANDOMDIR:
    fprintf(LOG,"random direction sampling\n"); break;
  }
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: r = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,HITRO_SET_R);
  fprintf(LOG,"\n%s: adaptive line sampling: %s",gen->genid,
	  (gen->variant&HITRO_VARFLAG_ADAPTLINE)?"on":"off");
  _unur_print_if_default(gen,HITRO_SET_ADAPTLINE);
  fprintf(LOG,"\n%s: use entire bounding rectangle: %s",gen->genid,
	  (gen->variant&HITRO_VARFLAG_BOUNDRECT)?"on":"off");
  _unur_print_if_default(gen,HITRO_SET_BOUNDRECT);
  fprintf(LOG,"\n%s: adaptive bounding rectangle: %s",gen->genid,
	  (gen->variant&HITRO_VARFLAG_ADAPTRECT)?"on":"off");
  _unur_print_if_default(gen,HITRO_SET_ADAPTRECT);
  if (gen->variant&HITRO_VARFLAG_ADAPTRECT) {
    fprintf(LOG,"\n%s:\tmultiplier = %g",gen->genid,GEN->adaptive_mult);
    _unur_print_if_default(gen,HITRO_SET_ADAPTMULT);
  }
  fprintf(LOG,"\n%s: use domain of distribution: %s\n",gen->genid,
	  (gen->variant&HITRO_VARFLAG_BOUNDDOMAIN)?"on":"off");
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  switch (gen->variant & HITRO_VARMASK_VARIANT) {
  case HITRO_VARIANT_COORD:
    fprintf(LOG,"%s: sampling routine = _unur_hitro_coord_sample_cvec()\n",gen->genid);
    break;
  case HITRO_VARIANT_RANDOMDIR:
    fprintf(LOG,"%s: sampling routine = _unur_hitro_randomdir_sample_cvec()\n",gen->genid);
    break;
  }

  fprintf(LOG,"%s: thinning = %d",gen->genid,GEN->thinning);
  _unur_print_if_default(gen,HITRO_SET_THINNING);
  fprintf(LOG,"\n%s: burn-in = %d",gen->genid,GEN->burnin);
  _unur_print_if_default(gen,HITRO_SET_BURNIN);
  fprintf(LOG,"\n%s:\n",gen->genid);
  _unur_matrix_print_vector( GEN->dim, GEN->x0, "starting point = ", LOG, gen->genid, "\t   ");

  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_hitro_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_hitro_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    fprintf(LOG,"%s: bounding rectangle%s:\n",gen->genid,
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? " [start for adaptive rectangle]" : "" );
    fprintf(LOG,"%s: vmax = %g\n",gen->genid, GEN->vumax[0]);
    _unur_matrix_print_vector( GEN->dim, GEN->vumin+1, "umin =", LOG, gen->genid, "\t   ");
    _unur_matrix_print_vector( GEN->dim, GEN->vumax+1, "umax =", LOG, gen->genid, "\t   ");
  }
  else {
    fprintf(LOG,"%s: upper bound vmax = %g %s\n",gen->genid, GEN->vumax[0],
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[start for adaptive bound]" : "" );
  }

  _unur_matrix_print_vector( GEN->dim+1, GEN->state, "starting state = ", LOG, gen->genid, "\t   ");
  
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);

  fflush(LOG);

} /* end of _unur_hitro_debug_init_finished() */

/*---------------------------------------------------------------------------*/

void
_unur_hitro_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HITRO_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->status == UNUR_SUCCESS) {
    fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
  }
  else {
    fprintf(LOG,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
  }
  fprintf(LOG,"%s:\n",gen->genid);

  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    fprintf(LOG,"%s: bounding rectangle%s:\n",gen->genid,
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? " [adaptive]" : "" );
    fprintf(LOG,"%s: vmax = %g\n",gen->genid, GEN->vumax[0]);
    _unur_matrix_print_vector( GEN->dim, GEN->vumin+1, "umin =", LOG, gen->genid, "\t   ");
    _unur_matrix_print_vector( GEN->dim, GEN->vumax+1, "umax =", LOG, gen->genid, "\t   ");
  }
  else {
    fprintf(LOG,"%s: upper bound vmax = %g %s\n",gen->genid, GEN->vumax[0],
	    (gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[adaptive]" : "" );
  }
  
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_hitro_debug_free() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_hitro_info( struct unur_gen *gen, int help )
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
  int i;
  double rc;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_distr_cvec_info_domain(gen);

  if ( distr->set & UNUR_DISTR_SET_MODE ) {
    _unur_string_append(info,"   mode      = ");
    _unur_distr_info_vector( gen, DISTR.mode, GEN->dim);
  }
  _unur_string_append(info,"\n");

  _unur_string_append(info,"   center    = ");
  _unur_distr_info_vector( gen, GEN->center, GEN->dim);
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
  _unur_string_append(info,"method: HITRO (HIT-and-run sampler with Ratio-Of-uniforms [MCMC])\n");
  _unur_string_append(info,"   variant = %s\n",
		      ((gen->variant & HITRO_VARMASK_VARIANT)==HITRO_VARIANT_COORD)
		      ? "coordinate sampling [default]" : "random direction sampling");

  _unur_string_append(info,"   r = %g\n", GEN->r);
  _unur_string_append(info,"   thinning = %d\n", GEN->thinning);

  _unur_string_append(info,"   adaptive line sampling = %s\n", 
		      (gen->variant&HITRO_VARFLAG_ADAPTLINE)?"on":"off");

  _unur_string_append(info,"   use entire bounding rectangle = %s\n",
		      (gen->variant&HITRO_VARFLAG_BOUNDRECT)?"on":"off");

  if (gen->variant&HITRO_VARFLAG_ADAPTRECT)
    _unur_string_append(info,"   adaptive bounding rectangle = on  [multiplier = %g]\n",
			GEN->adaptive_mult);
  else
    _unur_string_append(info,"   adaptive bounding rectangle = off\n");

  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");

  rc = unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize);

  if (gen->variant & HITRO_VARFLAG_BOUNDRECT) {
    _unur_string_append(info,"   bounding rectangle %s= ",
			(gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[adaptive] " : "" );
    for (i=0; i<GEN->dim; i++)
      _unur_string_append(info,"%s(%g,%g)", i?"x":"", GEN->vumin[i+1], GEN->vumax[i+1]);
    _unur_string_append(info," x (0,%g)\n", GEN->vumax[0]);
  }
  else {
    _unur_string_append(info,"   upper bound vmax = %g %s\n", GEN->vumax[0],
			(gen->variant & HITRO_VARFLAG_ADAPTRECT) ? "[adaptive]" : "" );
  }

  _unur_string_append(info,"   rejection constant =  %.2f  [approx.]\n", rc);
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");

    switch (gen->variant & HITRO_VARMASK_VARIANT) {
    case HITRO_VARIANT_COORD:
      _unur_string_append(info,"   variant_coordinate  [default]\n"); break;
    case HITRO_VARIANT_RANDOMDIR:
      _unur_string_append(info,"   variant_random_direction\n"); break;
    }

    _unur_string_append(info,"   r = %g  %s\n", GEN->r,
 			(gen->set & HITRO_SET_R) ? "" : "[default]");

    _unur_string_append(info,"   adaptiveline = %s  %s\n", 
			(gen->variant&HITRO_VARFLAG_ADAPTLINE)?"on":"off",
 			(gen->set & HITRO_SET_ADAPTLINE) ? "" : "[default]");

    _unur_string_append(info,"   boundingrectangle = %s  %s\n",
			(gen->variant&HITRO_VARFLAG_BOUNDRECT)?"on":"off",
 			(gen->set & HITRO_SET_BOUNDRECT) ? "" : "[default]");

    _unur_string_append(info,"   adaptiverectangle = %s  %s\n", 
			(gen->variant&HITRO_VARFLAG_ADAPTRECT)?"on":"off",
 			(gen->set & HITRO_SET_ADAPTRECT) ? "" : "[default]");

    if (gen->variant&HITRO_VARFLAG_ADAPTRECT)
      _unur_string_append(info,"   adaptive_multiplier = %g  %s\n", 
			  GEN->adaptive_mult,
			  (gen->set & HITRO_SET_ADAPTMULT) ? "" : "[default]");

   _unur_string_append(info,"   thinning = %d  %s\n", GEN->thinning,
 			(gen->set & HITRO_SET_THINNING) ? "" : "[default]");
   _unur_string_append(info,"   burnin = %d  %s\n", GEN->burnin,
 			(gen->set & HITRO_SET_THINNING) ? "" : "[default]");

    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_hitro_set_v( UNUR_PAR *parameters, double vmax );
       int unur_hitro_set_u( UNUR_PAR *parameters, const double *umin, const double *umax );
       int unur_hitro_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
    */
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_hitro_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
