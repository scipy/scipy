/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      nrou.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    naive ratio-of-uniforms method                               *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and (optionally) a bounding rectangle for the acceptance   *
 *      region.                                                              *
 *      Produce a value x consistent with its density                        *
 *      The bounding rectangle is computed numerically if it is not given.   * 
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *   OPTIONAL:                                                               *
 *      mode of the density                                                  *
 *      bounding rectangle of acceptance region                              *
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
 *   [1] Kinderman, A.J. and Monahan, F.J. (1977): Computer generation of    *
 *       random variables using the ratio of uniform deviates,               *
 *       ACM Trans. Math. Software 3(3), pp. 257--260.                       *
 *                                                                           *
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The ratio-of-uniforms method introduced in [1] is a flexible method that  *
 * is based on the following theorem:                                        *
 *                                                                           *
 * THEOREM:                                                                  *
 *    Let X be a random variable with density function f(x) = g(x) / G,      *
 *    where g(x) is a positive integrable function with support (x_0,x_1)    *
 *    not necessarily finite and G = integral g(x) dx.                       *
 *    If (U,V) is uniformly distributed in                                   *
 *       A = {(u,v): 0 < v <= sqrt(g(u/v)), x_0 < u/v < x_1},                *
 *    then X = V/U has probability density function f(x).                    *
 *                                                                           *
 * Generating point (U,V) uniformly distributed in A is done by rejection    *
 * from an enveloping region, usually from the minimal bounding rectangle.   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <utils/fmax_source.h>
#include <utils/unur_fp_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "nrou.h"
#include "nrou_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants:                                                                */

/* Scaling factor for the computed minimum bounding rectangle.               */
/* The computed rectangle  (0, vmax)x(umin[d], umax[d]) is scaled by this    */
/* factor, i.e. :                                                            */
/* vmax = vmax * ( 1+ NROU_RECT_SCALING)                                     */
/* umin = umin - (umax-umin)*NROU_RECT_SCALING/2.                            */
/* umax = umax + (umax-umin)*NROU_RECT_SCALING/2.                            */
#define NROU_RECT_SCALING (1.e-4)

/* The minimum bounding rectangle is computed by an numerical search         */
/* algorithm. In case of unbounded domain we need a very large bound         */
/* number as starting point:                                                 */
#define BD_MAX    (DBL_MAX/1000.)

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define NROU_VARFLAG_VERIFY   0x002u   /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define NROU_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NROU_SET_U       0x001u     /* set u values of bounding rectangle    */
#define NROU_SET_V       0x002u     /* set v values of bounding rectangle    */
#define NROU_SET_CENTER  0x004u     /* set center of distribution            */
#define NROU_SET_R       0x008u     /* set r-parameter                       */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NROU"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_nrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_nrou_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_nrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_nrou_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_nrou_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_nrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_nrou_sample( struct unur_gen *gen );
static double _unur_nrou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static double _unur_aux_bound_umax(double x, void *p);
static double _unur_aux_bound_umin(double x, void *p);
/*---------------------------------------------------------------------------*/
/* auxiliary functions to be used in the calculation of umin/umax            */
/*---------------------------------------------------------------------------*/

static int _unur_nrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_nrou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_nrou_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_nrou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_nrou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

#define _unur_nrou_getSAMPLE(gen) \
   ( ((gen)->variant & NROU_VARFLAG_VERIFY) \
     ? _unur_nrou_sample_check : _unur_nrou_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_nrou_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_nrou_par) );
  COOKIE_SET(par,CK_NROU_PAR);

  /* copy input */
  par->distr    = distr;          /* pointer to distribution object          */

  /* set default values */
  PAR->umin      = 0.;         /* u-boundary of bounding rectangle (unknown) */
  PAR->umax      = 0.;         /* u-boundary of bounding rectangle (unknown) */
  PAR->vmax      = 0.;         /* v-boundary of bounding rectangle (unknown) */
  PAR->center    = 0.;         /* center of distribution (default: 0)        */
  PAR->r         = 1.;         /* r-parameter of the generalized method      */
  
  par->method   = UNUR_METH_NROU;     /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_nrou_init;

  return par;

} /* end of unur_nrou_new() */

/*****************************************************************************/

int
unur_nrou_set_u( struct unur_par *par, double umin, double umax )
     /*----------------------------------------------------------------------*/
     /* Sets left and right u-boundary of bounding rectangle.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   umin ... left boundary of rectangle                                */
     /*   umax ... right boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* check new parameter for generator */
  if (!_unur_FP_greater(umax,umin)) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR->umin = umin;
  PAR->umax = umax;

  /* changelog */
  par->set |= NROU_SET_U;

  return UNUR_SUCCESS;

} /* end of unur_nrou_set_u() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_v( struct unur_par *par, double vmax )
     /*----------------------------------------------------------------------*/
     /* Sets upper v-boundary of bounding rectangle.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   vmax ... upper boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR->vmax = vmax;

  /* changelog */
  par->set |= NROU_SET_V;

  return UNUR_SUCCESS;

} /* end of unur_nrou_set_v() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_center( struct unur_par *par, double center )
     /*----------------------------------------------------------------------*/
     /* Set the center (approximate mode) of the PDF.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   center ... center of distribution                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* store data */
  PAR->center = center;

  /* changelog */
  par->set |= NROU_SET_CENTER;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_nrou_set_center() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* Sets r-parameter of the generalized ratio-of-uniforms method         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   r ... r-parameter                                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* check new parameter for generator */
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r<=0");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR->r = r;

  /* changelog */
  par->set |= NROU_SET_R;

  return UNUR_SUCCESS;

} /* end of unur_nrou_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, NROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | NROU_VARFLAG_VERIFY) : (par->variant & (~NROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_nrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_nrou_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, NROU, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= NROU_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~NROU_VARFLAG_VERIFY;

  SAMPLE = _unur_nrou_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_nrou_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_nrou_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params ... pointer to paramters for building generator object      */
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
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_NROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_nrou_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_nrou_check_par(gen) != UNUR_SUCCESS) {
    _unur_nrou_free(gen); return NULL;
  }

  /* compute bounding rectangle */
  if (_unur_nrou_rectangle(gen)!=UNUR_SUCCESS) {
    _unur_error(gen->genid , UNUR_ERR_GEN_CONDITION, "Cannot compute bounding rectangle");  
    _unur_nrou_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_nrou_debug_init(gen);
#endif

  return gen;

} /* end of _unur_nrou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_nrou_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int rcode;

  /* mark U and V as unknown */
  gen->set &= ~(NROU_SET_V | NROU_SET_U);

  /* check parameters */
  if ( (rcode = _unur_nrou_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* (re)set sampling routine */
  SAMPLE = _unur_nrou_getSAMPLE(gen);

  /* compute bounding rectangle */
  return _unur_nrou_rectangle(gen);
} /* end of _unur_nrou_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_nrou_create( struct unur_par *par )
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

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_nrou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_NROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_nrou_getSAMPLE(gen);
  gen->destroy = _unur_nrou_free;
  gen->clone = _unur_nrou_clone;
  gen->reinit = _unur_nrou_reinit;

  /* copy some parameters into generator object */
  GEN->umin  = PAR->umin;             /* left u-boundary of bounding rectangle */
  GEN->umax  = PAR->umax;             /* right u-boundary of bounding rectangle */
  GEN->vmax  = PAR->vmax;             /* upper v-boundary of bounding rectangle */
  GEN->center = PAR->center;          /* center of distribution */
  GEN->r = PAR->r;                    /* r-parameter of the generalized rou-method */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_nrou_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_nrou_create() */

/*---------------------------------------------------------------------------*/

int
_unur_nrou_check_par( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check parameters of given distribution and method                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{

  /* check for required data: center */
  if (!(gen->set & NROU_SET_CENTER))
    /* center not set via unur_nrou_set_center */
    /* use center of distribution instead.     */
    GEN->center = unur_distr_cont_get_center(gen->distr) ;

  return UNUR_SUCCESS;
} /* end of _unur_nrou_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_nrou_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_nrou_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_nrou_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_nrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_nrou_free() */

/*****************************************************************************/

double
_unur_nrou_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NROU_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    U = GEN->umin + _unur_call_urng(gen->urng) * (GEN->umax - GEN->umin);

    /* compute X */
    if (_unur_isone(GEN->r))
      X = U/V + GEN->center;
    else
      X = U/pow(V,GEN->r) + GEN->center;

    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;

    /* accept or reject */
    if (_unur_isone(GEN->r)) {
      /* normal rou-method with square-root */
      if (V*V <= PDF(X)) 
        return X;
    }
    else {
      /* generalized rou-method with pow-function */
      if (V <= pow(PDF(X), 1./(1.+GEN->r)) )
        return X;
    }
  }

} /* end of _unur_nrou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_nrou_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X,fx,sfx,xfx;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_NROU_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on rectangle */
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    U = GEN->umin + _unur_call_urng(gen->urng) * (GEN->umax - GEN->umin);
    
    /* compute x */
    if (_unur_isone(GEN->r))
      X = U/V + GEN->center;
    else
      X = U/pow(V,GEN->r) + GEN->center;
    
    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;
    
    /* evaluate PDF */
    fx = PDF(X);
    
    /* a point on the boundary of the region of acceptance
       has the coordinates ( (X-center) * (fx)^(r/(1+r)), (fx)^(1/(1+r)) ). */
    if (_unur_isone(GEN->r)) {
      /* normal rou-method with square-root */
      sfx = sqrt(fx);
      xfx = (X-GEN->center) * sfx;
    }
    else {
      /* generalized rou-method with pow-function */
      sfx = pow(fx, 1./(1.+GEN->r));
      xfx = (X-GEN->center) * pow(fx, GEN->r/(1.+GEN->r));
    }
    
    /* check hat */
    if ( ( sfx > (1.+DBL_EPSILON) * GEN->vmax )   /* avoid roundoff error with FP registers */
	 || (xfx < (1.+UNUR_EPSILON) * GEN->umin) 
	 || (xfx > (1.+UNUR_EPSILON) * GEN->umax) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
    
    /* accept or reject */
    if (_unur_isone(GEN->r)) {
      /* normal rou-method with square-root */
      if (V*V <= PDF(X))
        return X;
    }
    else {
      /* generalized rou-method with pow-function */
      if (V <= pow(PDF(X), 1./(1.+GEN->r)) )
        return X;
    }
  }

} /* end of _unur_nrou_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

double
_unur_aux_bound_umax(double x, void *p) 
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	
{
  struct unur_gen *gen;
  gen = p; /* typecast from void* to unur_gen* */
  
  if (_unur_isone(GEN->r)) 
    return (x-GEN->center) * sqrt( _unur_cont_PDF((x),(gen->distr)) );

  else
    return (x-GEN->center) * pow( _unur_cont_PDF((x),(gen->distr)),
				 GEN->r / (1.+ GEN->r) );
}

/*---------------------------------------------------------------------------*/

double
_unur_aux_bound_umin(double x, void *p)
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/	
{
  return (- _unur_aux_bound_umax(x,p)) ;
}

/*---------------------------------------------------------------------------*/

int
_unur_nrou_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_funct_generic faux; /* function to be minimized/maximized    */
  double mode; /* position of the distribution mode */
  double x, cx, sx, bx; /* parameters to be used in min/max search */ 

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_NROU_GEN, UNUR_ERR_COOKIE );

  /* boundary rectangle is already set */
  if ((gen->set & NROU_SET_U) && (gen->set & NROU_SET_V)) {
    return UNUR_SUCCESS;
  }

  /* starting point in min/max algorithm */
  cx=GEN->center;  
 
  /* --------------------------------------------------------------------- */
 
  /* calculation of vmax */
  if (!(gen->set & NROU_SET_V)) {
    /* user has not provided any upper bound for v */

    /* get (optional) mode if present or find it numerically */
    mode = unur_distr_cont_get_mode(gen->distr);

    if (!_unur_isfinite(mode))
      return UNUR_ERR_GENERIC;
      
    /* setting vmax to be (f(mode))^(1/(1+r)) */
    GEN->vmax = pow(PDF(mode), 1./(1.+GEN->r));

    /* additional scaling of boundary rectangle */
    GEN->vmax = GEN->vmax * ( 1. + NROU_RECT_SCALING);

    /* check for bounded rectangle */
    if (! _unur_isfinite(GEN->vmax)) {
      _unur_error(gen->genid , UNUR_ERR_GENERIC, "vmax not finite");  
      return UNUR_ERR_GENERIC;
    }
    
  }

  /* --------------------------------------------------------------------- */
  
  /* calculation of umin and umax */
  if (!(gen->set & NROU_SET_U)) {

    /* umin: */
    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umin;
    faux.params = gen;

    /* calculating start point for extremum search routine */
    sx = _unur_isfinite(DISTR.BD_LEFT) ? (cx+DISTR.BD_LEFT)/2. : (cx-1.); 
    bx = _unur_isfinite(DISTR.BD_LEFT) ? DISTR.BD_LEFT : (-BD_MAX);

    x = (_unur_FP_same(DISTR.BD_LEFT,cx)) 
      ? cx : _unur_util_find_max(faux, bx, cx, sx);
          
    while (!_unur_isfinite(x) && (fabs(bx) >= UNUR_EPSILON) ) { 
       /* _unur_util_find_max() could not yet find a suitable extremum */
       /* trying with a sequence of intervals with decreasing length   */
       bx = bx/10.; sx = bx/2.;  
       x = _unur_util_find_max(faux, bx, cx, sx);
    }
    /* umin */
    GEN->umin = -faux.f(x,faux.params);

    /* and now, an analogue calculation for umax */

    faux.f = (UNUR_FUNCT_GENERIC*) _unur_aux_bound_umax;
    faux.params = gen;

    /* calculating start point for extremum search routine */
    sx = _unur_isfinite(DISTR.BD_RIGHT) ? (cx+DISTR.BD_RIGHT)/2. : (cx+1.); 
    bx = _unur_isfinite(DISTR.BD_RIGHT) ? DISTR.BD_RIGHT : BD_MAX;

    x = (_unur_FP_same(DISTR.BD_RIGHT,cx)) 
      ? cx: _unur_util_find_max(faux, cx, bx, sx);
      
    while (!_unur_isfinite(x) && (fabs(bx) >= UNUR_EPSILON) ) { 
       /* _unur_util_find_max() could not yet find a suitable extremum */
       /* trying with a sequence of intervals with decreasing length   */
       bx = bx/10.; sx = bx/2.; 
       x = _unur_util_find_max(faux, cx, bx, sx);
    }
    /* umax */
    GEN->umax = faux.f(x,faux.params);

    /* additional scaling of boundary rectangle */
    GEN->umin = GEN->umin - (GEN->umax-GEN->umin)*NROU_RECT_SCALING/2.;
    GEN->umax = GEN->umax + (GEN->umax-GEN->umin)*NROU_RECT_SCALING/2.;

    /* check for bounded rectangle */
    if (! (_unur_isfinite(GEN->umin) && _unur_isfinite(GEN->umax))) {
       _unur_error(gen->genid , UNUR_ERR_GENERIC, "umin or umax not finite");  
       return UNUR_ERR_GENERIC;
    }

  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_nrou_rectangle() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_nrou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = nrou (naive ratio-of-uniforms)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_nrou_sample",gen->genid);
  if (gen->variant & NROU_VARFLAG_VERIFY) fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(LOG,"%s: r-parameter = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,NROU_SET_R);
  fprintf(LOG,"\n%s:\n",gen->genid);

  /* center */
  fprintf(LOG,"%s: center = %g\n",gen->genid,GEN->center);
  fprintf(LOG,"%s:\n",gen->genid);

  /* bounding rectangle */
  fprintf(LOG,"%s: Rectangle:\n",gen->genid);
  fprintf(LOG,"%s:    left  upper point = (%g,%g)\n",gen->genid,GEN->umin,GEN->vmax);
  fprintf(LOG,"%s:    right upper point = (%g,%g)\n",gen->genid,GEN->umax,GEN->vmax);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_nrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_nrou_info( struct unur_gen *gen, int help )
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
  double harea;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }

  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX ) 
      _unur_string_append(info,"\n[ Hint: %s\n\t%s ]\n",
			  "You may provide the \"mode\" or at least",
			  "the \"center\" (a point near the mode)."); 
  }
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: NROU (Naive Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"   r = %g\n\n", GEN->r);

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   bounding rectangle = (%g,%g) x (%g,%g)\n",
		      GEN->umin,GEN->umax, 0.,GEN->vmax);
  harea = (GEN->umax - GEN->umin) * GEN->vmax;
  _unur_string_append(info,"   area(hat) = %g\n", harea);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", 2. * harea / DISTR.area);
  else
    _unur_string_append(info,"= %.2f [approx.]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   r = %g  %s\n", GEN->r,
			(gen->set & NROU_SET_R) ? "" : "[default]");
    _unur_string_append(info,"   center = %g  %s\n",GEN->center,
			(gen->set & NROU_SET_CENTER) ? "" : "[default]");
    _unur_string_append(info,"   v = %g  %s\n", GEN->vmax,
			(gen->set & NROU_SET_V) ? "" : "[numeric.]");
    _unur_string_append(info,"   u = (%g, %g)  %s\n", GEN->umin,GEN->umax,
			(gen->set & NROU_SET_U) ? "" : "[numeric.]");
    if (gen->variant & NROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  if (help) {
    if ( !(gen->set & NROU_SET_V) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"v\" to avoid numerical estimate." );
    if ( !(gen->set & NROU_SET_U) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"u\" to avoid slow (and inexact) numerical estimates." );
    _unur_string_append(info,"\n");
  }

} /* end of _unur_nrou_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
