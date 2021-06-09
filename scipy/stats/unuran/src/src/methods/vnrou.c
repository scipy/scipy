/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vnrou.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
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
 *   [1] Wakefield J.C., Gelfand A.E., Smith A.F.M.                          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method.                                                             *
 *       Statistics and Computing (1991) 1, pp (129-133)                     *
 *                                                                           *
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h> 
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "vnrou.h"
#include "vnrou_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define VNROU_VARFLAG_VERIFY   0x002u   /* run verify mode                   */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define VNROU_DEBUG_REINIT   0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define VNROU_SET_U       0x001u     /* set u values of bounding rectangle   */
#define VNROU_SET_V       0x002u     /* set v values of bounding rectangle   */
#define VNROU_SET_R       0x008u     /* set r-parameter                      */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VNROU"         /* type of generator                         */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vnrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_vnrou_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vnrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vnrou_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_vnrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_vnrou_sample_cvec( struct unur_gen *gen, double *vec );
static int _unur_vnrou_sample_check( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static int _unur_vnrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute (minimal) bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_vnrou_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_vnrou_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_vnrou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_vnrou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec  /* data for distribution in generator object */
#define SAMPLE    gen->sample.cvec       /* pointer to sampling routine      */     
#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

#define _unur_vnrou_getSAMPLE(gen) \
   ( ((gen)->variant & VNROU_VARFLAG_VERIFY) \
     ? _unur_vnrou_sample_check : _unur_vnrou_sample_cvec )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_vnrou_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_vnrou_par) );
  COOKIE_SET(par,CK_VNROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR->r		= 1.; 	      /* r-parameter of the generalized method       */
  PAR->vmax      = 0.;         /* v-boundary of bounding rectangle (unknown)  */
  PAR->umin 	= NULL;       /* u-boundary of bounding rectangle (unknown)  */
  PAR->umax 	= NULL;       /* u-boundary of bounding rectangle (unknown)  */
  par->method   = UNUR_METH_VNROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */
  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_vnrou_init;

  return par;

} /* end of unur_vnrou_new() */

/*****************************************************************************/

int
unur_vnrou_set_u( struct unur_par *par, double *umin, double *umax )
     /*----------------------------------------------------------------------*/
     /* Sets left and right u-boundary of bounding rectangle.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   umin ... left boundary of rectangle                                */
     /*   umax ... right boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int d; /* index used in dimension loops (0 <= d < dim) */

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );

  /* check new parameter for generator */
  for (d=0; d<par->distr->dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
  }
  
  /* set values */
  PAR->umin = umin;
  PAR->umax = umax;
  
  /* changelog */
  par->set |= VNROU_SET_U;

  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_u() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_chg_u( struct unur_gen *gen, double *umin, double *umax )
     /*----------------------------------------------------------------------*/
     /* Sets left and right u-boundary of bounding rectangle.                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   umin ... left boundary of rectangle                                */
     /*   umax ... right boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int d; /* index used in dimension loops (0 <= d < dim) */

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( GENTYPE, umin, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, umax, UNUR_ERR_NULL );

  /* check new parameter for generator */
  for (d=0; d<GEN->dim; d++) {
    if (!_unur_FP_greater(umax[d],umin[d])) {
      _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"umax <= umin");
      return UNUR_ERR_PAR_SET;
    }
  }
  
  /* set values */
  memcpy(GEN->umin, umin, GEN->dim * sizeof(double));
  memcpy(GEN->umax, umax, GEN->dim * sizeof(double));
  
  /* changelog */
  gen->set |= VNROU_SET_U;

  return UNUR_SUCCESS;

} /* end of unur_vnrou_chg_u() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_set_v( struct unur_par *par, double vmax )
     /*----------------------------------------------------------------------*/
     /* Sets upper v-boundary of bounding rectangle.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par  ... pointer to parameter for building generator object        */
     /*   vmax ... upper boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  PAR->vmax = vmax;

  /* changelog */
  par->set |= VNROU_SET_V;

  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_v() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_chg_v( struct unur_gen *gen, double vmax )
     /*----------------------------------------------------------------------*/
     /* Sets upper v-boundary of bounding rectangle.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   vmax ... upper boundary of rectangle                               */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (vmax <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"vmax <= 0");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store values */
  GEN->vmax = vmax;

  /* changelog */
  gen->set |= VNROU_SET_V;

  return UNUR_SUCCESS;

} /* end of unur_vnrou_chg_v() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_set_r( struct unur_par *par, double r )
     /*----------------------------------------------------------------------*/
     /* Set the r-parameter for the generalized ratio-of-uniforms method.    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   r   ... r-parameter                                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, VNROU );

  /* check new parameter for generator */
  if (r <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"r<=0");
    return UNUR_ERR_PAR_SET;
  }

  /* store data */
  PAR->r = r;

  /* changelog */
  par->set |= VNROU_SET_R;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_r() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, VNROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | VNROU_VARFLAG_VERIFY) : (par->variant & (~VNROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vnrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_vnrou_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, VNROU, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cvec_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify bounding rectangle on */
    gen->variant |= VNROU_VARFLAG_VERIFY;
  else
    /* turn verify bounding rectangle off */
    gen->variant &= ~VNROU_VARFLAG_VERIFY;

  SAMPLE = _unur_vnrou_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vnrou_chg_verify() */

/*---------------------------------------------------------------------------*/

double 
unur_vnrou_get_volumehat( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /*   Get the volume of below the hat.                                   */
     /*   For normalized densities, i.e. when the volume below PDF is 1,     */
     /*   this value equals the rejection constant for the vnrou method.     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  double vol;
  int d;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, VNROU, INFINITY );

  /* compute volume of bounding rectangle */
  vol = GEN->vmax;  
  for (d=0; d<GEN->dim; d++) {
    vol *= (GEN->umax[d]-GEN->umin[d]);
  }
  /* compute volume of corresponding hat function */
  vol *= (GEN->r*GEN->dim+1);

  /* return result */
  return vol;
} /* end of unur_vnrou_get_volumehat() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_vnrou_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_VNROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VNROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_vnrou_create(par);
  _unur_par_free(par); 
  if (!gen) return NULL; 

  /* compute bounding rectangle */
  if (_unur_vnrou_rectangle(gen)!=UNUR_SUCCESS) {
    _unur_vnrou_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_vnrou_debug_init(gen);
#endif

  return gen;

} /* end of _unur_vnrou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_vnrou_reinit( struct unur_gen *gen )
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

  /* compute bounding rectangle */
  if ( (rcode = _unur_vnrou_rectangle(gen))!=UNUR_SUCCESS) {
    return rcode;
  }

  /* (re)set sampling routine */
  SAMPLE = _unur_vnrou_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug & VNROU_DEBUG_REINIT) _unur_vnrou_debug_init(gen);
#endif

  return UNUR_SUCCESS;
} /* end of _unur_vnrou_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vnrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VNROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_vnrou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_VNROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_vnrou_getSAMPLE(gen);
  gen->destroy = _unur_vnrou_free;
  gen->clone = _unur_vnrou_clone;
  gen->reinit = _unur_vnrou_reinit;

  /* copy parameters into generator object */
  GEN->dim   = gen->distr->dim;       /* dimension */
  GEN->r     = PAR->r;                /* r-parameter of the vnrou method */  
  GEN->vmax  = PAR->vmax;             /* upper v-boundary of bounding rectangle */
  
  /* allocate memory for u-boundary arrays */
  GEN->umin = _unur_xmalloc( GEN->dim * sizeof(double)); /* bounding rectangle */
  GEN->umax = _unur_xmalloc( GEN->dim * sizeof(double)); /* bounding rectangle */

  if (PAR->umin != NULL) memcpy(GEN->umin, PAR->umin, GEN->dim * sizeof(double));
  if (PAR->umax != NULL) memcpy(GEN->umax, PAR->umax, GEN->dim * sizeof(double));

  /* get center of the distribution */
  GEN->center = unur_distr_cvec_get_center(gen->distr);
 
#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_vnrou_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_vnrou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vnrou_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_vnrou_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VNROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate memory for u-arrays */
  CLONE->umin = _unur_xmalloc( GEN->dim * sizeof(double));
  CLONE->umax = _unur_xmalloc( GEN->dim * sizeof(double));
  
  /* copy parameters into clone object */
  memcpy(CLONE->umin, GEN->umin, GEN->dim * sizeof(double));
  memcpy(CLONE->umax, GEN->umax, GEN->dim * sizeof(double));

  /* copy data */
  CLONE->center = unur_distr_cvec_get_center(clone->distr);

  return clone;

#undef CLONE
} /* end of _unur_vnrou_clone() */

/*****************************************************************************/

int
_unur_vnrou_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{ 
  double U, V;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_VNROU_GEN,UNUR_ERR_COOKIE); 

  dim = GEN->dim;
 
  while (1) {

    /* generate point uniformly on rectangle */
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    for (d=0; d<dim; d++) {
      U = GEN->umin[d] + _unur_call_urng(gen->urng) * (GEN->umax[d] - GEN->umin[d]);
      vec[d] = U/pow(V,GEN->r) + GEN->center[d];
    }
    
    /* X[] inside domain ? */
    
    /* accept or reject */
    if (V <= pow(PDF(vec),1./(GEN->r * dim + 1.)))
      return UNUR_SUCCESS;
  }

} /* end of _unur_vnrou_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_vnrou_sample_check( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random sample vector (return)                              */
     /*----------------------------------------------------------------------*/
{ 
  double U, V;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  int hat_error;

  double fx,sfx,xfx;
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  
  COOKIE_CHECK(gen,CK_VNROU_GEN,UNUR_ERR_COOKIE); 

  dim = GEN->dim;
 
  while (1) {
    /* generate point uniformly on rectangle */
    while ( _unur_iszero(V = _unur_call_urng(gen->urng)) );
    V *= GEN->vmax;
    for (d=0; d<dim; d++) {
      U = GEN->umin[d] + _unur_call_urng(gen->urng) * (GEN->umax[d] - GEN->umin[d]);
      vec[d] = U/pow(V,GEN->r) + GEN->center[d];
    }
    
    /* X[] inside domain ? */

    /* evaluate PDF */
    fx = PDF(vec);
    
    /* a point on the boundary of the region of acceptance
       has the coordinates ( (vec[]-center[]) * (fx)^(r/r*dim+1)), fx^(1/r*dim+1) ). */
    sfx = pow( fx, 1./(GEN->r * dim+1.) );
    /* check hat */
    hat_error=0;
    if ( sfx > (1.+DBL_EPSILON) * GEN->vmax ) hat_error++;  
   
    sfx = pow( fx, GEN->r/(GEN->r * dim + 1.) );
    for (d=0; d<dim; d++) {
     xfx = (vec[d]-GEN->center[d]) * sfx;
     if ( (xfx < (1.+UNUR_EPSILON) * GEN->umin[d]) 
       || (xfx > (1.+UNUR_EPSILON) * GEN->umax[d]))
       hat_error++;
    }

    if (hat_error>0) _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");
 
    /* accept or reject */
    if (V <= pow(PDF(vec),1./( GEN->r * dim + 1.)))
      return UNUR_SUCCESS;
  }

} /* end of _unur_vnrou_sample_check() */

/*****************************************************************************/

void
_unur_vnrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_VNROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->umin) free(GEN->umin); 
  if (GEN->umax) free(GEN->umax);
  _unur_generic_free(gen);

} /* end of _unur_vnrou_free() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_vnrou_rectangle( struct unur_gen *gen )
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

  int d; /* index used in dimension loops (0 <= d < dim) */
  struct MROU_RECTANGLE *rr;
  int rectangle_compute;
  
  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_VNROU_GEN, UNUR_ERR_COOKIE );

  /* Boundary rectangle is already set */
  if ((gen->set & VNROU_SET_U) && (gen->set & VNROU_SET_V)) {
    return UNUR_SUCCESS;
  }

  /* Allocating and filling mrou_rectangle struct */
  rr = _unur_mrou_rectangle_new();

  rr->distr  = gen->distr;
  rr->dim    = GEN->dim;
  rr->umin   = GEN->umin;
  rr->umax   = GEN->umax;
  rr->r      = GEN->r;
  rr->center = GEN->center; 
  rr->genid  = gen->genid;
  
  /* calculate bounding rectangle */
  rectangle_compute = _unur_mrou_rectangle_compute(rr);
  
  if (!(gen->set & VNROU_SET_V)) {
     /* user has not provided any upper bound for v */
     GEN->vmax = rr->vmax;
  }

  if (!(gen->set & VNROU_SET_U)) {
    /* user has not provided any bounds for u */
    for (d=0; d<GEN->dim; d++) {
      GEN->umin[d] = rr->umin[d];
      GEN->umax[d] = rr->umax[d];
    }
  }

  free(rr);
  
  if (rectangle_compute != UNUR_SUCCESS)
    return UNUR_ERR_INF;
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_vnrou_rectangle() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_vnrou_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int d, dim; /* index used in dimension loops (0 <= d < dim) */
  double vol;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VNROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();
  dim = GEN->dim;

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = vnrou (naive ratio-of-uniforms)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  
  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_vnrou_sample",gen->genid);
  if (gen->variant & VNROU_VARFLAG_VERIFY) fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(LOG,"%s: r-parameter = %g",gen->genid, GEN->r);
  _unur_print_if_default(gen,VNROU_SET_R);
  fprintf(LOG,"\n%s:\n",gen->genid);

  /* print center */
  _unur_matrix_print_vector( GEN->dim, GEN->center, "center =", LOG, gen->genid, "\t   ");

  /* print bounding rectangle */
  fprintf(LOG,"%s: Rectangle:",gen->genid);
  if (!((gen->set & VNROU_SET_U) && (gen->set & VNROU_SET_V)))
    fprintf(LOG,"\t[computed]");
  else 
    fprintf(LOG,"\t[input]");
  fprintf(LOG,"\n");

  vol = GEN->vmax;
  fprintf(LOG,"%s:\tvmax = %g\n",gen->genid, GEN->vmax);
  for (d=0; d<dim; d++) {
    vol *= (GEN->umax[d]-GEN->umin[d]);
    fprintf(LOG,"%s:\tumin[%d],umax[%d] = (%g,%g)\n",gen->genid, 
	    d, d, GEN->umin[d], GEN->umax[d]);
  }
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s:\tvolume = %g\t(hat = %g)\n",gen->genid, vol, vol*(GEN->r*GEN->dim+1));
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_vnrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_vnrou_info( struct unur_gen *gen, int help )
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
  double hvol;

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
  _unur_string_append(info,"method: VNROU (Naive Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"   r = %g\n", GEN->r);
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");

  _unur_string_append(info,"   bounding rectangle = ");
  for (i=0; i<GEN->dim; i++)
    _unur_string_append(info,"%s(%g,%g)", i?"x":"", GEN->umin[i], GEN->umax[i]);
  _unur_string_append(info," x (0,%g)\n", GEN->vmax);

  hvol = GEN->vmax;
  for (i=0; i<GEN->dim; i++)
    hvol *= GEN->umax[i] - GEN->umin[i];
  _unur_string_append(info,"   volume(hat) = %g\n", hvol);

  _unur_string_append(info,"   rejection constant ");
  if ((distr->set & UNUR_DISTR_SET_PDFVOLUME) && _unur_isone(GEN->r))
    _unur_string_append(info,"= %g\n", (GEN->dim + 1.) * hvol / DISTR.volume);
  else
    _unur_string_append(info,"= %.2f  [approx.]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/((1.+GEN->dim)*samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");

    _unur_string_append(info,"   r = %g  %s\n", GEN->r,
			(gen->set & VNROU_SET_R) ? "" : "[default]");

    _unur_string_append(info,"   v = %g  %s\n", GEN->vmax,
 			(gen->set & VNROU_SET_V) ? "" : "[numeric.]");

    _unur_string_append(info,"   u = ");
    _unur_distr_info_vector( gen, GEN->umin, GEN->dim);
    _unur_string_append(info," -- ");
    _unur_distr_info_vector( gen, GEN->umax, GEN->dim);
    _unur_string_append(info,"%s\n",(gen->set & VNROU_SET_U) ? "" : "  [numeric.]"); 

    if (gen->variant & VNROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    _unur_string_append(info,"\n");
  }

  /* Hints */
  if (help) {
    if ( !(gen->set & VNROU_SET_V) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"v\" to avoid numerical estimate." );
    if ( !(gen->set & VNROU_SET_U) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"u\" to avoid slow (and inexact) numerical estimates." );
    _unur_string_append(info,"\n");
  }

} /* end of _unur_vnrou_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
