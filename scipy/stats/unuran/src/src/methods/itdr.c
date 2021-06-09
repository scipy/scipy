/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      itdr.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    inverse transformed density rejection                        *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF and pole.                                                  *
 *      Produce a value x consistent with its density.                       *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the density function                                      *
 *      location of pole                                                     *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      splitting point between pole and tail region                         *
 *      c-value for pole and tail region, repectively                        *
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
 *   [1] W. Hoermann, J. Leydold, and G. Derflinger (2007):                  *
 *       ACM Trans. Model. Comput. Simul. 17(4), pp.18.                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <utils/unur_fp_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "itdr.h"
#include "itdr_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants:                                                                */

/* maximum value for parameter c                                             */
/* (if this value is changed, then also change the value in the WARNING      */
/* message in the set calls unur_itdr_set_cp() and unur_itdr_set_cp()! )     */
#define C_MAX  (-0.1)

/* relative distance for computing derivatives numerically                   */
/* (1.e-6 yieldes better results than 1.e-8 when compared with exact values) */
#define DX (1.e-6)      

/* point near pole for estimating lim_{x->plole} ilc(x)                      */
/* the point x_i * NEAR_POLE is used                                         */
#define NEAR_POLE  (1.e-8)

/* precision when computing intersection points.                             */
/* since the performance is not very sensitive to the accuracy of the result */
/* this value need not be very small.                                        */
#define RESOLUTION_XI (1.e-5)

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define ITDR_VARFLAG_VERIFY   0x001u   /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define ITDR_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define ITDR_SET_XI      0x001u     /* set intersection point                */
#define ITDR_SET_CP      0x002u     /* set c-value for pole region           */
#define ITDR_SET_CT      0x004u     /* set c-value for tail region           */

/*---------------------------------------------------------------------------*/

#define GENTYPE "ITDR"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_itdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_itdr_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_itdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_itdr_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_itdr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_itdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_sample( struct unur_gen *gen );
static double _unur_itdr_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator.                                                    */
/*---------------------------------------------------------------------------*/

static int _unur_itdr_get_hat( struct unur_gen *gen );
static int _unur_itdr_get_hat_pole( struct unur_gen *gen );
static int _unur_itdr_get_hat_tail( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* construct hat function                                                    */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_lc( struct unur_gen *gen, double x );
/*---------------------------------------------------------------------------*/
/* compute local concavity at x                                              */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_ilc( struct unur_gen *gen, double x );
/*---------------------------------------------------------------------------*/
/* compute inverse local concavity at x                                      */
/*---------------------------------------------------------------------------*/

static double _unur_itdr_find_xt( struct unur_gen *gen, double b );
/*---------------------------------------------------------------------------*/
/* solves equation (x-b)*f'(x)+f(x)=0, where f is PDF of distribution        */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_itdr_debug_init( const struct unur_gen *gen, int error );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_itdr_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_itdr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_itdr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

/* the hat is computed for density PDF(sign*(x-pole)). Thus we have to       */
/* transform the internal x value back to the original scale before the PDF  */
/* is evaluated or the generated value is returned.                          */
/* likewise, all x values for the PDF given by the user has to transformed   */
/* into the internal scale.                                                  */
/* transformation: internal scale --> original scale */
#define I2O(x)    ( GEN->sign*(x) + GEN->pole )
/* transformation: original scale --> internal scale */
#define O2I(x)    ( GEN->sign * ((x)-GEN->pole) )

/* call to PDF and its derivative in internal scale */
#define PDF(x)    ( _unur_cont_PDF(I2O(x), gen->distr) )
#define dPDF(x)   ( GEN->sign * _unur_cont_dPDF(I2O(x), gen->distr) )

/* call to PDF and its derivative in original scale */
#define PDFo(x)   ( _unur_cont_PDF((x), gen->distr) )
#define dPDFo(x)  ( _unur_cont_dPDF((x), gen->distr) )



#define logPDF(x)   ( _unur_cont_logPDF(I2O(x), gen->distr) )   /* call to logPDF   */
#define dlogPDF(x)  ( GEN->sign * _unur_cont_dlogPDF(I2O(x), gen->distr) )  /* call to derivative of log PDF */

#define logPDFo(x)  ( _unur_cont_logPDF((x), gen->distr) )  /* call to logPDF   */
#define dlogPDFo(x) ( _unur_cont_dlogPDF((x), gen->distr) )  /* call to derivative of log PDF */


/*---------------------------------------------------------------------------*/
/* transformations */

/* T_c(x) */
#define T(c,x)   ( -pow((x), (c)) )
#define DT(c,x)  ( -(c)*pow((x), ((c)-1.)) )
#define TI(c,x)  ( pow(-(x), 1./(c)) )
#define FT(c,x)  ( -pow(-(x), ((c)+1.)/(c))*((c)/((c)+1.)) )
#define FTI(c,x) ( -pow(-(x)*((c)+1.)/(c), (c)/((c)+1.)) )

/* logarithm of inverse transformation */
#define logTI(c,x)  ( -log(-(x)) / (c) )

/* special case: T_{-1/2}(x)   ["square root transformation"] */
#define TsI(c,x)  ( 1./((x)*(x)) )
#define FTs(c,x)  ( -1./(x) )
#define FTsI(c,x) ( -1./(x) )

/*---------------------------------------------------------------------------*/

#define _unur_itdr_getSAMPLE(gen) \
   ( ((gen)->variant & ITDR_VARFLAG_VERIFY) \
     ? _unur_itdr_sample_check : _unur_itdr_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_itdr_new( const struct unur_distr *distr )
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

  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"dPDF");
    return NULL;
  }

  if (!(distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode (pole)");
    return NULL; 
  }

  if ( ! (_unur_isfinite(DISTR_IN.mode) &&
	  (_unur_FP_equal(DISTR_IN.mode,DISTR_IN.domain[0]) ||
	   _unur_FP_equal(DISTR_IN.mode,DISTR_IN.domain[1]))) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_PROP,"pole not on boundary of domain");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_itdr_par) );
  COOKIE_SET(par,CK_ITDR_PAR);

  /* copy input */
  par->distr    = distr;          /* pointer to distribution object          */

  /* set default values */
  PAR->xi = INFINITY;       /* intersection point lc(x)=ilc(x)               */
  PAR->cp = INFINITY;       /* c-value for pole region (unknown)             */
  PAR->ct = INFINITY;       /* c-value for tail region (unknown)             */
  
  par->method   = UNUR_METH_ITDR;     /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_itdr_init;

  return par;

} /* end of unur_itdr_new() */

/*****************************************************************************/

int
unur_itdr_set_xi( struct unur_par *par, double xi )
     /*----------------------------------------------------------------------*/
     /* Sets intersection point xi where lc(x) = ilc(x)                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   xi ... intersection point                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );

  /* check new parameter for generator */
  if (xi <= par->distr->data.cont.BD_LEFT || 
      xi >= par->distr->data.cont.BD_RIGHT) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"xi out of domain");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store value */
  PAR->xi = xi;

  /* changelog */
  par->set |= ITDR_SET_XI;

  return UNUR_SUCCESS;
} /* end of unur_itdr_set_bx() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_set_cp( struct unur_par *par, double cp )
     /*----------------------------------------------------------------------*/
     /* Sets c-value for transformation T for inverse density in pole region */
     /*                                                                      */
     /* parameters:                                                          */
     /*   cp ... c-value for pole region                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );

  /* check new parameter for generator */
  if ( cp > C_MAX || cp <= -1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"cp > -0.1 or <= -1");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store value */
  PAR->cp = cp;

  /* changelog */
  par->set |= ITDR_SET_CP;

  return UNUR_SUCCESS;
} /* end of unur_itdr_set_cp() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_set_ct( struct unur_par *par, double ct )
     /*----------------------------------------------------------------------*/
     /* Sets c-value for transformation T for density in tail region         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   ct ... c-value for tail region                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double range;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ITDR );

  /* check new parameter for generator */
  range = ( par->distr->data.cont.BD_RIGHT
	    - par->distr->data.cont.BD_LEFT );
  if ( ct > C_MAX || (ct <= -1. && !_unur_isfinite(range)) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ct > -0.1 or <= -1");
    return UNUR_ERR_PAR_SET;
  }
  
  /* store value */
  PAR->ct = ct;

  /* changelog */
  par->set |= ITDR_SET_CT;

  return UNUR_SUCCESS;
} /* end of unur_itdr_set_ct() */

/*---------------------------------------------------------------------------*/

double unur_itdr_get_xi( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get intersection point xi                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   xi       ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );

  return GEN->xi;
} /* end of unur_itdr_get_xi() */

/*---------------------------------------------------------------------------*/

double unur_itdr_get_cp( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get intersection c-value for pole region                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   cp       ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );

  return GEN->cp;
} /* end of unur_itdr_get_cp() */

/*---------------------------------------------------------------------------*/

double unur_itdr_get_ct( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get intersection c-value for tail region                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   ct       ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );

  return GEN->ct;
} /* end of unur_itdr_get_ct() */

/*---------------------------------------------------------------------------*/

double unur_itdr_get_area( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below hat                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   area     ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, ITDR, INFINITY );

  return GEN->Atot;
} /* end of unur_itdr_get_area() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, ITDR );

  /* we use a bit in variant */
  par->variant = (verify) 
    ? (par->variant | ITDR_VARFLAG_VERIFY) 
    : (par->variant & (~ITDR_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_itdr_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_itdr_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, ITDR, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= ITDR_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~ITDR_VARFLAG_VERIFY;

  SAMPLE = _unur_itdr_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_itdr_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_itdr_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_ITDR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_ITDR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_itdr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_itdr_check_par(gen) != UNUR_SUCCESS) {
    _unur_itdr_free(gen); return NULL;
  }

  /* create hat function */
  if (_unur_itdr_get_hat(gen) != UNUR_SUCCESS) {
#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_itdr_debug_init(gen,UNUR_FAILURE);
#endif
    _unur_itdr_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_itdr_debug_init(gen,UNUR_SUCCESS);
#endif

  return gen;

} /* end of _unur_itdr_init() */

/*---------------------------------------------------------------------------*/

int
_unur_itdr_reinit( struct unur_gen *gen )
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

  /* we do not use the given values for when we run reinit */
  gen->set &= ~(ITDR_SET_XI | ITDR_SET_CP | ITDR_SET_CT);

  /* check parameters */
  if ( (rcode = _unur_itdr_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* create hat function */
  rcode = _unur_itdr_get_hat(gen);

  /* (re)set sampling routine */
  SAMPLE = _unur_itdr_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug & ITDR_DEBUG_REINIT)
      _unur_itdr_debug_init(gen,rcode);
#endif

  return rcode;
} /* end of _unur_itdr_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_itdr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_ITDR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_itdr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_ITDR_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_itdr_getSAMPLE(gen);
  gen->destroy = _unur_itdr_free;
  gen->clone = _unur_itdr_clone;
  gen->reinit = _unur_itdr_reinit;

  /* copy data from distribution into generator object*/
  GEN->pole = DISTR.mode;  /* location of pole                    */

  /* copy some parameters into generator object */
  GEN->xi = PAR->xi;       /* intersection point lc(x)=ilc(x)     */
  GEN->cp = PAR->cp;       /* c-value for pole region             */
  GEN->ct = PAR->ct;       /* c-value for tail region             */

  /* initialize values */
  GEN->bx = INFINITY;      /* splitting point betw. pole and tail */
  GEN->xp = INFINITY;      /* design point in pole region         */
  GEN->xt = INFINITY;      /* design point in tail region         */
  GEN->alphap = INFINITY;  /* parameters for hat in pole region   */
  GEN->betap = INFINITY;
  GEN->Tfxt = INFINITY;    /* parameters for hat in tail region   */
  GEN->dTfxt = INFINITY;   /* parameters for hat in tail region   */
  GEN->by = INFINITY;      /* hat of pole region at bx            */
  GEN->Ap = INFINITY;      /* areas in upper pole region          */     
  GEN->Ac = INFINITY;      /* areas in central region             */     
  GEN->At = INFINITY;      /* areas in tail region                */     
  GEN->Atot = INFINITY;    /* total area below hat                */
  GEN->sy = 0.;            /* squeeze for central region          */
  GEN->sign = 1.;          /* region: +1 .. (-oo,0], -1 .. [0,oo) */
  GEN->bd_right = INFINITY; /* right boundary of shifted domain   */
  
#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_itdr_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_itdr_create() */

/*---------------------------------------------------------------------------*/

int
_unur_itdr_check_par( struct unur_gen *gen )
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

  /* copy data from distribution into generator object*/
  GEN->pole = DISTR.mode;  /* location of pole                    */

  /* estimate sign of region: +1 ... (-oo,0], -1 ... [0,oo) */
  do {
    if (_unur_isfinite(DISTR.BD_LEFT) && !_unur_isfinite(DISTR.BD_RIGHT)) {
      GEN->sign = 1.; 
      if (dPDFo(DISTR.BD_LEFT) <= 0.) break;
    }
    if (!_unur_isfinite(DISTR.BD_LEFT) && _unur_isfinite(DISTR.BD_RIGHT)) {
      GEN->sign = -1.;
      if (dPDFo(DISTR.BD_RIGHT) >= 0.) break;
    }
    if (_unur_isfinite(DISTR.BD_LEFT) && _unur_isfinite(DISTR.BD_RIGHT)) {
      GEN->sign = (PDFo(DISTR.BD_LEFT)>=PDFo(DISTR.BD_RIGHT)) ? 1. : -1.;
      if ( GEN->sign*dPDFo(DISTR.BD_LEFT) <= 0. &&
	   GEN->sign*dPDFo(DISTR.BD_RIGHT) <= 0. )
	break;
    }
    /* else */
    _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute sign of region");
    return UNUR_ERR_DISTR_PROP;
  } while (1);

  /* right boundary of shifted domain */
  GEN->bd_right = ( (GEN->sign > 0) 
		    ? DISTR.BD_RIGHT - GEN->pole
		    : GEN->pole - DISTR.BD_LEFT );

  return UNUR_SUCCESS;
} /* end of _unur_itdr_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_itdr_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_itdr_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ITDR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_itdr_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_itdr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_ITDR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_ITDR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_itdr_free() */

/*****************************************************************************/

double
_unur_itdr_sample( struct unur_gen *gen )
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
  double U, V, X, Y;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ITDR_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on (0,Atot) */
    U = _unur_call_urng(gen->urng) * GEN->Atot;

    /* generate pair (X,Y) below hat */

    if (U < GEN->Ap) {
      /* upper pole region */
      V = _unur_call_urng(gen->urng) * GEN->Ap;
      if (_unur_isfsame(GEN->cp, -0.5)) {
	/* square root transformation */
	Y = ( FTsI(GEN->cp, GEN->betap*V + FTs(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TsI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
      else {
	/* general T_c transformation */
	Y = ( FTI(GEN->cp, GEN->betap*V + FT(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
    }

    else if ((U -= GEN->Ap) < GEN->Ac) {
      /* central region */
      X = U * GEN->bx / GEN->Ac;
      Y = _unur_call_urng(gen->urng) * GEN->by;
      if (Y <= GEN->sy)
	/* squeeze acceptance */
	return (I2O(X));
    }

    else {
      /* tail region */
      U -= GEN->Ac;
      if (_unur_isfsame(GEN->ct, -0.5)) {
	/* square root transformation */
	X = GEN->xt + (FTsI(GEN->ct,
			   GEN->dTfxt*U
			   + FTs(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TsI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
      else {
	/* general T_c transformation */
	X = GEN->xt + (FTI(GEN->ct,
			   GEN->dTfxt*U
			   + FT(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
    }

    /* transform back into original scale */
    X = I2O(X);

    /* accept or reject */
    if (Y <= PDFo(X))
      return X;
  }

} /* end of _unur_itdr_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_sample_check( struct unur_gen *gen )
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
#define ht(x)  ( TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*((x)-GEN->xt)) )
#define hp(x)  ( (T(GEN->cp,(x)) - GEN->alphap) / GEN->betap )

  double U, V, X, Y;
  double fx, hx, sqx;  /* values of PDF, hat and squeeze at x */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ITDR_GEN,INFINITY);

  while (1) {
    /* generate point uniformly on (0,Atot) */
    U = _unur_call_urng(gen->urng) * GEN->Atot;

    /* generate pair (X,Y) below hat */

    if (U < GEN->Ap) {
      /* upper pole region */
      V = _unur_call_urng(gen->urng) * GEN->Ap;
      if (_unur_isfsame(GEN->cp, -0.5)) {
	/* square root transformation */
	Y = ( FTsI(GEN->cp, GEN->betap*V + FTs(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TsI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
      else {
	/* general T_c transformation */
	Y = ( FTI(GEN->cp, GEN->betap*V + FT(GEN->cp,GEN->alphap+GEN->betap*GEN->by))
	      - GEN->alphap ) / GEN->betap;
	X = U * TI(GEN->cp, GEN->alphap+GEN->betap*Y) / GEN->Ap;
      }
      hx = hp(X);
      sqx = 0.;
    }

    else if ((U -= GEN->Ap) < GEN->Ac) {
      /* central region */
      X = U * GEN->bx / GEN->Ac;
      Y = _unur_call_urng(gen->urng) * GEN->by;
      hx = hp(X);
      sqx = GEN->sy;
      /* no squeeze acceptance in verify mode */
      /* [ if (Y <= GEN->sy) return (I2O(X)); ] */
    }

    else {
      /* tail region */
      U -= GEN->Ac;
      if (_unur_isfsame(GEN->ct, -0.5)) {
	/* square root transformation */
	X = GEN->xt + (FTsI(GEN->ct,
			   GEN->dTfxt*U
			   + FTs(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TsI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
      else {
	/* general T_c transformation */
	X = GEN->xt + (FTI(GEN->ct,
			   GEN->dTfxt*U
			   + FT(GEN->ct,
				GEN->Tfxt + GEN->dTfxt*(GEN->bx-GEN->xt))
			   )
		       - GEN->Tfxt) / GEN->dTfxt;
	Y = ( _unur_call_urng(gen->urng)
	      * TI(GEN->ct, GEN->Tfxt + GEN->dTfxt*(X - GEN->xt)));
      }
      hx = ht(X);
      sqx = 0.;
    }

    /* transform back into original scale */
    X = I2O(X);

    /* compute PDF at x */
    fx = PDFo(X);

    /* verify hat and squeeze function */
    if ( (1.+UNUR_EPSILON) * hx < fx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");

    if ( (1.-UNUR_EPSILON) * sqx > fx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) < squeeze(x)");

    /* accept or reject */
    if (Y <= PDFo(X))
      return X;
  }

#undef ht
#undef hp
} /* end of _unur_itdr_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_itdr_get_hat( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* construct hat function                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ITDR_GEN,UNUR_ERR_COOKIE);

  /* Get candidate for bx */
  if (gen->set & ITDR_SET_XI) {
    /* bx set user; have to shift by pole */
    GEN->bx = O2I(GEN->xi);
  }
  else {
    /* compute intersection point of local concavity and inverse lc */
    GEN->bx = _unur_itdr_find_xt( gen, 0. );
    GEN->xi = I2O(GEN->bx);
    if (!_unur_isfinite(GEN->bx)) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute bx");
      return UNUR_ERR_DISTR_PROP;
    }
  }

  /* pole region          */
  if (_unur_itdr_get_hat_pole(gen) != UNUR_SUCCESS)
    return UNUR_ERR_DISTR_PROP;

  /* tail region          */
  if (_unur_FP_equal(GEN->bx, GEN->bd_right)) {
    GEN->At = 0.;
  }
  else {
    if (_unur_itdr_get_hat_tail(gen) != UNUR_SUCCESS)
      return UNUR_ERR_DISTR_PROP;
  }
  
  /* total area below hat */
  GEN->Atot = GEN->Ap + GEN->Ac + GEN->At;

  return UNUR_SUCCESS;
} /* end of _unur_itdr_get_hat() */

/*---------------------------------------------------------------------------*/

int
_unur_itdr_get_hat_pole( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* construct hat function in pole region                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define hp(x)  ( (T(cp,(x)) - GEN->alphap) / GEN->betap )

  double cp, xp;
  double pdf_bx;
  double near_pole, ilc_near_pole, pdf_near_pole, logpdf_near_pole;
  double ilc_bx = -INFINITY;

  /* get cp */
  if (gen->set & ITDR_SET_CP) {
    /* cp set by user */
    cp = GEN->cp;
  }
  else {
    ilc_bx = _unur_itdr_ilc(gen, GEN->bx);
    near_pole = GEN->bx*NEAR_POLE + fabs(GEN->pole)*DBL_EPSILON;
    ilc_near_pole = (DISTR.logpdf) 
      ? logPDF(near_pole) / log(near_pole)
      : log(PDF(near_pole)) / log(near_pole);
    cp = ilc_near_pole;
    if (cp > C_MAX) cp = C_MAX;
    if (cp <= -1.) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->cp = cp;
  }
  if (cp < -0.5)
    GEN->bx = _unur_min(2.*GEN->bx, GEN->bd_right);

  /* compute PDF at check points */
  pdf_bx = PDF(GEN->bx);
  near_pole = fabs(GEN->pole)*DBL_EPSILON;
  if (near_pole < 1.e-100) near_pole = 1.e-100;
  pdf_near_pole = logpdf_near_pole = INFINITY;
  while (1) {
    /* we have to search for a point with PDF(x) < INFINITY */
    if (DISTR.logpdf) {
      logpdf_near_pole = logPDF(near_pole);
      if (_unur_isfinite(logpdf_near_pole)) 
	break;
    }
    else {
      pdf_near_pole = PDF(near_pole);
      if (_unur_isfinite(pdf_near_pole)) 
	break;
    }
    near_pole *= 1000.;
    if (!_unur_isfinite(near_pole)) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
      return UNUR_ERR_DISTR_PROP;
    }
  }

  /* get design point xp */
  while (1) {
    xp = GEN->bx * pow(1.+cp, -1./cp);
    if ( !(xp > 0. && xp < GEN->bx) ) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: xp");
      return UNUR_ERR_DISTR_PROP;
    }

    /* parameters for hat */
    GEN->betap = DT(cp,xp) / dPDF(xp);
    GEN->alphap = T(cp,xp) - GEN->betap * PDF(xp);

    /* check hat */
    if ( hp(GEN->bx) < pdf_bx ||
	 (DISTR.logpdf && _unur_FP_less(log(hp(near_pole)), logpdf_near_pole)) ||
	 (DISTR.logpdf==NULL && _unur_FP_less(hp(near_pole), pdf_near_pole)) ) {
      if (gen->set & ITDR_SET_CP) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"inverse pdf not T_cp concave");
	return UNUR_ERR_DISTR_PROP;
      }
      /* try new value for cp */
      GEN->cp = cp = 0.9*cp-0.1;
      if (cp < ilc_bx) {
	/* we try ilc at bx first before we use an even smaller value for cp */
	GEN->cp = cp = ilc_bx; 
	ilc_bx = -INFINITY;
      }
      if (cp < -0.999) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for pole: cp");
	return UNUR_ERR_DISTR_PROP;
      }	
    }
    else break;
  }
  GEN->xp = xp;

  /* hat at intersection point */
  GEN->by = hp(GEN->bx);

  /* area below hat */
  GEN->Ap = -FT(cp, GEN->alphap + GEN->betap * GEN->by) / GEN->betap;
  GEN->Ac = GEN->by * GEN->bx;

  /* add squeeze for central region */
  GEN->sy = PDF(GEN->bx);

  return UNUR_SUCCESS;

#undef hp
} /* end of _unur_itdr_get_hat_pole() */

/*---------------------------------------------------------------------------*/

int
_unur_itdr_get_hat_tail( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* construct hat function in tail region                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define ht(x)  ( TI(ct, GEN->Tfxt + GEN->dTfxt*((x)-xt)) )
#define loght(x)  ( TI(ct, GEN->Tfxt + GEN->dTfxt*((x)-xt)) )

  double ct, xt;
  double lc_bx, lc_inf;
  double br;
  double bx = GEN->bx;

  /* get design point xt */
  GEN->xt = xt = _unur_itdr_find_xt( gen, bx );

  /* get ct */
  if (gen->set & ITDR_SET_CT) {
    /* ct set by user */
    ct = GEN->ct;
  }
  else {
    /* first try: use point between bx and xt */
    ct = _unur_itdr_lc(gen, 0.5*(bx + xt));
    /* check local concavity at right boundary */
    if ( _unur_isfinite(GEN->bd_right)) 
      lc_inf = _unur_itdr_lc(gen, GEN->bd_right);
    else { /* right boundary = infinity */
      if (DISTR.logpdf) {
	lc_inf = log(1.e100) / logPDF(1.e100);
	/* we need lim x->oo log(x) / log(f(x))    */
	/* however, this convergence is very slow. */
	/* so we add -0.1 to be on the save side.  */
	lc_inf += -0.01;
      }
      else {
	lc_inf = log(1.e10*bx) / log(PDF(1.e10*bx));
	/* we need lim x->oo log(x) / log(f(x))    */
	/* however, this convergence is very slow. */
	/* so we add -0.1 to be on the save side.  */
	lc_inf += -0.05;
      }
    }
    if (lc_inf < ct) ct = lc_inf;

    /* ct should not be too close to 0. */
    if (ct > C_MAX) ct = C_MAX;
    if (ct <= -1.) {
      _unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for tail: ct");
      return UNUR_ERR_DISTR_PROP;
    }
    GEN->ct = ct;
  }

  /* compute and check parameters for hat */
  lc_bx = _unur_itdr_lc(gen, bx);
  while (1) {
    /* parameters for hat */
    GEN->Tfxt = T(ct, PDF(xt));
    GEN->dTfxt = DT(ct, PDF(xt)) * dPDF(xt);

    /* check hat */
    br = 1000.*bx;  /* "very large x" */
    if (br > GEN->bd_right) br = GEN->bd_right;

    if ( ((GEN->Tfxt + GEN->dTfxt*(bx-xt)) >= 0.) ||
	 (DISTR.logpdf && (_unur_FP_less(loght(br),logPDF(br)) ||
			   _unur_FP_less(loght(bx), logPDF(bx)) )) ||
	 (DISTR.logpdf==NULL && (_unur_FP_less(ht(br),PDF(br)) ||
				 _unur_FP_less(ht(bx), PDF(bx)) )) ) {
      if (gen->set & ITDR_SET_CT) {
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"pdf not T_ct concave");
	return UNUR_ERR_DISTR_PROP;
      }
      /* try new value for ct */
      ct = 0.5*(ct + lc_bx);
      if (ct > GEN->ct || ct < -0.999 || _unur_FP_approx(ct,lc_bx)) {
	/* new ct value is even larger or too small or its time to stop */ 
	_unur_error(gen->genid,UNUR_ERR_DISTR_PROP,"cannot compute hat for tail: ct");
	return UNUR_ERR_DISTR_PROP;
      }
      GEN->ct = ct;
    }
    else
      break;
  }

  /* area below hat */
  GEN->At = (!_unur_isfinite(GEN->bd_right)) ? 0. 
    : FT(ct, GEN->Tfxt + GEN->dTfxt * (GEN->bd_right - xt)) / GEN->dTfxt;
  GEN->At += -FT(ct, GEN->Tfxt + GEN->dTfxt * (GEN->bx - xt)) / GEN->dTfxt;

  return UNUR_SUCCESS;

#undef ht
#undef loght
} /* end of _unur_itdr_get_hat_tail() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_lc( struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* compute local concavity at x                                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... point x                                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   local concavity                                                    */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double dx, f, df, ddf;

  if (DISTR.dlogpdf == NULL) {
    /* use PDF */

    f = PDF(x);
    df = dPDF(x);

    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;
    if (x+dx > GEN->bd_right)
      ddf = (dPDF(x)-dPDF(x-dx))/dx;
    else
      ddf = (dPDF(x+dx)-dPDF(x-dx))/(2.*dx);
    
    return 1. - ddf*f/(df*df); 
  }

  else {
    /* use logarithm of PDF */

    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;

    if (x+dx > GEN->bd_right)
      return (1./dlogPDF(x) - 1./dlogPDF(x-dx))/dx;
    else
      return (1./dlogPDF(x+dx) - 1./dlogPDF(x-dx))/(2.*dx);
  }

} /* end of _unur_itdr_lc() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_ilc( struct unur_gen *gen, double x )
     /*----------------------------------------------------------------------*/
     /* compute inverse local concavity at x                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... point x                                                    */
     /*                                                                      */
     /* return:                                                              */
     /*   inverse local concavity                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{

  if (DISTR.dlogpdf == NULL) {
    /* use PDF */
    double dx, df, ddf;

    df = dPDF(x);

    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;

    if (x+dx > GEN->bd_right)
      ddf = (dPDF(x)-dPDF(x-dx))/dx;
    else
      ddf = (dPDF(x+dx)-dPDF(x-dx))/(2.*dx);

    return 1.+x*ddf/(df); 
  }

  else {
    /* use logarithm of PDF */
    double dx, dlf, ddlf;

    dlf = dlogPDF(x);

    dx = x * DX + fabs(GEN->pole) * UNUR_SQRT_DBL_EPSILON;
    if (x-dx <= 0.) dx = x;
    
    if (x+dx > GEN->bd_right)
      ddlf = (dlogPDF(x)-dlogPDF(x-dx))/dx;
    else
      ddlf = (dlogPDF(x+dx)-dlogPDF(x-dx))/(2.*dx);
    
    return 1.+x*(dlf + ddlf/dlf); 
  }

} /* end of _unur_itdr_ilc() */

/*---------------------------------------------------------------------------*/

double
_unur_itdr_find_xt( struct unur_gen *gen, double b )
     /*----------------------------------------------------------------------*/
     /* solves equation (x-b)*f'(x)+f(x)=0, where f is PDF of distribution   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   startx ... starting point for finding x                            */
     /*                                                                      */
     /* return:                                                              */
     /*   solution xi                                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  /* function for finding root */
#define FKT(x) ( DISTR.dlogpdf \
                 ? (1./((x)-b) + dlogPDF(x)) \
                 : (((x)-b)*dPDF(x) + PDF(x)) )

  double xl, xu;  /* lower and upper boundary of bracket */
  double xn;      /* new guess for root */

  /* check parameter */
  if (b < 0.) return INFINITY;

  /* find appropriate starting value */
  xl = b + _unur_max(1., (fabs(GEN->pole)+b)*UNUR_SQRT_DBL_EPSILON);
  if (xl > GEN->bd_right) xl = GEN->bd_right;
  while (!_unur_isfinite(FKT(xl)) || _unur_iszero(PDF(xl)) ) {
    xl = 0.5*(xl + b);
    if (!_unur_isfinite(xl) || _unur_FP_same(xl,b)) return INFINITY;
  }
  xu = xl;

  /* check upper value xu */
  if (_unur_FP_greater(xu,GEN->bd_right)) return GEN->bd_right;

  /* find bracket for root */
  if (FKT(xl)>0.) {
    do {
      xl = xu;
      xu += xu - b;
      if (!_unur_isfinite(xu) || xu < (1.+2.*DBL_EPSILON)*xl)
	/* unable to proceed --> break to avoid infinite loop */
	return INFINITY;
      if (xu >= GEN->bd_right) 
	/* we have reached right boundary */
	return GEN->bd_right;
    } while(FKT(xu) > 0.);
  }
  else { /* FKT(xl)<=0. */
    do {
      xu = xl;
      xl = 0.5*(xl + b);
      if (!_unur_isfinite(xl)) return INFINITY;
    } while(FKT(xl) < 0.);
  }

  /* use bisection to find root */
  while(xu > (1.+RESOLUTION_XI)*xl) {
    xn = 0.5*(xl+xu);
    if(FKT(xn)>0.) 
      xl = xn;
    else 
      xu = xn;
  }

  /* return point */
  return 0.5*(xl+xu);

#undef FKT
} /* end of _unur_itdr_find_xt() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_itdr_debug_init( const struct unur_gen *gen, int error )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   error ... error code of hat generation                             */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ITDR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = itdr (inverse transformed density rejection)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_itdr_sample",gen->genid);
  if (gen->variant & ITDR_VARFLAG_VERIFY) fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);

  /* parameters */
  fprintf(LOG,"%s: sign = %g\n",gen->genid, GEN->sign);
  fprintf(LOG,"%s: pole = %g\n",gen->genid, GEN->pole);
  fprintf(LOG,"%s: bd_right = %g\n",gen->genid, GEN->bd_right);
  fprintf(LOG,"%s: xi = %g",gen->genid, GEN->xi);
  fprintf(LOG,"%s\n", (gen->set & ITDR_SET_XI) ? "" : " [computed]");
  fprintf(LOG,"%s: bx = %g\n",gen->genid, GEN->bx);

  fprintf(LOG,"%s: pole region:\n",gen->genid);
  fprintf(LOG,"%s:\tcp = %g",gen->genid, GEN->cp);
  fprintf(LOG,"%s\n", (gen->set & ITDR_SET_CP) ? "" : " [computed]");
  fprintf(LOG,"%s:\txp = %g\n",gen->genid, GEN->xp);
  fprintf(LOG,"%s:\talphap = %g, betap = %g\n",gen->genid, GEN->alphap, GEN->betap);
  fprintf(LOG,"%s:\tby = %g\n",gen->genid, GEN->by);
  fprintf(LOG,"%s:\tsy = %g\n",gen->genid, GEN->sy);

  fprintf(LOG,"%s: tail region:\n",gen->genid);
  fprintf(LOG,"%s:\tct = %g",gen->genid, GEN->ct);
  fprintf(LOG,"%s\n", (gen->set & ITDR_SET_CT) ? "" : " [computed]");
  fprintf(LOG,"%s:\txt = %g\n",gen->genid, GEN->xt);
  fprintf(LOG,"%s:\tTfxt = %g, dTfxt = %g\n",gen->genid, GEN->Tfxt, GEN->dTfxt);

  fprintf(LOG,"%s: Area = %g + %g + %g = %g\n",gen->genid,
	  GEN->Ap, GEN->Ac, GEN->At, GEN->Atot);

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: **** INIT %s ***\n",gen->genid,
	  (error==UNUR_SUCCESS) ? "successful" : "failed" );   
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_itdr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_itdr_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   pole/mode = %g\n", DISTR.mode);
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: ITDR (Inverse Transformed Density Rejection -- 2 point method)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g  [ = %g + %g + %g ]\n",
		      GEN->Atot, GEN->Ap, GEN->Ac, GEN->At);
  _unur_string_append(info,"   rejection constant = ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"%g\n", GEN->Atot/DISTR.area);
  else
    _unur_string_append(info,"%.2f  [approx. ]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   cp = %g  %s\n", GEN->cp,
			(gen->set & ITDR_SET_CP) ? "" : " [computed]");
    _unur_string_append(info,"   ct = %g  %s\n", GEN->cp,
			(gen->set & ITDR_SET_CT) ? "" : " [computed]");
    _unur_string_append(info,"   xi = %g  %s\n", GEN->xi,
			(gen->set & ITDR_SET_XI) ? "" : " [computed]");
    if (gen->variant & ITDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_itdr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
