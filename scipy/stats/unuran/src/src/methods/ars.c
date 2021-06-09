/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ars.c                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    adaptive rejection sampling (Gilks & Wild)                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given logPDF of a log-concave distribution                           *
 *      produce a value x consistent with its density                        *
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
 *   REFERENCES and DESCRIPTION:                                             *
 *                                                                           *
 *   This is special implementation of the TDR method that requires the      *
 *   logPDF of a distribution and uses solely adaptive rejection sampling    *
 *   for finding constrution points as proposed by Gilks & Wild (1992).      *
 *                                                                           *
 *   See tdr.c for further information.                                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Store intersection point in structure when defined.                       */

#ifdef DEBUG_STORE_IP 
#  undef DEBUG_STORE_IP
#endif

/* #define DEBUG_STORE_IP 1 */

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "ars.h"
#include "ars_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define ARS_VARFLAG_VERIFY    0x0100u   /* flag for verifying mode           */
#define ARS_VARFLAG_PEDANTIC  0x0800u   /* whether pedantic checking is used */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define ARS_DEBUG_REINIT    0x00000002u  /* print parameters after reinit    */
#define ARS_DEBUG_IV        0x00000010u
#define ARS_DEBUG_SPLIT     0x00010000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define ARS_SET_CPOINTS        0x001u
#define ARS_SET_N_CPOINTS      0x002u
#define ARS_SET_PERCENTILES    0x004u
#define ARS_SET_N_PERCENTILES  0x008u
#define ARS_SET_RETRY_NCPOINTS 0x010u
#define ARS_SET_MAX_IVS        0x020u
#define ARS_SET_MAX_ITER       0x040u   /* maximum number of iterations      */

/*---------------------------------------------------------------------------*/

#define GENTYPE "ARS"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ars_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_ars_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ars_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ars_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_ars_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_ars_sample( struct unur_gen *generator );
static double _unur_ars_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_ars_starting_cpoints( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create list of construction points for starting segments.                 */
/* if user has not provided such points compute these by means of the        */
/* "equi-angle rule".                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_ars_starting_intervals( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute intervals from given starting construction points.                */
/*---------------------------------------------------------------------------*/

static int _unur_ars_interval_parameter( struct unur_gen *gen, struct unur_ars_interval *iv );
/*---------------------------------------------------------------------------*/
/* compute all necessary data for interval.                                  */
/*---------------------------------------------------------------------------*/

static struct unur_ars_interval *_unur_ars_interval_new( struct unur_gen *gen,
							 double x, double logfx );
/*---------------------------------------------------------------------------*/
/* make a new interval with left construction point x.                       */
/*---------------------------------------------------------------------------*/

static int _unur_ars_tangent_intersection_point( struct unur_gen *gen,
						 struct unur_ars_interval *iv, double *ipt );
/*---------------------------------------------------------------------------*/
/* compute cutting point of interval into left and right part.               */
/*---------------------------------------------------------------------------*/

static double _unur_ars_interval_logarea( struct unur_gen *gen, struct unur_ars_interval *iv,
					  double slope, double x );
/*---------------------------------------------------------------------------*/
/* compute log of area below piece of hat or squeeze in interval.            */
/*---------------------------------------------------------------------------*/

static int _unur_ars_interval_split( struct unur_gen *gen,
				     struct unur_ars_interval *iv_old, double x, double logfx );
/*---------------------------------------------------------------------------*/
/* split am interval point x. return 0 if not successful.                    */                                           
/*---------------------------------------------------------------------------*/

static int _unur_ars_improve_hat( struct unur_gen *gen, struct unur_ars_interval *iv,
				  double x, double logfx);
/*---------------------------------------------------------------------------*/
/* improve hat function by splitting interval                                */
/*---------------------------------------------------------------------------*/

static int _unur_ars_make_area_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make table of areas and compute largest area for rescaling.               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_init_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after (almost empty generator) object has been created.             */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_init_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized.                               */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_reinit_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before reinitialization of generator starts.                        */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_reinit_retry( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before second trial of reinitialization of generator starts.        */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_reinit_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been reinitialized.                             */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals                                                  */
/*---------------------------------------------------------------------------*/

static void _unur_ars_debug_split_start( const struct unur_gen *gen,
					 const struct unur_ars_interval *iv,
					 double x, double logfx );
static void _unur_ars_debug_split_stop( const struct unur_gen *gen,
					const struct unur_ars_interval *iv_left,
					const struct unur_ars_interval *iv_right );
/*---------------------------------------------------------------------------*/
/* print before and after an interval has been split (not / successfully).   */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_ars_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_ars_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_ars_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define logPDF(x)  _unur_cont_logPDF((x),(gen->distr))   /* call to logPDF   */
#define dlogPDF(x) _unur_cont_dlogPDF((x),(gen->distr))  /* call to derivative of log PDF */

/* areas in intervaled relative to maximal area */
#define scaled_logarea(iv)  ((iv)->logAhat - GEN->logAmax)
#define scaled_area(iv)     (exp(scaled_logarea(iv)))

/* log of function rescaled relative to maximal area */
#define rescaled_logf(logf) ((logf) - GEN->logAmax)

/*---------------------------------------------------------------------------*/

#define _unur_ars_getSAMPLE(gen) \
   ( ((gen)->variant & ARS_VARFLAG_VERIFY) \
     ? _unur_ars_sample_check : _unur_ars_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_ars_new( const struct unur_distr* distr )
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

  if (DISTR_IN.logpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"logPDF"); return NULL; }
  if (DISTR_IN.dlogpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of logPDF"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_ars_par) );
  COOKIE_SET(par,CK_ARS_PAR);

  /* copy input */
  par->distr              = distr;  /* pointer to distribution object        */

  /* set default values */
  PAR->starting_cpoints    = NULL;   /* pointer to array of starting points  */
  PAR->n_starting_cpoints  = 2;      /* number of starting points            */
  PAR->percentiles         = NULL;   /* pointer to array of percentiles      */
  PAR->n_percentiles       = 2;      /* number of percentiles                */
  PAR->retry_ncpoints      = 30;     /* number of cpoints for second trial of reinit */
  PAR->max_ivs             = 200;    /* maximum number of intervals          */
  PAR->max_iter            = 10000;  /* maximum number of iterations         */
 
  par->method   = UNUR_METH_ARS;     /* method                               */
  par->variant  = 0u;                /* default variant                      */

  par->set      = 0u;               /* inidicate default parameters          */
  par->urng     = unur_get_default_urng(); /* use default URNG               */
  par->urng_aux = par->urng;               /* no special auxilliary URNG     */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_ars_init;

  return par;

} /* end of unur_ars_new() */

/*****************************************************************************/

int
unur_ars_set_max_intervals( struct unur_par *par, int max_ivs )
     /*----------------------------------------------------------------------*/
     /* set maximum number of intervals                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ivs   ... maximum number of intervals                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );

  /* check new parameter for generator */
  if (max_ivs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of intervals < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ivs = max_ivs;

  /* changelog */
  par->set |= ARS_SET_MAX_IVS;

  return UNUR_SUCCESS;

} /* end of unur_ars_set_max_intervals() */

/*---------------------------------------------------------------------------*/

int 
unur_ars_set_cpoints( struct unur_par *par, int n_cpoints, const double *cpoints )
     /*----------------------------------------------------------------------*/
     /* set construction points for hat function                             */
     /* and/or its number for initialization                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   n_cpoints ... number of starting construction points               */
     /*   cpoints   ... pointer to array of starting construction points     */
     /*                 (NULL for changing only the number of default points)*/
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );

  /* check starting construction points */
  if (n_cpoints < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 2. using defaults");
    n_cpoints = 2;
    cpoints = NULL;
  }

  if (cpoints)
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_cpoints; i++ )
      if (cpoints[i] <= cpoints[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }

  /* store date */
  PAR->starting_cpoints = cpoints;
  PAR->n_starting_cpoints = n_cpoints;

  /* changelog */
  par->set |= ARS_SET_N_CPOINTS | ((cpoints) ? ARS_SET_CPOINTS : 0);

  return UNUR_SUCCESS;

} /* end of unur_ars_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_ars_set_reinit_percentiles( struct unur_par *par, int n_percentiles, const double *percentiles )
     /*----------------------------------------------------------------------*/
     /* set percentiles for construction points for hat function             */
     /* and/or its number for re-initialization                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par           ... pointer to parameter for building generator      */
     /*   n_percentiles ... number of percentiles                            */
     /*   percentiles   ... pointer to array of percentiles                  */
     /*                     (NULL for using a rule of thumb)                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );

  /* check given percentiles */
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }

  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
    
  if (percentiles) {
    /* percentiles must be strictly monontonically increasing */
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }

  /* store date */
  PAR->percentiles = percentiles;
  PAR->n_percentiles = n_percentiles;

  /* changelog */
  par->set |= ARS_SET_N_PERCENTILES | ((percentiles) ? ARS_SET_PERCENTILES : 0);

  return UNUR_SUCCESS;

} /* end of unur_ars_set_reinit_percentiles() */

/*---------------------------------------------------------------------------*/

int
unur_ars_chg_reinit_percentiles( struct unur_gen *gen, int n_percentiles, const double *percentiles )
     /*----------------------------------------------------------------------*/
     /* change percentiles for construction points for hat function          */
     /* and/or its number for re-initialization                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen           ... pointer to generator object                      */
     /*   n_percentiles ... number of percentiles                            */
     /*   percentiles   ... pointer to array of percentiles                  */
     /*                     (NULL for using a rule of thumb)                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );

  /* check given percentiles */
  if (n_percentiles < 2 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles < 2. using defaults");
    n_percentiles = 2;
    percentiles = NULL;
  }

  if (n_percentiles > 100 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of percentiles > 100. using 100");
    n_percentiles = 100;
  }
    
  if (percentiles) {
    /* percentiles must be strictly monontonically increasing */
    for( i=1; i<n_percentiles; i++ ) {
      if (percentiles[i] <= percentiles[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }
      if (percentiles[i] < 0.01 || percentiles[i] > 0.99) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"percentiles out of range");
	return UNUR_ERR_PAR_SET;
      }
    }
  }

  /* store date */
  GEN->n_percentiles = n_percentiles;
  GEN->percentiles = _unur_xrealloc( GEN->percentiles, n_percentiles * sizeof(double) );
  if (percentiles) {
    memcpy( GEN->percentiles, percentiles, n_percentiles * sizeof(double) );
  }
  else {
    if (n_percentiles == 2) {
      GEN->percentiles[0] = 0.25;
      GEN->percentiles[1] = 0.75;
    }
    else {
      for (i=0; i<n_percentiles; i++ )
	GEN->percentiles[i] = (i + 1.) / (n_percentiles + 1.);
    }
  }

  /* changelog */
  gen->set |= ARS_SET_N_PERCENTILES | ((percentiles) ? ARS_SET_PERCENTILES : 0);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_ars_chg_reinit_percentiles() */

/*---------------------------------------------------------------------------*/

int
unur_ars_set_reinit_ncpoints( struct unur_par *par, int ncpoints )
     /*----------------------------------------------------------------------*/
     /* set number of construction points for second trial of reinit         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator           */
     /*   ncpoints ... number of construction points                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );

  /* check number */
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->retry_ncpoints = ncpoints;

  /* changelog */
  par->set |= ARS_SET_RETRY_NCPOINTS; 

  return UNUR_SUCCESS;

} /* end of unur_ars_set_reinit_ncpoints() */

/*---------------------------------------------------------------------------*/

int
unur_ars_chg_reinit_ncpoints( struct unur_gen *gen, int ncpoints )
     /*----------------------------------------------------------------------*/
     /* change number of construction points for second trial of reinit      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen           ... pointer to generator object                      */
     /*   ncpoints ... number of construction points                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );

  /* check number */
  if (ncpoints < 10 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of construction points < 10");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  GEN->retry_ncpoints = ncpoints;

  /* changelog */
  gen->set |= ARS_SET_RETRY_NCPOINTS; 

  return UNUR_SUCCESS;

} /* end of unur_ars_chg_reinit_ncpoints() */

/*---------------------------------------------------------------------------*/

int
unur_ars_set_max_iter( struct unur_par *par, int max_iter )
     /*----------------------------------------------------------------------*/
     /* set maximum number of iterations                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   max_iter ... maximum number of iterations                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );

  /* check new parameter for generator */
  if (max_iter < 1) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of iterations");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_iter = max_iter;

  /* changelog */
  par->set |= ARS_SET_MAX_ITER;

  return UNUR_SUCCESS;

} /* end of unur_ars_set_max_iter() */

/*---------------------------------------------------------------------------*/

int
unur_ars_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, ARS );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | ARS_VARFLAG_VERIFY) : (par->variant & (~ARS_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_ars_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_ars_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  /* we use a bit in variant */
  gen->variant = (verify) 
    ? (gen->variant | ARS_VARFLAG_VERIFY) 
    : (gen->variant & (~ARS_VARFLAG_VERIFY));

  /* sampling routines */
  SAMPLE = _unur_ars_getSAMPLE(gen);
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_ars_chg_verify() */

/*---------------------------------------------------------------------------*/

int
unur_ars_set_pedantic( struct unur_par *par, int pedantic )
     /*----------------------------------------------------------------------*/
     /* turn pedantic mode on/off                                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   pedantic ... 0 = no pedantic mode, !0 = use pedantic mode          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   pedantic is the default                                            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, ARS );

  /* we use a bit in variant */
  par->variant = (pedantic) ? (par->variant | ARS_VARFLAG_PEDANTIC) : (par->variant & (~ARS_VARFLAG_PEDANTIC));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_ars_set_pedantic() */

/*---------------------------------------------------------------------------*/

double
unur_ars_get_loghatarea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get log of area below hat                                            */
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
  _unur_check_gen_object( gen, ARS, INFINITY );

  return log(GEN->Atotal) + GEN->logAmax;

} /* end of unur_ars_get_loghatarea() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_ars_init( struct unur_par *par )
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
/*   int i,k; */

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_ARS ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_ARS_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_ars_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ars_debug_init_start(gen);
#endif

  /* get starting points */
  if (_unur_ars_starting_cpoints(gen)!=UNUR_SUCCESS) {
    _unur_ars_free(gen); return NULL;
  }

  /* compute intervals for given starting points */
  if (_unur_ars_starting_intervals(gen)!=UNUR_SUCCESS) {
    _unur_ars_free(gen); return NULL;
  }

  /* update maximal number of intervals */
  if (GEN->n_ivs > GEN->max_ivs) {
    GEN->max_ivs = GEN->n_ivs;
  }

  /* make initial table of areas */
  _unur_ars_make_area_table(gen);

  /* is there any hat at all ? */
  if (GEN->Atotal <= 0. || !_unur_isfinite(GEN->Atotal)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points.");
    _unur_ars_free(gen);
    return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ars_debug_init_finished(gen);
#endif
  
  /* creation of generator object successfull */
  gen->status = UNUR_SUCCESS;

  /* o.k. */
  return gen;

} /* end of _unur_ars_init() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_reinit( struct unur_gen *gen )
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
  struct unur_ars_interval *iv,*next;
  double *bak_cpoints;
  int bak_n_cpoints;
  int i;
  int n_trials;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, ARS, UNUR_ERR_GEN_INVALID );

  /* first trial */
  n_trials = 1;

  /* which construction points should be used ? */
  if (gen->set & ARS_SET_N_PERCENTILES) {
    if (GEN->starting_cpoints==NULL || (GEN->n_starting_cpoints != GEN->n_percentiles)) {
      GEN->n_starting_cpoints = GEN->n_percentiles;
      GEN->starting_cpoints = _unur_xrealloc( GEN->starting_cpoints, GEN->n_percentiles * sizeof(double));
    }
    for (i=0; i<GEN->n_percentiles; i++) {
      GEN->starting_cpoints[i] = unur_ars_eval_invcdfhat( gen, GEN->percentiles[i] );
      if (!_unur_isfinite(GEN->starting_cpoints[i])) 
	/* we cannot use these starting points --> skip to second trial immediately */
	n_trials = 2;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & ARS_DEBUG_REINIT)
    _unur_ars_debug_reinit_start(gen);
#endif

  /* make backup of cpoints */
  bak_n_cpoints = GEN->n_starting_cpoints;
  bak_cpoints = GEN->starting_cpoints;

  for (;; ++n_trials) {
    /* free linked list of intervals */
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
    GEN->iv = NULL;
    GEN->n_ivs = 0;
    GEN->Atotal = 0.;
    GEN->logAmax = 0.;

    if (n_trials > 2) {
      /* we have done our best */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points for reinit");
      GEN->n_starting_cpoints = bak_n_cpoints;
      GEN->starting_cpoints = bak_cpoints;
      return UNUR_FAILURE;
    }

    if (n_trials > 1) {
      /* second trial */
      GEN->n_starting_cpoints = GEN->retry_ncpoints;
      GEN->starting_cpoints = NULL;
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (gen->debug & ARS_DEBUG_REINIT)
	_unur_ars_debug_reinit_retry(gen);
#endif
    }

    /* get starting points */
    if (_unur_ars_starting_cpoints(gen)!=UNUR_SUCCESS)
      continue;

    /* compute intervals for given starting points */
    if (_unur_ars_starting_intervals(gen)!=UNUR_SUCCESS)
      continue;

    /* update maximal number of intervals */
    if (GEN->n_ivs > GEN->max_ivs)
      GEN->max_ivs = GEN->n_ivs;
    
    /* make table of areas */
    _unur_ars_make_area_table(gen);
    
    /* is there any hat at all ? */
    if (GEN->Atotal <= 0.)
      continue;

    /* reinit successful */
    break;

  }

  /* clean up */
  if (n_trials > 1) {
    GEN->n_starting_cpoints = bak_n_cpoints;
    GEN->starting_cpoints = bak_cpoints;
  }

  /* (re)set sampling routine */
  SAMPLE = _unur_ars_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & ARS_DEBUG_REINIT)
    _unur_ars_debug_reinit_finished(gen);
#endif

  return UNUR_SUCCESS;
} /* end of _unur_ars_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ars_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_ARS_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_ars_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_ARS_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_ars_getSAMPLE(gen);
  gen->destroy = _unur_ars_free;
  gen->clone = _unur_ars_clone;
  gen->reinit = _unur_ars_reinit;

  /* set all pointers to NULL */
  GEN->iv          = NULL;
  GEN->n_ivs       = 0;
  GEN->percentiles = NULL;
  GEN->Atotal      = 0.;
  GEN->logAmax     = 0.;

  /* copy starting points */
  GEN->n_starting_cpoints = PAR->n_starting_cpoints;
  if (PAR->starting_cpoints) {
    GEN->starting_cpoints = _unur_xmalloc( PAR->n_starting_cpoints * sizeof(double) );
    memcpy( GEN->starting_cpoints, PAR->starting_cpoints, PAR->n_starting_cpoints * sizeof(double) );
  }
  else {
    GEN->starting_cpoints = NULL;
  }

  /* copy percentiles */
  if (gen->set & ARS_SET_N_PERCENTILES)
    unur_ars_chg_reinit_percentiles( gen, PAR->n_percentiles, PAR->percentiles );

  /* copy all other parameters */
  GEN->retry_ncpoints = PAR->retry_ncpoints;   /* number of cpoints for second trial of reinit */

  /* bounds for adding construction points  */
  GEN->max_ivs = _unur_max(2*PAR->n_starting_cpoints,PAR->max_ivs);  /* maximum number of intervals */

  /* bound for number of repetitions of rejection loop */
  GEN->max_iter = PAR->max_iter;

  /* copy variant */
  gen->variant = par->variant;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_ars_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_ars_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_ars_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_ars_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_ars_interval *iv, *clone_iv, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of intervals */
  clone_iv = NULL;
  clone_prev = NULL;
  for (iv = GEN->iv; iv != NULL; iv = iv->next) {
    /* copy segment */
    clone_iv = _unur_xmalloc( sizeof(struct unur_ars_interval) );
    memcpy( clone_iv, iv, sizeof(struct unur_ars_interval) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->iv = clone_iv;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_iv;
    }
    /* next step */
    clone_prev = clone_iv;
  }
  /* terminate linked list */
  if (clone_iv) clone_iv->next = NULL;

  /* copy starting points */
  if (GEN->starting_cpoints) {
    CLONE->starting_cpoints = _unur_xmalloc( GEN->n_starting_cpoints * sizeof(double) );
    memcpy( CLONE->starting_cpoints, GEN->starting_cpoints, GEN->n_starting_cpoints * sizeof(double) );
  }

  /* copy percentiles */
  if (GEN->percentiles) {
    CLONE->percentiles = _unur_xmalloc( GEN->n_percentiles * sizeof(double) );
    memcpy( CLONE->percentiles, GEN->percentiles, GEN->n_percentiles * sizeof(double) );
  }

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_ars_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_ARS ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_ars_debug_free(gen);
#endif

  /* free linked list of intervals */
  {
    struct unur_ars_interval *iv,*next;
    for (iv = GEN->iv; iv != NULL; iv = next) {
      next = iv->next;
      free(iv);
    }
  }

  /* free list of starting points */
  if (GEN->starting_cpoints) 
    free (GEN->starting_cpoints);

  /* free list of percentiles */
  if (GEN->percentiles) 
    free (GEN->percentiles);

  /* free other memory not stored in list */
  _unur_generic_free(gen);

} /* end of _unur_ars_free() */

/*****************************************************************************/

double
_unur_ars_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (original variant by Gilks & Wild)             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*   squeeze(x) = f(x0) * exp(sq * (x-x0))                              */
     /*                                                                      */
     /*   left hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                   */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = f(x1) * exp( (Tf)'(x1) *  (x-x1) )                  */
     /*   generation:                                                        */
     /*      X = x1 + 1/(Tf)'(x1) * \log( (Tf)'(x1)/f(x1) * U + 1 )          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_ars_interval *iv, *cp;
  double U, logV;                   /* (log of) uniform random number        */
  double X;                         /* generated point                       */
  double logfx, logsqx, loghx;      /* log of density, squeeze, and hat at X */
  double x0, logfx0, dlogfx0, fx0;  /* construction point and logPDF at x0   */
  int n_trials;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ARS_GEN,INFINITY);

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 

  for (n_trials=0; n_trials<GEN->max_iter; ++n_trials) {

    /* sample from U(0,1) */
    U = _unur_call_urng(gen->urng);

    /* evaluate inverse of hat CDF */

    /* find interval by sequential search */
    /* Remark: there is no need for a guide table as we only generate one point! */
    iv =  GEN->iv;
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* l.h.s. or r.h.s. of hat */
    if (-U < (scaled_area(iv) * iv->Ahatr_fract)) { /* right */
      cp = iv->next;
      /* U unchanged */
    }
    else {                /* left */
      cp = iv;
      U += scaled_area(iv);
    }

    /* PDF at x0 */
    x0 = cp->x;
    logfx0 = cp->logfx;
    dlogfx0 = cp->dlogfx;
    fx0 = exp(rescaled_logf(logfx0));

    /* random variate */
    if (_unur_iszero(dlogfx0))
      X = x0 + U / fx0;
    else {
      double t = dlogfx0 * U / fx0;
      if (fabs(t) > 1.e-6)
	X = x0 + log(t + 1.) * U / (fx0 * t);
      /* x = x0 + log(t + 1.) / dlogfx0; is cheaper but numerical unstable */
      else if (fabs(t) > 1.e-8)
	/* use Taylor series */
	X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
      else
	X = x0 + U / fx0 * (1 - t/2.);
    }
    
    /* log of hat at x */
    loghx = rescaled_logf(logfx0) + dlogfx0*(X - x0);

    /* log of a random point between 0 and hat at x */
    logV = log(_unur_call_urng(gen->urng)) + loghx;
    
    /* log of spueeze at x */
    logsqx = rescaled_logf(iv->logfx) + iv->sq*(X - iv->x);
 
    /* below squeeze ? */
    if (logV <= logsqx)
      return X;
    
    /* log of PDF at x */
    logfx = logPDF(X);

    /* below PDF ? */
    if (logV <= rescaled_logf(logfx))
      return X;

    /* being above PDF is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      /* first check for valid values of X and logf(X) */
      if (! (_unur_isfinite(X) && _unur_isfinite(logfx)) ) {
	X = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */
	logfx = logPDF(X);
      }
      if ( (_unur_ars_improve_hat( gen, iv, X, logfx) != UNUR_SUCCESS)
	   && (gen->variant & ARS_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    /* else reject and try again */

  }

  /* number of trials exceeded */
  _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"max number of iterations exceeded");
  return UNUR_INFINITY;

} /* end of _unur_ars_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_ars_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results                             */
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
  struct unur_ars_interval *iv, *cp;
  double U, logV;                   /* (log of) uniform random number        */
  double X;                         /* generated point                       */
  double logfx, logsqx, loghx;      /* log of density, squeeze, and hat at X */
  double x0, logfx0, dlogfx0, fx0;  /* construction point and logPDF at x0   */
  int n_trials;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_ARS_GEN,INFINITY);

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"empty generator object");
    return INFINITY;
  } 

  for (n_trials=0; n_trials<GEN->max_iter; ++n_trials) {

    /* sample from U(0,1) */
    U = _unur_call_urng(gen->urng);

    /* evaluate inverse of hat CDF */

    /* find interval by sequential search */
    /* Remark: there is no need for a guide table as we only generate one point! */
    iv =  GEN->iv;
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* l.h.s. or r.h.s. of hat */
    if (-U < (scaled_area(iv) * iv->Ahatr_fract)) { /* right */
      cp = iv->next;
      /* U unchanged */
    }
    else {                /* left */
      cp = iv;
      U += scaled_area(iv);
    }

    /* PDF at x0 */
    x0 = cp->x;
    logfx0 = cp->logfx;
    dlogfx0 = cp->dlogfx;
    fx0 = exp(rescaled_logf(logfx0));

    /* random variate */
    if (_unur_iszero(dlogfx0))
      X = x0 + U / fx0;
    else {
      double t = dlogfx0 * U / fx0;
      if (fabs(t) > 1.e-6)
	X = x0 + log(t + 1.) * U / (fx0 * t);
      /* x = x0 + log(t + 1.) / dlogfx0; is cheaper but numerical unstable */
      else if (fabs(t) > 1.e-8)
	/* use Taylor series */
	X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
      else
	X = x0 + U / fx0 * (1 - t/2.);
    }
    
    /* log of hat at x */
    loghx = rescaled_logf(logfx0) + dlogfx0*(X - x0);

    /* log of spueeze at x */
    logsqx = rescaled_logf(iv->logfx) + iv->sq*(X - iv->x);
 
    /* log of PDF at x */
    logfx = logPDF(X);

    /* check result */
    if (X < DISTR.BD_LEFT || X > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
    }
    if (_unur_FP_greater(rescaled_logf(logfx), loghx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not log-concave!");
    }
    if (_unur_FP_less(rescaled_logf(logfx), logsqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not log-concave!");
    }

    /* log of a random point between 0 and hat at x */
    logV = log(_unur_call_urng(gen->urng)) + loghx;

    /* below squeeze ? */
    if (logV <= logsqx)
      return X;

    /* below PDF ? */
    if (logV <= rescaled_logf(logfx))
      return X;

    /* being above PDF is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      /* first check for valid values of X and logf(X) */
      if (! (_unur_isfinite(X) && _unur_isfinite(logfx)) ) {
	X = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */
	logfx = logPDF(X);
      }
      if ( (_unur_ars_improve_hat( gen, iv, X, logfx) != UNUR_SUCCESS)
	   && (gen->variant & ARS_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    /* reject and try again */

  }

  /* number of trials exceeded */
  _unur_warning(gen->genid,UNUR_ERR_GEN_SAMPLING,"max number of iterations exceeded");
  return UNUR_INFINITY;

} /* end of _unur_ars_sample_check() */

/*---------------------------------------------------------------------------*/

double
unur_ars_eval_invcdfhat( const struct unur_gen *gen, double U )
     /*----------------------------------------------------------------------*/
     /* evaluate the inverse of the hat CDF at u                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   U   ... argument for inverse CDF (0<=U<=1, no validation!)         */
     /*                                                                      */
     /* return:                                                              */
     /*   inverse of hat CDF.                                                */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_ars_interval *iv, *cp;
  double X;                         /* generated point                       */
  double x0, logfx0, dlogfx0, fx0;  /* construction point and logPDF at x0   */


  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_ARS ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY; 
  }
  COOKIE_CHECK(gen,CK_ARS_GEN,INFINITY);

  if ( U<0. || U>1.) {
    _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"argument u not in [0,1]");
  }

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 

  /* validate argument */
  if (U<=0.) return DISTR.domain[0];
  if (U>=1.) return DISTR.domain[1];

  /* find interval by sequential search */
  /* Remark: there is no need for a guide table as we only generate one point! */
  iv =  GEN->iv;
  U *= GEN->Atotal;
  while (iv->Acum < U) {
    iv = iv->next;
  }
  
  /* rescale U: U in (-A_hat, 0) */
  U -= iv->Acum;
  
  /* l.h.s. or r.h.s. of hat */
  if (-U < (scaled_area(iv) * iv->Ahatr_fract)) { /* right */
    cp = iv->next;
    /* U unchanged */
  }
  else {                /* left */
    cp = iv;
    U += scaled_area(iv);
  }
  
  /* PDF at x0 */
  x0 = cp->x;
  logfx0 = cp->logfx;
  dlogfx0 = cp->dlogfx;
  fx0 = exp(rescaled_logf(logfx0));
  
  /* X = H^{-1}(U) in interval*/
  if (_unur_iszero(dlogfx0))
    X = x0 + U / fx0;
  else {
    double t = dlogfx0 * U / fx0;
    if (fabs(t) > 1.e-6)
      X = x0 + log(t + 1.) * U / (fx0 * t);
    /* x = x0 + log(t + 1.) / dlogfx0; is cheaper but numerical unstable */
    else if (fabs(t) > 1.e-8)
      /* use Taylor series */
      X = x0 + U / fx0 * (1 - t/2. + t*t/3.);
    else
      X = x0 + U / fx0 * (1 - t/2.);
  }
  
  return X;
  
} /* end of unur_ars_eval_invcdfhat() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_ars_improve_hat( struct unur_gen *gen, struct unur_ars_interval *iv,
			  double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* improve hat function by splitting interval                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   iv         ... pointer to interval that has to be split            */
     /*   x          ... splitting point                                     */
     /*   logfx      ... value of logPDF at splitting point                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... improving hat successful                       */
     /*   others          ... error: PDF not monotone in interval            */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* add construction point */
  result = _unur_ars_interval_split(gen, iv, x, logfx);
  if (result!=UNUR_SUCCESS && result!=UNUR_ERR_SILENT) {
    /* condition for PDF is violated! */
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
    if (gen->variant & ARS_VARFLAG_PEDANTIC) {
      /* replace sampling routine by dummy routine that just returns INFINITY */
      SAMPLE = _unur_sample_cont_error;
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /* splitting successful --> update table of areas */
  _unur_ars_make_area_table(gen);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_ars_improve_hat() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_starting_cpoints( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting intervals.                  */
     /* if not provided as arguments compute these                           */
     /* by means of the "equiangular rule" from AROU.                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_ars_interval *iv;
  double left_angle, right_angle, diff_angle, angle;
  double x, logfx, logfx_last;
  int is_increasing;
  int i;
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  
  /* reset counter of intervals */
  GEN->n_ivs = 0;

  /* prepare for computing construction points */
  if (!GEN->starting_cpoints) {
    /* angles of boundary of domain */
    left_angle =  _unur_FP_is_minus_infinity(DISTR.BD_LEFT) ? -M_PI/2. : atan(DISTR.BD_LEFT);
    right_angle = _unur_FP_is_infinity(DISTR.BD_RIGHT)      ? M_PI/2.  : atan(DISTR.BD_RIGHT);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (GEN->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = DISTR.BD_LEFT;
  is_increasing = TRUE;
    
  logfx = logfx_last = _unur_isfinite(x) ? logPDF(x) : -INFINITY;
  iv = GEN->iv = _unur_ars_interval_new( gen, x, logfx );
  if (iv == NULL) return UNUR_ERR_GEN_DATA;  /* logPDF(x) overflow */

  /* now all the other points */
  for( i=0; i<=GEN->n_starting_cpoints; i++ ) {

    /* construction point */
    if (i < GEN->n_starting_cpoints) {
      if (GEN->starting_cpoints) {
	/* construction points provided by user */
	x = GEN->starting_cpoints[i];
	/* check starting point */
	if (x < DISTR.BD_LEFT || x > DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equiangular rule" */
	angle += diff_angle;
	x = tan( angle );
      }
    }
    else {
      /* the very last interval. it is rather a "virtual" interval to store
	 the right vertex of the last interval, i.e., the right boundary point. */
      x = DISTR.BD_RIGHT;
    }

    /** TODO: check if two construction points are too close ??
	check if a point is too close to mode ??  */

    /* value of PDF at starting point */
    logfx = _unur_isfinite(x) ? logPDF(x) : -INFINITY;

    /* check value of PDF at starting point */
    if (!is_increasing && logfx > logfx_last * (1.+DBL_EPSILON)) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal!");
      return UNUR_ERR_GEN_CONDITION;
    }

    /* check whether we are outside of support of PDF */
    if (!_unur_isfinite(logfx) && !_unur_isfinite(logfx_last) ) {
      /* we do not need two such point */
      if (is_increasing) {
	/* PDF is still increasing, i.e., constant 0 til now */
	if (i<GEN->n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the PDF is constant 0 on all construction points.
	     then we need both boundary points. */
	  iv->x = x;  /* we only have to change x, everything else remains unchanged */
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with logPDF(x) > -INFINITY */
	break;
    }
    
    /* need a new interval */
    iv->next = _unur_ars_interval_new( gen, x, logfx );
    if (iv->next == NULL) return UNUR_ERR_GEN_DATA;    /* logPDF(x) overflow */
    
    /* skip pointer to current interval */
    iv = iv->next;

    /* PDF still increasing ? */
    if (is_increasing && logfx < logfx_last)
      is_increasing = FALSE;

    /* store last computed values */
    logfx_last = logfx;

  }

  /* we have left the loop with the right boundary of the support of PDF
     make shure that we will never use iv for sampling. */
  iv->logAhat = -INFINITY;
  iv->Ahatr_fract = iv->sq = 0.;
  iv->Acum = INFINITY;
#ifdef DEBUG_STORE_IP 
  iv->ip = iv->x;
#endif
  iv->next = NULL;         /* terminate list */
  --(GEN->n_ivs);          /* we do not count the terminating interval */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ars_starting_cpoints() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_starting_intervals( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute intervals for starting points                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_ars_interval *iv, *iv_new, *iv_tmp;
  double x, logfx;              /* construction point, value of logPDF at x  */
  
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(GEN->iv,UNUR_ERR_NULL);  COOKIE_CHECK(GEN->iv,CK_ARS_IV,UNUR_ERR_COOKIE);
  
  /* compute paramters for all intervals */
  for( iv=GEN->iv; iv->next != NULL; ) {

    /* compute parameters for interval */
    switch (_unur_ars_interval_parameter(gen, iv)) {
    case UNUR_SUCCESS:      /* computation of parameters for interval successful */
      /* skip to next interval */
      iv = iv->next;
      continue;
    case UNUR_ERR_INF:      /* interval unbounded */
      /* split interval */
      break;
    case UNUR_ERR_SILENT:   /* construction points too close */
      /* we have to remove this last interval from list */
      /* (the last construction point in the list is a boundary point.
	 thus we might change the domain of the distribution.
	 however, we only cut off a piece that is beyond the precesion
	 of the floating point arithmetic.)  */
      iv_tmp = iv->next;
      iv->next = iv->next->next;
      free(iv_tmp);
      --(GEN->n_ivs);
      
      if (iv->next==NULL) {
	/* last (virtuel) interval in list.
	   make sure that we will never use this segment */
	iv->logAhat = -INFINITY;
	iv->Ahatr_fract = iv->sq = 0.;
	iv->Acum = INFINITY;
      }
      continue;
    default:     /* PDF not T-concave */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area below hat infinite.
       insert new construction point. */
    x = _unur_arcmean(iv->x,iv->next->x);  /* use mean point in interval */

    /* value of logPDF at x */
    logfx = logPDF(x);

    /* add a new interval, but check if we had to used too many intervals */
    if (GEN->n_ivs >= GEN->max_ivs) {
      /* we do not want to create too many intervals */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded hat!");
      return UNUR_ERR_GEN_CONDITION;
    }
    iv_new = _unur_ars_interval_new( gen, x, logfx );
    if (iv_new == NULL) return UNUR_ERR_GEN_DATA;  /* logPDF(x) overflow */


    /* if fx is 0, then we can cut off the tail of the distribution
       (since it must be T-concave)  */
    if (!_unur_isfinite(logfx) ) {

      if (!_unur_isfinite(iv->logfx) ) {
	/* cut off left tail */
	iv_new->next = iv->next;
	free(iv);
	--(GEN->n_ivs);
	GEN->iv = iv_new;
	iv = iv_new;
      }

      else if (!_unur_isfinite(iv->next->logfx) ) {
	/* cut off right tail */
	free(iv->next);
	--(GEN->n_ivs);
	iv->next = iv_new;
      }

      else {
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave!");
	free(iv_new);
	return UNUR_ERR_GEN_CONDITION;
      }
    }

    else {
      /* insert new interval into linked list */
      iv_new->next = iv->next;
      iv->next = iv_new;
    }

  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ars_starting_intervals() */

/*---------------------------------------------------------------------------*/

struct unur_ars_interval *
_unur_ars_interval_new( struct unur_gen *gen, double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* get new interval and compute left construction point at x.           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   x       ... left point of new interval                             */
     /*   logfx   ... value of logPDF at x                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new interval                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_ars_interval *iv;
/*   double dfx; */

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,NULL);

  /* first check logfx */
  if (!(logfx < INFINITY)) {
    /* overflow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"logPDF(x) overflow");
    return NULL;
  }

  /* we need a new segment */
  iv = _unur_xmalloc( sizeof(struct unur_ars_interval) );
  iv->next = NULL; /* add eol marker */
  ++(GEN->n_ivs);   /* increment counter for intervals */
  COOKIE_SET(iv,CK_ARS_IV);

  /* avoid uninitialized variables */
  iv->logAhat = -INFINITY;
  iv->Acum = iv->Ahatr_fract = 0.;
  iv->sq = 0.;
#ifdef DEBUG_STORE_IP 
  iv->ip = 0.;
#endif

  /* make left construction point in interval */
  iv->x = x;              /* point x */
  iv->logfx = logfx;      /* value of logPDF at x */

  /* derivative of transformed density */
  iv->dlogfx = _unur_isfinite(logfx) ? dlogPDF(x) : INFINITY;
  
  /* the program requires dlogPDF > -INFINITY */
  if ( !(iv->dlogfx > -INFINITY))
    iv->dlogfx = INFINITY;

  return iv;

} /* end of _unur_ars_interval_new() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_interval_parameter( struct unur_gen *gen, struct unur_ars_interval *iv )
     /*----------------------------------------------------------------------*/
     /* compute intersection point of tangents and                           */
     /* the area below the hat  (Gilks & Wild variant)                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   iv   ... pointer to interval                                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... do not add this construction point             */
     /*   UNUR_ERR_INF    ... area = INFINITY                                */
     /*   others          ... error (PDF not T-concave)                      */
     /*----------------------------------------------------------------------*/
{
  double logAhatl, logAhatr;  /* log of areas below hat at l.h.s. and r.h.s. 
				 of intersection point, resp.                */

  double ip = 0.;             /* intersection point of tangents              */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE); 

  /* check interval on the right side of iv */
  CHECK_NULL(iv->next,UNUR_ERR_NULL);  COOKIE_CHECK(iv->next,CK_ARS_IV,UNUR_ERR_COOKIE); 

  /* get intersection point of tangents.
     used to partition interval into left hand part (construction point of tangent
     on the left hand boundary) and right hand part (construction point of tangent
     on the left hand boundary). */
  if ( _unur_ars_tangent_intersection_point(gen,iv,&ip)!=UNUR_SUCCESS )
    return UNUR_ERR_GEN_CONDITION;

#ifdef DEBUG_STORE_IP 
  iv->ip = ip;
#endif
 
  /* squeeze and area below squeeze */
  if (_unur_isfinite(iv->logfx) && _unur_isfinite(iv->next->dlogfx) ) {

    /* we do not compute the slope when the construction points
       are too close. at least 8 significant digits should remain. */
    if (_unur_FP_approx(iv->x, iv->next->x) )
      return UNUR_ERR_SILENT;   /* construction points too close */

    /* slope of transformed squeeze */
    iv->sq = (iv->next->logfx - iv->logfx) / (iv->next->x - iv->x);

    /* check squeeze */
    /* we have to take care about round off error.
       the following accepts PDFs with might be a little bit not T_concave */
    if ( ( (iv->sq > iv->dlogfx      && (!_unur_FP_approx(iv->sq,iv->dlogfx)) ) ||
	   (iv->sq < iv->next->dlogfx && (!_unur_FP_approx(iv->sq,iv->next->dlogfx)) ) )
	 && iv->next->dlogfx < INFINITY ) {   

      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"Squeeze too steep/flat. PDF not T-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  else {  /* no squeeze */
    iv->sq = -INFINITY;
  }

  /* volume below hat */
  logAhatl = _unur_ars_interval_logarea( gen, iv, iv->dlogfx, ip);
  logAhatr = _unur_ars_interval_logarea( gen, iv->next, iv->next->dlogfx, ip);

  /* areas below head unbounded ? */
  if (! (logAhatl < INFINITY && logAhatr < INFINITY) )
    return UNUR_ERR_INF;

  /* total area */
  iv->logAhat = (logAhatl > logAhatr) 
    ? logAhatl+log(1+exp(logAhatr-logAhatl)) 
    : logAhatr+log(1+exp(logAhatl-logAhatr)) ;       /* = log( Ahatr + Ahatl )     */

  iv->Ahatr_fract = 1./(1.+exp(logAhatl-logAhatr));  /*  = Ahatr / (Ahatr + Ahatl) */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ars_interval_parameter() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_interval_split( struct unur_gen *gen, struct unur_ars_interval *iv_oldl, double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* split interval iv_oldl into two intervals at point x                 */
     /*   old interval -> left hand side                                     */
     /*   new interval -> right hand side                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   iv_oldl ... pointer to interval                                    */
     /*   x       ... left point of new segment                              */
     /*   logfx   ... value of logPDF at x                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... if successful                                  */
     /*   UNUR_ERR_SILENT ... if no intervals are splitted                   */
     /*   others          ... error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_ars_interval *iv_newr;  /* pointer to new interval */
  struct unur_ars_interval iv_bak;    /* space for backing up data of interval */
  int success, success_r;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv_oldl,UNUR_ERR_NULL);  COOKIE_CHECK(iv_oldl,CK_ARS_IV,UNUR_ERR_COOKIE);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & ARS_DEBUG_SPLIT)
    _unur_ars_debug_split_start( gen,iv_oldl,x,logfx );
#endif

  /* the splitting point must be inside the interval */
  if (x < iv_oldl->x || x > iv_oldl->next->x) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"splitting point not in interval!");
    return UNUR_ERR_SILENT;
  }

  /* back up data */
  memcpy(&iv_bak, iv_oldl, sizeof(struct unur_ars_interval));

  /* check if the new interval is completely outside the support of PDF */
  if (!_unur_isfinite(logfx)) {
    
    /* one of the two boundary points must be 0, too! */
    if (!_unur_isfinite(iv_oldl->logfx)) {
      /* chop off left part (it's out of support) */
      iv_oldl->x = x;
    }
    else if (!_unur_isfinite(iv_oldl->next->logfx)) {
      /* chop off right part (it's out of support) */
      iv_oldl->next->x = x;
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not log-concave");
      return UNUR_ERR_GEN_CONDITION;
    }
    
    /* compute parameters for chopped interval */
    success = _unur_ars_interval_parameter(gen, iv_oldl);
    
    /* we did not add a new interval */
    iv_newr = NULL;
  }

  else {
    
    /* we need a new interval */
    iv_newr = _unur_ars_interval_new( gen, x, logfx );
    if (iv_newr == NULL) {
      /* logPDF(x) overflow */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    
    /* insert into linked list */
    iv_newr->next = iv_oldl->next;
    iv_oldl->next = iv_newr;
    
    /* compute parameters for interval */
    success   = _unur_ars_interval_parameter(gen, iv_oldl);
    success_r = _unur_ars_interval_parameter(gen, iv_newr);
    
    /* worst of success and success_r */
    if (success_r!=UNUR_SUCCESS)
      if ((success_r!=UNUR_ERR_SILENT&&success_r!=UNUR_ERR_INF) ||
	  (success==UNUR_SUCCESS||success==UNUR_ERR_SILENT||success==UNUR_ERR_INF))
	success = success_r;
  }

  /* successfull ? */
  if (success!=UNUR_SUCCESS) {
    /* cannot split interval at given point */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split interval at given point.");
    if (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not log-concave");

    /* the case of unbounded hat is treated as round-off error for
       very steep tangents. so we simply do not add this construction point. */

    /* restore old interval */
    memcpy(iv_oldl, &iv_bak, sizeof(struct unur_ars_interval));

    /* decrement counter for intervals and free unused interval */
    if (iv_newr) {
      --(GEN->n_ivs);
      free( iv_newr );
    }

  return ( (success!=UNUR_ERR_SILENT && success!=UNUR_ERR_INF)
	   ? UNUR_ERR_GEN_CONDITION : UNUR_SUCCESS );
  }

  /* successful */

#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & ARS_DEBUG_SPLIT) {
    /* update total area below hat and squeeze */
    GEN->Atotal = ( GEN->Atotal - scaled_area(&iv_bak) + scaled_area(iv_oldl) + ((iv_newr) ? scaled_area(iv_newr) : 0.) );
    /* write info into LOG file */
    _unur_ars_debug_split_stop( gen,iv_oldl,iv_newr );
  }
#endif

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ars_interval_split() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_tangent_intersection_point( struct unur_gen *gen, struct unur_ars_interval *iv, double *ipt )
     /*----------------------------------------------------------------------*/
     /* compute cutting point of interval into left and right part.          */
     /* (1) use intersection point of tangents of transformed hat.           */
     /* (2) use mean point if (1) is unstable due to roundoff errors.        */
     /* (3) use boundary point which is closer to the mode. this is          */
     /*     important when the transformed tagents are extremely steep.      */
     /*     (This might cause a serious roundoff error while sampling.)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   iv  ... pointer to interval                                        */
     /*   ipt ... pointer to intersection point                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(iv,UNUR_ERR_NULL);   COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE); 

  /*
     case: there is no tangent at one of the boundary points of the interval
           (then the slope is INFINITY)
     or
     case: the tangents are too steep  (--> case (3))
  */
  if ( iv->dlogfx > 1.e+140 ) {
    *ipt = iv->x;        /* intersection point = left boundary of interval */
    return UNUR_SUCCESS;
  }
  if ( iv->next->dlogfx < -1.e+140 || _unur_FP_is_infinity(iv->next->dlogfx)) {
    *ipt = iv->next->x;   /* intersection point = right boundary of interval */
    return UNUR_SUCCESS;
  }
  /** TODO: 1.e+140 (= sqrt(DBL_MAX) / 1.e15) is arbitrary  **/

  /* test for T-concavity */
  if ( _unur_FP_less( iv->dlogfx, iv->next->dlogfx ) ) {

    /* it might happen because of round-off errors
       that iv->next->dTfx is almost zero although it should be large.
       thus we ignore this case. */
    if ( fabs(iv->dlogfx) < DBL_EPSILON * fabs(iv->next->dlogfx) ) {
      *ipt = iv->x;        /* intersection point = left boundary of interval */
      iv->dlogfx = INFINITY;
      return UNUR_SUCCESS;
    }
    else if ( fabs(iv->next->dlogfx) < DBL_EPSILON * fabs(iv->dlogfx) ) {
      *ipt = iv->next->x;   /* intersection point = right boundary of interval */
      iv->next->dlogfx = INFINITY;
      return UNUR_SUCCESS;
    }
    else {
      if (_unur_FP_approx(iv->dlogfx, iv->next->dlogfx)) {
        /* use mean point */
        *ipt = 0.5 * (iv->x + iv->next->x);
        return UNUR_SUCCESS;
      }
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"dTfx0 < dTfx1 (x0<x1). PDF not log-concave!");
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /** TODO: the following test is too sensitve to roundoff errors **/
  /*    if (iv->next->logfx > iv->x + iv->dlogfx*(iv->next->x - iv->x)) { */
  /*      _unur_warning(gen->genid,UNUR_ERR_INIT,"tangent below PDF, not log-concave!"); */
  /*      return UNUR_ERR_INIT; */
  /*    } */
  
  /* case (2): computing intersection of tangents is unstable */
  if (_unur_FP_approx(iv->dlogfx, iv->next->dlogfx)) {
    /* use mean point */
    *ipt = 0.5 * (iv->x + iv->next->x);
    return UNUR_SUCCESS;
  }

  /* case (1): compute intersection point of tangents (regular case) */
  *ipt = ( (iv->next->logfx - iv->logfx - iv->next->dlogfx * iv->next->x + iv->dlogfx * iv->x) /
	   (iv->dlogfx - iv->next->dlogfx) );

  /* check position of intersection point */
  if (_unur_FP_less(*ipt, iv->x) || _unur_FP_greater(*ipt, iv->next->x))
    /* intersection point of tangents not in interval.
       This is mostly the case for numerical reasons.
       Thus we use the center of the interval instead.
       if the PDF not T-concave, it will catched at a later
       point when we compare slope of tangents and squeeze. */
    *ipt = 0.5 * (iv->x + iv->next->x);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_ars_tangent_intersection_point() */

/*---------------------------------------------------------------------------*/

double
_unur_ars_interval_logarea( struct unur_gen *gen ATTRIBUTE__UNUSED, 
			      struct unur_ars_interval *iv, double slope, double x )
     /*---------------------------------------------------------------------------*/
     /* compute log of area below piece of hat or squeeze in                      */
     /* interval [iv->x,x] or [x,iv->x]                                           */
     /*                                                                           */
     /* parameters:                                                               */
     /*   gen   ... pointer to generator object                                   */
     /*   iv    ... pointer to interval that stores construction point of tangent */
     /*   slope ... slope of tangent or secant of transformed PDF                 */
     /*   x     ... boundary of integration domain                                */
     /*                                                                           */
     /* return:                                                                   */
     /*   log of area                                                             */
     /*                                                                           */
     /* error:                                                                    */
     /*   return INFINITY                                                         */
     /*                                                                           */
     /* comment:                                                                  */
     /*   x0    ... construction point of tangent (= iv->x)                       */
     /*                                                                           */
     /*   area = | \int_{x0}^x \exp(Tf(x0) + slope*(t-x0)) dt |                   */
     /*        = f(x0) * |x - x0|                              if slope = 0       */
     /*        = | f(x0)/slope * (\exp(slope*(x-x0))-1) |      if slope != 0      */
     /*                                                                           */
     /*---------------------------------------------------------------------------*/
{
  double x0, logfx0;
  double logxdiff;
  double t, logt;

  /* check arguments */
  CHECK_NULL(iv,INFINITY);   COOKIE_CHECK(iv,CK_ARS_IV,INFINITY); 

  /* length of interval > 0 ? */
  if (_unur_FP_same(x, iv->x))
    return -INFINITY;

  /* if the construction point is at infinity, we cannot compute an area.
     (in this case we should have x == iv->x == INFINITY). */
  if (!_unur_isfinite(iv->x)) 
    return INFINITY;

  /* unbounded? */
  if ( !_unur_isfinite(slope)    ||
       (_unur_FP_is_minus_infinity(x) && slope<=0.) ||
       (_unur_FP_is_infinity(x)       && slope>=0.)  )   /* we have set (Tf)'(x) = INFINITY, if f(x)=0 */
    return INFINITY;

  /* construction point x0 of tangent and log of PDF at x0 */
  x0 = iv->x;
  logfx0 = iv->logfx;

  /* log of |x - x=| */
  logxdiff = log(fabs(x - x0));

  /* case: hat/squeeze constant --> area = f(x0) * |x - x0| */
  if (_unur_iszero(slope))
    return (_unur_isfinite(x) ? logfx0 + logxdiff : INFINITY);

  /* case: domain unbounded --> area = f(x0) / |slope| */
  if (!_unur_isfinite(x))
    return (logfx0 - log(fabs(slope)));

  /* case bounded domain --> area = | f(x0)/slope * (\exp(slope*(x-x0))-1) | */
  /* have to deal with numerical problems when x \approx x0                  */
  t = slope * (x - x0);
  logt = log(fabs(slope)) + logxdiff;

  if (fabs(t) > 1.e-6) {
    if (t > MAXLOG / 10.)
      return ( logfx0 + logxdiff + t - logt );
    else
      return ( logfx0 + logxdiff + log( fabs(exp(t) - 1.) ) - log(fabs(t)) );
  }

  else  /* use Taylor series */
    return (logfx0 + logxdiff + log1p(t/2. + t*t/6.));

} /* end of _unur_ars_interval_logarea() */

/*---------------------------------------------------------------------------*/

int
_unur_ars_make_area_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make table of areas and compute largest area for rescaling           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_ars_interval *iv;
  double Acum;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_ARS_GEN,UNUR_ERR_COOKIE);

  /* first we need the maximum area in intervals as scaling factor */
  GEN->logAmax = -INFINITY;
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE);
    if (GEN->logAmax < iv->logAhat)
      GEN->logAmax = iv->logAhat;
  }

  /* cumulated areas in intervals */
  Acum = 0.;            /* area below hat */
  for (iv = GEN->iv; iv != NULL; iv = iv->next ) {
    COOKIE_CHECK(iv,CK_ARS_IV,UNUR_ERR_COOKIE);
    Acum += scaled_area(iv);
    iv->Acum = Acum;
  }

  /* total area below hat */
  GEN->Atotal = Acum;

  return UNUR_SUCCESS;
} /* end of _unur_ars_make_area_table() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_init_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print after (almost empty generator) object has been created.        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ARS (Adaptive Rejection Sampling)\n",gen->genid);
  fprintf(LOG,"%s: transformation T_c(x) = log(x)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  if (gen->distr_is_privatecopy)
    fprintf(LOG,"%s: use private copy of distribution object\n",gen->genid);
  else
    fprintf(LOG,"%s: use pointer to external distribution object (dangerous!)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_ars_sample",gen->genid);
  if (gen->variant & ARS_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: maximum number of intervals        = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(gen,ARS_SET_MAX_IVS);
  fprintf(LOG,"\n");

  fprintf(LOG,"%s: maximum number of iterations       = %d",gen->genid,GEN->max_iter);
  _unur_print_if_default(gen,ARS_SET_MAX_ITER);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fprintf(LOG,"%s: number of starting points = %d",gen->genid,GEN->n_starting_cpoints);
  _unur_print_if_default(gen,ARS_SET_N_CPOINTS);
  fprintf(LOG,"\n%s: starting points:",gen->genid);
  if (gen->set & ARS_SET_CPOINTS)
    for (i=0; i<GEN->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",GEN->starting_cpoints[i]);
    }
  else
    fprintf(LOG," use \"equidistribution\" rule [default]");
  fprintf(LOG,"\n%s:\n",gen->genid);
  
  fflush(LOG);

} /* end of _unur_ars_debug_init_start() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_init_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into LOG file                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  _unur_ars_debug_intervals(gen,"INIT completed",TRUE);

  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_ars_debug_init_finished() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_reinit_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before reinitialization into LOG file     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int i;
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: *** Re-Initialize generator object ***\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  if (gen->set & ARS_SET_N_PERCENTILES) {
    fprintf(LOG,"%s: use percentiles of old hat as starting points for new hat:",gen->genid);
    for (i=0; i<GEN->n_percentiles; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",GEN->percentiles[i]);
    }
    fprintf(LOG,"\n%s: starting points:",gen->genid);
    for (i=0; i<GEN->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",GEN->starting_cpoints[i]);
    }
    fprintf(LOG,"\n");
  }
  else {
    fprintf(LOG,"%s: use starting points given at init\n",gen->genid);
  }
  
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_ars_debug_reinit_start() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_reinit_retry( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before second trial of reinitialization   */
     /* into LOG file.                                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: *** Re-Initialize failed  -->  second trial ***\n",gen->genid);
  fprintf(LOG,"%s: use equal-area-rule with %d points\n",gen->genid,GEN->retry_ncpoints);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_ars_debug_reinit_retry() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_reinit_finished( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after reinitialization into LOG file      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  _unur_ars_debug_intervals(gen," *** Generator reinitialized ***",TRUE);

  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_ars_debug_reinit_finished() */

/*---------------------------------------------------------------------------*/

void 
_unur_ars_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas )
     /*----------------------------------------------------------------------*/
     /* write list of intervals into LOG file (orig. variant by Gilks & Wild)*/
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen         ... pointer to generator object                        */
     /*   header      ... header for table                                   */
     /*   print_areas ... whether table of areas should be printed           */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  struct unur_ars_interval *iv;
  double Ahat, Ahatl, Ahatr;
  double sAhatl, sAhatr, Atotal, logAmax;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  if (header) fprintf(LOG,"%s:%s\n",gen->genid,header);

  fprintf(LOG,"%s:Intervals: %d\n",gen->genid,GEN->n_ivs);
  if (GEN->iv) {
    if (gen->debug & ARS_DEBUG_IV) {
#ifdef DEBUG_STORE_IP 
      fprintf(LOG,"%s: Nr.            tp            ip       logf(tp)     dlogf(tp)       squeeze\n",gen->genid);
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
	fprintf(LOG,"%s:[%3d]: %#12.6g  %#12.6g  %#12.6g  %#12.6g  %#12.6g\n", gen->genid, i,
		iv->x, iv->ip, iv->logfx, iv->dlogfx, iv->sq);
      }
#else
      fprintf(LOG,"%s: Nr.            tp       logf(tp)     dlogf(tp)       squeeze\n",gen->genid);
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
	fprintf(LOG,"%s:[%3d]: %#12.6g  %#12.6g  %#12.6g  %#12.6g\n", gen->genid, i,
		iv->x, iv->logfx, iv->dlogfx, iv->sq);
      }
#endif
      COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
      fprintf(LOG,"%s:[...]: %#12.6g                %#12.6g  %#12.6g\n", gen->genid,
	      iv->x, iv->logfx, iv->dlogfx);
    }
    fprintf(LOG,"%s:\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: No intervals !\n",gen->genid);

  if (!print_areas || GEN->Atotal <= 0.) return;

  /* print and sum areas below squeeze and hat */
  Atotal = GEN->Atotal;
  logAmax = GEN->logAmax;
  if (gen->debug & ARS_DEBUG_IV) {
    fprintf(LOG,"%s:Areas in intervals relative to maximum:\t[ log(A_max) = %g ]\n",gen->genid, logAmax);
    fprintf(LOG,"%s: Nr.\tbelow hat (left and right)\t\t   cumulated\n",gen->genid);
    sAhatl = sAhatr = 0.;
    if (GEN->iv) {
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID); 
	Ahat = scaled_area(iv);
	sAhatr += Ahatr = Ahat * iv->Ahatr_fract;
	sAhatl += Ahatl = Ahat - Ahatr;
	fprintf(LOG,"%s:[%3d]: %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		Ahatl, Ahatr, Ahat * 100. / Atotal,
		iv->Acum, iv->Acum * 100. / Atotal);
      }
      fprintf(LOG,"%s:       ------------------------  ---------  +\n",gen->genid);
      fprintf(LOG,"%s: Sum :        %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / Atotal);
      fprintf(LOG,"%s:\n",gen->genid);
    }
  }

  /* summary of areas */
  fprintf(LOG,"%s: A(total) = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_ars_debug_intervals() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into LOG file           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  if (gen->status == UNUR_SUCCESS) {
    fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    _unur_ars_debug_intervals(gen,NULL,TRUE);
  }
  else {
    fprintf(LOG,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
    _unur_ars_debug_intervals(gen,"Intervals after failure:",FALSE);
  }
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_ars_debug_free() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_split_start( const struct unur_gen *gen, 
			     const struct unur_ars_interval *iv,
			     double x, double logfx )
     /*----------------------------------------------------------------------*/
     /* write info about splitting interval                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   iv    ... pointer to interval                                      */
     /*   x     ... split at this point                                      */
     /*   logfx ... value of logPDF at x                                     */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  double Ahat;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  CHECK_NULL(iv,RETURN_VOID);   COOKIE_CHECK(iv,CK_ARS_IV,RETURN_VOID);

  LOG = unur_get_stream();

  Ahat = scaled_area(iv);

  fprintf(LOG,"%s: split interval at x = %g \t\tlogf(x) = %g\n",gen->genid,x,logfx);
  fprintf(LOG,"%s: old interval:\n",gen->genid);
  fprintf(LOG,"%s:   left  construction point = %-12.6g\tlogf(x) = %-12.6g\n",gen->genid,iv->x,iv->logfx);
  fprintf(LOG,"%s:   right construction point = %-12.6g\tlogf(x) = %-12.6g\n",gen->genid,iv->next->x,iv->next->logfx);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\t[ relative to A_max ]\n",gen->genid,
	  Ahat * (1.-iv->Ahatr_fract), Ahat * iv->Ahatr_fract, Ahat*100./GEN->Atotal);

  fflush(LOG);

} /* end of _unur_ars_debug_split_start() */

/*---------------------------------------------------------------------------*/

void
_unur_ars_debug_split_stop( const struct unur_gen *gen, 
			    const struct unur_ars_interval *iv_left, 
			    const struct unur_ars_interval *iv_right )
     /*----------------------------------------------------------------------*/
     /* write info about new splitted intervals                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iv_left  ... pointer to new left hand interval                     */
     /*   iv_right ... pointer to new right hand interval                    */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  double Ahat;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);       COOKIE_CHECK(gen,CK_ARS_GEN,RETURN_VOID);
  CHECK_NULL(iv_left,RETURN_VOID);   COOKIE_CHECK(iv_left,CK_ARS_IV,RETURN_VOID);

  if (iv_right == NULL) iv_right = iv_left;

  LOG = unur_get_stream();

  fprintf(LOG,"%s: inserted point:\n",gen->genid);
  fprintf(LOG,"%s: x = %g, logf(x) = %g, dlogf(x) = %g, squeeze = %g:\n",
	  gen->genid, iv_right->x, iv_right->logfx, iv_right->dlogfx, iv_right->sq);
  fprintf(LOG,"%s: new intervals:\n",gen->genid);
  fprintf(LOG,"%s:   left   construction point = %g\n",gen->genid, iv_left->x);
  if (iv_left != iv_right)
    fprintf(LOG,"%s:   middle construction point = %g\n",gen->genid, iv_right->x);
  fprintf(LOG,"%s:   right  construction point = %g\n",gen->genid, iv_right->next->x);

  fprintf(LOG,"%s: left interval:\n",gen->genid);
  Ahat = scaled_area(iv_left);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\t[ relative to A_max ]\n",gen->genid,
	  Ahat * (1.-iv_left->Ahatr_fract), Ahat * iv_left->Ahatr_fract,
	  Ahat * 100./GEN->Atotal);

  if (iv_left == iv_right)
    fprintf(LOG,"%s: interval chopped.\n",gen->genid);
  else {
    fprintf(LOG,"%s: right interval:\n",gen->genid);
    Ahat = scaled_area(iv_right);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\t[ relative to A_max ]\n",gen->genid,
	    Ahat * (1.-iv_right->Ahatr_fract), Ahat * iv_right->Ahatr_fract,
	    Ahat * 100./GEN->Atotal);
  }

  fprintf(LOG,"%s: total areas:\n",gen->genid);
  fprintf(LOG,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_ars_debug_split_stop() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_ars_info( struct unur_gen *gen, int help )
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
  int n_ivs_bak;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);

  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = logPDF dlogPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");
      
  /* method */
  _unur_string_append(info,"method: ARS (Transformed Density Rejection -- Gilks&Wild variant)\n");
  _unur_string_append(info,"   T_c(x) = log(x)  ... c = 0\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g  [ log = %g ]\n", 
		      GEN->Atotal*exp(GEN->logAmax), log(GEN->Atotal) + GEN->logAmax);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", GEN->Atotal*exp(GEN->logAmax)/DISTR.area);
  else {
    n_ivs_bak = GEN->n_ivs;
    GEN->n_ivs = GEN->max_ivs+1;
    _unur_string_append(info,"= %.3f  [approx.]\n",
			unur_test_count_urn(gen,samplesize,0,NULL)/(2.*samplesize));
    GEN->n_ivs = n_ivs_bak;
  }
  _unur_string_append(info,"   # intervals = %d\n", GEN->n_ivs);
  _unur_string_append(info,"\n");
  
  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   cpoints = %d  %s\n", GEN->n_starting_cpoints,
			(gen->set & ARS_SET_N_CPOINTS) ? "" : "[default]");

    if (gen->variant & ARS_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    if (gen->variant & ARS_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");
    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_ars_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
       int unur_ars_set_cpoints( UNUR_PAR *parameters, int n_cpoints, const double *cpoints );
       int unur_ars_set_reinit_percentiles( UNUR_PAR *parameters, int n_percentiles, const double *percentiles );
       int unur_ars_set_reinit_ncpoints( UNUR_PAR *parameters, int ncpoints );
    */
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_ars_info() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
