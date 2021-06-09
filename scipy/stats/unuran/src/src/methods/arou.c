/*****************************************************************************
 *                                                                           *
 *          unuran -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      arou.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    ratio-of-uniforms with enveloping polygon                    *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a T-concave distribution;                               *
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
 *   REFERENCES:                                                             *
 *   [1] Leydold J. (2000): Automatic Sampling with the ratio-of-uniforms    *
 *       method, ACM TOMS, forthcoming                                       *
 *                                                                           *
 *   [2] Kinderman, A.J. and Monahan, F.J. (1977): Computer generation of    *
 *       random variables using the ratio of uniform deviates,               *
 *       ACM Trans. Math. Software 3(3), pp. 257--260.                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The ratio-of-uniforms method introduced in [2] is a flexible method that  *
 * is based on the following theorem:                                        *
 *                                                                           *
 * THEOREM:                                                                  *
 *    Let X be a random variable with density function f(x) = g(x) / G,      *
 *    where g(x) is a positive integrable function with support (x_0,x_1)    *
 *    not necessarily finite and G = integral g(x) dx.                       *
 *    If (V,U) is uniformly distributed in                                   *
 *       A = {(v,u): 0 < u <= sqrt(g(v/u)), x_0 < v/u < x_1},                *
 *    then X = V/U has probability density function f(x).                    *
 *                                                                           *
 * Generating point (V,U) uniformly distributed in A is done by rejection    *
 * from an enveloping region, usually from the minimal bounding rectangle.   *
 *                                                                           *
 * The implemented algorithm uses the fact, that for many distribtions,      *
 * A is convex. Then we easily can construct an enveloping polygon by means  *
 * of tangent lines and a squeeze region by means of secants.                *
 * The resulting algorithm is very fast, since we can sample from the        *
 * squeeze region with immedate acceptance (thus only one uniform random     *
 * number is necessary). The region between envelope and squeeze consists    *
 * of triangles and can be made arbitrarily small (see [1] for details).     *
 *                                                                           *
 * Distributions with a convex set A are characterized by the following      *
 * theorem that shows a connection to transformed density rejection TDR.     *
 *                                                                           *
 * THEOREM:                                                                  *
 *    A is convex if and only if g is T-concave with transformation          *
 *    T(x) = -1/sqrt(x), i.e., -1/sqrt(g(x)) is a concave function.          *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * The region A is divided into segments, that contain the origin. Each      *
 * segments consist of an inner triangle (the squeeze region) and the        *
 * outer triangle (the region between envelope and squeeze). We have to      *
 * compute the areas of these triangles.                                     *
 *                                                                           *
 * To generate from the distribution we have to sample from a discrete       *
 * random variate with probability vector proportional to these areas to get *
 * one of these triangles (inner or outer). We use indexed search (or guide  *
 * tables) to perform this task ([1], see also description of DIS).          *
 * When we have an inner triangle (squeeze), we reuse the uniform random     *
 * variate to get a point (v,u) uniformly distributed on the edge opposite   *
 * to the origin and return the ratio x = v/u (Thus generation from the      *
 * squeeze region is equivalent to the inversion method.)                    *
 * When whe have an outer triangle, we have to sample a point (v,u)          *
 * uniformly distributed in this triangle. If u <= g(v/u) then return the    *
 * ratio x = v/u, otherwise reject.                                          *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Algorithm AROU                                                            *
 *                                                                           *
 * [Required]                                                                *
 * PDF f(x), construction points c_1,...,c_n                                 *
 *                                                                           *
 * [Setup]                                                                   *
 *  1: Construct inner triangles S_i^s and outer triangles S_i^o             *
 *  2: Foreach triangle Do                                                   *
 *  3:    Compute respective areas A_i^s and A_i^o.                          *
 *                                                                           *
 * [Generate]                                                                *
 *  4: Generate I proportional to (A_1^s, A_1^o; A_2^s, A_2^o; ...).         *
 *  5: If inner triangle S_i^s Then                                          *
 *  6:    Compute point (v,u) on edge (Reuse u.r.n. from 4).                 *
 *  7:    Return v/u.                                                        *
 *  8: Else                                                                  *
 *  9:    Generate point (V,U) uniformly distributed in outer triangle S_i^o.*
 * 10:    If u <= f(v/u) Return V/U.                                         *
 * 11:    Goto 4.                                                            *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "arou.h"
#include "arou_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define AROU_VARFLAG_VERIFY     0x01u   /* flag for verifying mode           */
#define AROU_VARFLAG_USECENTER  0x02u   /* flag whether center is used as cpoint or not */
#define AROU_VARFLAG_PEDANTIC   0x04u   /* whether pedantic checking is used */
#define AROU_VARFLAG_USEDARS    0x10u   /* whether DARS is used in setup or not */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define AROU_DEBUG_SEGMENTS     0x00000010u   /* print list of segments      */
#define AROU_DEBUG_SPLIT        0x00010000u   /* trace splitting of segments */
#define AROU_DEBUG_DARS         0x00020000u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define AROU_SET_CENTER         0x001u
#define AROU_SET_STP            0x002u
#define AROU_SET_N_STP          0x004u
#define AROU_SET_GUIDEFACTOR    0x010u
#define AROU_SET_MAX_SQHRATIO   0x020u
#define AROU_SET_MAX_SEGS       0x040u
#define AROU_SET_USE_DARS       0x100u
#define AROU_SET_DARS_FACTOR    0x200u

/*---------------------------------------------------------------------------*/

#define GENTYPE "AROU"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_arou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_arou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_arou_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_arou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_arou_sample( struct unur_gen *gen );
static double _unur_arou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_arou_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create list of construction points for starting segments.                 */
/* if user has not provided such points compute these by means of the        */
/* "equi-angle rule".                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_arou_get_starting_segments( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute segments from given starting construction points.                 */
/*---------------------------------------------------------------------------*/

static double _unur_arou_compute_x( double v, double u );
/*---------------------------------------------------------------------------*/
/* compute point x from (v,u) tuple.                                         */
/*---------------------------------------------------------------------------*/

static int _unur_arou_run_dars( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* run derandomized adaptive rejection sampling.                             */
/*---------------------------------------------------------------------------*/

static struct unur_arou_segment *_unur_arou_segment_new( struct unur_gen *gen, double x, double fx );
/*---------------------------------------------------------------------------*/
/* make a new segment with left construction point x = v/u.                  */
/*---------------------------------------------------------------------------*/

static int _unur_arou_segment_parameter( struct unur_gen *gen, struct unur_arou_segment *seg );
/*---------------------------------------------------------------------------*/
/* compute all necessary data for segment.                                   */
/* return:                                                                   */
/*   UNUR_SUCCESS     ... on success                                         */
/*   UNUR_ERR_SILENT  ... do not add this construction point                 */
/*   UNUR_ERR_INF     ... if area = INFINITY                                 */
/*   other error code ... error (PDF not T-concave)                          */
/*---------------------------------------------------------------------------*/

static int _unur_arou_segment_split( struct unur_gen *gen, struct unur_arou_segment *seg_old, double x, double fx );
/*---------------------------------------------------------------------------*/
/* split a segment at point (direction) x.                                   */
/*---------------------------------------------------------------------------*/

static int _unur_arou_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static double _unur_arou_segment_arcmean( struct unur_arou_segment *seg );
/*---------------------------------------------------------------------------*/
/* compute the "arcmean" of the two construction points of a segement.       */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_init( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_dars_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print header before runniung derandomized adaptive rejection sampling.    */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_dars( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has run derandomized adaptive rejection sampling.   */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_segments( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for segments.                                                  */
/*---------------------------------------------------------------------------*/

static void _unur_arou_debug_split_start( const struct unur_gen *gen, 
					  const struct unur_arou_segment *seg,
					  double x, double fx );
static void _unur_arou_debug_split_stop( const struct unur_gen *gen, 
					 const struct unur_arou_segment *seg_left,
					 const struct unur_arou_segment *seg_right );
/*---------------------------------------------------------------------------*/
/* print before and after a segment has been split (not / successfully).     */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_arou_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_arou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_arou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))  /* call to PDF           */
#define dPDF(x)   _unur_cont_dPDF((x),(gen->distr)) /* call to derivative of PDF */

/*---------------------------------------------------------------------------*/

#define _unur_arou_getSAMPLE(gen) \
   ( ((gen)->variant & AROU_VARFLAG_VERIFY) \
     ? _unur_arou_sample_check : _unur_arou_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_arou_new( const struct unur_distr *distr )
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
  _unur_check_NULL(GENTYPE,distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.pdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PDF"); return NULL;
  }
  if (DISTR_IN.dpdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"derivative of PDF"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_arou_par) );
  COOKIE_SET(par,CK_AROU_PAR);

  /* copy input */
  par->distr              = distr;  /* pointer to distribution object        */

  /* set default values */
  PAR->guide_factor        = 2.;     /* size of guide table / number of intervals */
  PAR->darsfactor          = 0.99;   /* factor for (derandomized) ARS.
                                       do not add a new construction point in a segment, 
				       where abiguous region is too small, i.e. if
				       area / (|S^e\S^s|/number of segments) < darsfactor */

  PAR->starting_cpoints    = NULL;   /* pointer to array of starting points   */
  PAR->n_starting_cpoints  = 30;     /* number of starting points             */
  PAR->max_segs            = 100;    /* maximum number of segments            */
  PAR->max_ratio           = 0.99;   /* do not add construction points if
				       ratio r_n = |P^s| / |P^e| > max_ratio */

  par->method   = UNUR_METH_AROU;             /* method                      */
  par->variant  = ( AROU_VARFLAG_USECENTER |
		    AROU_VARFLAG_USEDARS );   /* default variant             */
  par->set      = 0u;                      /* inidicate default parameters   */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = par->urng;               /* no special auxilliary URNG     */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_arou_init;

  return par;

} /* end of unur_arou_new() */

/*****************************************************************************/

int
unur_arou_set_usedars( struct unur_par *par, int usedars )
     /*----------------------------------------------------------------------*/
     /* set flag for using DARS (derandomized adaptive rejection sampling).  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usedars   ... 0 = do not use,  !0 = use DARS                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using not using DARS is the default                                */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* we use a bit in variant */
  par->variant = (usedars) ? (par->variant | AROU_VARFLAG_USEDARS) : (par->variant & (~AROU_VARFLAG_USEDARS));

  /* changelog */
  par->set |= AROU_SET_USE_DARS;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_arou_set_usedars() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_darsfactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for derandomized adaptive rejection sampling              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... parameter for DARS                                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* check new parameter for generator */
  if (factor < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"DARS factor < 0");
    return UNUR_ERR_PAR_SET;
  }
    
  /* store date */
  PAR->darsfactor = factor;

  /* changelog */
  par->set |= AROU_SET_DARS_FACTOR;

  return UNUR_SUCCESS;

} /* end of unur_arou_set_darsfactor() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_cpoints( struct unur_par *par, int n_stp, const double *stp )
     /*----------------------------------------------------------------------*/
     /* set construction points for envelope                                 */
     /* and/or its number for initialization                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   n_stp  ... number of starting points                               */
     /*   stp    ... pointer to array of starting points                     */
     /*              (NULL for changing only the number of default points)   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* check starting construction points */
  /* we always use the boundary points as additional starting points,
     so we do not count these here! */
  if (n_stp < 0 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"number of starting points < 0");
    return UNUR_ERR_PAR_SET;
  }

  if (stp) 
    /* starting points must be strictly monontonically increasing */
    for( i=1; i<n_stp; i++ )
      if (stp[i] <= stp[i-1]) {
	_unur_warning(GENTYPE,UNUR_ERR_PAR_SET,
		      "starting points not strictly monotonically increasing");
	return UNUR_ERR_PAR_SET;
      }

  /* store date */
  PAR->starting_cpoints = stp;
  PAR->n_starting_cpoints = n_stp;

  /* changelog */
  par->set |= AROU_SET_N_STP | ((stp) ? AROU_SET_STP : 0);

  return UNUR_SUCCESS;

} /* end of unur_arou_set_cpoints() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"guide table size < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->guide_factor = factor;

  /* changelog */
  par->set |= AROU_SET_GUIDEFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_arou_set_guidefactor() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_max_sqhratio( struct unur_par *par, double max_ratio )
     /*----------------------------------------------------------------------*/
     /* set bound for ratio A(squeeze) / A(hat)                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_ratio ... upper bound for ratio to add a new construction point*/
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* check new parameter for generator */
  if (max_ratio < 0. || max_ratio > 1. ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"ratio A(squeeze)/A(hat) not in [0,1]");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_ratio = max_ratio;

  /* changelog */
  par->set |= AROU_SET_MAX_SQHRATIO;

  return UNUR_SUCCESS;

} /* end of unur_arou_set_max_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_arou_get_sqhratio( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get ratio A(squeeze) / A(hat)                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   ratio    ... on success                                            */
     /*   INFINITY ... on error                                              */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  _unur_check_gen_object( gen, AROU, INFINITY );

  return (GEN->Asqueeze / GEN->Atotal);

} /* end of unur_arou_get_sqhratio() */

/*---------------------------------------------------------------------------*/

double
unur_arou_get_hatarea( const struct unur_gen *gen )
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
  _unur_check_gen_object( gen, AROU, INFINITY );

  return GEN->Atotal;

} /* end of unur_arou_get_hatarea() */

/*---------------------------------------------------------------------------*/

double
unur_arou_get_squeezearea( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* get area below squeeze                                               */
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
  _unur_check_gen_object( gen, AROU, INFINITY );

  return GEN->Asqueeze;

} /* end of unur_arou_get_squeezearea() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_max_segments( struct unur_par *par, int max_segs )
     /*----------------------------------------------------------------------*/
     /* set maximum number of segments                                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   max_segs  ... maximum number of segments                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* check new parameter for generator */
  if (max_segs < 1 ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"maximum number of segments < 1");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->max_segs = max_segs;

  /* changelog */
  par->set |= AROU_SET_MAX_SEGS;

  return UNUR_SUCCESS;

} /* end of unur_arou_set_max_segments() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_usecenter( struct unur_par *par, int usecenter )
     /*----------------------------------------------------------------------*/
     /* set flag for using center as construction point                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   usecenter ... 0 = do not use,  !0 = use                            */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   using center as construction point is the default                  */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, AROU );

  /* we use a bit in variant */
  par->variant = (usecenter) ? (par->variant | AROU_VARFLAG_USECENTER) : (par->variant & (~AROU_VARFLAG_USECENTER));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_arou_set_usecenter() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_verify( struct unur_par *par, int verify )
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

  /* check input */
  _unur_check_par_object( par, AROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | AROU_VARFLAG_VERIFY) : (par->variant & (~AROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_arou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_arou_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, AROU, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= AROU_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~AROU_VARFLAG_VERIFY;

  SAMPLE = _unur_arou_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_arou_chg_verify() */

/*---------------------------------------------------------------------------*/

int
unur_arou_set_pedantic( struct unur_par *par, int pedantic )
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

  /* check input */
  _unur_check_par_object( par, AROU );

  /* we use a bit in variant */
  par->variant = (pedantic) ? (par->variant | AROU_VARFLAG_PEDANTIC) : (par->variant & (~AROU_VARFLAG_PEDANTIC));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_arou_set_pedantic() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_arou_init( struct unur_par *par )
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
  int i,k;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_AROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_AROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_arou_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* get starting points */
  if (_unur_arou_get_starting_cpoints(par,gen)!=UNUR_SUCCESS ) {
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_arou_debug_init(par,gen);
#endif
    _unur_par_free(par); _unur_arou_free(gen);
    return NULL;
  }

  /* compute segments for given starting points */
  if ( _unur_arou_get_starting_segments(gen)!=UNUR_SUCCESS ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug) _unur_arou_debug_init(par,gen);
#endif
    _unur_par_free(par); _unur_arou_free(gen);
    return NULL;
  }

  /* we have to update the maximal number of segments,
     if the user wants more starting points. */
  if (GEN->n_segs > GEN->max_segs) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"maximal number of segments too small. increase.");
    GEN->max_segs = GEN->n_segs;
  }

  if (gen->variant & AROU_VARFLAG_USEDARS) {
    /* run derandomized adaptive rejection sampling (DARS) */

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & AROU_DEBUG_DARS) {
      /* make initial guide table (only necessary for writing debug info) */
      _unur_arou_make_guide_table(gen);
      /* write info into LOG file */
      _unur_arou_debug_init(par,gen);
      _unur_arou_debug_dars_start(gen);
    }
#endif

    for (i=0; i<3; i++) {
      /* we make several tries */

      /* run DARS */
      if ( _unur_arou_run_dars(gen)!=UNUR_SUCCESS ) {
	_unur_par_free(par); _unur_arou_free(gen);
	return NULL;
      }
   
      /* make initial guide table */
      _unur_arou_make_guide_table(gen);

      /* check if DARS was completed */
      if (GEN->n_segs < GEN->max_segs) {
  	/* ran ARS instead */
	for (k=0; k<5; k++)
	  _unur_sample_cont(gen);
      }
      else
	break;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) {
      if (gen->debug & AROU_DEBUG_DARS)
	_unur_arou_debug_dars(gen);
      else 
  	_unur_arou_debug_init(par,gen);
    }
#endif
  }

  else { /* do not run DARS */
    /* make initial guide table */
    _unur_arou_make_guide_table(gen);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_arou_debug_init(par,gen);
#endif

  }

  /* free parameters */
  _unur_par_free(par);

  /* is there any envelope at all ? */
  if (GEN->Atotal <= 0. || !_unur_isfinite(GEN->Atotal)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"bad construction points");
    _unur_arou_free(gen);
    return NULL;
  }

  /* o.k. */
  return gen;

} /* end of _unur_arou_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_arou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_AROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_arou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_AROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_arou_getSAMPLE(gen);
  gen->destroy = _unur_arou_free;
  gen->clone = _unur_arou_clone;

  /* set all pointers to NULL */
  GEN->seg         = NULL;
  GEN->n_segs      = 0;
  GEN->guide       = NULL;
  GEN->guide_size  = 0;
  GEN->Atotal      = 0.;
  GEN->Asqueeze    = 0.;

  /* copy some parameters into generator object */
  GEN->guide_factor = PAR->guide_factor; /* relative size of guide tables      */

  /* bounds for adding construction points  */
  GEN->max_segs = PAR->max_segs;      /* maximum number of segments            */
#ifdef UNUR_ENABLE_INFO
  GEN->max_segs_info = PAR->max_segs;   /* ... for info string */
#endif
  GEN->max_ratio = PAR->max_ratio;    
  GEN->darsfactor = PAR->darsfactor;

  /* get center */
  if ( (gen->distr->set & UNUR_DISTR_SET_CENTER) ||
       (gen->distr->set & UNUR_DISTR_SET_MODE) ) {
    GEN->center = unur_distr_cont_get_center(gen->distr);
    /* center must be in domain */
    GEN->center = _unur_max(GEN->center,DISTR.BD_LEFT);
    GEN->center = _unur_min(GEN->center,DISTR.BD_RIGHT);
    gen->set |= AROU_SET_CENTER;
  }
  else {
    GEN->center = 0.;
    /* we cannot use the center as construction point */
    gen->variant &= ~AROU_VARFLAG_USECENTER;
  }

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_arou_info;
#endif

  /* return pointer to (almost empty) generator object */
  return(gen);

} /* end of _unur_arou_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_arou_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_arou_gen*)clone->datap)

  struct unur_gen *clone;
  struct unur_arou_segment *seg,*next, *clone_seg, *clone_prev;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy linked list of segments */
  clone_seg = NULL;
  clone_prev = NULL;
  for (seg = GEN->seg; seg != NULL; seg = next) {
    /* copy segment */
    clone_seg = _unur_xmalloc( sizeof(struct unur_arou_segment) );
    memcpy( clone_seg, seg, sizeof(struct unur_arou_segment) );
    if (clone_prev == NULL) {
      /* starting point of linked list */
      CLONE->seg = clone_seg;
    }
    else {
      /* insert into linked list */
      clone_prev->next = clone_seg;
      clone_prev->rtp  = clone_seg->ltp;
      clone_prev->drtp = clone_seg->dltp;
    }
    /* next step */
    next = seg->next;
    clone_prev = clone_seg;
  }
  /* terminate linked list */
  if (clone_seg) clone_seg->next = NULL;

  /* make new guide table */
  CLONE->guide = NULL;
  _unur_arou_make_guide_table(clone);

  /* finished clone */
  return clone;

#undef CLONE
} /* end of _unur_arou_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_arou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_AROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* write info into LOG file */
#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug) _unur_arou_debug_free(gen);
#endif

  /* free linked list of segments */
  {
    struct unur_arou_segment *seg,*next;
    for (seg = GEN->seg; seg != NULL; seg = next) {
      next = seg->next;
      free(seg);
    }
  }

  /* free other memory not stored in list */
  if (GEN->guide) free(GEN->guide);

  _unur_generic_free(gen);

} /* end of _unur_arou_free() */

/*****************************************************************************/

double
_unur_arou_sample( struct unur_gen *gen )
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
  /** TODO: check uniform random number: u != 0 and u != 1 ??  **/

  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_arou_segment *seg;
  int result_split;
  double R,R1,R2,R3,tmp,x,fx,u;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_AROU_GEN,INFINITY);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    R = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    seg =  GEN->guide[(int) (R * GEN->guide_size)];
    R *= GEN->Atotal;
    while (seg->Acum < R) {
      seg = seg->next;
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,INFINITY);

    /* reuse of uniform random number */
    R = seg->Acum - R;

    /* inside or outside squeeze */
    if (R < seg->Ain) {
      /* inside */
      /* reuse of random number.                       */
      /* We can avoid R = (seg->Ain - R) / seg->Ain    */
      return( ( seg->Ain * seg->rtp[0] + R * (seg->ltp[0] - seg->rtp[0]) ) /
	      ( seg->Ain * seg->rtp[1] + R * (seg->ltp[1] - seg->rtp[1]) ) );
    }

    else {
      /* outside */

      /* from now on we use the auxilliary generator
	 (it can be the same as the main generator) */
      urng = gen->urng_aux;

      /* three uniform random numbers with R1 + R2 + R3 = 1 */
      R1 = (R - seg->Ain) / seg->Aout;  /* reuse of random number (good ?? ) */
      R2 = _unur_call_urng(urng);
      if (R1>R2) { tmp = R1; R1=R2; R2=tmp; }  /* swap */
      R3 = 1.-R2;
      R2 -= R1;

      /* point (v,u) and ratio x = v/u */
      u = seg->ltp[1]*R1 + seg->rtp[1]*R2 + seg->mid[1]*R3;
      x = (seg->ltp[0]*R1 + seg->rtp[0]*R2 + seg->mid[0]*R3) / u;

      /* density at x */
      fx = PDF(x);

      /* being outside the squeeze is bad. improve the situation! */
      if (GEN->n_segs < GEN->max_segs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  result_split = _unur_arou_segment_split(gen,seg,x,fx);
	  if ( !(result_split == UNUR_SUCCESS || result_split == UNUR_ERR_SILENT) ) {
	    /* condition for PDF is violated! */
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	    if (gen->variant & AROU_VARFLAG_PEDANTIC) {
	      /* replace sampling routine by dummy routine that just returns INFINITY */
	      SAMPLE = _unur_sample_cont_error;
	      return INFINITY;
	    }
	  }
	  else {
	    /* splitting successful --> update guide table */ 
	    _unur_arou_make_guide_table(gen);
	  }
	}
	else 
	  /* no more construction points (avoid too many second if statements above) */
	  GEN->max_segs = GEN->n_segs;
      }

      /* if inside region of acceptance, return ratio x */
      if (u*u <= fx) 
	return x;
    }
  }
} /* end of _unur_arou_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_arou_sample_check( struct unur_gen *gen )
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
  /** TODO: check uniform random number: u != 0 and u != 1 ??  **/

  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_arou_segment *seg;
  int result_split;
  double R,R1,R2,R3,tmp,x,fx,u,sqx,a;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_AROU_GEN,INFINITY);

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U(0,1) */
    R = _unur_call_urng(urng);

    /* look up in guide table and search for segment */
    seg =  GEN->guide[(int) (R * GEN->guide_size)];
    R *= GEN->Atotal;
    while (seg->Acum < R) {
      seg = seg->next;
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,INFINITY);

    /* reuse of uniform random number */
    R = seg->Acum - R;

    /* inside or outside squeeze */
    if (R < seg->Ain) {
      /* inside */
      /* reuse of random number.                       */
      /* We can avoid R = (seg->Ain - R) / seg->Ain    */
      x = ( ( seg->Ain * seg->rtp[0] + R * (seg->ltp[0] - seg->rtp[0]) ) /
	    ( seg->Ain * seg->rtp[1] + R * (seg->ltp[1] - seg->rtp[1]) ) );

      /* density at x */
      fx = PDF(x);

      /* compute value of squeeze at x, i.e., we have to solve 
	 a*ltp[0] + (1-a)*rtp[0] == a*ltp[1] + (1-a)*rtp[1] */
      a = ( (seg->rtp[0] - x * seg->rtp[1]) / 
	    (seg->rtp[0] - seg->ltp[0] + x * (seg->ltp[1] - seg->rtp[1])) );
      sqx = a * seg->ltp[1] + (1.-a) * seg->rtp[1];

      /* test for T-concavity */
      if (sqx*sqx > fx * (1.+UNUR_EPSILON))
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave.");

      return x;
    }

    else {
      /* outside */

      /* from now on we use the auxilliary generator
	 (it can be the same as the main generator) */
      urng = gen->urng_aux;

      /* three uniform random numbers with R1 + R2 + R3 = 1 */
      R1 = (R - seg->Ain) / seg->Aout;  /* reuse of random number (good ?? ) */
      R2 = _unur_call_urng(urng);
      if (R1>R2) { tmp = R1; R1=R2; R2=tmp; }  /* swap */
      R3 = 1.-R2;
      R2 -= R1;

      /* point (v,u) and ratio x = v/u */
      u = seg->ltp[1]*R1 + seg->rtp[1]*R2 + seg->mid[1]*R3;
      x = (seg->ltp[0]*R1 + seg->rtp[0]*R2 + seg->mid[0]*R3) / u;

      /* density at x */
      fx = PDF(x);

      /* compute value of squeeze at x, i.e., we have to solve 
	 a*ltp[0] + (1-a)*rtp[0] == a*ltp[1] + (1-a)*rtp[1] */
      a = ( (seg->rtp[0] - x * seg->rtp[1]) / 
	    (seg->rtp[0] - seg->ltp[0] + x * (seg->ltp[1] - seg->rtp[1])) );
      sqx = a * seg->ltp[1] + (1.-a) * seg->rtp[1];

      /* test for T-concavity */
      if (sqx*sqx > fx)
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave.");

      /* being outside the squeeze is bad. improve the situation! */
      if (GEN->n_segs < GEN->max_segs) {
	if (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) {
	  result_split = _unur_arou_segment_split(gen,seg,x,fx);
	  if ( !(result_split == UNUR_SUCCESS || result_split == UNUR_ERR_SILENT) ) {
	    /* condition for PDF is violated! */
	    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	    if (gen->variant & AROU_VARFLAG_PEDANTIC) {
	      /* replace sampling routine by dummy routine that just returns INFINITY */
	      SAMPLE = _unur_sample_cont_error;
	      return INFINITY;
	    }
	  }
	  else {
	    /* splitting successful --> update guide table */ 
	    _unur_arou_make_guide_table(gen);
	  }
	}
	else 
	  /* no more construction points (avoid to many second if statement above */
	  GEN->max_segs = GEN->n_segs;
      }

      /* if inside region of acceptance, return ratio x */
      if (u*u <= fx) 
	return x;
    }
  }

} /* end of _unur_arou_sample_check() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_arou_get_starting_cpoints( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* list of construction points for starting segments.                   */
     /* if not provided as arguments compute these                           */
     /* by means of the "equi-angle rule".                                   */
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
  struct unur_arou_segment *seg, *seg_new;
  double left_angle, right_angle, diff_angle, angle;
  double x, x_last, fx, fx_last;
  int i, use_center, is_center, is_increasing;

  /* check arguments */
  CHECK_NULL(par,UNUR_ERR_NULL);  COOKIE_CHECK(par,CK_AROU_PAR,UNUR_ERR_COOKIE);
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);

  /* initialize boolean */
  is_center = FALSE;

  /* use center as construction point ? */
  use_center = (gen->variant & AROU_VARFLAG_USECENTER) ? TRUE : FALSE;

  /* reset counter of segments */
  GEN->n_segs = 0;

  /* prepare for computing construction points */
  if (!PAR->starting_cpoints) {
    /* move center into  x = 0 */
    /* angles of boundary of domain */
    left_angle =  ( DISTR.BD_LEFT  <= -INFINITY ) ? -M_PI/2. : atan(DISTR.BD_LEFT  - GEN->center);  
    right_angle = ( DISTR.BD_RIGHT >= INFINITY )  ? M_PI/2.  : atan(DISTR.BD_RIGHT - GEN->center);
    /* we use equal distances between the angles of the cpoints   */
    /* and the boundary points                                    */
    diff_angle = (right_angle-left_angle) / (PAR->n_starting_cpoints + 1);
    angle = left_angle;
  }
  else
    diff_angle = angle = 0.;   /* we do not need these variables in this case */

  /* the left boundary point */
  x = x_last = DISTR.BD_LEFT;
  fx = fx_last = (x <= -INFINITY) ? 0. : PDF(x);
  seg = GEN->seg = _unur_arou_segment_new( gen, x, fx );
  if (seg == NULL) return UNUR_ERR_GEN_CONDITION;  /* case of error */
  is_increasing = 1; /* assume PDF(x) is increasing for the first construction points */

  /* now all the other points */
  for( i=0; i<=PAR->n_starting_cpoints; i++ ) {

    /* starting point */
    if (i < PAR->n_starting_cpoints) {
      if (PAR->starting_cpoints) {   
	/* construction points provided by user */
	x = PAR->starting_cpoints[i];
	/* check starting point */
	if (x <= DISTR.BD_LEFT || x >= DISTR.BD_RIGHT) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting point out of domain");
	  continue;
	}
	if (x<=x_last) {
	  _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"starting points not increasing -> skip");
	  continue;
	}
      }
      else {
	/* compute construction points by means of "equidistance" rule */
	angle += diff_angle;
	x = tan( angle ) + GEN->center;
      }
    }
    else {
      /* the very last segment. it is rather a "virtual" segment to store 
	 the right vertex of the last segment, i.e., the right boundary point. */
      x = DISTR.BD_RIGHT;
    }

    /* insert center ? */
    if (use_center && x >= GEN->center) {
      use_center = FALSE;   /* we use the center only once (of course) */
      is_center = TRUE;     /* the next construction point is the center */
      if (x>GEN->center) {
	x = GEN->center;   /* use the center now ... */
	--i;              /* and push the orignal starting point back on stack */
	if (!PAR->starting_cpoints)
	  angle -= diff_angle; /* we have to compute the starting point in this case */
      }
      /* else: x == GEN->center --> nothing to do */
    }
    else
      is_center = FALSE;

    /* value of PDF at starting point */
    fx = (x >= INFINITY) ? 0. : PDF(x);

    /* check value of PDF at starting point */
    if (!is_increasing && fx > fx_last) {
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not unimodal");
      return UNUR_ERR_GEN_CONDITION;
    }

    if (fx <= 0. && fx_last <= 0.) {
      /* we do not need two such points */
      if (is_increasing) {
	/* PDF is still increasing, i.e., constant 0 til now */
	if (i<PAR->n_starting_cpoints) {
	  /* and it is not the right boundary.
	     otherwise the PDF is constant 0 on all construction points.
	     then we need both boundary points. */
	  /* we only have to change tangent line v/u = x,
	     everything else remains unchanged */
	  seg->dltp[0] = -1.;
	  seg->dltp[1] = x;
	  /* seg->dltp[0] = -1; seg->dltp[2] = 0.;  not changed */
	  x_last = x;
	  continue;   /* next construction point */
	}
      }
      else
	/* there should be no more points with PDF(x) > 0 */
	break;
    }

    /* need a new segment */
    seg_new = _unur_arou_segment_new( gen, x, fx );
    if (seg_new == NULL) {
      /* case of error */
      seg->next = NULL;  /* derminate list (for listing in debugging mode) */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* append to linked list */
    seg->next =seg_new;
    seg->rtp = seg_new->ltp;
    seg->drtp = seg_new->dltp;
    
    /* next step */
    seg = seg_new;

    /* PDF still increasing ? */
    if (is_increasing && fx < fx_last)
      is_increasing = 0;

    /* store last computed values */
    x_last = x;
    fx_last = fx;
  }

  /* we have left the loop with the right boundary of the support of PDF
     make shure that we will never use seg for sampling. */
  seg->Ain = seg->Aout = 0.;
  seg->Acum = INFINITY;
  seg->next = NULL;         /* terminate list */
  --(GEN->n_segs);           /* we do not count this segment */

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_arou_get_starting_cpoints() */

/*****************************************************************************/

int
_unur_arou_get_starting_segments( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute segments for starting points                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define MAX_IT   (1000)      /* maximal number of iterations to avoid 
				infinite loop in case of numerical errors    */

  struct unur_arou_segment *seg, *seg_new, *seg_tmp; 
  double x,fx;              /* construction point, value of PDF at x         */
  int n_it = 0;             /* counter for number of iterations              */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);

  /* compute paramters for all segments */
  for( seg=GEN->seg; seg->next != NULL; ) {

    /* compute parameters for semgent */
    switch (_unur_arou_segment_parameter(gen, seg)) {
    case UNUR_SUCCESS:     /* computation of parameters for segment successful */
      /* skip to next segment. */
      seg = seg->next;
      continue;
    case UNUR_ERR_SILENT:    /* construction points too close */
      /* we have to remove this last segment from list */
      /* (the last construction point in the list is a boundary point.
         thus we might change the domain of the distribution.
         however, we only cut off a piece that is beyond the precesion
         of the floating point arithmetic.)  */
      if (seg->next != NULL) {
	seg_tmp = seg->next;
	seg->next = seg->next->next;
	seg->rtp = seg->next->ltp;
	seg->drtp = seg->next->dltp;
	free(seg_tmp);
	--(GEN->n_segs);
      }
      else { /* seg->next==NULL */
        /* last (virtuel) interval in list.
           make shure that we will never use this segment */
	seg->Ain = seg->Aout = 0.;
	seg->Acum = INFINITY;
      }
      continue;
    case UNUR_ERR_INF:    /* segment unbounded */
      /* split segment */
      break;
    default:     /* PDF not T-concave */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* next iteration step */
    ++n_it;

    if (n_it > MAX_IT) {
      /* maximal number of iterations exceeded */
      /* assume PDF not T-concave              */
      return UNUR_ERR_GEN_CONDITION;
    }

    /* area in segment infinite. insert new construction point. */
    x = _unur_arou_segment_arcmean(seg);  /* use mean point in segment */

    /* value of PDF at x */
    fx = PDF(x);
    if (fx < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF < 0");
      return UNUR_ERR_GEN_DATA;
    }

    /* add a new segment, but check if we had to used too many segments */
    if (GEN->n_segs >= GEN->max_segs) {
      /* we do not want to create too many segments */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot create bounded envelope!");
      return UNUR_ERR_GEN_CONDITION;
    }
    seg_new = _unur_arou_segment_new( gen, x, fx );
    if (seg_new == NULL) return UNUR_ERR_GEN_CONDITION;  /* case of error */

    /* insert into linked list */
    seg_new->next = seg->next;
    seg->next = seg_new;

    /* right vertices */
    seg_new->rtp = seg->rtp;
    seg_new->drtp = seg->drtp;
    seg->rtp = seg_new->ltp;
    seg->drtp = seg_new->dltp;

  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_arou_get_starting_segments() */

/*****************************************************************************/

struct unur_arou_segment *
_unur_arou_segment_new( struct unur_gen *gen, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* get new segment and compute left construction point at x.            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... left point of new segment                                  */
     /*   fx  ... value of PDF at x                                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to new segment                                             */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg;
  double u,v,dfx;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,NULL);

  /* first check fx */
  if (fx<0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.");
    return NULL;
  }
  if (_unur_FP_is_infinity(fx)) {
    /* over flow */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) overflow");
    return NULL;
  }

  /* we need a new segment */
  seg = _unur_xmalloc( sizeof(struct unur_arou_segment) );
  seg->next = NULL; /* add eol marker */
  ++(GEN->n_segs);   /* increment counter for segments */
  COOKIE_SET(seg,CK_AROU_SEG);

  /* initialize some entries in seg */
  seg->Ain = seg->Aout = seg->Acum = 0.;
  seg->mid[0] = seg->mid[1] = 0.;

  /* make left construction point in segment */

  /* case: x out of support */
  if ( _unur_iszero(fx) ) {
    seg->ltp[0] = 0.;   /* vertex == origin */
    seg->ltp[1] = 0.;
    if (x <= -INFINITY || x >= INFINITY ) {
      /* tangent line == line u = 0 (i.e., v-axis) */
      seg->dltp[0] = 0.;   /* dv */
      seg->dltp[1] = 1.;   /* du */
      seg->dltp[2] = 0.;   /* v * dv + u * du */
    }
    else {
      /* tangent line == line v/u = x */
      seg->dltp[0] = -1.;  /* dv */
      seg->dltp[1] = x;    /* du */
      seg->dltp[2] = 0.;   /* v * dv + u * du */
    }
    return seg;
  }

  /* case: x in support */
  /* boundary point */
  u = sqrt( fx );
  v = x * u;
  seg->ltp[0] = v;
  seg->ltp[1] = u; 

  /* tangent line at tp */

  /* compute derivative of PDF at tp x */
  dfx = dPDF(x);

  /* subcase: derivative bounded 
     use derivative for tangent line */
  if ( dfx > -INFINITY && dfx < INFINITY ) {
    seg->dltp[0] = -dfx / u;             /* dv */    /** TODO: possible overflow **/
    seg->dltp[1] = 2 * u + dfx * x / u;  /* du */    /** TODO: possible overflow **/
    seg->dltp[2] = seg->dltp[0] * v + seg->dltp[1] * u;
    return seg;
  }

  /* subcase: derivative unbounded.
     use straight line through origin and vertex */
  seg->dltp[0] = -u;   /* dv */
  seg->dltp[1] = v;    /* du */
  seg->dltp[2] = 0.;

  return seg;

} /* end of _unur_arou_segment_new() */

/*****************************************************************************/

/* maximal distance of intersection point from origin 
   compared to distance of construction points to origin */
#define MAX_NORM_OF_INTERSECTION_POINT  1.e6

/*---------------------------------------------------------------------------*/

int
_unur_arou_segment_parameter( struct unur_gen *gen, struct unur_arou_segment *seg )
     /*----------------------------------------------------------------------*/
     /* compute all parameters for a segment.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   seg ... pointer to segment                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS     ... on success                                    */
     /*   UNUR_ERR_SILENT  ... do not add this construction point            */
     /*   UNUR_ERR_INF     ... if area = INFINITY                            */
     /*   other error code ... error (PDF not T-concave)                     */
     /*----------------------------------------------------------------------*/
{
  double coeff_det, cramer_det[2];
  double norm_vertex;      /* sum of 1-norms of vertices */
  double diff_tangents;    /* difference between coefficients of tangents */
  double det_bound;        /* bound for determinant for Cramer's rule */
  double tmp_a, tmp_b;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(seg,UNUR_ERR_NULL);  COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);

  /* sum of 1-norms of vertices */
  norm_vertex = fabs(seg->ltp[0]) + fabs(seg->ltp[1]) + fabs(seg->rtp[0]) + fabs(seg->rtp[1]);

  /* area inside the squeeze */
  seg->Ain = (seg->ltp[1] * seg->rtp[0] - seg->ltp[0] * seg->rtp[1]) / 2.;

  /* due to our ordering of construction points, seg->Ain must be >= 0 ! */
  if( seg->Ain < 0. ) {
    /* this could happen by round-off errors when the segment has angle extremely
       close to 0. Check this.  */
    if (fabs(seg->Ain) < 1.e-8 * norm_vertex) {
      seg->Ain = seg->Aout = 0.;
    }
    else {
      /* This should not happen:
	 non-ascending ordering of construction points */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    }
    return UNUR_ERR_SILENT;
  }

  /* we use Cramer's rule to compute intersection of tangent lines.
     (well, we could save one multiplication otherwise)             */
  coeff_det     = seg->dltp[0] * seg->drtp[1] - seg->dltp[1] * seg->drtp[0];
  cramer_det[0] = seg->dltp[2] * seg->drtp[1] - seg->dltp[1] * seg->drtp[2];
  cramer_det[1] = seg->dltp[0] * seg->drtp[2] - seg->dltp[2] * seg->drtp[0];

  /* we there are two possibilities for singular coefficent matrix:

     either the outer triangle is unlimited.
            then the two tangents are distinct but parallel and 
	    the corresponding linear equation has no solution, i.e.
	    coeff_det == 0 but cramer_det[0] != 0 or cramer_det[1] != 0.

     or     the outer triangle degenerates to a line segment and has area 0.
            then the two tangents are equal and
	    the corresponding linear equation has a nonunique solution, i.e.
	    coeff_det == cramer_det[0] == cramer_det[1] == 0.
  */

  /* we to not allow that the outer triangles becomes too large.
     so if the 1-norm of intersection point is too large compared
     to norm_vertex we assume that this triangle is unbounded.
     we thus avoid numerical errors.
     (we use the 1-norm here since it much easier to handle.)

     However, this might also happen due to roundoff errors,
     when the real position is extremely close to the secant.
     But at least we are on the save side.
     We only make an exception when coeff_det == 0, since otherwise there
     might be some problems with distributions like U(1,2).
  */
  det_bound = fabs(coeff_det) * norm_vertex * MAX_NORM_OF_INTERSECTION_POINT;

  /* difference between coefficents of tangents */
  diff_tangents = ( fabs(seg->dltp[0] - seg->drtp[0]) + fabs(seg->dltp[1] - seg->drtp[1])
		    + fabs(seg->dltp[2] - seg->drtp[2]) );

  /* case: intersection point exists and is unique */
  if (!_unur_iszero(coeff_det) && !_unur_iszero(diff_tangents)) {

    /* first check whether we can compute the intersection point */
    if ( fabs(cramer_det[0]) > det_bound || fabs(cramer_det[1]) > det_bound ) {
      /* case: triangle is assumed to be unbounded */	     
      /*      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"outer triangle assumed unbounded  1"); */
      seg->Aout = INFINITY;
      return UNUR_ERR_INF;
    }

    /* compute intersection point */
    seg->mid[0] = cramer_det[0] / coeff_det;
    seg->mid[1] = cramer_det[1] / coeff_det;

    /* area outside the squeeze */
    seg->Aout = ( (seg->ltp[0] - seg->mid[0]) * (seg->rtp[1] - seg->mid[1])
		  - (seg->ltp[1] - seg->mid[1]) * (seg->rtp[0] - seg->mid[0])) / 2.;

    /* due to our ordering of construction points, seg->Aout must be >= 0
       for a regular triangle.
       Thus if seg->Aout < 0, then the intersection point of tangents is on
       the WRONG side of the secant through vertices of segment,
       i.e. the "triangle" outside the squeeze region is unbounded.

       However this might also happen due to roundoff errors when 
       the real position is extremely close to the secant.
       We can distinguish between these two case by means of the u-coordinate
       of the intersection point. If it is on the wrong side of the secant,
       then we have seg->mid[1] < 0 (Otherwise we have either a round-off error
       or the PDF is not T-concave.) */
    if( seg->mid[1] < 0. ) {
      /* _unur_warning(gen->genid,UNUR_ERR_GENERIC,"outer triangle unbounded  2"); */
      seg->Aout = INFINITY;
      return UNUR_ERR_INF;
    }

    /* at last check result.
       we must have:
         (*) seg->Aout > 0
         (*) intersection point right of left construction point
         (*) intersection point left of right construction point */
    if ( seg->Aout > 0. ) {
      tmp_a = seg->mid[0] * seg->ltp[1];
      tmp_b = seg->ltp[0] * seg->mid[1];
      if ( ! _unur_FP_less(tmp_a, tmp_b) ) {
	tmp_a = seg->mid[0] * seg->rtp[1];
	tmp_b = seg->rtp[0] * seg->mid[1];
	if ( ! _unur_FP_greater(tmp_a, tmp_b) )
	  /* everything o.k. */
	  return UNUR_SUCCESS;
      }
    }

    /* there are two cases why the above check failed:
       (1) the PDF is not T-concave
       (2) small roundoff errors.
    */

    if ( !_unur_iszero(seg->ltp[1]) && !_unur_iszero(seg->rtp[1]) ) {
      tmp_a = seg->ltp[0] * seg->rtp[1];
      tmp_b = seg->rtp[0] * seg->ltp[1];
      if ( _unur_FP_equal(tmp_a, tmp_b) ) {
	/* we assume that left construction point = right construction point */
	/* (construction points are too close) */
	seg->Ain = seg->Aout = 0.;
	return UNUR_ERR_SILENT;
      }
    }

    /* _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF not T-concave"); */

    /* The outer area is unbounded (this might happen if the given PDF
       is not T-concave or due to round off errors.
       We assume round off errors. If the PDF is not T-concave we
       will exceed any bound for the number of segments.    */

    /* However due to round off errors, Aout might have become < 0
       when the boundary of region is (almost) a straight line
       and the area outside squeeze should be (almost) zero. */

    if (!(fabs(seg->Aout) < fabs(seg->Ain) * UNUR_EPSILON)) {
      seg->Aout = INFINITY;
      return UNUR_ERR_INF;
    }
    /* else we assume round-off errors and set Aout = 0.
       i.e. go to the remaining case below */
  }

  /* remaining case: triangle degenerates to a line segment, i.e.
                     intersection point exists but is not unique */

  /* boundary of region is (almost) a straight line
     and area outside squeeze is (almost) zero.
     use middle point as intersection point and 
     set area outside the squeeze to 0.
  */
  /*    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"outer triangle is line"); */
  seg->mid[0] =  0.5 * (seg->ltp[0] + seg->rtp[0]);
  seg->mid[1] =  0.5 * (seg->ltp[1] + seg->rtp[1]);
  seg->Aout = 0.;
  
  /* now it should be o.k. */
  return UNUR_SUCCESS;
  
} /* end if _unur_arou_segment_parameter() */

/*---------------------------------------------------------------------------*/

#undef MAX_NORM_INTERSECTION

/*****************************************************************************/

int
_unur_arou_segment_split( struct unur_gen *gen, struct unur_arou_segment *seg_oldl, double x, double fx )
     /*----------------------------------------------------------------------*/
     /* insert new segment                                                   */
     /*   old segment -> left hand side                                      */
     /*   new segment -> right hand side                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   seg_oldl ... pointer to segment                                    */
     /*   x        ... left point of new segment                             */
     /*   fx       ... value of PDF at x                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS      ... on success                                   */
     /*   UNUR_ERR_SILENT   ... if no intervals are splitted                 */
     /*   other error codes ... on error                                     */
     /*   1  ... if successful                                               */
     /*   0  ... if */
     /*  -1  ... error                                                       */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg_newr;    /* pointer to newly created segment */
  struct unur_arou_segment seg_bak;      /* space for saving data of segment */
  double Adiff;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);      COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(seg_oldl,UNUR_ERR_NULL); COOKIE_CHECK(seg_oldl,CK_AROU_SEG,UNUR_ERR_COOKIE);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug & AROU_DEBUG_SPLIT) 
      _unur_arou_debug_split_start( gen,seg_oldl,x,fx );
#endif

  /* we only add a new construction point, if the relative area is large enough */
  if (GEN->n_segs * seg_oldl->Aout / (GEN->Atotal - GEN->Asqueeze) < GEN->darsfactor )
    return UNUR_SUCCESS;

  /* check for data error */
  if (fx < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PDF(x) < 0.!");
    return UNUR_ERR_GEN_DATA;
  }

  /* backup data */
  memcpy(&seg_bak, seg_oldl, sizeof(struct unur_arou_segment));

  /* PDF at x is 0. */
  if (fx <= 0.) {
    if (seg_oldl->rtp[1] <= 0. && seg_oldl->rtp[0] <= 0. ) {
      /* just chop off the right part of segment */
      /* we only have to change tangent line v/u = x
	 at the right hand vertex */
      seg_oldl->drtp[1] = x;    /* du */
    }
    else if (seg_oldl->ltp[1] <= 0. && seg_oldl->ltp[0] <= 0. ) {
      /* just chop off the left part of segment */
      /* we only have to change tangent line v/u = x
	 at the left hand vertex */
      seg_oldl->dltp[1] = x;    /* du */
    }
    else {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return UNUR_ERR_SHOULD_NOT_HAPPEN;
    }
    
    /* parameters of new segment */
    if( _unur_arou_segment_parameter(gen,seg_oldl)!=UNUR_SUCCESS ) {
      /* PDF not T-concave or area in segment not bounded */
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot chop segment at given point");
      /* error --> restore the old segment */
      memcpy(seg_oldl, &seg_bak, sizeof(struct unur_arou_segment));
      return UNUR_ERR_SILENT;
    }

    /* for _unur_arou_debug_split_stop only */
    seg_newr = seg_oldl;
  }

  else {  /* fx > 0 */

    /* need new segment */
    seg_newr = _unur_arou_segment_new(gen,x,fx);
    if (seg_newr == NULL) return UNUR_ERR_GEN_DATA;  /* case of error */
    
    /* link into list */
    seg_newr->next = seg_oldl->next;
    seg_oldl->next = seg_newr;
    
    /* right vertices */
    seg_newr->rtp = seg_oldl->rtp;
    seg_newr->drtp = seg_oldl->drtp;
    seg_oldl->rtp = seg_newr->ltp;
    seg_oldl->drtp = seg_newr->dltp;
    
    /* parameters of new segments */
    if( _unur_arou_segment_parameter(gen,seg_oldl)!=UNUR_SUCCESS || 
	/* PDF not T-concave or */
	_unur_arou_segment_parameter(gen,seg_newr)!=UNUR_SUCCESS 
	/* area in segment not bounded */ ) {
    
      /* new construction point not suitable --> do not add */
      _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"Cannot split segment at given point.");
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (gen->debug & AROU_DEBUG_SPLIT)
	_unur_arou_debug_split_stop( gen,seg_oldl,seg_newr );
#endif
      
      /* restore the old segment */
      memcpy(seg_oldl, &seg_bak, sizeof(struct unur_arou_segment));

      /* decrement counter for segments and free unused segment */
      if (seg_newr) {
	--(GEN->n_segs); 
	free( seg_newr );
      }

      return UNUR_ERR_SILENT;
    }
  }

  /* successful */

  /* update total area below hat and squeeze */
  Adiff =  - seg_bak.Ain  + seg_oldl->Ain  + ((seg_newr!=seg_oldl) ? seg_newr->Ain : 0. );
  GEN->Asqueeze += Adiff;
  Adiff += - seg_bak.Aout + seg_oldl->Aout + ((seg_newr!=seg_oldl) ? seg_newr->Aout : 0. );
  GEN->Atotal += Adiff;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & AROU_DEBUG_SPLIT) 
    _unur_arou_debug_split_stop( gen,seg_oldl,seg_newr );
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_arou_segment_split() */


/*****************************************************************************/

double
_unur_arou_compute_x( double v, double u )
     /*----------------------------------------------------------------------*/
     /* compute point x from (v,u) tuple.                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   v      ... numerator                                               */
     /*   u      ... denominator                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   point x of (x,y) tuple                                             */
     /*                                                                      */
     /* remark:                                                              */
     /*   if (v,u)=(0,0) then x is set to INFINITY                           */
     /*----------------------------------------------------------------------*/
{
  if (!_unur_iszero(u)) return v/u;
  else if (v<0.)        return -INFINITY;
  else                  return INFINITY;
} /* end of _unur_arou_compute_x() */

/*---------------------------------------------------------------------------*/

int
_unur_arou_run_dars( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* run derandomized adaptive rejection sampling.                         */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg, *seg_next;
  double Atot, Asqueezetot;    /* total area below hat and squeeze, resp. */
  double Alimit;               /* threshhold value for splitting interval */
  int n_splitted = 1;          /* count splitted intervals */
  int splitted;                /* result of splitting routine */
  double xl, xr;               /* boundary of interval */
  double xsp, fxsp;            /* splitting point in interval */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);     COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);

  /* there is no need to run DARS when the DARS factor is INFINITY */
  if (_unur_FP_is_infinity(GEN->darsfactor))
    return UNUR_SUCCESS;

  /* first we need the total areas below hat and squeeze.
     (This is only necessary, when _unur_arou_make_guide_table() has not been
     called!)                                                                */
  Atot = 0.;            /* area below hat */
  Asqueezetot = 0.;     /* area below squeeze */
  for (seg = GEN->seg; seg != NULL; seg = seg->next ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
    Asqueezetot += seg->Ain;
    Atot += seg->Ain + seg->Aout;
  }
  GEN->Atotal = Atot;
  GEN->Asqueeze = Asqueezetot;

  /* now split intervals */
  while ( (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) &&
	  (GEN->n_segs < GEN->max_segs) ) {

    /* compute threshhold value. every interval with area between
       hat and squeeze greater than this value will be splitted.  */
    if (GEN->n_segs > 1)
      Alimit = GEN->darsfactor * ( (GEN->Atotal - GEN->Asqueeze) / GEN->n_segs );
    else
      /* we split every interval if there are only one interval */
      Alimit = 0.; 

    /* reset counter for splitted intervals */
    n_splitted = 0;

    /* for all intervals do ... */
    for (seg = GEN->seg; seg->next != NULL; seg = seg->next ) {
      COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);

      /* do not exceed the maximum number of intervals */
      if (GEN->n_segs >= GEN->max_segs)
	break;

      /* we skip over all intervals where the area between hat and
	 squeeze does not exceed the threshhold value.             */
      if (seg->Aout <= Alimit) 
	continue;  /* goto next interval */

      /* store pointer to next interval */
      seg_next = seg->next;

      /* boundary of interval */
      xl = _unur_arou_compute_x(seg->ltp[0],seg->ltp[1]);
      xr = _unur_arou_compute_x(seg->rtp[0],seg->rtp[1]);
      if (xl>xr) xl = -INFINITY;   /* (0,0) is mapped to INF instead of -INF */

      /* However ... */
      if ( _unur_FP_is_minus_infinity(xl)
	   && _unur_FP_same(seg->dltp[0],-1.) && _unur_iszero(seg->dltp[2]) )
  	/* boundary of domain given by tangent */
	xl = seg->dltp[1];

      if ( _unur_FP_is_infinity(xr)
	   && _unur_FP_same(seg->drtp[0],-1.) && _unur_iszero(seg->drtp[2]) )
  	/* boundary of domain given by tangent */
	xr = seg->drtp[1];

      /* get splitting point (arc-mean rule) */
      xsp = _unur_arcmean(xl,xr);

      /* value of PDF at splitting point */
      fxsp = PDF(xsp);

      /* now split interval at given point */
      splitted = _unur_arou_segment_split(gen, seg, xsp, fxsp);

      if (splitted == UNUR_SUCCESS) {
      	/* splitting successful */
      	++n_splitted;

	/* If we have chopped off the left or right tail,
	   we can simply continue with seg->next.
	   Otherwise, if we have split the segment, then
	   seg->next points to the newly created segment
	   and there is no need to split this in this loop;
	   thus we skip seg->next and continue we seg->next->next.
	   we can distinguish between these two cases by comparing
	   the stored pointer seg_next with seg->next. */
	if (seg->next != seg_next)
	  seg = seg->next;
      }
      else if (splitted != UNUR_ERR_SILENT) {
	/* some serious error occurred */
	/* condition for PDF is violated! */
	_unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
	return UNUR_ERR_GEN_CONDITION;
      }
      /* else: could not split construction points: too close (?) */
    }

    if (n_splitted == 0) {
      /* we are not successful in splitting any inteval.
	 abort to avoid endless loop */
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: no intervals could be splitted.");
      break;
    }
  }

  /* ratio between squeeze and hat o.k. ? */
  if ( GEN->max_ratio * GEN->Atotal > GEN->Asqueeze ) {
    if ( GEN->n_segs >= GEN->max_segs )
      _unur_warning(gen->genid,UNUR_ERR_GENERIC,"DARS aborted: maximum number of intervals exceeded.");
    _unur_warning(gen->genid,UNUR_ERR_GENERIC,"hat/squeeze ratio too small.");
  }
  else {
    /* no more construction points */
    GEN->max_segs = GEN->n_segs;
  }

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_arou_run_dars() */

/*****************************************************************************/

int
_unur_arou_make_guide_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make a guide table for indexed search                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  struct unur_arou_segment *seg;
  double Acum, Aincum, Astep;
  int max_guide_size;
  int j;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_AROU_GEN,UNUR_ERR_COOKIE);

  /* allocate blocks for guide table (if necessary).
     (we allocate blocks for maximal guide table.) */
  if (!GEN->guide) {
    max_guide_size = (GEN->guide_factor > 0.) ? ((int)(GEN->max_segs * GEN->guide_factor)) : 1;
    if (max_guide_size <= 0) max_guide_size = 1;   /* protect against overflow */
    GEN->guide = _unur_xmalloc( max_guide_size * sizeof(struct unur_arou_segment*) );
  }

  /* first we need cumulated areas in segments */
  Acum = 0.;       /* area in enveloping polygon */
  Aincum = 0.;     /* area in squeeze */
  for (seg = GEN->seg; seg != NULL; seg = seg->next ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
    Acum += seg->Ain + seg->Aout;
    Aincum += seg->Ain;
    seg->Acum = Acum;
  }

  /* total area below hat */
  GEN->Atotal = Acum;
  GEN->Asqueeze = Aincum;

  /* actual size of guide table */
  GEN->guide_size = (int)(GEN->n_segs * GEN->guide_factor);
  /* we do not vary the relative size of the guide table,
     since it has very little influence on speed */

  /* make table (use variant 2; see dgt.c) */
  Astep = GEN->Atotal / GEN->guide_size;
  Acum=0.;
  for( j=0, seg=GEN->seg; j < GEN->guide_size; j++ ) {
    COOKIE_CHECK(seg,CK_AROU_SEG,UNUR_ERR_COOKIE);
    while( seg->Acum < Acum )
      if( seg->next != NULL )    /* skip to next segment if it exists */
        seg = seg->next;
      else {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
    GEN->guide[j] = seg;
    Acum += Astep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide[j] = seg;

  return UNUR_SUCCESS;
} /* end of _unur_arou_make_guide_table() */

/*****************************************************************************/

double
_unur_arou_segment_arcmean( struct unur_arou_segment *seg )
     /*----------------------------------------------------------------------*/
     /* compute "arctan mean" of two numbers expressed as v/u, u>=0          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   seg ... pointer to segment                                         */
     /*                                                                      */
     /* return:                                                              */
     /*   mean                                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*                                                                      */
     /* comment:                                                             */
     /*   "arctan mean" = tan(0.5*(arctan(x0)+arctan(x1)))                   */
     /*                                                                      */
     /*   a combination of arithmetical mean (for x0 and x1 close to 0)      */
     /*   and the harmonic mean (for |x0| and |x1| large).                   */
     /*----------------------------------------------------------------------*/
{
  double xl, xr;

  /* check arguments */
  CHECK_NULL(seg,INFINITY);  COOKIE_CHECK(seg,CK_AROU_SEG,INFINITY);

  /* if u != 0 ... x is stored in tp (= v/u)   */
  /* else      ... x is stored in tangent dltp */
  xl = (seg->ltp[1] > 0.) ? (seg->ltp[0] / seg->ltp[1]) :
    ( _unur_iszero(seg->dltp[0]) ? -INFINITY : (seg->dltp[1]) );

  xr = (seg->rtp[1] > 0.) ? (seg->rtp[0] / seg->rtp[1]) :
    ( _unur_iszero(seg->drtp[0]) ? INFINITY : (seg->drtp[1]) );

  return _unur_arcmean(xl,xr);

} /* end of _unur_arou_segment_arcmean() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_arou_debug_init( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator after setup into LOG file                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_AROU_PAR,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = ratio-of-uniforms method with enveloping polygon\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_arou_sample",gen->genid);
  if (gen->variant & AROU_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: center = %g",gen->genid,GEN->center);
  _unur_print_if_default(gen,AROU_SET_CENTER);
  if (gen->variant & AROU_VARFLAG_USECENTER)
    fprintf(LOG,"\n%s: use center as construction point",gen->genid);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fprintf(LOG,"%s: maximum number of segments         = %d",gen->genid,GEN->max_segs);
  _unur_print_if_default(gen,AROU_SET_MAX_SEGS);
  fprintf(LOG,"\n%s: bound for ratio  Asqueeze / Atotal = %g%%",gen->genid,GEN->max_ratio*100.);
  _unur_print_if_default(gen,AROU_SET_MAX_SQHRATIO);
  fprintf(LOG,"\n%s:\n",gen->genid);

  if (gen->variant & AROU_VARFLAG_USEDARS) {
    fprintf(LOG,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(gen,AROU_SET_USE_DARS);
    fprintf(LOG,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(gen,AROU_SET_DARS_FACTOR);
  }
  else {
    fprintf(LOG,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(gen,AROU_SET_USE_DARS);
  }
  fprintf(LOG,"\n%s:\n",gen->genid);

  fprintf(LOG,"%s: sampling from list of segments: indexed search (guide table method)\n",gen->genid);
  fprintf(LOG,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN->guide_factor);
  _unur_print_if_default(gen,AROU_SET_GUIDEFACTOR);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fprintf(LOG,"%s: number of starting points = %d",gen->genid,PAR->n_starting_cpoints);
  _unur_print_if_default(gen,AROU_SET_N_STP);
  fprintf(LOG,"\n%s: starting points:",gen->genid);
  if (gen->set & AROU_SET_STP)
    for (i=0; i<PAR->n_starting_cpoints; i++) {
      if (i%5==0) fprintf(LOG,"\n%s:\t",gen->genid);
      fprintf(LOG,"   %#g,",PAR->starting_cpoints[i]);
    }
  else
    fprintf(LOG," use \"equidistribution\" rule [default]");
  fprintf(LOG,"\n%s:\n",gen->genid);
  
  _unur_arou_debug_segments(gen);

  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_arou_debug_init() */

/*****************************************************************************/

void 
_unur_arou_debug_dars_start( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print header before runniung DARS into LOG file                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: DARS started **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(gen,AROU_SET_DARS_FACTOR);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fflush(LOG);
} /* end of _unur_arou_debug_dars_start() */

/*---------------------------------------------------------------------------*/

void
_unur_arou_debug_dars( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print infor after generator has run DARS into LOG file               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS finished **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_arou_debug_segments(gen);
  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: DARS completed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);
} /* end of _unur_arou_debug_dars() */

/*****************************************************************************/

void
_unur_arou_debug_free( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator before destroying into LOG file           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: GENERATOR destroyed **********************\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_arou_debug_segments(gen);
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_arou_debug_free() */

/*****************************************************************************/

void
_unur_arou_debug_segments( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write list of segments into LOGfile                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  struct unur_arou_segment *seg;
  double sAin, sAout, Atotal;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:Segments: %d\n",gen->genid,GEN->n_segs);
  if ((gen->debug & AROU_DEBUG_SEGMENTS) && GEN->seg != NULL) {
    fprintf(LOG,"%s: Nr.\t    left touching point\t\t   intersection point\t\t tangent at left touching point\n",gen->genid);
    for (seg = GEN->seg, i=0; seg->next!=NULL; seg=seg->next, i++) {
      COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID); 
      fprintf(LOG,"%s:[%3d]: (%-12.6g,%-12.6g)   (%-12.6g,%-12.6g)   (%-12.6g,%-12.6g,%-12.6g)\n", gen->genid, i,
	      seg->ltp[0],seg->ltp[1],
	      seg->mid[0],seg->mid[1],
	      seg->dltp[0],seg->dltp[1],seg->dltp[2]);
    }
    COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID); 
    fprintf(LOG,"%s:[...]: (%-12.6g,%-12.6g)\n", gen->genid,seg->ltp[0],seg->ltp[1]);
  }
  fprintf(LOG,"%s:\n",gen->genid);

  if (GEN->Atotal <= 0.) {
    fprintf(LOG,"%s: Construction of enveloping polygon not successful\n",gen->genid);
    fprintf(LOG,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!!!!!!!\n",gen->genid);
    fprintf(LOG,"%s:\n",gen->genid);
    Atotal = -1.;   /* to avoid floating point exceptions */
  }
  else {
    Atotal = GEN->Atotal;
  }

  /* print and sum areas inside and outside of squeeze */
  if ((gen->debug & AROU_DEBUG_SEGMENTS) && GEN->seg != NULL) {
    fprintf(LOG,"%s:Areas in segments:\n",gen->genid);
    fprintf(LOG,"%s: Nr.\t inside squeeze\t\t   outside squeeze\t     total segment\t\tcumulated\n",gen->genid);
    sAin = sAout = 0.;
    for (seg = GEN->seg, i=0; seg->next!=NULL; seg=seg->next, i++) {
      COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID); 
      sAin += seg->Ain;
      sAout += seg->Aout;
      fprintf(LOG,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      seg->Ain, seg->Ain * 100. / Atotal,
	      seg->Aout, seg->Aout * 100. / Atotal,
	      seg->Ain + seg->Aout, (seg->Ain + seg->Aout) * 100. / Atotal,
	      seg->Acum, seg->Acum * 100. / Atotal);
    }
    fprintf(LOG,"%s:\t----------  ---------  |  ----------  ---------  |  ----------  ---------  +\n",gen->genid);
    fprintf(LOG,"%s: Sum : %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-11.6g(%6.3f%%)\n",
	    gen->genid,
	    sAin, sAin * 100./Atotal,
	    sAout, sAout * 100./Atotal,
	    sAin + sAout, (sAin + sAout) * 100./Atotal);
    fprintf(LOG,"%s:\n",gen->genid);
  }

  /* summary of areas */
  fprintf(LOG,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(LOG,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(LOG,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_arou_debug_segments() */

/*****************************************************************************/

void
_unur_arou_debug_split_start( const struct unur_gen *gen,
			      const struct unur_arou_segment *seg, 
			      double x, double fx )
     /*----------------------------------------------------------------------*/
     /* write info about splitting segment                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   seg ... pointer to segment                                         */
     /*   x   ... split at this point                                        */
     /*   fx  ... value of PDF at x                                          */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  CHECK_NULL(seg,RETURN_VOID);  COOKIE_CHECK(seg,CK_AROU_SEG,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: split segment at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(LOG,"%s: old segment:\n",gen->genid);

  fprintf(LOG,"%s:   left  construction point = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	  gen->genid, seg->ltp[0], seg->ltp[1], seg->ltp[0]/seg->ltp[1], sqrt(seg->ltp[1]) );

  fprintf(LOG,"%s:   intersection point       = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\n",
	  gen->genid, seg->mid[0], seg->mid[1], seg->mid[0]/seg->mid[1]);

  fprintf(LOG,"%s:   right construction point = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	  gen->genid, seg->rtp[0], seg->rtp[1], seg->rtp[0]/seg->rtp[1], sqrt(seg->rtp[1]) );

  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg->Ain, seg->Ain * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  seg->Aout, seg->Aout * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  (seg->Ain + seg->Aout), (seg->Ain +seg->Aout) * 100./GEN->Atotal);

  fflush(LOG);

} /* end of _unur_arou_debug_split_start() */

/*****************************************************************************/

void
_unur_arou_debug_split_stop( const struct unur_gen *gen, 
			     const struct unur_arou_segment *seg_left,
			     const struct unur_arou_segment *seg_right )
     /*----------------------------------------------------------------------*/
     /* write info about new splitted segments                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   iv_left  ... pointer to new left hand segment                      */
     /*   iv_right ... pointer to new right hand segment                     */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int chopped;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);        COOKIE_CHECK(gen,CK_AROU_GEN,RETURN_VOID);
  CHECK_NULL(seg_left,RETURN_VOID);   COOKIE_CHECK(seg_left,CK_AROU_SEG,RETURN_VOID);
  CHECK_NULL(seg_right,RETURN_VOID);  COOKIE_CHECK(seg_right,CK_AROU_SEG,RETURN_VOID);

  LOG = unur_get_stream();

  /* whether segment was split or just chopped */
  chopped = (seg_left==seg_right) ? 1 : 0;

  if (chopped)
    fprintf(LOG,"%s: new segment (chopped):\n",gen->genid);
  else
    fprintf(LOG,"%s: new segments:\n",gen->genid);

  fprintf(LOG,"%s:   left  construction point  = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	  gen->genid, seg_left->ltp[0], seg_left->ltp[1], seg_left->ltp[0]/seg_left->ltp[1], sqrt(seg_left->ltp[1]) );
  
  fprintf(LOG,"%s:   intersection point        = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\n",
	  gen->genid, seg_left->mid[0], seg_left->mid[1], seg_left->mid[0]/seg_left->mid[1] );
    
  if (chopped) {
    fprintf(LOG,"%s:   right construction point  = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	    gen->genid, seg_left->rtp[0], seg_left->rtp[1], seg_left->rtp[0]/seg_left->rtp[1], sqrt(seg_left->rtp[1]) );
  }
  else {
    fprintf(LOG,"%s:   middle construction point = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	    gen->genid, seg_left->rtp[0], seg_left->rtp[1], seg_left->rtp[0]/seg_left->rtp[1], sqrt(seg_left->rtp[1]) );
    
    fprintf(LOG,"%s:   intersection point        = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\n",
	    gen->genid, seg_right->mid[0], seg_right->mid[1], seg_right->mid[0]/seg_right->mid[1] );
    
    fprintf(LOG,"%s:   right construction point  = (%-12.6g,%-12.6g)\t x = v/u = %-12.6g\tf(x) = %-12.6g\n",
	    gen->genid, seg_right->rtp[0], seg_right->rtp[1], seg_right->rtp[0]/seg_right->rtp[1], sqrt(seg_right->rtp[1]) );
  }
  
  if (!chopped) {
    fprintf(LOG,"%s: left segment:\n",gen->genid);
    fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_left->Ain, seg_left->Ain * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_left->Aout, seg_left->Aout * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    (seg_left->Ain + seg_left->Aout), (seg_left->Ain +seg_left->Aout) * 100./GEN->Atotal);
    
    fprintf(LOG,"%s: right segment:\n",gen->genid);
    fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_right->Ain, seg_right->Ain * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    seg_right->Aout, seg_right->Aout * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g\t(%6.3f%%)\n",gen->genid,
	    (seg_right->Ain + seg_right->Aout), (seg_right->Ain +seg_right->Aout) * 100./GEN->Atotal);
  }

  fprintf(LOG,"%s: total areas:\n",gen->genid);
  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t(%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (GEN->Atotal - GEN->Asqueeze) * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);

  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_arou_debug_split_stop() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_arou_info( struct unur_gen *gen, int help )
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

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF dPDF\n");
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
    if ( !(distr->set & (UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_MODE )) ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
			  "You may provide a point near the mode as \"center\"."); 
  }
  _unur_string_append(info,"\n");
      
  /* method */
  _unur_string_append(info,"method: AROU (Automatic Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g\n", GEN->Atotal);

  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", GEN->Atotal/(0.5*DISTR.area));
  else
    _unur_string_append(info,"<= %g\n", GEN->Atotal/GEN->Asqueeze);

  _unur_string_append(info,"   area ratio squeeze/hat = %g\n",
 		      GEN->Asqueeze/GEN->Atotal);

  _unur_string_append(info,"   # segments = %d\n", GEN->n_segs);
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");

    _unur_string_append(info,"   max_sqhratio = %g  %s\n", GEN->max_ratio,
			(gen->set & AROU_SET_MAX_SQHRATIO) ? "" : "[default]");
    _unur_string_append(info,"   max_segments = %d  %s\n", GEN->max_segs_info,
			(gen->set & AROU_SET_MAX_SEGS) ? "" : "[default]");

    if (gen->variant & AROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    if (gen->variant & AROU_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");

    /* Not displayed:
       int unur_arou_set_usedars( UNUR_PAR *parameters, int usedars );
       int unur_arou_set_darsfactor( UNUR_PAR *parameters, double factor );
       int unur_arou_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
       int unur_arou_set_usecenter( UNUR_PAR *parameters, int usecenter );
       int unur_arou_set_guidefactor( UNUR_PAR *parameters, double factor );
    */
    _unur_string_append(info,"\n");
  }

  /* Hints */
  if (help) {
    if ( !(gen->set & AROU_SET_MAX_SQHRATIO) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"max_sqhratio\" closer to 1 to decrease rejection constant." );
    if (GEN->Asqueeze/GEN->Atotal < GEN->max_ratio) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You should increase \"max_segments\" to obtain the desired rejection constant." );
    _unur_string_append(info,"\n");
  }

} /* end of _unur_arou_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
