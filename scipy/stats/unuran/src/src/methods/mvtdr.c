/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    multivariate transformed density rejection                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given (logarithm of the) PDF of a log-concave distribution;          *
 *      produce a value x consistent with its density.                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold             *
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
 *   [1] Leydold, J. (1998): A Rejection Technique for Sampling from         *
 *       Log-Concave Multivariate Distributions,                             *
 *       ACM TOMACS 8(3), pp. 254-280.                                       *
 *                                                                           *
 *   [2] Hoermann, W., J. Leydold, and G. Derflinger (2004):                 *
 *       Automatic Nonuniform Random Variate Generation, Springer, Berlin.   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <distr/cvec.h>
#include <distributions/unur_distributions.h>
#include <utils/fmax_source.h>
#include <utils/matrix_source.h>
#include <specfunct/unur_specfunct_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "tdr.h"
#include "mvtdr.h"
#include "mvtdr_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* relative size of guide table compared to number of cones */
#define GUIDE_TABLE_SIZE    1

/* find proper touching point using Brent's algorithm: */
/* acceptable tolerance for Brent's algorithm          */
#define FIND_TP_TOL         0.001   

/* fine tuning of generator: */
/* a number is considered to be zero if abs is below                         */
/*   TOLERANCE * PDF(center) / dim                                           */                                
#define TOLERANCE           (1.e-8)

/* gamma variate generator (using method TDR) */
#define MVTDR_TDR_SQH_RATIO (0.95)

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define MVTDR_VARFLAG_VERIFY     0x01u   /* flag for verifying hat           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define MVTDR_DEBUG_VERTEX      0x00000010u   /* print list of vertices      */
#define MVTDR_DEBUG_CONE        0x00000020u   /* print list of conesces      */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define MVTDR_SET_STEPSMIN        0x001u   /* min number of triangulation steps */
#define MVTDR_SET_MAXCONES        0x002u   /* max number of cones            */
#define MVTDR_SET_BOUNDSPLITTING  0x004u   /* bound for splitting cones      */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MVTDR"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_mvtdr_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_mvtdr_simplex_sample( const struct unur_gen *gen, double *U );
/*---------------------------------------------------------------------------*/
/* sample point uniformly on standard simplex.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvtdr_gammagen( struct unur_gen *gen, double alpha );
/*---------------------------------------------------------------------------*/
/* create a gamma random variate generator with shape parameter alpha.       */
/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/* Hat.                                                                      */
/*****************************************************************************/

static int _unur_mvtdr_create_hat( struct unur_gen *gen );
/* compute cones and hat function */


/*****************************************************************************/
/* CONES.                                                                    */
/*****************************************************************************/

static int _unur_mvtdr_initial_cones( struct unur_gen *gen );
/* get initial cones */

static CONE *_unur_mvtdr_cone_new( struct unur_gen *gen );
/* get new (empty) vertex object */

static int _unur_mvtdr_cone_center( struct unur_gen *gen, CONE *c );
/* computer center of cone */

static int _unur_mvtdr_cone_params( struct unur_gen *gen, CONE *c );
/* compute parameters for hat for a cone   (expect touching point and volume below hat) */

static double _unur_mvtdr_cone_logH( struct unur_gen *gen, CONE *c );
/* calculate log of volume below hat for given touching point */

static int _unur_mvtdr_cone_split( struct unur_gen *gen, CONE *c, int step );
/* split a cone */

static int _unur_mvtdr_triangulate( struct unur_gen *gen, int step, int all);
/* make one triangulation step */

static int _unur_mvtdr_cone_height( struct unur_gen *gen, CONE *c );
/* calculate height of pyramid (cone) */

static int _unur_mvtdr_max_gamma( struct unur_gen *gen );
/* compute upper bound for gamma variates */

/*****************************************************************************/
/* optimal distance for touching points                                      */
/*****************************************************************************/

static double _unur_mvtdr_tp_min_aux(double t, void *p);
/* auxiliary function to be used with _unur_util_brent(). */

static double _unur_mvtdr_tp_min( double t, void *p );
/* wrapper for _unur_mvtdr_cone_hatvolume();
   sets cone->tp;
   funtion that must be minimized for optimal touching point */

static int _unur_mvtdr_tp_find( struct unur_gen *gen, CONE *c );
/* find optimal touching point for cone */

static int _unur_mvtdr_tp_search( struct unur_gen *gen, TP_ARG *a );
/* search for proper touching point */

static int _unur_mvtdr_tp_bracket( struct unur_gen *gen, TP_ARG *a );
/* search for proper bracket of minimum of tp_f2min() */

/*****************************************************************************/
/* VERTICES.                                                                 */
/*****************************************************************************/

static int _unur_mvtdr_initial_vertices( struct unur_gen *gen );
/* get vertices of initial cones */

static VERTEX *_unur_mvtdr_vertex_new( struct unur_gen *gen );
/* get new (empty) vertex object */

static int _unur_mvtdr_number_vertices( struct unur_gen *gen, int level );
/* number of vertices in given triangulation level */

static VERTEX *_unur_mvtdr_vertex_on_edge( struct unur_gen *gen, VERTEX **vl );
/* compute new vertex on edge */


/*****************************************************************************/
/* hash table for storing EDGES.                                             */
/*****************************************************************************/

static int _unur_mvtdr_etable_new( struct unur_gen *gen, int size );
/* make new hash table */

static void _unur_mvtdr_etable_free( struct unur_gen *gen );
/* free hash table */

static VERTEX *_unur_mvtdr_etable_find_or_insert( struct unur_gen *gen, VERTEX **vidx );
/* find or insert entry in hash table, return pointer to vertex */

/*****************************************************************************/

static int _unur_mvtdr_make_guide_table( struct unur_gen *gen );
/* create guide table */


/*****************************************************************************/
/* Debug.                                                                    */
/*****************************************************************************/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_debug_init_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after (almost empty generator) object has been created.             */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_debug_init_finished( const struct unur_gen *gen, int successful );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized.                               */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_debug_vertices( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print list of vertices.                                                   */
/*---------------------------------------------------------------------------*/

static void _unur_mvtdr_debug_cones( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print list of cones.                                                      */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_mvtdr_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_mvtdr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mvtdr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cvec_PDF((x),(gen->distr))         /* call to PDF    */
#define logPDF(x) _unur_cvec_logPDF((x),(gen->distr))      /* call to logPDF */
#define dPDF(r,x) _unur_cvec_dPDF((r),(x),(gen->distr))    /* call to dPDF   */
#define dlogPDF(r,x) _unur_cvec_dlogPDF((r),(x),(gen->distr))   /* call to dlogPDF */

/* an auxiliary generator for gamma variates */
#define GEN_GAMMA  gen->gen_aux

/*---------------------------------------------------------------------------*/
/* Transformation                                                            */

/*  this version supports T(x) = log(x) only                                 */
#define T(x)       (log(x))       /* transformation function                 */
#define T_deriv(x) (1./(x))       /* first derivative of transform funtion   */
#define T_inv(x)   (exp(x))       /* inverse of transform function           */

/*---------------------------------------------------------------------------*/
/* return codes for _unur_mvtdr_tp_bracket()                                 */

#define TP_LEFT    1              /* minimum in left point                   */
#define TP_MIDDLE  2              /* minimum in middle point                 */
#define TP_RIGHT   3              /* minimum in right point                  */
#define TP_BRACKET 4              /* bracket found                           */
 
/*---------------------------------------------------------------------------*/
/* since there is only file scope or program code, we abuse the              */
/* #include directive.                                                       */

/**  Public: User Interface (API)                                           **/
#include "mvtdr_newset.ch"

/**  Private                                                                **/
#include "mvtdr_init.ch"
#include "mvtdr_sample.ch"
#include "mvtdr_debug.ch"
#include "mvtdr_info.ch"

/*---------------------------------------------------------------------------*/
