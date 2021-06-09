/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      ninv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    numerical inversion of cumulative distribution function      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the CDF                                                   *
 *      newton's method: additional pointer to the PDF                       *
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
 *   [1] Neumaier A. (to be published):                                      *
 *       Introduction to numerical analysis,                                 *
 *       Cambridge University Press.                                         *
 *                                                                           *
 *   [2] Dahlquist, G. and Bj{\"o}rck, {\AA} (2008):                         *
 *       Numerical methods in scientific computing, Vol 1.,                  *
 *       SIAM, Philadelphia, PA.                                             *
 *                                                                           *
 *   [3] Brent, R. (1973):                                                   *
 *       Algorithms for Minimizing without Derivatives,                      *
 *       Prentice-Hall, Englewood Cliffs, NJ.                                *
 *       Republished by Dover Publications, Mineola, NY, 2002.               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Numerical inversion is a method for generating random variables          *
 *  using the CDF (and in case of newton's method the PDF).                  *
 *                                                                           *
 *  THEOREM:                                                                 *
 *     Let X be a random variable with CDF F(x).                             *
 *     Then the F(X) are  uniformly distributed.                             *
 *                                                                           *
 *  COROLLARY:                                                               *
 *     Starting with uniformly distributed random variables U,               *
 *     the F^(-1)(U) have F(x) as CDF.                                       *
 *                                                                           *
 *  Starting with an U, the task is to find a X fulfilling:                  *
 *    F(X) - U = 0.                                                          *
 *                                                                           *
 *  Numerical algorithms to find zeros that are used in NINV are variants of * 
 *  newton's method (damped newton to guarantee improvement) and             *
 *  the regula falsi ( stabilized regula falsi preserving sign change; at    *
 *  first an interval with sign change is determined).                       *
 *                                                                           *
 *  In both cases it is possible to specify the maximal number of            *
 *  iterations, a desired accuracy in X and starting values for the          *
 *  algorithms.                                                              *
 *  Instead of starting values it is also possible to use a table            *
 *  containing suitable starting values.                                     *
 *  If neither the table nor explicit starting values are used,              *
 *  NINV chooses as starting values:                                         *
 *     newton's method:  x:     CDF(x) = 0.5                                 *
 *     regula falsi:     x1,x2: CDF(x1) = 1 - CDF(x2) = 0.05                 *
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
#include "ninv.h"
#include "ninv_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* Starting interval for Regula Falsi includes this percentage of all        */
/* univariate random numbers (must be > 0. and < 1.)                         */
#define INTERVAL_COVERS  (0.5)

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

#define NINV_VARFLAG_NEWTON   0x1u   /* use Newton's method                  */
#define NINV_VARFLAG_REGULA   0x2u   /* use regula falsi [default]           */
#define NINV_VARFLAG_BISECT   0x4u   /* use bisection method                 */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define NINV_DEBUG_REINIT    0x00000002u   /* print parameters after reinit  */
#define NINV_DEBUG_TABLE     0x00000010u   /* print table                    */
#define NINV_DEBUG_CHG       0x00001000u   /* print changed parameters       */
#define NINV_DEBUG_SAMPLE    0x01000000u   /* trace sampling                 */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define NINV_SET_MAX_ITER     0x001u   /* number of maximal iterations       */
#define NINV_SET_X_RESOLUTION 0x002u   /* maximal tolerated relative x-error */
#define NINV_SET_U_RESOLUTION 0x004u   /* maximal tolerated (abs.) u-error   */
#define NINV_SET_START        0x008u   /* intervals at start (left/right)    */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

/*........................*/
/*  file: ninv_newset.ch  */
/*........................*/

/* See ninv.h for 'new', 'set', and 'get' calls. */


/*......................*/
/*  file: ninv_init.ch  */
/*......................*/

static struct unur_gen *_unur_ninv_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_ninv_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_create_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create the table with starting points                                     */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_compute_start( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* get starting points for numerical inversion                               */
/*---------------------------------------------------------------------------*/


/*........................*/
/*  file: ninv_sample.ch  */
/*........................*/

static double _unur_ninv_sample_newton( struct unur_gen *gen );
static double _unur_ninv_sample_regula( struct unur_gen *gen );
static double _unur_ninv_sample_bisect( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/


/*........................*/
/*  file: ninv_newton.ch  */
/*........................*/

static double _unur_ninv_newton( const struct unur_gen *gen, double u);
/*---------------------------------------------------------------------------*/
/* algorithm: newton method                                                  */
/*---------------------------------------------------------------------------*/


/*........................*/
/*  file: ninv_regula.ch  */
/*........................*/

static double _unur_ninv_regula( const struct unur_gen *gen, double u );
/*---------------------------------------------------------------------------*/
/* algorithm: regula falsi                                                   */
/*---------------------------------------------------------------------------*/

static double _unur_ninv_bisect( const struct unur_gen *gen, double u );
/*---------------------------------------------------------------------------*/
/* algorithm: bisection method                                               */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_bracket( const struct unur_gen *gen, double u, 
			       double *xl, double *fl, double *xu, double *fu );
/*---------------------------------------------------------------------------*/
/* find a bracket (enclosing interval) for root of CDF(x)-u.                 */
/*---------------------------------------------------------------------------*/

static int _unur_ninv_accuracy( const struct unur_gen *gen,
				double x_resol, double u_resol,
				double x0, double f0, double x1, double f1 );
/*---------------------------------------------------------------------------*/
/* check accuracy goal for approximate root.                                 */
/*---------------------------------------------------------------------------*/


/*.......................*/
/*  file: ninv_debug.ch  */
/*.......................*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into LOG file               */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_sample( const struct unur_gen *gen, 
				     double u, double x, double fx, int iter );
/*---------------------------------------------------------------------------*/
/* trace sampling.                                                           */
/*---------------------------------------------------------------------------*/

static void _unur_ninv_debug_chg_truncated( const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* trace changes of the truncated domain.                                    */
/*---------------------------------------------------------------------------*/
#endif


/*......................*/
/*  file: pinv_info.ch  */
/*......................*/

#ifdef UNUR_ENABLE_INFO
static void _unur_ninv_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif


/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_ninv_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_ninv_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */
#define CDF(x)    _unur_cont_CDF((x),(gen->distr))    /* call to CDF         */

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CONT *
_unur_ninv_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant) {
  case NINV_VARFLAG_NEWTON:
    return _unur_ninv_sample_newton;
  case NINV_VARFLAG_BISECT:
    return _unur_ninv_sample_bisect;
  case NINV_VARFLAG_REGULA:
  default:
    return _unur_ninv_sample_regula;
  }
} /* end of _unur_ninv_getSAMPLE() */

/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
/* since there is only file scope or program code, we abuse the              */
/* #include directive.                                                       */

/**  Public: User Interface (API)                                           **/
#include "ninv_newset.ch"

/**  Private                                                                **/
#include "ninv_init.ch"
#include "ninv_sample.ch"
#include "ninv_newton.ch"
#include "ninv_regula.ch"
#include "ninv_debug.ch"
#include "ninv_info.ch"

/*---------------------------------------------------------------------------*/
