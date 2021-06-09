/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tabl.h                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    rejection form piecewise constant hat                        *
 *              (Ahren's table method)                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a unimodal distribution                                 *
 *      produce random variate X with its density                            *
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
 *   [1] Ahrens J. H. (1993): Sampling from general distributions by         *
 *       suboptimal division of domains,                                     *
 *       Grazer Math. Berichte 319, 30pp.                                    *
 *                                                                           *
 *   [2] Ahrens J. H. (1995): An one-table method for sampling from          *
 *       continuous and discrete distributions,                              *
 *       Computing 54(2), pp. 127-146                                        *
 *                                                                           *
 *   [3] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *   SEE ALSO:                                                               *
 *   [4] Gilks, W. R. and Wild,  P. (1992):                                  *
 *       Adaptive rejection sampling for Gibbs sampling,                     *
 *       Applied Statistics 41, pp. 337-348                                  *
 *                                                                           *
 *   [5] Zaman, A. (1996), Generation of Random Numbers from an Arbitrary    *
 *       Unimodal Density by Cutting Corners, unpublished manuskript         *
 *       available at http://chenab.lums.edu.pk/~arifz/                      *
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
#include "tabl.h"
#include "tabl_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

#define TABL_DEFAULT_COMPUTATION_LIMIT  1.e20
/*
   The domain has to be finite for method TABL.
   Thus the domain has to be truncated where the regions (tails) that are
   chopped off must not of "computational relevance".
   By default the points for chopping of the tails are +/- this figure.
*/

#define TABL_N_RETRY_DARS  5
#define TABL_N_RUN_ARS     10
/* 
   Sometimes DARS fails then we need random splitting points for the intervals,
   i.e. we have to run ARS a few times an continue (retry) with DARS.
   These figures give the number of retries and the number of sampling using
   ARS, respectively.
*/

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

/* how to sample from genertor object */
#define TABL_VARMASK_VARIANT      0x000fu   /* indicates which variant       */
#define TABL_VARIANT_IA           0x0001u   /* use immediate acceptance      */

#define TABL_VARIANT_FAST         0x0002u   /* use single array for data 
					       (not implemented)             */

/* indicate how to split interval */
#define TABL_VARMASK_SPLIT        0x00f0u  /* split at        computation     convergence of hat */
#define TABL_VARFLAG_SPLIT_POINT  0x0010u  /* sampled point    none            slowest          */
#define TABL_VARFLAG_SPLIT_MEAN   0x0020u  /* mean point       slower          better           */
#define TABL_VARFLAG_SPLIT_ARC    0x0040u  /* "arcmean"        very slow       very good for almost unbounded domain */

/* indicate if starting intervals have to be split */
#define TABL_VARFLAG_USEEAR       0x0100u  /* use equal area rule (SPLIT A in [1])   */
#define TABL_VARFLAG_USEDARS      0x0200u  /* use main subdivisions (SPLIT B in [1]) 
                                             (= derandomized ARS)                   */
#define TABL_VARFLAG_PEDANTIC     0x0400u  /* whether pedantic checking is used */
#define TABL_VARFLAG_VERIFY       0x0800u  /* flag for verifying mode        */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define TABL_DEBUG_IV        0x00000100u /* show intervals                   */
#define TABL_DEBUG_IV_START  0x00000200u /* show starting intervals          */
#define TABL_DEBUG_EAR       0x00000400u /* show intervals after EAR         */
#define TABL_DEBUG_DARS      0x00000800u /* indicate DARS in LOG file        */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TABL_SET_GUIDEFACTOR      0x0001u
#define TABL_SET_SLOPES           0x0004u
#define TABL_SET_AREAFRACTION     0x0008u
#define TABL_SET_MAX_IVS          0x0010u
#define TABL_SET_MAX_SQHRATIO     0x0020u
#define TABL_SET_N_STP            0x0040u
#define TABL_SET_STP              0x0080u
#define TABL_SET_BOUNDARY         0x0100u
#define TABL_SET_USE_EAR          0x0200u
#define TABL_SET_USE_DARS         0x0400u
#define TABL_SET_DARS_FACTOR      0x0800u

/*---------------------------------------------------------------------------*/

#define GENTYPE "TABL"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tabl_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_tabl_rh_sample( struct unur_gen *gen );
static double _unur_tabl_rh_sample_check( struct unur_gen *gen );
static double _unur_tabl_ia_sample( struct unur_gen *gen );
static double _unur_tabl_ia_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_get_intervals_from_slopes( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals from slopes                                    */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_get_intervals_from_cpoints( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute starting intervals from given cpoints                             */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_compute_intervals( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute all intervals (by splitting starting intervals/slopes)            */
/*---------------------------------------------------------------------------*/

static struct unur_tabl_interval *
_unur_tabl_run_equalarearule( struct unur_par *par, struct unur_gen *gen, struct unur_tabl_interval *iv_slope );
/*---------------------------------------------------------------------------*/
/* split starting intervals according to [1]                                 */
/* SPLIT A (equal areas rule)                                                */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_run_dars( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* run derandomized adaptive rejection sampling.                             */
/* (split starting intervals according to [1] SPLIT B,                       */
/* but instead of the iteration in [1] use "arcmean".                        */
/*---------------------------------------------------------------------------*/

static int
_unur_tabl_split_interval( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			   double x, double fx, unsigned split_mode );
/*---------------------------------------------------------------------------*/
/* split interval (replace old one by two new ones in same place)            */
/*---------------------------------------------------------------------------*/

static int
_unur_tabl_improve_hat( struct unur_gen *gen, struct unur_tabl_interval *iv, 
			double x, double fx );
/*---------------------------------------------------------------------------*/
/* improve hat function by splitting interval                                */
/*---------------------------------------------------------------------------*/

static int _unur_tabl_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static double _unur_tabl_eval_cdfhat( struct unur_gen *gen, double x );
/*---------------------------------------------------------------------------*/
/* evaluate CDF of hat at x.                                                 */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_init_start( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after (almost empty generator) object has been created.             */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_init_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print header before runniung derandomized adaptive rejection sampling.    */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tabl_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals.                                                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_tabl_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_tabl_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_tabl_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)    _unur_cont_PDF((x),(gen->distr))    /* call to PDF         */

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CONT *
_unur_tabl_getSAMPLE( struct unur_gen *gen )
{
  if (gen->variant & TABL_VARIANT_IA)
    /* immediate acceptance */
    return (gen->variant & TABL_VARFLAG_VERIFY) 
      ? _unur_tabl_ia_sample_check 
      : _unur_tabl_ia_sample;
  else
    /* "classical" acceptance/rejection method */
    return (gen->variant & TABL_VARFLAG_VERIFY) 
      ? _unur_tabl_rh_sample_check 
      : _unur_tabl_rh_sample;
} /* end of _unur_tabl_getSAMPLE() */

/*---------------------------------------------------------------------------*/
/* since there is only file scope or program code, we abuse the              */
/* #include directive.                                                       */

/**  Public: User Interface (API)                                           **/
#include "tabl_newset.ch"

/**  Private                                                                **/
#include "tabl_init.ch"
#include "tabl_sample.ch"
#include "tabl_debug.ch"
#include "tabl_info.ch"

/* #include "tabl_fast.ch" ... not implemented (only experimental code)      */

/*---------------------------------------------------------------------------*/
