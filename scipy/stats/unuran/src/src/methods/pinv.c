/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      pinv.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    Polynomial interpolation based INVersion of CDF              *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Compute values of CDF incrementally and interpolate resulting points *
 *      by polynomials.                                                      *
 *                                                                           *
 *      Integration:   adaptive Gauss-Lobatto integration with 5 points      *
 *      Interpolation: Newton recursion for interpolating polynomial         *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the PDF, center of distribution                           *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      pointer to CDF, pointer to derivative of PDF                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2008-2010 Wolfgang Hoermann and Josef Leydold             *
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
 *   [1] Derflinger, Gerhard, H{\"o}rmann, Wolfgang, and Leydold, Josef:     *
 *       Random Variate Generation by Numerical Inversion when only the      *
 *       Density Is Known. ACM Trans. Model. Comput. Simul. (2010),          *
 *       to appear.                                                          *
 *       See also Research Report Series of the Department of                *
 *       Statistics and Mathematics at WU Vienna, No. 90, June 2009.         *
 *       http://epub.wu.ac.at/, oai:epub.wu-wien.ac.at:epub-wu-01_f41        *
 *                                                                           *
 *   [2] H{\"o}rmann, Wolfgang and Leydold, Josef:                           *
 *       Continuous Random Variate Generation by Fast Numerical Inversion,   *
 *       ACM Trans. Model. Comput. Simul. 13(4), pp. 347--362, 2003.         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Method PINV combines numerical integration with interpolation of the    *
 *   inverse CDF.                                                            *
 *                                                                           *
 *   1.  Preprocessing:                                                      *
 *                                                                           *
 *   1a.   Estimate computationally relevant domain (support) of PDF         *
 *         (finite interval where PDF is above some threshold value).        *
 *            _unur_pinv_relevant_support()                                  *
 *            _unur_pinv_searchborder()                                      *
 *                                                                           *
 *   1b.   Compute area below PDF over relevant domain approximately.        *
 *            _unur_pinv_approx_pdfarea()                                    *
 *                                                                           *
 *   1c.   Compute computational domain where inverse CDF is approximated    *
 *         (interval where we safely can compute coefficients of             *
 *         interpolating polynomial).                                        *
 *            _unur_pinv_computational_domain()                              *
 *            _unur_pinv_cut()                                               *
 *            _unur_pinv_cut_bisect()                                        *
 *            _unur_pinv_cut_CDF()                                           *
 *                                                                           *
 *   1d.   Compute area below PDF over relevant domain with requested        *
 *         accuracy and store subinterval boundaries and corresponding       *
 *         from adaptive integration.                                        *
 *            _unur_pinv_pdfarea()                                           *
 *                                                                           *
 *   2.  Interpolation:                                                      *
 *                                                                           *
 *   2a.   Compute coefficients for interpolating polynomial for             *
 *         fixed (sub-) interval.                                            *
 *            _unur_pinv_newton_create()                                     *
 *            _unur_pinv_chebyshev_points()                                  *
 *            _unur_pinv_newton_cpoints()                                    *
 *                                                                           *
 *   2b.   Evaluate interpolating polynomial.                                *
 *            _unur_pinv_newton_eval()                                       *
 *                                                                           *
 *   2c.   Estimate approximation error for given interpolating polynomial.  *
 *            _unur_pinv_newton_maxerror()                                   *
 *            _unur_pinv_newton_testpoints()                                 *
 *            _unur_pinv_cubic_hermite_monotone()                            *
 *                                                                           *
 *                                                                           *
 *   Currently the following methods are implemented:                        *
 *                                                                           *
 *   Quadrature (Integration):                                               *
 *     Adaptive Gauss-Lobatto integration with 5 points.                     *
 *        see utils/lobatto.c                                                *
 *                                                                           *
 *   Interpolation:                                                          *
 *     Newton recursion for coefficients of polynomial                       *
 *     ("Newton interpolation").                                             *
 *                                                                           *
 *     Hermite interpolation [2] (as the limit case with double and          *
 *     tripple nodes).                                                       *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

/* enable additional code for development */
/* #define PINV_DEVEL */

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <urng/urng.h>
#include <tests/unuran_tests.h>
#include <utils/lobatto_source.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "pinv.h"
#include "pinv_struct.h"


/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* -- Global parameters                                                      */
#define MAX_ORDER   (17)
/* Maximum order of Newton interpolation polynomial */

#define PINV_UERROR_CORRECTION  (0.9)
/* PINV tries to create an approximation of the inverse CDF where the        */
/* U-error is bounded by the given u-resolution. However, the error caused   */
/* by cutting off the tails introduces some additional errors which must be  */
/* corrected. Thus the upper bound for the pure approximation error of the   */
/* interpolation is set to PINV_UERROR_CORRECTION * u-resolution.            */

#define PINV_DEFAULT_MAX_IVS  (10000)
/* Default for maximum number of subintervals */

/* -- Gauss-Lobatto integration                                              */

#define PINV_MAX_LOBATTO_IVS  (20001)
/* Maximum  number of subintervals for adaptive Gauss-Lobatto integration.   */
/* We keep this number fixed (independent of the maximum number of           */
/* subintervals for Newton interpolation), since the number of these         */
/* subintervals does neither depend on the order of the interpolating        */
/* polynomial nor on the requested u-resolution.                             */
/*                                                                           */
/* Remark: This parameter MUST NOT be smaller than 2 !!                      */

/* -- 1. Preprocessing                                                       */

#define PINV_PDFLLIM    (1.e-13)
/* Threshold value used for finding the boundary of the computational        */
/* domain. When starting the search at some point x0 then the search stops   */
/* when a point x is found where PDF(x) approx. PDF(x0) * PINV_PDFLLIM.      */

#define PINV_UERROR_AREA_APPROX  (1.e-5)
/* Tolerated relative area when computing the area below the PDF             */
/* approximately in Step 1b.                                                 */

#define PINV_TAILCUTOFF_FACTOR   (0.05)
#define PINV_TAILCUTOFF_MAX      (1.e-10) 
/* For unbounded domains the tails has to be cut off. We use the given       */
/* u-resolution for finding the cut points. (The probability for each of the */
/* chopped regions should be less than                                       */
/* HINV_TAILCUTOFF_FACTOR * u-resolution.)                                   */
/* However, the tail probabilities should not be greater than some threshold */
/* value, given by PINV_TAILCUTOFF_MAX which reflects the precision of the   */
/* used stream of uniform pseudo-random numbers (typically about 2^32).      */
/* However, for computational reasons we use a value that is at least twice  */
/* the machine epsilon for the right hand boundary.                          */

/* -- 2. Newton interpolation                                                */

#define PINV_UTOL_CORRECTION  (0.05)
/* We use a smaller tolerance when computing the Gauss-Lobatto integral for  */
/* the PDF between construction points of the Newton polynomial.             */

#define PINV_MAX_ITER_IVS    (10 * GEN->max_ivs)
/* maximum number of iterations for computing intervals for Newtwon          */
/* interpolation polynomials. Obviously it should be larger than the maximum */
/* number of intervals (or the latter can be reduced).                       */

#define PINV_GUIDE_FACTOR  (1)
/* relative size of guide table for finding the subinterval corresponding    */
/* to the given U-value.                                                     */


/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define PINV_VARIANT_PDF      0x0010u   /* use PDF and Lobatto integration
	                                   [ if not present, use CDF ]       */
#define PINV_VARIANT_UPOINTS  0x0040u   /* use Chebyshev points in u scale   */
#define PINV_VARIANT_KEEPCDF  0x0080u   /* keep table for integration        */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define PINV_DEBUG_REINIT    0x00000002u   /* print parameters after reinit  */
#define PINV_DEBUG_TABLE     0x00000010u   /* print table                    */
#define PINV_DEBUG_SEARCHBD  0x00010000u   /* trace search boundary          */
#define PINV_DEBUG_ITABLE    0x00020000u   /* print table of integral values */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define PINV_SET_ORDER          0x0001u  /* order of polynomial              */
#define PINV_SET_ORDER_COR      0x1000u  /* ... corrected                    */
#define PINV_SET_SMOOTH         0x0002u  /* smoothness of interpolant        */
#define PINV_SET_SMOOTH_COR     0x2000u  /* ... corrected                    */
#define PINV_SET_U_RESOLUTION   0x0004u  /* maximal error in u               */
#define PINV_SET_UPOINTS        0x0008u  /* use Chebyshev points in u scale  */
#define PINV_SET_BOUNDARY       0x0010u  /* boundary of computational region */
#define PINV_SET_SEARCHBOUNDARY 0x0020u  /* search for boundary              */
#define PINV_SET_VARIANT        0x0040u  /* variant of algorithm             */
#define PINV_SET_MAX_IVS        0x0080u  /* maximum number of subintervals   */
#define PINV_SET_KEEPCDF        0x0100u  /* keep table for integration       */

/*---------------------------------------------------------------------------*/

#define GENTYPE "PINV"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

/*........................*/
/*  file: pinv_newset.ch  */
/*........................*/

/* See pinv.h for 'new', 'set', and 'get' calls. */


/*......................*/
/*  file: pinv_init.ch  */
/*......................*/

static struct unur_gen *_unur_pinv_init (struct unur_par *par);
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

/* static int _unur_pinv_reinit (struct unur_gen *gen); */
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_create (struct unur_par *par);
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_check_par (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_pinv_clone (const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_free (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_make_guide_table (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_eval_PDF (double x, struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* call to PDF.                                                              */
/*---------------------------------------------------------------------------*/


/*........................*/
/*  file: pinv_sample.ch  */
/*........................*/

static double _unur_pinv_sample (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_eval_approxinvcdf (const struct unur_gen *gen, double u);
/*---------------------------------------------------------------------------*/
/* evaluate interpolation of inverse CDF at u.                               */
/*---------------------------------------------------------------------------*/

/* declared in pinv.h: */
/*   double unur_pinv_eval_approxinvcdf (const struct unur_gen *gen, double u); */
/*   int unur_pinv_estimate_error (const UNUR_GEN *gen, int samplesize, double *max_error, double *MAE); */


/*......................*/
/*  file: pinv_prep.ch  */
/*......................*/

static int _unur_pinv_preprocessing (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1. Find computational domain and compute PDF area.                        */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_relevant_support (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1a. Estimate computationally relevant domain (support) of PDF             */
/*     (finite interval where PDF is above some threshold value).            */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_searchborder (struct unur_gen *gen, double x0, double bound,
				       double *dom, int *search);
/*---------------------------------------------------------------------------*/
/* [1a.] find left or right hand border of relevant domain.                  */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_approx_pdfarea (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1b. Compute area below PDF over relevant domain approximately.            */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_pdfarea (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1d. Compute area below PDF with requested accuracy and                    */
/*     store intermediate results from adaptive integration.                 */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_computational_domain (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1c. Compute computational domain where inverse CDF is approximated        */
/*     (interval where we safely can compute coefficients of                 */
/*     interpolating polynomial).                                            */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_cut (struct unur_gen *gen, double w, double dw, double crit);
/*---------------------------------------------------------------------------*/
/* [1c.] calculate cut-off points for computational domain of distribution.  */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_cut_bisect (struct unur_gen *gen, double x0, double x1);
/*---------------------------------------------------------------------------*/
/* [1c.] calculate cut-off points as boudary of bounded support of PDF.      */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_computational_domain_CDF (struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* 1c. Compute computational domain where inverse CDF is approximated        */
/*     (interval where we safely can compute coefficients of                 */
/*     interpolating polynomial).                                            */
/*     Use CDF.                                                              */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_cut_CDF( struct unur_gen *gen, double dom, double x0, double ul, double uu );
/*---------------------------------------------------------------------------*/
/* [1c.] calculate cut-off points for computational domain of distribution   */
/*       using CDF.                                                          */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_Udiff (struct unur_gen *gen, double x, double h, double *fx);
/*---------------------------------------------------------------------------*/
/* compute difference CDF(x+h)-CDF(x) (approximately), where CDF is the      */
/* integral of the given (quasi-) density.                                   */
/*---------------------------------------------------------------------------*/


/*........................*/
/*  file: pinv_newton.ch  */
/*........................*/

static int _unur_pinv_create_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create table for Newton interpolation                                     */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_chebyshev_points (double *pt, int order, int smooth);
/*---------------------------------------------------------------------------*/
/* [2a.] Compute Chebyshev points.                                           */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_newton_cpoints (double *xval, int order, struct unur_pinv_interval *iv, 
				      double h, double *chebyshev, int smooth, int use_upoints);
/*---------------------------------------------------------------------------*/
/* [2a.] Compute points for construct interpolating polynomial.              */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_newton_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
				     double *xval);
/*---------------------------------------------------------------------------*/
/* 2a. Compute coefficients for Newton interpolation.                        */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_linear_create (struct unur_gen *gen, struct unur_pinv_interval *iv, 
				     double *xval);
/*---------------------------------------------------------------------------*/
/* [2a.] Compute coefficients for linear interpolation.                      */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_newton_eval (double q, double *ui, double *zi, int order);
/*---------------------------------------------------------------------------*/
/* 2b. evaluate Newton polynomial.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_pinv_newton_maxerror (struct unur_gen *gen, struct unur_pinv_interval *iv, 
					  double *xval, int use_linear);
/*---------------------------------------------------------------------------*/
/* 2c. estimate maximal error of Newton interpolation in subinterval         */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_newton_testpoints (double *utest, double ui[], int order);
/*---------------------------------------------------------------------------*/
/* [2c.] calculate the local maxima of the interpolation polynomial          */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_linear_testpoints (double *utest, double *ui, int order);
/*---------------------------------------------------------------------------*/
/* [2c.] create table of test points for linear interpolation                */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_cubic_hermite_is_monotone();
/*---------------------------------------------------------------------------*/
/* [2c.] check monotonicity of cubic Hermite interpolation                   */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_interval( struct unur_gen *gen, int i, double x, double cdfx );
/*---------------------------------------------------------------------------*/
/* make a new interval i with left boundary point x and CDF(x).              */
/*---------------------------------------------------------------------------*/

static int _unur_pinv_lastinterval( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* update size of array and set all uninitialized values to 0.               */
/*---------------------------------------------------------------------------*/


/*.......................*/
/*  file: pinv_debug.ch  */
/*.......................*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_init_start (const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* print before setup starts.                                                */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_init (const struct unur_gen *gen, int ok);
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_relevant_support (const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* print relevant domain                                                     */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_pdfarea (const struct unur_gen *gen, int approx);
/*---------------------------------------------------------------------------*/
/* print estimated area below PDF                                            */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_computational_domain (const struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* print computational domain                                                */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_intervals (const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print starting points or table for algorithms into LOG file.              */
/*---------------------------------------------------------------------------*/

static void _unur_pinv_debug_create_table (const struct unur_gen *gen,
					   int iter, int n_incr_h, int n_decr_h,
					   int n_use_linear);
/*---------------------------------------------------------------------------*/
/* print data that have been collected while creating polynomials.           */
/*---------------------------------------------------------------------------*/
#endif


/*......................*/
/*  file: pinv_info.ch  */
/*......................*/

#ifdef UNUR_ENABLE_INFO
static void _unur_pinv_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_pinv_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_pinv_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define PDF(x)  (_unur_pinv_eval_PDF((x),(gen)))      /* call to PDF         */
#define dPDF(x) (_unur_cont_dPDF((x),(gen->distr)))   /* call to derivate of PDF */
#define CDF(x)  (_unur_cont_CDF((x),(gen->distr)))    /* call to CDF         */

/*---------------------------------------------------------------------------*/

#define _unur_pinv_getSAMPLE(gen)  (_unur_pinv_sample)

/*---------------------------------------------------------------------------*/



/*---------------------------------------------------------------------------*/
/* since there is only file scope or program code, we abuse the              */
/* #include directive.                                                       */

/**  Public: User Interface (API)                                           **/
#include "pinv_newset.ch"

/**  Private                                                                **/
#include "pinv_init.ch"
#include "pinv_sample.ch"
#include "pinv_prep.ch"
#include "pinv_newton.ch"
#include "pinv_debug.ch"
#include "pinv_info.ch"

/*---------------------------------------------------------------------------*/
