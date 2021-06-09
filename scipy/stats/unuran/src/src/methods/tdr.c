/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr.c                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PDF of a T-concave distribution                                *
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
 *   [1] Hoermann W. (1995): A rejection technique for sampling from         *
 *       T-concave distributions, ACM TOMS 21, p. 182-193                    *
 *                                                                           *
 *   [2] Chen, H. C. and Asau, Y. (1974): On generating random variates      *
 *       from an empirical distribution, AIIE Trans. 6, pp. 163-166          *
 *                                                                           *
 *   [3] Gilks, W. R. and Wild,  P. (1992): Adaptive rejection sampling      *
 *       for Gibbs sampling, Applied Statistics 41, pp. 337-348              *
 *                                                                           *
 *   [4] Evans, M. and Swartz, T. (1998): Random variable generation         *
 *       using concavity properties of transformed densities,                *
 *       Journal of Computational and Graphical Statistics 7(4), pp. 514-528 *
 *                                                                           *
 *   [5] Derflinger, G. and Hoermann, W. (1999): The optimal selection of    *
 *       hat functions for rejection algorithms, preprint                    *
 *                                                                           *
 *   [6] Leydold, J. (1999): Automatic Sampling with the ratio-of-uniforms   *
 *       method, preprint, 15pp.                                             *
 *                                                                           *
 *   [7] Leydold, J and Hoermann, W. (2001): ?????, preprint                 *
 *                                                                           *
 *****************************************************************************

 *****************************************************************************
 * Variant GW (Gilks & Wild)                                                 *
 *****************************************************************************
 *                                                                           *
 * Transformed density rejection (see [1,3,4]) is an acceptance/rejection    *
 * technique that uses the fact, that the probability density function f(x)  *
 * for many distribution is T-concave, to construct a hat function.          *
 * That is, there exists a transformation T(x), such that T(f(x)) is         *
 * concave. Then it is easy to construct a majorizing or hat function Th(x)  *
 * for the transformed density T(F(x)) by the pointwise minima of several    *
 * tangents. Transforming this back into the original scale gives the        *
 * hat h(x) = T^(-1)(Th(x)). Squeezes can be constructed by secants.         *
 *                                                                           *
 * However, transformation T(x) has to satisfy (see [1]):                    *
 *                                                                           *
 *  (1) T(f(x)) is concave.                                                  *
 *  (2) T is differentiable.                                                 *
 *  (3) T is strictly monotonically increasing (i.e., T'(x) > 0),            *
 *      which implies that T^(-1) exists.                                    *
 *  (4) the area below the hat is finite.                                    *
 *  (5) it is easy to generate from the hat distribution.                    *
 *                                                                           *
 * We use a family T_c(x) of distribution defined by (see [1]):              *
 *                                                                           *
 *           c    |   T_c(x)                                                 *
 *       ---------+---------------                                           *
 *           0    |   log(x)                                                 *
 *         -1/2   |   -1/sqrt(x)                                             *
 *        -1<c<0  |   -x^c                                                   *
 *                                                                           *
 * Figures 1 and 2 show the situation for the standard normal distribution   *
 * with T(x) = log(x) and two construction points at -1 and 0.3 for the      *
 * tangents.                                                                 *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *           figure 1: transformed density with hat and squeeze              *
 *                     T(x)   = log(x)                                       *
 *                     PDF(x) = exp(-x^2)                                    *
 *                                                                           *
 *                                .                                          *
 *     -3        -2        -1     ....0         1         2         3        *
 *    --+---------+---------X----.--*****-X-----+---------+---------+--->    *
 *                              .**     / **.                                *
 *                             *     /       *...                            *
 *                            *   /           *  ....                        *
 *             construction  * / squeeze       *     ...                     *
 *                   point  X                   *       ...                  *
 *                         *                     *         ....              *
 *                        *                       *            ...           *
 *                       *                         *              ..         *
 *                      .                                           .        *
 *                     .*                           *                        *
 *                    .*                             *                       *
 *                   .                                                       *
 *                  . *                               *                      *
 *                 . *                                 *                     *
 *                .                                                          *
 *               .  * transformed PDF                   *                    *
 *  transformed .                                                            *
 *         hat .   *                                     *                   *
 *            .                                                              *
 *           .    *                                       *                  *
 *          .                                                                *
 *         .     *                                         *                 *
 *        .                                                                  *
 *       .      *                                           *                *
 *      .                                                                    *
 *             *                                             *               *
 *                                                                           *
 *                                                                           *
 *            *                                               *              *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *                   figure 2: density with hat and squeeze                  *
 *                             T(x)   = log(x)                               *
 *                             PDF(x) = exp(-x^2)                            *
 *                                                                           *
 *                                ..                                         *
 *                                . .                                        *
 *                                .  .                                       *
 *                               .    .                                      *
 *                               .   ***.                                    *
 *                               . **   **                                   *
 *                              . *       X                                  *
 *                              .*       /|*                                 *
 *                              .*      / |* .                               *
 *                             .*      /  | * .                              *
 *                             .*      /  | *  .                             *
 *                             *      /   |  *  ..                           *
 *                             *     /    |  *    .                          *
 *                            *    /      |   *    ..                        *
 *                           *   /squeeze |    *     .                       *
 *         construction      * /          |    *      ..                     *
 *                point --- X/            |     *       ..                   *
 *                         *|             |      *        ..                 *
 *                         *|             |      *          ..               *
 *                        * |             |       *           ..             *
 *                        * |             |       *             ...          *
 *                       *  |             |        *               ..        *
 *                      .*  |             |        *                         *
 *                     .*   |             |         *                        *
 *                    .*    |             |          *                       *
 *                   . *    |             |           *                      *
 *                  . *     |             |           *                      *
 *           hat  .. *      |             |            *                     *
 *               .  *       |             |             *                    *
 *             ..  * PDF    |             |              *                   *
 *          ...   *         |             |               *                  *
 *       ...    **          |             |                **                *
 *      .     **            |             |                  **              *
 *        ****              |             |                    ****          *
 *    --**--------+---------X---------+---X-----+---------+--------**--->    *
 *     -3        -2        -1         0         1         2         3        *
 *                                                                           *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * To generate from the hat distribution we have to partition the domain     *
 * of the PDF by means of the construction points of the tangents into       *
 * several intervals. Each interval has to be divided into two parts         *
 * according to the tangent line that defines the hat function.              *
 * Although the interval can be divided at any point in the interval,        *
 * the best possible splitting point obviously is the intersection point     *
 * of respective the tangents at the left and right construction point of    *
 * the interval. However there are two cases when this choice is not         *
 * appropriate: (1) The computation of this point is unstable. Then we use   *
 * the mean point of the interval. (2) One of the two construction points    *
 * is very much greater (or smaller) than the other one. Since we generate   *
 * points by adding or subtracting a random number to the construction       *
 * point, it then happens, that we have the difference of two large numbers, *
 * which results in serious roundoff error. Then it is better to choose the  *
 * intersection point such that we have only one interval with the point     *
 * nearer to the origin as the only construction point.                      *
 * (This procedure is a little bit different from [1] or [3], but the author *
 * of this program finds the implemented version more convenient.) Then we   *
 * have to compute the area below the hat for each of the two parts of the   *
 * interval i (A_l^i for the left hand part, A_r^i for the right hand part)  *
 * for each interval. Then we have to sample from a discrete random variate  *
 * with probability vector proportional to (A_l^0,A_r^0; A_l^1,A_r^1; ...;   *
 * A_l^n,A_r^n) to get one of the intervals. We use indexed search (or       *
 * guide tables) to perform this task ([2], see also description of DIS).    *
 * Then we can sample from the hat distribution in this part of the          *
 * interval by inversion. Notice that we can use reuse the uniform random    *
 * number from the first step (i.e., sampling an interval) for this second   *
 * step without danger. The whole generation can be seen as inversion from   *
 * the hat distribution.                                                     *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Algorithm TDR                                                             *
 *                                                                           *
 * [Required]                                                                *
 * PDF f(x), transformation T(x), construction points c_1,...,c_n            *
 *                                                                           *
 * [Setup]                                                                   *
 *  1: Construct hat h(x) and squeeze s(x).                                  *
 *  2: Foreach Interval (c_i,c_{i+1}) Do                                     *
 *  3:    Compute areas A_l^i, A_r^i, A_squeeze.                             *
 *                                                                           *
 * [Generate]                                                                *
 *  4: Generate I proportional to (A_l^1,A_r^1; ...).                        *
 *  5: Generate X with PDF proportional to h|I.                              *
 *  6: Generate U ~ U(0,1).                                                  *
 *  7: If U * h(x) <= min(f(c_i),f(c_{i+1})) Return X.                       *
 *  8: If U * h(x) <= s(X) Return X.                                         *
 *  9: If U * h(x) <= f(X) Return X.                                         *
 * 10: Goto 4.                                                               *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * There are several alternatives for the choice of the construction points  *
 *                                                                           *
 * (1) Adaptive rejection sampling:                                          *
 *     Use two points on both sides of the mode of the PDF (provided         *
 *     that the hat has finite area).                                        *
 *     Whenever the PDF has to be evaluated in step 8, add a new             *
 *     construction point at X.                                              *
 *                                                                           *
 * (2) Optimal construction points:                                          *
 *     [5] have show an asymptotic formula for optimal construction points,  *
 *     i.e., with minimiza the ration Ahat/Asqueeze for a given number of    *
 *     points. Although this formula is only valid for infinitely many       *
 *     construction points, it works very well even for 3 or 4 points.       *
 *                                                                           *
 *     [1] gives a formula for three optimal construction points.            *
 *                                                                           *
 * (3) Empirical formulas:                                                   *
 *     [6] uses points c_i, such that arctan(c_i) are equidistributed.       *
 *     ("equiangular rule").                                                 *
 *     In most cases it has turned out as an acceptable good choice.         *
 *     (We made one change to [6]: If the mode is known, all points are      *
 *     moved, such the mode becomes the center of these points.)             *
 *                                                                           *
 * This implementation uses (3) to get good points for the start and adds    *
 * additional construction points by means of method (1).                    *
 * (The number of starting points and the maximum number of points can       *
 * be set for each generator.)                                               *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 * Implementation details:                                                   *
 *                                                                           *
 * (1) If possible the boundary points of the domain are used as             *
 *     construction points. (So do not include these in a list of starting   *
 *     points.) As a consequence 0 starting points are allowed.              *
 *                                                                           *
 * (2) To make use of the automatic generation of starting points call       *
 *     `unur_set_cpoints() with the NULL pointer as the last argument        *
 *     and the number of construction points (besides the boundary points    *
 *     and the mode) as its third argument.                                  *
 *                                                                           *
 * (3) If the mode is given, we use a tangent with slope 0 at this point;    *
 *     even when this does not result in the best possible hat (e.g. for     *
 *     the exponential distribution).                                        *
 *                                                                           *
 * (4) Starting points and the mode outside the domain are ignored.          *
 *                                                                           *
 * (5) If the given starting points does not result in a hat with finite     *
 *     area, the program tries to find some proper additional construction   *
 *     points by splitting interval with infinite area. (Here we use the     *
 *     "equiangular rule" again.)                                            *
 *                                                                           *
 * (6) The use of the additional "fast" squeeze in step 7 is due to a        *
 *     suggestion of G. Derflinger. It speeds up the generation about 20 %   *
 *     when many construction points are used.                               *
 *                                                                           *
 *****************************************************************************

 *****************************************************************************
 * Variant PS (Proportional Squeezes)                                        *
 *****************************************************************************
 *                                                                           *
 * The original algorithm can be modified in the following way:              *
 * Use squeezes that are proportional to the hat function in each of the     *
 * intervals. This decreases the marginal generation times;                  *
 * see [7] for details.                                                      *
 *                                                                           *
 * Technical remarks:                                                        *
 * The boundary points of the intervals are now the intersection points      *
 * of the tangents to the transformed density. However we have used the      *
 * same structure for an interval as in variant GW. To have a similar        *
 * code we have stored the left boundary of an interval together with the    *
 * construction point. Thus the right boundary point is stored in the        *
 * interval structure in the list of intervals.                              *
 * We also have added an empty interval structure in this list that just     *
 * marks the right boundary of the domain of the PDF.                        *
 *                                                                           *
 *****************************************************************************

 *****************************************************************************
 * Variant IA (Immedate Acceptance)                                          *
 *****************************************************************************
 *                                                                           *
 * IA is a composition method that uses squeezes that are proportional to    *
 * the hat function in each of the intervals (as in PS). We immediately      *
 * accept all points below the squeeze and use an acceptance/rejection       *
 * technique only we fall in the region between hat and squeeze. It works    *
 * very well since it is easy to generate points in this region (in          *
 * opposition to variant GW).                                                *
 *                                                                           *
 * For technical details see variant PS.                                     *
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
#include "tdr.h"
#include "tdr_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define TDR_VARMASK_T          0x000fu   /* indicates transformation         */
#define TDR_VAR_T_SQRT         0x0001u   /* T(x) = -1/sqrt(x)                */
#define TDR_VAR_T_LOG          0x0002u   /* T(x) = log(x)                    */
#define TDR_VAR_T_POW          0x0003u   /* T(x) = -x^c                      */

#define TDR_VARMASK_VARIANT    0x00f0u   /* indicates which variant          */
#define TDR_VARIANT_GW         0x0010u   /* original variant (Gilks&Wild)    */
#define TDR_VARIANT_PS         0x0020u   /* use proportional squeeze         */
#define TDR_VARIANT_IA         0x0030u   /* use immediate acceptance
					    (requires prop. squeeze)         */

#define TDR_VARFLAG_VERIFY     0x0100u   /* flag for verifying mode          */
#define TDR_VARFLAG_USECENTER  0x0200u   /* whether center is used as cpoint or not */
#define TDR_VARFLAG_USEMODE    0x0400u   /* whether mode is used as cpoint or not */
#define TDR_VARFLAG_PEDANTIC   0x0800u   /* whether pedantic checking is used */
#define TDR_VARFLAG_USEDARS    0x1000u   /* whether DARS is used in setup or not */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define TDR_DEBUG_IV           0x00000010u
#define TDR_DEBUG_SPLIT        0x00010000u
#define TDR_DEBUG_DARS         0x00020000u
#define TDR_DEBUG_SAMPLE       0x01000000u
#define TDR_DEBUG_REINIT       0x00000020u  /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define TDR_SET_CENTER         0x0002u
#define TDR_SET_STP            0x0001u
#define TDR_SET_N_STP          0x0002u
#define TDR_SET_PERCENTILES    0x0004u
#define TDR_SET_N_PERCENTILES  0x0008u
#define TDR_SET_RETRY_NCPOINTS 0x0010u
#define TDR_SET_GUIDEFACTOR    0x0020u
#define TDR_SET_C              0x0040u
#define TDR_SET_MAX_SQHRATIO   0x0080u
#define TDR_SET_MAX_IVS        0x0100u
#define TDR_SET_USE_DARS       0x0200u
#define TDR_SET_DARS_FACTOR    0x0400u

/*---------------------------------------------------------------------------*/

#define GENTYPE "TDR"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdr_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_make_gen( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Make generator object.                                                    */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_gw_sample( struct unur_gen *generator );
static double _unur_tdr_gw_sample_check( struct unur_gen *generator );
static double _unur_tdr_ps_sample( struct unur_gen *generator );
static double _unur_tdr_ps_sample_check( struct unur_gen *generator );
static double _unur_tdr_ia_sample( struct unur_gen *generator );
static double _unur_tdr_ia_sample_check( struct unur_gen *generator );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_gw_eval_invcdfhat( const struct unur_gen *generator, double u,
					   double *hx, double *fx, double *sqx,
					   struct unur_tdr_interval **iv,
					   struct unur_tdr_interval **cpt );
static double _unur_tdr_ps_eval_invcdfhat( const struct unur_gen *generator, double u,
					   double *hx, double *fx, double *sqx,
					   struct unur_tdr_interval **iv );
/*---------------------------------------------------------------------------*/
/* auxiliary routines to evaluate the inverse of the hat CDF.                */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_tdr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_starting_cpoints( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create list of construction points for starting segments.                 */
/* if user has not provided such points compute these by means of the        */
/* "equi-angle rule".                                                        */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_starting_intervals( struct unur_gen *gen );
static int _unur_tdr_gw_starting_intervals( struct unur_gen *gen );
static int _unur_tdr_ps_starting_intervals( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute intervals from given starting construction points.                */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_run_dars( struct unur_gen *gen );
static int _unur_tdr_gw_dars( struct unur_gen *gen );
static int _unur_tdr_ps_dars( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* run derandomized adaptive rejection sampling.                             */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_gw_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv );
static int _unur_tdr_ps_interval_parameter( struct unur_gen *gen, struct unur_tdr_interval *iv );
/*---------------------------------------------------------------------------*/
/* compute all necessary data for interval.                                  */
/*---------------------------------------------------------------------------*/

static struct unur_tdr_interval *_unur_tdr_interval_new( struct unur_gen *gen, 
							 double x, double fx, int is_mode );
/*---------------------------------------------------------------------------*/
/* make a new segment with left construction point x.                        */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_tangent_intersection_point( struct unur_gen *gen,
						 struct unur_tdr_interval *iv, double *ipt );
/*---------------------------------------------------------------------------*/
/* compute cutting point of interval into left and right part.               */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_interval_area( struct unur_gen *gen, struct unur_tdr_interval *iv,
				       double slope, double x );
/*---------------------------------------------------------------------------*/
/* compute area below piece of hat or squeeze in interval.                   */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_interval_xxarea( struct unur_gen *gen, struct unur_tdr_interval *iv,
					 double slope, double x );
/*---------------------------------------------------------------------------*/
/* compute the interal of x times hat or squeeze ("expected value").         */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_eval_intervalhat( struct unur_gen *gen,
					  struct unur_tdr_interval *iv, double x );
/*---------------------------------------------------------------------------*/
/* evaluate hat at x in interval.                                            */
/*---------------------------------------------------------------------------*/

static double _unur_tdr_eval_cdfhat( struct unur_gen *gen, double x );
/*---------------------------------------------------------------------------*/
/* evaluate CDF of hat at x.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_gw_interval_split( struct unur_gen *gen, 
					struct unur_tdr_interval *iv_old, double x, double fx );
static int _unur_tdr_ps_interval_split( struct unur_gen *gen, 
					struct unur_tdr_interval *iv_old, double x, double fx );
/*---------------------------------------------------------------------------*/
/* split am interval point x. return 0 if not successful.                    */                                           
/*---------------------------------------------------------------------------*/

static int _unur_tdr_gw_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv, 
				     double x, double fx);
static int _unur_tdr_ps_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv, 
				     double x, double fx);
/*---------------------------------------------------------------------------*/
/* improve hat function by splitting interval                                */
/*---------------------------------------------------------------------------*/

static int _unur_tdr_make_guide_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make a guide table for indexed search.                                    */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_init_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after (almost empty generator) object has been created.             */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_init_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_dars_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print header before runniung derandomized adaptive rejection sampling.    */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_dars_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has run derandomized adaptive rejection sampling.   */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_reinit_start( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before reinitialization of generator starts.                        */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_reinit_retry( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before second trial of reinitialization of generator starts.        */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_reinit_finished( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been reinitialized.                             */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_free( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print before generater is destroyed.                                      */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas );
static void _unur_tdr_gw_debug_intervals( const struct unur_gen *gen, int print_areas );
static void _unur_tdr_ps_debug_intervals( const struct unur_gen *gen, int print_areas );
/*---------------------------------------------------------------------------*/
/* print data for intervals                                                  */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_gw_debug_sample( const struct unur_gen *gen,
				       const struct unur_tdr_interval *iv,
				       const struct unur_tdr_interval *pt, 
				       double x, double fx, double hx, double sqx );
static void _unur_tdr_ps_debug_sample( const struct unur_gen *gen, 
				       const struct unur_tdr_interval *iv,
				       double x, double fx, double hx, double sqx );
/*---------------------------------------------------------------------------*/
/* print data while sampling from generators.                                */
/*---------------------------------------------------------------------------*/

static void _unur_tdr_gw_debug_split_start( const struct unur_gen *gen, 
					    const struct unur_tdr_interval *iv,
					    double x, double fx );
static void _unur_tdr_gw_debug_split_stop( const struct unur_gen *gen, 
					   const struct unur_tdr_interval *iv_left,
					   const struct unur_tdr_interval *iv_right );
static void _unur_tdr_ps_debug_split_start( const struct unur_gen *gen, 
					    const struct unur_tdr_interval *iv_left,
					    const struct unur_tdr_interval *iv_right,
					    double x, double fx );
static void _unur_tdr_ps_debug_split_stop( const struct unur_gen *gen, 
					   const struct unur_tdr_interval *iv_left,
					   const struct unur_tdr_interval *iv_middle,
					   const struct unur_tdr_interval *iv_right );
/*---------------------------------------------------------------------------*/
/* print before and after an interval has been split (not / successfully).   */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_tdr_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif


/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_tdr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_tdr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

#define PDF(x)     _unur_cont_PDF((x),(gen->distr))      /* call to PDF      */
#define dPDF(x)    _unur_cont_dPDF((x),(gen->distr))     /* call to derivative of PDF */
#define logPDF(x)  _unur_cont_logPDF((x),(gen->distr))   /* call to logPDF   */
#define dlogPDF(x) _unur_cont_dlogPDF((x),(gen->distr))  /* call to derivative of log PDF */

/*---------------------------------------------------------------------------*/

static UNUR_SAMPLING_ROUTINE_CONT *
_unur_tdr_getSAMPLE( struct unur_gen *gen )
{
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:    /* original variant (Gilks&Wild) */
    return (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_gw_sample_check : _unur_tdr_gw_sample;
  case TDR_VARIANT_IA:    /* immediate acceptance */
    return (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ia_sample_check : _unur_tdr_ia_sample;
  case TDR_VARIANT_PS:    /* proportional squeeze */
  default:
    return (gen->variant & TDR_VARFLAG_VERIFY) ? _unur_tdr_ps_sample_check : _unur_tdr_ps_sample;
  }
} /* end of _unur_tdr_getSAMPLE() */

/*---------------------------------------------------------------------------*/
/* since there is only file scope or program code, we abuse the              */
/* #include directive.                                                       */

/**  Public: User Interface (API)                                           **/
#include "tdr_newset.ch"

/**  Private                                                                **/
#include "tdr_init.ch"
#include "tdr_sample.ch"
#include "tdr_debug.ch"
#include "tdr_info.ch"

/*---------------------------------------------------------------------------*/
