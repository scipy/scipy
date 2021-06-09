/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: tdr.h                                                             *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method TDR                                *
 *         (Transformed Density Rejection)                                   *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
 *****************************************************************************/

/* 
   =METHOD  TDR   Transformed Density Rejection

   =UP  Methods_for_CONT

   =REQUIRED T-concave PDF, dPDF

   =OPTIONAL mode

   =SPEED Set-up: slow, Sampling: fast

   =REINIT supported

   =REF  [GWa92] [HWa95] [HLD04: Cha.4]

   =DESCRIPTION
      TDR is an acceptance/rejection method that uses the concavity of a
      transformed density to construct hat function and squeezes
      automatically. Such PDFs are called T-concave. Currently the
      following transformations are implemented and can be selected by
      setting their @code{c}-values by a unur_tdr_set_c() call:

      @table @code
      @item c = 0
      T(x) = log(x)
      @item c = -0.5
      T(x) = -1/sqrt(x) @ @ @ @ @ (Default)
      @end table
   
      In future releases the transformations T(x) = -(x)^c will be
      available for any c with 0 > c > -1.
      Notice that if a PDF is T-concave for a c then it also T-concave
      for every c'<c. However the performance decreases when c' is
      smaller than c. For computational reasons we suggest the usage of 
      c = -0.5 (this is the default). 
      For c <= -1 the hat is not bounded any more if the domain of the
      PDF is unbounded. But in the case of a bounded domain using
      method TABL is preferred to a TDR with c < -1 (except in a few
      special cases).
      
      We offer three variants of the algorithm. 

      @table @code
      @item GW
      squeezes between construction points
      @item PS
      squeezes proportional to hat function  @ @ @ @ @ (Default)
      @item IA
      same as variant PS but uses a compositon method with
      ``immediate acceptance'' in the region below the squeeze.
      @end table

      @code{GW} has a slightly faster setup but higher marginal generation
      times.
      @code{PS} is faster than @code{GW}. @code{IA} uses less uniform
      random numbers and is therefore faster than @code{PS}.

      It is also possible to evaluate the inverse of the CDF of the hat distribution
      directly using the unur_tdr_eval_invcdfhat() call.
      
      There are lots of parameters for these methods, see below.
      
      It is possible to use this method for correlation induction by
      setting an auxiliary uniform random number generator via the
      unur_set_urng_aux() call. (Notice that this must be done after a
      possible unur_set_urng() call.)
      When an auxiliary generator is used then the number of
      uniform random numbers from the first URNG that are used for one
      generated random variate is constant and given in the following table:

      @table @code
      @item GW ... 2
      @item PS ... 2
      @item IA ... 1
      @end table
      
      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_tdr_set_verify() and unur_tdr_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.

      For densities with modes not close to 0 it is suggested to set
      either the mode or the center of the distribution by the
      unur_distr_cont_set_mode() or unur_distr_cont_set_center() call.
      The latter is the approximate location of the mode or the mean
      of the distribution. This location provides some information
      about the main part of the PDF and is used to avoid numerical
      problems.

      It is possible to use this method for generating from truncated
      distributions. It even can be changed for an existing generator
      object by an unur_tdr_chg_truncated() call.

      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.

      @emph{Important:} The ratio between the area below the hat and
      the area below the squeeze changes when the sampling region is
      restricted. Especially it becomes (very) small when sampling
      from the (far) tail of the distribution. Then it is better to
      create a new generator object for the tail of the distribution
      only.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_tdr_new( const UNUR_DISTR* distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_tdr_set_c( UNUR_PAR *parameters, double c );
/* 
   Set parameter @var{c} for transformation T. 
   Currently only values between 0 and -0.5 are allowed.
   If @code{c} is between 0 and -0.5 it is set to -0.5.

   Default is @code{-0.5}.
*/


int unur_tdr_set_variant_gw( UNUR_PAR *parameters );
/* 
   Use original version with squeezes between construction points as
   proposed by Gilks & Wild  (1992).
*/

int unur_tdr_set_variant_ps( UNUR_PAR *parameters );
/* 
   Use squeezes proportional to the hat function. This is faster than
   the original version.
   This is the default.
*/

int unur_tdr_set_variant_ia( UNUR_PAR *parameters );
/* 
   Use squeezes proportional to the hat function together with a
   composition method that required less uniform random numbers.
*/

int unur_tdr_set_usedars( UNUR_PAR *parameters, int usedars );
/* 
   If @var{usedars} is set to TRUE, ``derandomized adaptive rejection
   sampling'' (DARS) is used in setup.
   Intervals where the area between hat and squeeze is too
   large compared to the average area between hat and squeeze
   over all intervals are split.
   This procedure is repeated until the ratio between area below squeeze
   and area below hat exceeds the bound given by 
   unur_tdr_set_max_sqhratio() call or the maximum number of intervals is 
   reached. Moreover, it also aborts when no more intervals can be
   found for splitting.

   For finding splitting points the following rules are used (in
   this order, i.e., is if the first rule cannot be applied, the next
   one is used):
   @enumerate
   @item
     Use the expected value of adaptive rejection sampling.
   @item
     Use the arc-mean rule (a mixture of arithmetic mean and harmonic
     mean).
   @item
     Use the arithmetic mean of the interval boundaries.
   @end enumerate
   Notice, however, that for unbounded intervals neither rule 1 nor rule
   3 can be used.

   As an additional feature, it is possible to choose amoung these
   rules. 
   If @var{usedars} is set to @code{1} or TRUE the expected point
   (rule 1) is used (it switches to rule 2 for a particular 
   interval if rule 1 cannot be applied).
   If it is set to @code{2} the arc-mean rule is used.
   If it is set to @code{3} the mean is used.
   Notice that rule 3 can only be used if the domain of the
   distribution is bounded. It is faster than the other two methods
   but for heavy-tailed distribution and large domain the hat
   converges extremely slowly.

   The default depends on the given construction points.
   If the user has provided such points via a unur_tdr_set_cpoints()
   call, then @var{usedars} is set to FALSE by default, i.e.,
   there is no further splitting.
   If the user has only given the number of construction points (or
   only uses the default number), then @var{usedars} is set to TRUE
   (i.e., use rule 1).
*/

int unur_tdr_set_darsfactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for ``derandomized adaptive rejection sampling''.
   This factor is used to determine the intervals that are ``too
   large'', that is, all intervals where the area between squeeze and
   hat is larger than @var{factor} times the average area over all
   intervals between squeeze and hat.
   Notice that all intervals are split when @var{factor} is set to
   @code{0.}, and that there is no splitting at all when @var{factor}
   is set to @code{UNUR_INFINITY}.

   Default is @code{0.99}. There is no need to change this parameter.
*/

int unur_tdr_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
/* 
   Set construction points for the hat function. If @var{stp} is NULL
   than a heuristic rule of thumb is used to get @var{n_stp}
   construction points. This is the default behavior. 

   The default number of construction points is 30.
*/

int unur_tdr_set_reinit_percentiles( UNUR_PAR *parameters, int n_percentiles, const double *percentiles );
/* */ 

int unur_tdr_chg_reinit_percentiles( UNUR_GEN *generator, int n_percentiles, const double *percentiles );
/* 
   By default, when the @var{generator} object is reinitialized, it
   used the same construction points as for the initialization
   procedure.
   Often the underlying distribution object has been changed only
   moderately. For example, the full conditional distribution of a
   multivariate distribution. 
   In this case it might be more appropriate to use
   percentilesm of the hat function for the last (unchanged)
   distribution. @var{percentiles} must then be a pointer to an
   ordered array of numbers between @code{0.01} and @code{0.99}.
   If @var{percentiles} is NULL, then a heuristic rule of thumb is
   used to get @var{n_percentiles} values for these percentiles.
   Notice that @var{n_percentiles} must be at least @code{2},
   otherwise defaults are used.
   (Then the first and third quartiles are used by default.) 
*/

int unur_tdr_set_reinit_ncpoints( UNUR_PAR *parameters, int ncpoints );
/* */ 

int unur_tdr_chg_reinit_ncpoints( UNUR_GEN *generator, int ncpoints );
/* 
   When reinit fails with the given construction points or the percentiles 
   of the old hat function, another trial is undertaken with @var{ncpoints}
   construction points. @var{ncpoints} must be at least @code{10}.

   Default: @code{50}
 */

int unur_tdr_chg_truncated(UNUR_GEN *gen, double left, double right);
/*
   Change the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call. The hat function will not be changed and
   there is no need to run unur_reinit().

   @emph{Important:}
   The ratio between the area below the hat and the area below the
   squeeze changes when the sampling region is restricted. In particular
   it becomes (very) large when sampling from the (far) tail of the
   distribution. Then it is better to create a generator object for the
   tail of distribution only.

   @emph{Important:}
   This call does not work for variant @code{IA} (immediate
   acceptance). In this case UNU.RAN switches @emph{automatically} to 
   variant @code{PS}.

   @emph{Important:}
   It is not a good idea to use adaptave rejection sampling while 
   sampling from a domain that is a strict subset of the domain that
   has been used to construct the hat.
   For that reason adaptive adding of construction points is
   @emph{automatically disabled} by this call.

   @emph{Important:} If the CDF of the hat is (almost) the same 
   for @var{left} and @var{right} and (almost) equal to @code{0} or
   @code{1}, then the truncated domain is not changed and the call
   returns an error code.
*/


int unur_tdr_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area below squeeze) / (area below hat).
   It must be a number between 0 and 1.
   When the ratio exceeds the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use 0 if no construction points should be added after the setup.
   Use 1 if added new construction points should not be stopped until
   the maximum number of construction points is reached.

   Default is @code{0.99}.
*/

double unur_tdr_get_sqhratio( const UNUR_GEN *generator );
/* 
   Get the current ratio (area below squeeze) / (area below hat)
   for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_tdr_get_hatarea( const UNUR_GEN *generator );
/* 
   Get the area below the hat for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_tdr_get_squeezearea( const UNUR_GEN *generator );
/* 
   Get the area below the squeeze for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

int unur_tdr_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals.
   No construction points are added after the setup when the number of
   intervals suceeds @var{max_ivs}. 
   It is increased automatically to twice the number of construction
   points if this is larger.

   Default is @code{100}.
*/

int unur_tdr_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* 
   Use the center as construction point. Default is TRUE.
*/

int unur_tdr_set_usemode( UNUR_PAR *parameters, int usemode );
/* 
   Use the (exact!) mode as construction point.
   Notice that the behavior of the algorithm is different to simply
   adding the mode in the list of construction points via a
   unur_tdr_set_cpoints() call. In the latter case the mode is treated
   just like any other point. However, when @code{usemode} is TRUE, the
   tangent in the mode is always set to 0. Then the hat of the
   transformed density can never cut the x-axis which must never
   happen if c < 0, since otherwise the hat would not be bounded.

   Default is TRUE.
*/

int unur_tdr_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is 2.
*/

int unur_tdr_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_tdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_tdr_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* 
   Sometimes it might happen that unur_init() has been executed
   successfully. But when additional construction points are added by
   adaptive rejection sampling, the algorithm detects that the
   PDF is not T-concave. 

   With @var{pedantic} being TRUE, the
   sampling routine is exchanged by a routine that simply returns
   @code{UNUR_INFINITY}. Otherwise the new point is not added to the
   list of construction points. At least the hat function remains
   T-concave.

   Setting @var{pedantic} to FALSE allows sampling from a
   distribution which is ``almost'' T-concave and small errors are
   tolerated. However it might happen that the hat function cannot be
   improved significantly. When the hat functions that has been
   constructed by the unur_init() call is extremely large then it
   might happen that the generation times are extremely high
   (even hours are possible in extremely rare cases).

   Default is FALSE.
*/

double unur_tdr_eval_invcdfhat( const UNUR_GEN *generator, double u, 
				double *hx, double *fx, double *sqx );
/* 
   Evaluate the inverse of the CDF of the hat distribution at @var{u}.
   As a side effect the values of the hat, the density, and the squeeze
   at the computed point @i{x} are stored in @var{hx}, @var{fx}, and
   @var{sqx}, respectively. However, these computations are suppressed
   if the corresponding variable is set to NULL.

   If @var{u} is out of the domain [0,1] then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).

   @emph{Important:}
   This call does not work for variant @code{IA} (immediate
   acceptance). In this case the hat CDF is evaluated as if
   variant @code{PS} is used.

   @emph{Notice}: This function always evaluates the inverse CDF of
   the hat distribution. A call to unur_tdr_chg_truncated() call
   has no effect.
*/


/* =END */
/*---------------------------------------------------------------------------*/
/* Internal routines                                                         */

int _unur_tdr_is_ARS_running( const UNUR_GEN *generator );
/* 
   Check whether more points will be added by adaptive rejection sampling.
*/

/*---------------------------------------------------------------------------*/
