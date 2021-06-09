/*****************************************************************************
 *                                                                           *
 *          UNU.RAN -- Universal Non-Uniform Random number generator         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: pinv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method PINV                               *
 *         (Polynomial interpolation based INVersion of CDF)                 *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
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
 *****************************************************************************/

/* 
   =METHOD PINV  Polynomial interpolation based INVersion of CDF

   =UP  Methods_for_CONT

   =REQUIRED PDF

   =OPTIONAL domain, center, CDF, derivative of PDF

   =REF [DHLa08]

   =SPEED Set-up: (very) slow, Sampling: (very) fast

   =REINIT not implemented

   =DESCRIPTION
      PINV is a variant of numerical inversion, where the inverse CDF
      is approximated using Newton's interpolating formula. 
      The interval [0,1] is split into several subintervals. In each
      of these the inverse CDF is constructed at nodes 
      @unurmath{(CDF(x),x)} for some points @i{x} in this subinterval.
      If the PDF is given, then the CDF is computed numerically
      from the given PDF using adaptive Gauss-Lobatto
      integration with 5 points. Subintervals are split until the
      requested accuracy goal is reached.

      The method is not exact, as it only produces random variates of
      the approximated distribution. Nevertheless, the maximal
      tolerated approximation error can be set to be the resolution
      (but of course is bounded by the machine precision).
      We use the u-error @unurmath{|U-CDF(X)|} to measure the error
      for @i{X} = "approximate inverse CDF"(@i{U}).
      Notice that very small values of the u-resolution are possible
      but increase the cost for the setup step.
      We call the maximal tolerated u-error the @emph{u-resolution} of
      the algorithm in the sequel.

      Both the order of the interpolating polynomial and the
      u-resolution can be selected.

      The interpolating polynomials have to be computed in a setup
      step. However, it only works for distributions with bounded
      domain; for distributions with unbounded domain the tails are
      cut off such that the probability for the tail regions is
      small compared to the given u-resolution. 

      The construction of the interpolation polynomial only works when
      the PDF is unimodal or when the PDF does not vanish between two
      modes. 

      There are some restrictions for the given distribution:
      @itemize
      @item 
      The support of the distribution (i.e., the region where the PDF
      is strictly positive) must be connected. In practice this means,
      that the region where PDF is "not too small" must be connected.
      Unimodal densities satisfy this condition.
      If this condition is violated then the domain of the
      distribution might be truncated.
      @item
      When the PDF is integrated numerically, then the given PDF must
      be continuous and should be smooth.
      @item
      The PDF must be bounded.
      @item
      The algorithm has problems when the distribution has heavy tails
      (as then the inverse CDF becomes very steep at 0 or 1)
      and the requested u-resolution is very small.
      E.g., the Cauchy distribution is likely to show this problem
      when the requested u-resolution is less then @code{1.e-12}.
      @end itemize
      Regions with very small PDF values or heavy tails might lead to
      an abortion of the set-up or (even worse) the approximation
      error might become larger than requested, since the (computation of the)
      interpolating polynomial becomes numerically unstable.
      
      @emph{Remark:}
      We also have implemented experimental variants.
      However, we observed that these variants are more sensitive to
      round-off errors, especially in the right hand tail and we
      @emph{do not recommend} their usage unless there are severe
      reasons.

      @itemize @minus
      @item
      Use a function that implements the CDF instead of numerical
      integration of the PDF. 
      
      @item
      Use Hermite interpolation instead of Newton interpolation.
      Thus the first (and second) derivative of the interpolating
      polynomial coincides with that of the inverse CDF.
      Consequently the interpolant is also (twice) differentiable even
      at the interval boundaries.
      This variant can be seen as limiting case of Newton
      interpolation with double (or triple) points as nodes.

      We have used a @emph{smoothness} parameter to control this
      feature. However, besides numerical problems we observed that
      this variant requires more intervals and thus larger setup times
      and higher memory consumptions.
      @end itemize


   =HOWTOUSE
      PINV works for continuous univariate distribution objects with
      given PDF. The corresponding distribution object should contain a
      typical point of the distribution, i.e., a point where the PDF
      is not too small, e.g., (a point near) the mode.
      However, it is important that the center is @strong{not} the
      pole of the distribution (or a point too close to the pole).
      It can be set using a unur_distr_cont_set_center() or
      a unur_distr_cont_set_mode() call. If neither is set, or if the 
      given center cannot be used, then a simple search routine tries
      to find an appropriate point for the center. 

      It is recommended that the domain of the distribution with
      bounded domain is specified using a unur_distr_cont_set_domain()
      call. Otherwise, the boundary is searched numerically which
      might be rather expensive, especially when this boundary point
      is @code{0}.
      
      When sampling from truncated distributions with extreme
      truncation points, it is recommended to provide the log-density 
      using unur_distr_cont_set_logpdf() and the mode.
      Then the PDF is rescaled such that the PDF at the mode is 1.
      Thus the algorithm is numerically more stable.

      The inverse CDF is interpolated using Newton polynomials.
      The order of this polynomial can be set by means of a
      unur_pinv_set_order() call.
      
      The smoothness of the interpolant at interval boundaries can be
      controlled using a unur_pinv_set_smoothness() call. 
      Then Hermite interpolation instead of Newton interpolation is
      used. (The former can be seen as a limiting case of Newton
      interpolation with double (or triple) points.)
      However, using higher smoothness is @emph{not recommended}
      unless differentiability at the interval boundaries is
      important.

      For distributions with unbounded domains the tails are cut
      off such that the probability for the tail regions is small
      compared to the given u-resolution. For finding these cut points
      the algorithm starts with the region @code{[-1.e100,1.e100]}. For
      the exceptional case where this does not work these starting
      points can be changed via a unur_pinv_set_boundary() call.

      This method is not exact, as it only produces random variates of 
      the approximated distribution. Nevertheless, the numerical error
      in "u-direction" (i.e., |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be controlled
      by means of unur_pinv_set_u_resolution().
      However, the maximal error of this approximation is only
      estimated. For very small u-resolutions the actual approximation
      error might be (slightly) larger than the requested u-resolution.
      (Of course the size of this value depends on the given PDF.)
      If this error is crucial for an application we recommend to
      compute this error using unur_pinv_estimate_error() which runs a
      small Monte Carlo simulation.
      See also the documentation for function
      unur_pinv_set_u_resolution() and the remark given there.

      The number of required subintervals heavily depends on the order
      of the interpolating polynomial and the requested u-resolution:
      it increases when order or u-resolution are decreased.
      It can be checked using a unur_pinv_get_n_intervals() call.
      The maximum number of such subintervals is fixed but can be
      increased using a unur_pinv_set_max_intervals() call.
      If this maximum number is too small then the set-up aborts with
      a corresponding error message.

      It is also possible to use the CDF of the distribution instead
      of the PDF. Then the distribution object must contain a pointer
      to the CDF. Moreover, this variant of the algorithm has to be
      switched on using an unur_pinv_set_usecdf() call.
      Notice, however, that the setup for this variant is numerically
      less stable than using integration of the PDF (the default
      variant). Thus using the CDF is @emph{not recommended}.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_pinv_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_pinv_set_order( UNUR_PAR *parameters, int order);
/* 
   Set order of interpolation. Valid orders are between @code{3} and
   @code{17}. Higher orders result in fewer intervals for the
   approximations. 

   Default: @code{5}.
*/

int unur_pinv_set_smoothness( UNUR_PAR *parameters, int smoothness);
/* 
   Set smoothness of interpolant. By construction the interpolant is
   piecewise polynomial and thus smooth on each of the intervals
   where these polynomials are constructed. At the interval
   boundaries, however, it usually not be differentiable.
   Method PINV also implements variants of Newton interpolation where
   the first (or second) derivative of the interpolating
   polynomial coincides with the respective derivative of the inverse
   CDF at the nodes. The the interpolant is (twice) differentiable
   even at the interval boundaries.
   These variants can be seen as limiting case of Newton interpolation
   with double (or triple) points as nodes and are known as Hermite
   interpolation.

   Possible values for @var{smoothness}:

   @multitable @columnfractions .1 .25 .60
   @headitem Value @tab Effect @tab Requirements
   @item @code{0} 
   @tab continuous
   @tab requires PDF (or CDF)

   @item @code{1} 
   @tab differentiable
   @tab requires PDF (optional: CDF), @*
   order of polynomial must be odd

   @item @code{2} 
   @tab twice differentiable
   @tab requires PDF and its derivative (optional: CDF), @*
   order must be 5, 8, 11, 14 or 17
   @end multitable

   If the order of the polynomial does not satisfy the given
   condition, then it is increased to the next larger possible value.

   @emph{Remark:}
   A higher smoothness parameter usually results in a higher number of
   intervals and thus a higher setup time and memory consumption.
   We also observed that higher smoothness parameters make the
   algorithm more sensible for round-off error. Then the setup fails.

   @emph{Remark:}
   If the interpolating polynomial cannot be constructed for the
   requested smoothness on a particular interval, 
   then the smoothness parameter is reduced for that interval.

   @emph{Remark:}
   For order @code{3} and smoothness @code{1} (cubic Hermite
   interpolation) monotonicity is guaranteed by checking a simple
   monotonicity condition for the coefficients of the polynomials.

   @emph{Remark:}
   Using @var{smoothness} larger than @code{0} is 
   @emph{not recommended} unless differentiability at the interval
   boundaries is important for ones application.

   Default: @code{0}.
*/

int unur_pinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* 
   Set maximal tolerated u-error. Values of @var{u_resolution} must
   at least @code{1.e-15} and @code{1.e-5} at most.
   Notice that the resolution of most uniform random number sources is
   @unurmath{2^{-32}} = @code{2.3e-10}. Thus a value of @code{1.e-10} 
   leads to an inversion algorithm that could be called exact. For most
   simulations slightly bigger values for the maximal error are enough
   as well. 

   Smaller values for @var{u_resolution} increase the number of
   subinterval that are necessary for the approximation of the inverse
   CDF. For very small values (less then @code{1.e-12}) this number
   might exceed the maximum number of such intervals. However, this
   number can be increased using a unur_pinv_set_max_intervals() call.

   @emph{Remark:}
   We ran many experiments and found that the observed u-error was
   always smaller than the given @var{u_resolution} whenever this
   value was @code{1.e-12}. For values smaller than @code{1e-13} the
   maximal observed u-error was slightly larger. One use @code{1.e-15}
   if best approximation is required. However, then the actual u-error
   can be as large as @code{1.e-14}.

   @strong{Warning!}
   These figures are based on our experiments (with some tolerance
   added to be on the safe side). There is no guarantee for these error
   estimates for a particular distribution.

   Default is @code{1.e-10}.
*/

int unur_pinv_set_use_upoints( UNUR_PAR *parameters, int use_upoints );
/* 
   If @var{use_upoints} is TRUE, then the nodes of the interpolating
   polynomial are constructed by means of Chebyshev points in u-scale
   not in x-scale. This results is a better approximation but almost
   doubles the number of PDF or CDF evaluations during the setup.
   (This is an experimental feature.)
   
   Default: FALSE
*/

int unur_pinv_set_usepdf( UNUR_PAR *parameters );
/* 
   Use PDF (if available) to compute approximate inverse CDF.

   This is the default.
*/

int unur_pinv_set_usecdf( UNUR_PAR *parameters );
/* 
   Use CDF (if available) to compute approximate inverse CDF.
   This variant is intend for running experiments with method PINV.

   @emph{Remark:}
   We ran many experiments and found that for small values of the
   given @var{u_resolution} (less than @code{1.e-12}) the setup fails
   for distributions with heavy tails. We found that using the PDF
   (instead of the CDF) is numerically more stable.
   This is especially the case when the smoothness parameter is set
   to @code{1} or @code{2}.

   Using the CDF is @strong{not recommended}.
*/

int unur_pinv_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set @var{left} and @var{right} point for finding the cut-off points
   for the "computational domain", i.e., the domain that covers the
   essential part of the distribution.
   The cut-off points are computed such that the tail probabilities
   are smaller than given by unur_pinv_set_u_resolution().
   It is usually safe to use a large interval.
   However, @code{+/- UNUR_INFINITY} is not allowed.

   @emph{Important}: This call does not change the domain of the
   given distribution itself. But it restricts the domain for the
   resulting random variates.

   Default: intersection of @code{[-1.e100,+1.e100]} and the given
   domain of the distribution.
*/

int unur_pinv_set_searchboundary( UNUR_PAR *parameters, int left, int right );
/* 
   If @var{left} or @var{right} is set to FALSE then the respective
   boundary as given by a unur_pinv_set_boundary() call is used
   without any further computations.
   However, these boundary points might cause numerical problems
   during the setup when PDF returns @code{0} ``almost everywhere''.
   If set to TRUE (the default) then the computational interval is
   shortened to a more sensible region by means of a search algorithm.
   Switching off this search is useful, e.g., for the Gamma(2)
   distribution where the left border @code{0} is fixed and finite.

   @emph{Remark:}
   The searching algorithm assumes that the support of the distribution
   is connected.

   @emph{Remark:}
   Do not set this parameter to FALSE except when searching for
   cut-off points fails and one wants to try with precomputed values.

   Default: TRUE.
*/

int unur_pinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals. @var{max_ivs} must be at least
   @code{100} and at most @code{1000000}.

   Default is @code{10000}.
*/

int unur_pinv_get_n_intervals( const UNUR_GEN *generator ); 
/* 
   Get number of intervals used for interpolation in 
   the generator object.
   It returns @code{0} in case of an error.
*/

int unur_pinv_set_keepcdf( UNUR_PAR *parameters, int keepcdf);
/* 
   If the PDF is given, then the CDF is computed numerically
   from the given PDF using adaptive Gauss-Lobatto integration.
   Thus a table of CDF points is stored to keep the number of 
   evaluations of the PDF minimal. Usually this table is discarded
   when the setup is completed.
   If @var{keepcdf} is TRUE, then this table is kept and can be used
   to compute the CDF of the underlying distribution by means of 
   function unur_pinv_eval_approxcdf().
   This option is ignored when unur_pinv_set_usecdf() is called.
   
   Default: FALSE
*/

double unur_pinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
/*
   Evaluate interpolation of inverse CDF at @var{u}.
   If @var{u} is out of the domain (0,1) then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).
*/

double unur_pinv_eval_approxcdf( const UNUR_GEN *generator, double x );
/*
   Evaluate (approximate) CDF at @var{x}. If the PDF of the
   distribution is given, then adaptive Gauss-Lobatto integration is
   used to compute the CDF.
   If the PDF is used to create the generator object, then the
   table of integral values must not removed at the end of setup and thus 
   unur_pinv_set_keepcdf() must be called.
*/

int unur_pinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
/*
   Estimate maximal u-error and mean absolute error (MAE) for @var{generator}
   by means of Monte-Carlo simulation with sample size @var{samplesize}.
   The results are stored in @var{max_error} and @var{MAE}, respectively.

   It returns @code{UNUR_SUCCESS} if successful. 
*/

/* =END */
/*---------------------------------------------------------------------------*/

/*
   Remark:
   PINV needs a first guess for the area below the PDF. This value can
   be set using the unur_distr_cont_set_pdfarea() call. Otherwise
   @code{1} is used.
   There is no necessity to set this area unless it differs from 1 by
   several orders of magnitude. 
*/

/*---------------------------------------------------------------------------*/
