/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hinv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HINV                               *
 *         (Hermite interpolation based INVersion of CDF)                    *
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
   =METHOD HINV  Hermite interpolation based INVersion of CDF

   =UP  Methods_for_CONT

   =REQUIRED CDF

   =OPTIONAL PDF, dPDF

   =REF [HLa03] [HLD04: Sect.7.2, Alg.7.1]

   =SPEED Set-up: (very) slow, Sampling: (very) fast

   =REINIT supported

   =DESCRIPTION
      HINV is a variant of numerical inversion, where the inverse CDF
      is approximated using Hermite interpolation, i.e., the interval 
      [0,1] is split into several intervals and in each interval the
      inverse CDF is approximated by polynomials constructed by means
      of values of the CDF and PDF at interval boundaries. This makes
      it possible to improve the accuracy by splitting a particular
      interval without recomputations in unaffected intervals. Three
      types of splines are implemented: linear, cubic, and quintic
      interpolation. For linear interpolation only the CDF is
      required. Cubic interpolation also requires PDF and quintic
      interpolation PDF and its derivative. 

      These splines have to be computed in a setup step. However, it
      only works for distributions with bounded domain; for
      distributions with unbounded domain the tails are chopped off
      such that the probability for the tail regions is small compared
      to the given u-resolution. 

      The method is not exact, as it only produces random variates of
      the approximated distribution. Nevertheless, the maximal
      numerical error in "u-direction" (i.e. |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be set to the
      required resolution (within machine precision).  
      Notice that very small values of the u-resolution are possible
      but may increase the cost for the setup step.

      As the possible maximal error is only estimated in the setup it
      may be necessary to set some special design points for computing
      the Hermite interpolation to guarantee that the maximal u-error
      can not be bigger than desired. Such points are points where the
      density is not differentiable or has a local extremum. Notice
      that there is no necessity to do so. However, if you do not
      provide these points to the algorithm there might be a small
      chance that the approximation error is larger than the given
      u-resolution, or that the required number of intervals is larger
      than necessary.

   =HOWTOUSE
      HINV works for continuous univariate distribution objects with
      given CDF and (optional) PDF. It uses Hermite interpolation of
      order 1, 3 [default] or 5. The order can be set by means of
      unur_hinv_set_order().

      For distributions with unbounded domains the tails are chopped
      off such that the probability for the tail regions is small
      compared to the given u-resulution. For finding these cut points
      the algorithm starts with the region @code{[-1.e20,1.e20]}. For
      the exceptional case where this might be too small (or one knows
      this region and wants to avoid this search heuristics) it can be
      directly set via a unur_hinv_set_boundary() call.

      It is possible to use this method for generating from truncated
      distributions. It even can be changed for an existing generator
      object by an unur_hinv_chg_truncated() call.

      This method is not exact, as it only produces random variates of 
      the approximated distribution. Nevertheless, the numerical error
      in "u-direction" (i.e. |U-CDF(X)|, for 
      X = "approximate inverse CDF"(U) |U-CDF(X)|) can be controlled
      by means of unur_hinv_set_u_resolution().

      The possible maximal error is only estimated in the setup. Thus
      it might be necessary to set some special design points for
      computing the Hermite interpolation to guarantee that the
      maximal u-error can not be bigger than desired. Such points
      (e.g. extremal points of the PDF, points with infinite
      derivative) can be set using using the
      unur_hinv_set_cpoints() call. 
      If the mode for a unimodal distribution is set in the distribution
      object this mode is automatically used as design-point if the
      unur_hinv_set_cpoints() call is not used.

      As already mentioned the maximal error of this approximation is 
      only estimated. If this error is crucial for an application we
      recommend to compute this error using unur_hinv_estimate_error()
      which runs a small Monte Carlo simulation.
      
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.
      The values given by the last unur_hinv_chg_truncated() call will be 
      then changed to the values of the domain of the underlying distribution
      object. Moreover, starting construction points (nodes) that are given by
      a unur_hinv_set_cpoints() call are ignored when unur_reinit() is
      called.
      It is important to note that for a distribution from the 
      UNU.RAN library of standard distributions
      (@pxref{Stddist,,Standard distributions})
      the normalization constant has to be updated using the 
      unur_distr_cont_upd_pdfarea() call whenever its parameters have been
      changed by means of a unur_distr_cont_set_pdfparams() call.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_hinv_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hinv_set_order( UNUR_PAR *parameters, int order);
/* 
   Set order of Hermite interpolation. Valid orders are
   @code{1}, @code{3}, and @code{5}.
   Notice that @var{order} greater than @code{1} requires the density 
   of the distribution, and @var{order} greater than @code{3} even
   requires the derivative of the density. Using @var{order} @code{1}
   results for most distributions in a huge number of intervals
   and is therefore not recommended. If the maximal error in
   u-direction is very small (say smaller than @code{1.e-10}),
   @var{order} @code{5} is recommended as it leads to considerably 
   fewer design points, as long there are no poles or heavy tails.

   @emph{Remark:} When the target distribution has poles or (very) heavy
   tails @var{order} @code{5} (i.e., quintic interpolation) is 
   numerically less stable and more sensitive to round-off errors than
   @var{order} @code{3} (i.e., cubic interpolation).

   Default is @code{3} if the density is given and @code{1} otherwise.
*/

int unur_hinv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* 
   Set maximal error in u-direction. However, the given u-error must not
   be smaller than machine epsilon (@code{DBL_EPSILON}) and should not be
   too close to this value. As the resolution of most uniform random
   number sources is 2^(-32) = @code{2.3e-10}, a value of @code{1.e-10}
   leads to an inversion algorithm that could be called exact. For most
   simulations slightly bigger values for the maximal error are enough
   as well. 

   Remark: The u-error might become larger than @var{u_resolution} due
   to rescaling of floating point numbers when the domain of the
   distribution is truncated by a unur_hinv_chg_truncated() call.

   Default is @code{1.e-10}.
*/

int unur_hinv_set_cpoints( UNUR_PAR *parameters, const double *stp, int n_stp );
/* 
   Set starting construction points (nodes) for Hermite interpolation. 
      
   As the possible maximal error is only estimated in the setup
   it may be necessary to set some special design points for
   computing the Hermite interpolation to guarantee that the
   maximal u-error can not be bigger than desired. We suggest to 
   include as special design points all local extrema of the density,
   all points where the density is not differentiable, and isolated
   points inside of the domain with density 0. 
   If there is an interval with density constant equal to 0 inside of
   the given domain of the density, both endpoints of this interval 
   should be included as special design points. Notice that there is no
   necessity to do so. However, if these points are not provided to
   the algorithm the approximation error might be larger than the
   given u-resolution, or the required number of intervals could be
   larger than necessary.

   @emph{Important}: Notice that the given points must be in
   increasing order and they must be disjoint. 

   @emph{Important}: The boundary point of the computational region
   must not be given in this list!
   Points outside the boundary of the computational region are ignored.

   Default is for unimodal densities - if known - the mode of the 
   density, if it is not equal to the border of the domain. 
*/

int unur_hinv_set_boundary( UNUR_PAR *parameters, double left, double right );
/* 
   Set the left and right boundary of the computational interval.
   Of course @code{+/- UNUR_INFINITY} is not allowed.
   If the CDF at @var{left} and @var{right} is not close to the
   respective values @code{0.} and @code{1.} then this interval is
   increased by a (rather slow) search algorithm.

   @emph{Important}: This call does not change the domain of the
   given distribution itself. But it restricts the domain for the
   resulting random variates.

   Default is @code{1.e20}.
*/

int unur_hinv_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is @code{1}.
*/

int unur_hinv_set_max_intervals( UNUR_PAR *parameters, int max_ivs );
/* 
   Set maximum number of intervals. No generator object is created if
   the necessary number of intervals for the Hermite interpolation 
   exceeds @var{max_ivs}. It is used to prevent the algorithm to eat up
   all memory for very badly shaped CDFs.

   Default is @code{1000000} (1.e6).
*/

int unur_hinv_get_n_intervals( const UNUR_GEN *generator );
/* 
   Get number of nodes (design points) used for Hermite interpolation in 
   the generator object. The number of intervals is the number of
   nodes minus 1.
   It returns an error code in case of an error.
*/

double unur_hinv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
/*
   Evaluate Hermite interpolation of inverse CDF at @var{u}.
   If @var{u} is out of the domain [0,1] then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).

   @emph{Notice}: When the domain has been truncated by a  
   unur_hinv_chg_truncated() call then the inverse CDF of the
   truncated distribution is returned.
*/

int unur_hinv_chg_truncated( UNUR_GEN *generator, double left, double right );
/*
   Changes the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call. The tables of splines are not recomputed.
   Thus it might happen that the relative error for the generated variates
   from the truncated distribution is greater than the bound for the
   non-truncated distribution. This call also fails when the CDF values
   of the boundary points are too close, i.e. when only a few different
   floating point numbers would be computed due to round-off errors
   with floating point arithmetic.

   Remark: The u-error might become larger than the @var{u_resolution}
   given by a unur_hinv_set_u_resolution() call due to rescaling of
   floating point numbers when the domain of the distribution is
   truncated.

   When failed an error code is returned.

   @emph{Important}: Always check the return code since the domain is
   not changed in case of an error.
*/

int unur_hinv_estimate_error( const UNUR_GEN *generator, int samplesize, double *max_error, double *MAE );
/*
   Estimate maximal u-error and mean absolute error (MAE) for
   @var{generator} by means of a (quasi-) Monte-Carlo simulation with
   sample size @var{samplesize}.
   The results are stored in @var{max_error} and @var{MAE}, respectively.

   It returns @code{UNUR_SUCCESS} if successful. 
*/

/* =END */
/*---------------------------------------------------------------------------*/
