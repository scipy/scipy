/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: discr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  DISCR  (discrete univariate distribution)                   *
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

/*---------------------------------------------------------------------------*/

/* 
   =NODEX   DISCR   Discrete univariate distributions

   =UP Distribution_objects [50]

   =DESCRIPTION
      The calls in this section can be applied to discrete
      univariate distributions.

      @itemize @minus 
      @item Create a @command{new} instance of a discrete univariate
      distribution.

      @item Handle and evaluate 
      distribution function (CDF, @command{cdf}) and 
      probability mass function (PMF, @command{pmf}).
      The following is important:
      @itemize .
      @item @command{pmf} need not be normalized, i.e.,
      any summable nonnegative function on the set of intergers can be
      used. 
      @item @command{cdf} must be a distribution function, i.e. it
      must be monotonically increasing with range [0,1].
      @item If @command{cdf} and @command{pdf} are used together for a
      pariticular generation method, then @command{pmf} must be
      normalized, i.e. it must sum to 1.
      @end itemize

      @item Alternatively, @command{cdf} and @command{pdf} can be
      provided as @command{str}ings instead of function pointers.

      @item Some generation methods require a (finite) probability
      vector (PV, @command{pv}), i.e. an array of @code{double}s.
      It can be automatically computed if the @command{pmf} is
      given but @command{pv} is not.

      @item Set (and change) parameters (@command{pmfparams}) and the
      total sum (@command{pmfsum}) of the given PMF or PV.

      @item Set the @command{mode} of the distribution. 

      @item Set the @command{domain} of the distribution. 

      @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate discrete distributions (DISCR).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_discr_new( void );
/* 
   Create a new (empty) object for a univariate discrete distribution.
*/


/* ==DOC
   @subsubheading Essential parameters

   There are two interfaces for discrete univariate distributions:
   Either provide a (finite) probability vector (PV).
   Or provide a probability mass function (PMF). For the latter
   case there are also a couple of derived parameters that are not
   required when a PV is given.

   It is not possible to set both a PMF and a PV directly. However, the
   PV can be computed from the PMF (or the CDF if no PMF is available)
   by means of a unur_distr_discr_make_pv() call.
   If both the PV and the PMF are given in the distribution object it
   depends on the generation method which of these is used.
*/

int unur_distr_discr_set_pv( UNUR_DISTR *distribution, const double *pv, int n_pv );
/* 
   Set finite probability vector (PV) for the @var{distribution}. It is not
   necessary that the entries in the given PV sum to 1.
   @var{n_pv} must be positive. However, there is no testing
   whether all entries in @var{pv} are non-negative. 

   If no domain has been set, then the left boundary is set to
   @code{0}, by default. If @var{n_pv} is too large, e.g. because
   left boundary + @var{n_pv} exceeds the range of integers, 
   then the call fails. 

   Notice that it is not possible to set both a PV and a PMF or CDF.
   If the PMF or CDF is set first one cannot set the PV.
   If the PMF or CDF is set first after a PV is set, the latter is 
   removed (and recomputed using unur_distr_discr_make_pv() when required).
*/

int unur_distr_discr_make_pv( UNUR_DISTR *distribution );
/* 
   Compute a PV when a PMF or CDF is given. However, when the
   domain is not given or is too large and the sum over the PMF is given
   then the (right) tail of the @var{distribution} is chopped off such that
   the probability for the tail region is less than 1.e-8.
   If the sum over the PMF is not given a PV of maximal length is
   computed.

   The maximal size of the created PV is bounded by the macro
   @code{UNUR_MAX_AUTO_PV} that is defined in @file{unuran_config.h}.

   If successful, the length of the generated PV is returned.
   If the sum over the PMF on the chopped tail is not neglible small
   (i.e. greater than 1.e-8 or unknown) than the 
   negative of the length of the PV is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.

   Notice that the left boundary of the PV is set to @code{0} by
   default when a discrete distribution object is created from
   scratch.

   If computing a PV fails for some reasons, an error code is returned and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
*/

int unur_distr_discr_get_pv( const UNUR_DISTR *distribution, const double **pv );
/* 
   Get length of PV of the @var{distribution} and set pointer
   @var{pv} to array of probabilities. If no PV is given,
   an error code is returned and @var{pv} is set to NULL.@*
   (It does not call unur_distr_discr_make_pv()!)
*/

int unur_distr_discr_set_pmf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *pmf );
/* */

int unur_distr_discr_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_DISCR *cdf );
/* 
   Set respective pointer to the PMF and the CDF of the @var{distribution}.
   These functions must be of type
   @code{double funct(int k, const UNUR_DISTR *distr)}.

   It is important to note that all these functions must return a
   result for all integers @var{k}. E.g., if the domain of a given
   PMF is the interval @{1,2,3,@dots{},100@}, than the given function
   must return @code{0.0} for all points outside this interval.

   The default domain for the PMF or CDF is [@code{0}, @code{INT_MAX}].
   The domain can be changed using a unur_distr_discr_set_domain() call.

   It is not possible to change such a function. Once the PMF or
   CDF is set it cannot be overwritten. A new distribution object
   has to be used instead.

   Notice that it is not possible to set both a PV and a PMF or CDF.
   If the PMF or CDF is set first one cannot set the PV.
   If the PMF or CDF is set first after a PV is set, the latter is 
   removed (and recomputed using unur_distr_discr_make_pv() when required).
*/

int unur_distr_discr_set_invcdf( UNUR_DISTR *distribution, UNUR_IFUNCT_DISCR *invcdf );
/* 
   Set inverse CDF of the @var{distribution}.
   @var{invcdf} must be a pointer must be of type
   @code{int funct(double x, const UNUR_DISTR *distr)},
   i.e., it should return a @code{double}.
*/


UNUR_FUNCT_DISCR *unur_distr_discr_get_pmf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_DISCR *unur_distr_discr_get_cdf( const UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the PMF and the CDF of the 
   @var{distribution}. The pointer is of type
   @code{double funct(int k, const UNUR_DISTR *distr)}.
   If the corresponding function is not available for the @var{distribution},
   the NULL pointer is returned.
*/

UNUR_IFUNCT_DISCR *unur_distr_discr_get_invcdf( const UNUR_DISTR *distribution );
/* 
   Get pointer to the inverse CDF of the @var{distribution}. 
   The pointer is of type
   @code{int funct(double x, const UNUR_DISTR *distr)}.
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/


double unur_distr_discr_eval_pv(int k, const UNUR_DISTR *distribution );
/* */

double unur_distr_discr_eval_pmf( int k, const UNUR_DISTR *distribution );
/* */

double unur_distr_discr_eval_cdf( int k, const UNUR_DISTR *distribution );
/* 
   Evaluate the PV, PMF, and the CDF, respectively, at k.
   Notice that @var{distribution} must not be the NULL pointer.
   If no PV is set for the @var{distribution}, then
   unur_distr_discr_eval_pv() behaves like unur_distr_discr_eval_pmf().
   If the corresponding function is not available for the @var{distribution},
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.

   @emph{IMPORTANT:}
   In the case of a truncated standard distribution these calls always
   return the respective values of the @emph{untruncated} distribution!
*/

int unur_distr_discr_eval_invcdf( double u, const UNUR_DISTR *distribution );
/* 
   Evaluate the inverse CDF at @var{u}.
   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{INT_MAX} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.

   @emph{IMPORTANT:}
   In the case of a truncated standard distribution these calls always
   return the respective values of the @emph{untruncated} distribution!
*/


int unur_distr_discr_set_pmfstr( UNUR_DISTR *distribution, const char *pmfstr );
/* 
   This function provides an alternative way to set a PMF of the 
   @var{distribution}.
   @var{pmfstr} is a character string that contains the formula
   for the PMF, see @ref{StringFunct,,Function String}, for details.
   See also the remarks for the unur_distr_discr_set_pmf() call.

   It is not possible to call this funtion twice or to call this
   function after a unur_distr_discr_set_pmf() call.
*/

int unur_distr_discr_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
/* 
   This function provides an alternative way to set a CDF; analogously
   to the unur_distr_discr_set_pmfstr() call.
*/

char *unur_distr_discr_get_pmfstr( const UNUR_DISTR *distribution );
/* */

char *unur_distr_discr_get_cdfstr( const UNUR_DISTR *distribution );
/* 
   Get pointer to respective string for PMF and CDF
   of @var{distribution} that is given via the string interface.
   This call allocates memory to produce this string. It should be
   freed when it is not used any more.
*/

int unur_distr_discr_set_pmfparams( UNUR_DISTR *distribution, const double *params, int n_params );
/* 
   Set array of parameters for @var{distribution}. There is an upper limit
   for the number of parameters @var{n_params}. It is given by the
   macro @code{UNUR_DISTR_MAXPARAMS} in @file{unuran_config.h}. (It is set to
   5 but can be changed to any appropriate nonnegative number.)
   If @var{n_params} is negative or exceeds this limit no parameters
   are copied into the @var{distribution} object and @code{unur_errno}
   is set to @code{UNUR_ERR_DISTR_NPARAMS}. 


   For standard distributions from the UNU.RAN library the parameters
   are checked. Moreover, the domain is updated automatically unless it
   has been changed before by a unur_distr_discr_set_domain() call.
   If the given parameters are invalid for the standard distribution,
   then no parameters are set and an error code is returned.
   Notice that the given parameter list for such a distribution is
   handled in the same way as in the corresponding @command{new}
   calls, i.e. optional parameters for the PDF that are not present in
   the given list are (re-)set to their default values.

   @emph{Important:} Integer parameter must be given as @code{double}s.
*/

int unur_distr_discr_get_pmfparams( const UNUR_DISTR *distribution, const double **params );
/* 
   Get number of parameters of the PMF and set pointer
   @var{params} to array of parameters. If no parameters are stored
   in the object, an error code is returned and @code{params} is set to
   NULL.
*/

int unur_distr_discr_set_domain( UNUR_DISTR *distribution, int left, int right );
/* 
   Set the left and right borders of the domain of the
   @var{distribution}. This can also be used to truncate an existing
   distribution. For setting the boundary to @unurmath{\pm\infty} use
   @code{INT_MIN} and @code{INT_MAX}, respectively.
   If @var{right} is not strictly greater than @var{left} no domain
   is set and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
   It is allowed to use this call to increase the domain.
   If the PV of the discrete distribution is used,
   than the right boudary is ignored (and internally set to 
   @var{left} + size of PV @math{- 1}).
   Notice that @code{INT_MIN} and @code{INT_MAX} are interpreted as
   (minus/plus) infinity.

   Default: [@code{0}, @code{INT_MAX}].
*/

int unur_distr_discr_get_domain( const UNUR_DISTR *distribution, int *left, int *right );
/* 
   Get the left and right borders of the domain of the
   @var{distribution}. If the domain is not set explicitly 
   the interval [@code{INT_MIN}, @code{INT_MAX}] is assumed and returned.
   When a PV is given then the domain is set automatically to 
   [@code{0},size of PV @math{- 1}].
*/


/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the essential
   parameters has been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_discr_set_mode( UNUR_DISTR *distribution, int mode );
/* 
   Set mode of @var{distribution}.
*/

int unur_distr_discr_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the @var{distribution}. This call works properly
   for distribution objects from the 
   UNU.RAN library of standard distributions 
   when the corresponding function is available.
   Otherwise a (slow) numerical mode finder is used. It only works properly
   for unimodal probability mass functions. If it failes
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_discr_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of @var{distribution}. If the mode is not marked as known, 
   unur_distr_discr_upd_mode() is called to compute the mode. If this
   is not successful @code{INT_MAX} is returned and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the distribution at all.)
*/


int unur_distr_discr_set_pmfsum( UNUR_DISTR *distribution, double sum );
/* 
   Set the sum over the PMF. If @code{sum} is non-positive, no
   sum is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 

   For a distribution object created by the 
   UNU.RAN library of standard distributions you always should use
   the unur_distr_discr_upd_pmfsum(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_discr_upd_pmfsum( UNUR_DISTR *distribution );
/*
   Recompute the sum over the PMF of the @var{distribution}. 
   In most cases the normalization constant is recomputed and thus the
   sum is 1. This call works for distribution objects from the UNU.RAN
   library of standard distributions when the corresponding function
   is available. When a PV, a PMF with finite domain, or a CDF is
   given, a simple generic function which uses a naive summation loop
   is used. If this computation is not possible, an error code is
   returned and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}. 

   The call does not work for distributions from the 
   UNU.RAN library of standard distributions with truncated
   domain when the CDF is not available.
*/

double unur_distr_discr_get_pmfsum( UNUR_DISTR *distribution );
/* 
   Get the sum over the PMF of the @var{distribution}. If this sum is
   not known, unur_distr_discr_upd_pmfsum() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/* =END */

/*---------------------------------------------------------------------------*/
