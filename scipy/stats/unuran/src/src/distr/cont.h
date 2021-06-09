/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cont.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CONT  (continuous univariate distribution)                  *
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
   =NODEX   CONT   Continuous univariate distributions

   =UP Distribution_objects [10]

   =DESCRIPTION
      The calls in this section can be applied to continuous
      univariate distributions.

      @itemize @minus 
      @item Create a @command{new} instance of a continuous univariate
      distribution.

      @item Handle and evaluate 
      distribution function (CDF, @command{cdf}), 
      probability density function (PDF, @command{pdf}) and the 
      derivative of the density function (@command{dpdf}).
      The following is important:
      @itemize .
      @item @command{pdf} need not be normalized, i.e.,
      any integrable nonnegative function can be used. 
      @item @command{dpdf} must the derivate of the function provided
      as @command{pdf}.
      @item @command{cdf} must be a distribution function, i.e. it
      must be monotonically increasing with range [0,1].
      @item If @command{cdf} and @command{pdf} are used together for a
      pariticular generation method, then @command{pdf} must be the
      derivate of the @command{cdf}, i.e., it must be normalized.
      @end itemize

      @item Handle and evaluate 
      the logarithm of the probability density function (logPDF,
      @command{logpdf}) and the derivative of the logarithm of the
      density function (@command{dlogpdf}).

      Some methods use the logarithm of the density if available.

      @item Set (and change) parameters (@command{pdfparams}) and the
      area below the graph (@command{pdfarea}) of the given density.

      @item Set the @command{mode} (or pole) of the distribution. 

      @item Set the @command{center} of the distribution. 
      It is used by some generation methods to adjust the parameters
      of the generation algorithms to gain better performance. It can
      be seens as the location of the ``central part'' of the
      distribution. 

      @item Some generation methods require the hazard rate
      (@command{hr}) of the distribution instead of its @command{pdf}.

      @item Alternatively, @command{cdf}, @command{pdf}, @command{dpdf},
      and @command{hr} can be provided as @command{str}ings instead of
      function pointers.

      @item Set the @command{domain} of the distribution.  Notice that
      the library also can handle truncated distributions, i.e.,
      distributions that are derived from (standard) distributions by
      simply restricting its domain to a subset. However, there is a
      subtle difference between changing the domain of a distribution
      object by a unur_distr_cont_set_domain() call and changing the
      (truncated) domain for an existing generator object. The domain
      of the distribution object is used to create the generator
      object with hats, squeezes, tables, etc.  Whereas truncating the
      domain of an existing generator object need not necessarily
      require a recomputation of these data.  Thus by a
      @command{unur_<method>_chg_truncated} call (if available) the
      sampling region is restricted to the subset of the domain of the
      given distribution object. However, generation methods that
      require a recreation of the generator object when the domain is
      changed have a @command{unur_<method>_chg_domain} call instead.
      For these calls there are of course no restrictions on the given
      domain (i.e., it is possible to increase the domain of the
      distribution) (@pxref{Methods}, for details).

      @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/
/* 
   Routines for handling univariate continuous distributions (CONT).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cont_new( void );
/* 
   Create a new (empty) object for univariate continuous distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cont_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *pdf );
/* */

int unur_distr_cont_set_dpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *dpdf );
/* */

int unur_distr_cont_set_cdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *cdf );
/* */

int unur_distr_cont_set_invcdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *invcdf );
/* 
   Set respective pointer to the probability density function (PDF),
   the derivative of the probability density function (dPDF), the
   cumulative distribution function (CDF), and the inverse CDF of the
   @var{distribution}.
   Each of these function pointers must be of type
   @code{double funct(double x, const UNUR_DISTR *distr)}.

   Due to the fact that some of the methods do not require a
   normalized PDF the following is important:

   @itemize @minus
   @item
   The given CDF must be the cumulative distribution function of
   the (non-truncated) distribution. If a distribution from the 
   UNU.RAN library of standard distributions
   (@pxref{Stddist,,Standard distributions})
   is truncated, there is no need to change the CDF.

   @item
   If both the CDF and the PDF are used (for a method or for order
   statistics), the PDF must be the derivative of the CDF.
   If a truncated distribution for one of the standard distributions
   from the UNU.RAN library of standard distributions is used,
   there is no need to change the PDF.

   @item
   If the area below the PDF is required for a given distribution
   it must be given by the unur_distr_cont_set_pdfarea() call.
   For a truncated distribution this must be of course the integral of
   the PDF in the given truncated domain.
   For distributions from the UNU.RAN library of standard
   distributions this is done automatically by the
   unur_distr_cont_upd_pdfarea() call.

   @end itemize

   It is important to note that all these functions must return a
   result for all values of @var{x}. Eg., if the domain of a given
   PDF is the interval [-1,1], then the given function must return
   @code{0.0} for all points outside this interval.
   In case of an overflow the PDF should return 
   @code{UNUR_INFINITY}.

   It is not possible to change such a function. Once the PDF or
   CDF is set it cannot be overwritten. This also holds when the 
   logPDF is given or when the PDF
   is given by the unur_distr_cont_set_pdfstr() or
   unur_distr_cont_set_logpdfstr() call.
   A new distribution object has to be used instead.
*/


UNUR_FUNCT_CONT *unur_distr_cont_get_pdf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_CONT *unur_distr_cont_get_dpdf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_CONT *unur_distr_cont_get_cdf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_CONT *unur_distr_cont_get_invcdf( const UNUR_DISTR *distribution );
/* 
   Get the respective pointer to the PDF, the derivative of the 
   PDF, the CDF, and the inverse CDF of the @var{distribution}. The
   pointer is of type @code{double funct(double x, const UNUR_DISTR *distr)}.
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_cont_eval_pdf( double x, const UNUR_DISTR *distribution );
/* */

double unur_distr_cont_eval_dpdf( double x, const UNUR_DISTR *distribution );
/* */

double unur_distr_cont_eval_cdf( double x, const UNUR_DISTR *distribution );
/* */

double unur_distr_cont_eval_invcdf( double u, const UNUR_DISTR *distribution );
/* 
   Evaluate the PDF, derivative of the PDF, the CDF, and the inverse
   CDF at @var{x} and @var{u},respectively. 
   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.

   @emph{IMPORTANT:}
   In the case of a truncated standard distribution these calls always
   return the respective values of the @emph{untruncated} distribution!
*/

int unur_distr_cont_set_logpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *logpdf );
/* */

int unur_distr_cont_set_dlogpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *dlogpdf );
/* */

int unur_distr_cont_set_logcdf( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *logcdf );
/* */

UNUR_FUNCT_CONT *unur_distr_cont_get_logpdf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_CONT *unur_distr_cont_get_dlogpdf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCT_CONT *unur_distr_cont_get_logcdf( const UNUR_DISTR *distribution );
/* */

double unur_distr_cont_eval_logpdf( double x, const UNUR_DISTR *distribution );
/* */

double unur_distr_cont_eval_dlogpdf( double x, const UNUR_DISTR *distribution );
/* */

double unur_distr_cont_eval_logcdf( double x, const UNUR_DISTR *distribution );
/* 
   Analogous calls for the logarithm of the density distribution functions.
*/


int unur_distr_cont_set_pdfstr( UNUR_DISTR *distribution, const char *pdfstr );
/* 
   This function provides an alternative way to set a PDF and its
   derivative of the @var{distribution}.
   @var{pdfstr} is a character string that contains the formula
   for the PDF, see @ref{StringFunct,,Function String}, for details.
   The derivative of the given PDF is computed automatically.
   See also the remarks for the unur_distr_cont_set_pdf() call.

   It is not possible to call this funtion twice or to call this
   function after a unur_distr_cont_set_pdf() call.
*/

int unur_distr_cont_set_cdfstr( UNUR_DISTR *distribution, const char *cdfstr );
/* 
   This function provides an alternative way to set a CDF; analogously
   to the unur_distr_cont_set_pdfstr() call.
   The PDF and its derivative of the given CDF are computed automatically.
*/

char *unur_distr_cont_get_pdfstr( const UNUR_DISTR *distribution );
/* */

char *unur_distr_cont_get_dpdfstr( const UNUR_DISTR *distribution );
/* */

char *unur_distr_cont_get_cdfstr( const UNUR_DISTR *distribution );
/* 
   Get pointer to respective string for PDF, derivate of PDF, and CDF
   of @var{distribution} that is given as string (instead of a
   function pointer).
   This call allocates memory to produce this string. It should be
   freed when it is not used any more.
*/

int unur_distr_cont_set_pdfparams( UNUR_DISTR *distribution, const double *params, int n_params );
/* 
   Sets array of parameters for @var{distribution}. There is an upper limit
   for the number of parameters @code{n_params}. It is given by the
   macro @code{UNUR_DISTR_MAXPARAMS} in @file{unuran_config.h}. (It is set to
   5 by default but can be changed to any appropriate nonnegative number.)
   If @var{n_params} is negative or exceeds this limit no parameters
   are copied into the distribution object and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.

   For standard distributions from the UNU.RAN library the parameters 
   are checked. Moreover, the domain is updated automatically unless it
   has been changed before by a unur_distr_cont_set_domain() call.
   If the given parameters are invalid for the standard distribution,
   then no parameters are set and an error code is returned.
   Notice, that the given parameter list for such a distribution is
   handled in the same way as in the corresponding @command{new}
   calls, i.e. optional parameters for the PDF that are not present in
   the given list are (re-)set to their default values.

   @strong{Important:} If the parameters of a distribution from the 
   UNU.RAN library of standard distributions
   (@pxref{Stddist,,Standard distributions})
   are changed, then neither its mode nor the normalization 
   constant are updated. Please use the respective calls
   unur_distr_cont_upd_mode() and unur_distr_cont_upd_pdfarea().
   Moreover, if the domain has been changed by a
   unur_distr_cont_set_domain() it is not automatically updated, either.
   Updating the normalization constant is in particular very important,
   when the CDF of the distribution is used.
*/

int unur_distr_cont_get_pdfparams( const UNUR_DISTR *distribution, const double **params );
/* 
   Get number of parameters of the PDF and set pointer @var{params} to
   array of parameters. If no parameters are stored in the object, an
   error code is returned and @code{params} is set to NULL.
   
   @emph{Important:} Do @strong{not} change the entries in @var{params}!
*/

int unur_distr_cont_set_pdfparams_vec( UNUR_DISTR *distribution, int par, const double *param_vec, int n_param_vec );
/* 
   This function provides an interface for additional vector parameters for a
   continuous @var{distribution}.

   It sets the parameter with number @var{par}. 
   @var{par} indicates directly which of the parameters is set and
   must be a number between @code{0} and @code{UNUR_DISTR_MAXPARAMS}-1
   (the upper limit of possible parameters defined in
   @file{unuran_config.h}; it is set to 5 but can be changed to any 
   appropriate nonnegative number.)

   The entries of a this parameter are given by the array @var{param_vec}
   of size @var{n_param_vec}. 

   If @var{param_vec} is NULL then the corresponding entry is cleared.

   If an error occurs no parameters are copied into the parameter
   object @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cont_get_pdfparams_vec( const UNUR_DISTR *distribution, int par, const double **param_vecs );
/* 
   Get parameter of the PDF with number @var{par}.
   The pointer to the parameter array is stored in @var{param_vecs}, its
   size is returned by the function.
   If the requested parameter is not set, then an error code is returned
   and @code{params} is set to NULL.

   @emph{Important:} Do @strong{not} change the entries in @var{param_vecs}!
*/

int unur_distr_cont_set_logpdfstr( UNUR_DISTR *distribution, const char *logpdfstr );
/* */

char *unur_distr_cont_get_logpdfstr( const UNUR_DISTR *distribution );
/* */

char *unur_distr_cont_get_dlogpdfstr( const UNUR_DISTR *distribution );
/* */

int unur_distr_cont_set_logcdfstr( UNUR_DISTR *distribution, const char *logcdfstr );
/* */

char *unur_distr_cont_get_logcdfstr( const UNUR_DISTR *distribution );
/* 
   Analogous calls for the logarithm of the density and distribution functions.
*/


int unur_distr_cont_set_domain( UNUR_DISTR *distribution, double left, double right );
/* 
   Set the left and right borders of the domain of the
   distribution. This can also be used to truncate an existing
   distribution. For setting the boundary to @unurmath{\pm\infty}
   use @code{+/- UNUR_INFINITY}.
   If @var{right} is not strictly greater than @var{left} no domain
   is set and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.

   @emph{Important:} For some technical reasons it is assumed that the density 
   is unimodal and thus monotone on either side of the mode! This is used in
   the case when the given mode is outside of the original domain. Then the
   mode is set to the corresponding boundary of the new domain.
   If this result is not the desired it must be changed by using a 
   unur_distr_cont_set_mode() call (or a unur_distr_cont_upd_mode()
   call). The same holds for the center of the distribution.
*/

int unur_distr_cont_get_domain( const UNUR_DISTR *distribution, double *left, double *right );
/* 
   Get the left and right borders of the domain of the
   distribution. If the domain is not set @code{+/- UNUR_INFINITY} is
   assumed and returned. No error is reported in this case.
*/

int unur_distr_cont_get_truncated( const UNUR_DISTR *distribution, double *left, double *right );
/* 
   Get the left and right borders of the (truncated) domain of the
   distribution. For non-truncated distribution this call is
   equivalent to the unur_distr_cont_get_domain() call.

   This call is only useful in connection with a unur_get_distr() call
   to get the boundaries of the sampling region of a generator object.
*/

int unur_distr_cont_set_hr( UNUR_DISTR *distribution, UNUR_FUNCT_CONT *hazard );
/* 
   Set pointer to the hazard rate (HR) of the @var{distribution}.

   The @emph{hazard rate} (or failure rate) is a mathematical way of
   describing aging. If the lifetime @i{X} is a random variable with
   density @i{f(x)} and CDF @i{F(x)} the hazard rate @i{h(x)} 
   is defined as @i{h(x) = f(x) / (1-F(x))}.
   In other words, @i{h(x)} represents the (conditional) rate of
   failure of a unit that has survived up to time @i{x} with
   probability @i{1-F(x)}. 
   The key distribution is the exponential distribution as it has
   constant hazard rate of value 1. Hazard rates tending to infinity
   describe distributions with sub-exponential tails whereas
   distributions with hazard rates tending to zero have heavier tails
   than the exponential distribution. 

   It is important to note that all these functions must return a
   result for all floats @i{x}. In case of an overflow the PDF should
   return @code{UNUR_INFINITY}.

   @strong{Important}: Do not simply use @i{f(x) / (1-F(x))}, since
   this is numerically very unstable and results in numerical noise
   if @i{F(x)} is (very) close to 1. Moreover, if the density @i{f(x)}
   is known a generation method that uses the density is more 
   appropriate.

   It is not possible to change such a function. Once the HR is set it
   cannot be overwritten. This also holds when the HR is given by the
   unur_distr_cont_set_hrstr() call. A new distribution object has to
   be used instead.
*/

UNUR_FUNCT_CONT *unur_distr_cont_get_hr( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the hazard rate of the @var{distribution}. The
   pointer is of type 
   @code{double funct(double x, const UNUR_DISTR *distr)}.
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
*/

double unur_distr_cont_eval_hr( double x, const UNUR_DISTR *distribution );
/* 
   Evaluate the hazard rate at @var{x}. 
   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cont_set_hrstr( UNUR_DISTR *distribution, const char *hrstr );
/* 
   This function provides an alternative way to set a hazard rate and its
   derivative of the @var{distribution}.
   @var{hrstr} is a character string that contains the formula
   for the HR, see @ref{StringFunct,,Function String}, for details.
   See also the remarks for the unur_distr_cont_set_hr() call.

   It is not possible to call this funtion twice or to call this
   function after a unur_distr_cont_set_hr() call.
*/

char *unur_distr_cont_get_hrstr( const UNUR_DISTR *distribution );
/* 
   Get pointer to string for HR of @var{distribution} that is given
   via the string interface. This call allocates memory to produce
   this string. It should be freed when it is not used any more.
*/


/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the essential
   parameters has been set or changed (and the parameter is required
   for the chosen method).
*/

int unur_distr_cont_set_mode( UNUR_DISTR *distribution, double mode );
/* 
   Set mode of @var{distribution}. The @var{mode} must be contained in
   the domain of @var{distribution}. Otherwise the mode is not set and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.
   For distributions with unbounded density, this call is used to set
   the pole of the PDF. Notice that the PDF should then return
   UNUR_INFINITY at the pole.
   Notice that the mode is adjusted when the domain is set, see the
   remark for the unur_distr_cont_set_domain() call.
*/

int unur_distr_cont_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the @var{distribution}. This call works
   properly for distribution objects from the UNU.RAN library of
   standard distributions when the corresponding function is
   available.  Otherwise a (slow) numerical mode finder based on
   Brent's algorithm is used. If it failes @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/

double unur_distr_cont_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of @var{distribution}. If the mode is not marked as known, 
   unur_distr_cont_upd_mode() is called to compute the mode. If this
   is not successful @code{UNUR_INFINITY} is returned and 
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the distribution at all.)
*/

int unur_distr_cont_set_center( UNUR_DISTR *distribution, double center );
/* 
   Set center of the @var{distribution}. The center is used by some
   methods to shift the distribution in order to decrease numerical
   round-off error. If not given explicitly a default is used.

   @emph{Important:} This call does not check whether the center is
   contained in the given domain. 

   Default: The mode, if set by a unur_distr_cont_set_mode() or 
   unur_distr_cont_upd_mode() call; otherwise @code{0}.
*/

double unur_distr_cont_get_center( const UNUR_DISTR *distribution );
/* 
   Get center of the @var{distribution}. It always returns some point
   as there always exists a default for the center, see
   unur_distr_cont_set_center().
*/

int unur_distr_cont_set_pdfarea( UNUR_DISTR *distribution, double area );
/* 
   Set the area below the PDF. If @code{area} is non-positive, no
   area is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 
   
   For a distribution object created by the 
   UNU.RAN library of standard distributions you always should use
   the unur_distr_cont_upd_pdfarea(). Otherwise there might be
   ambiguous side-effects.
*/

int unur_distr_cont_upd_pdfarea( UNUR_DISTR *distribution );
/*
   Recompute the area below the PDF of the distribution. 
   It only works for distribution objects from the
   UNU.RAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 

   This call also sets the normalization constant such that the given
   PDF is the derivative of a given CDF, i.e. the area is 1.
   However, for truncated distributions the area is smaller than 1.

   The call does not work for distributions from the 
   UNU.RAN library of standard distributions with truncated
   domain when the CDF is not available.
*/

double unur_distr_cont_get_pdfarea( UNUR_DISTR *distribution );
/* 
   Get the area below the PDF of the distribution. If this area is
   not known,@* unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/* =END */

/*---------------------------------------------------------------------------*/

