/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cvec.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         type  CVEC  (continuous multivariate distribution)                *
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
   =NODEX   CVEC   Continuous multivariate distributions

   =UP Distribution_objects [30]

   =DESCRIPTION
      The following calls handle multivariate distributions.
      However, the requirements of particular generation methods is not
      as unique as for univariate distributions. Moreover, random vector
      generation methods are still under development.
      The below functions are a first attempt to handle this situation. 
      
      Notice that some of the parameters -- when given carelessly -- might
      contradict to others. For example: Some methods require the
      marginal distribution and some methods need a standardized form of
      the marginal distributions, where the actual mean and variance is
      stored in the mean vector and the covariance matrix, respectively. 
      
      We also have to mention that some methods might abuse some of the
      parameters. Please read the discription of the chosen sampling
      method carfully.
      
      The following kind of calls exists:
      
      @itemize @minus 
      @item Create a @command{new} instance of a continuous multivariate
      distribution; 

      @item Handle and evaluate 
      probability density function (PDF, @command{pdf}) and the 
      gradient of the density function (@command{dpdf}).
      The following is important:
      @itemize .
      @item @command{pdf} need not be normalized, i.e.,
      any integrable nonnegative function can be used. 
      @item @command{dpdf} must the derivate of the function provided
      as @command{pdf}.
      @end itemize

      @item Handle and evaluate 
      the logarithm of the probability density function (logPDF,
      @command{logpdf}) and the gradient of the logarithm of the
      density function (@command{dlogpdf}).

      Some methods use the logarithm of the density if available.

      @item Set (and change) parameters (@command{pdfparams}) and the
      volume below the graph (@command{pdfvol}) of the given density.

      @item Set @command{mode} and @command{mean} of the distribution. 

      @item Set the @command{center} of the distribution. 
      It is used by some generation methods to adjust the parameters
      of the generation algorithms to gain better performance. It can
      be seens as the location of the ``central part'' of the
      distribution. 

      @item Handle the @command{covar}iance matrix of the distribution and
      its @command{cholesky} and @command{inv}verse matrices.

      @item Set the @command{rankcorr}elation matrix of the distribution.

      @item Deal with @command{marginal} distributions.

      @item Set domain of the distribution.

   @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling multivariate continuous distributions (CVEC).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cvec_new( int dim );
/* 
   Create a new (empty) object for multivariate continuous
   distribution. @var{dim} is the number of components of the random
   vector (i.e. its dimension). It is also possible to use dimension 1.
   Notice, however, that this is treated as a distribution of random 
   vectors with only one component and not as a distribution of 
   real numbers. For the latter unur_distr_cont_new() should be used
   to create an object for a univariate distribution.
*/

/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cvec_set_pdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *pdf );
/* 
   Set respective pointer to the PDF of the @var{distribution}.
   This function must be of type 
   @code{double funct(const double *x, UNUR_DISTR *distr)},
   where @var{x} must be a pointer to a double array of appropriate
   size (i.e. of the same size as given to the unur_distr_cvec_new()
   call).

   It is not necessary that the given PDF is normalized, i.e. the
   integral need not be 1. 
   Nevertheless the volume below the PDF can be provided by a
   unur_distr_cvec_set_pdfvol() call.

   It is not possible to change the PDF. Once the PDF is set it cannot
   be overwritten. This also holds when the logPDF is given.
   A new distribution object has to be used instead.
*/

int unur_distr_cvec_set_dpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dpdf );
/* 
   Set pointer to the gradient of the PDF. The type of this function must be
   @code{int funct(double *result, const double *x, UNUR_DISTR *distr)},
   where @var{result} and @var{x} must be pointers to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).
   The gradient of the PDF is stored in the array @var{result}.
   The function should return an error code in case of an error and must
   return @code{UNUR_SUCCESS} otherwise.

   The given function must be the gradient of the function
   given by a unur_distr_cvec_set_pdf() call.

   It is not possible to change the gradient of the PDF. Once the dPDF
   is set it cannot be overwritten. This also holds when the gradient
   of the logPDF is given.
   A new distribution object has to be used instead.
*/

int unur_distr_cvec_set_pdpdf( UNUR_DISTR *distribution, UNUR_FUNCTD_CVEC *pdpdf );
/* 
   Set pointer to partial derivatives of the PDF. The type of this function must be
   @code{double funct(const double *x, int coord, UNUR_DISTR *distr)},
   where @var{x} must be a pointer to a double array of appropriate
   size (i.e. of the same size as given to the unur_distr_cvec_new()
   call). @var{coord} is the coordinate for which the partial dervative should be
   computed.

   Notice that @var{coord} must be an integer from @{0,@dots{},dim-1@}. 
   
   It is not possible to change the partial derivative of the PDF. Once the pdPDF
   is set it cannot be overwritten. This also holds when the partial derivative 
   of the logPDF is given.
   A new distribution object has to be used instead.
*/

UNUR_FUNCT_CVEC *unur_distr_cvec_get_pdf( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the PDF of the @var{distribution}. The
   pointer is of type 
   @code{double funct(const double *x, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the
   @var{distribution}, the NULL pointer is returned.
*/

UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dpdf( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the gradient of the PDF of the
   @var{distribution}. The pointer is of type 
   @code{int double funct(double *result, const double *x, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the
   @var{distribution}, the NULL pointer is returned.
*/

UNUR_FUNCTD_CVEC *unur_distr_cvec_get_pdpdf( const UNUR_DISTR *distribution );
/* 
   Get the pointer to the partial derivative of the PDF of the @var{distribution}.
   The pointer is of type 
   @code{double funct(const double *x, int coord, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the
   @var{distribution}, the NULL pointer is returned.
*/

double unur_distr_cvec_eval_pdf( const double *x, UNUR_DISTR *distribution );
/* 
   Evaluate the PDF of the @var{distribution} at @var{x}.
   @var{x} must be a pointer to a double array of appropriate size
   (i.e. of the same size as given to the unur_distr_cvec_new() call)
   that contains the vector for which the function has to be evaluated.

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   @var{distribution}, @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_eval_dpdf( double *result, const double *x, UNUR_DISTR *distribution );
/* 
   Evaluate the gradient of the PDF of the @var{distribution} at
   @var{x}. 
   The result is stored in the double array @var{result}.
   Both @var{result} and @var{x} must be pointer to double arrays of
   appropriate size (i.e. of the same size as given to the
   unur_distr_cvec_new() call).

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   @var{distribution}, an error code is returned and @code{unur_errno}
   is set to @code{UNUR_ERR_DISTR_DATA} (@var{result} is left unmodified).
*/

double unur_distr_cvec_eval_pdpdf( const double *x, int coord, UNUR_DISTR *distribution );
/* 
   Evaluate the partial derivative of the PDF of the @var{distribution} 
   at @var{x} for the coordinate @var{coord}.
   @var{x} must be a pointer to a double array of appropriate size
   (i.e. of the same size as given to the unur_distr_cvec_new() call)
   that contains the vector for which the function has to be evaluated.

   Notice that @var{coord} must be an integer from @{0,@dots{},dim-1@}. 

   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the
   @var{distribution}, @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_set_logpdf( UNUR_DISTR *distribution, UNUR_FUNCT_CVEC *logpdf );
/* */

int unur_distr_cvec_set_dlogpdf( UNUR_DISTR *distribution, UNUR_VFUNCT_CVEC *dlogpdf );
/* */

int unur_distr_cvec_set_pdlogpdf( UNUR_DISTR *distribution, UNUR_FUNCTD_CVEC *pdlogpdf );
/* */

UNUR_FUNCT_CVEC *unur_distr_cvec_get_logpdf( const UNUR_DISTR *distribution );
/* */

UNUR_VFUNCT_CVEC *unur_distr_cvec_get_dlogpdf( const UNUR_DISTR *distribution );
/* */

UNUR_FUNCTD_CVEC *unur_distr_cvec_get_pdlogpdf( const UNUR_DISTR *distribution );
/* */

double unur_distr_cvec_eval_logpdf( const double *x, UNUR_DISTR *distribution );
/* */

int unur_distr_cvec_eval_dlogpdf( double *result, const double *x, UNUR_DISTR *distribution );
/* */

double unur_distr_cvec_eval_pdlogpdf( const double *x, int coord, UNUR_DISTR *distribution );
/* 
   Analogous calls for the logarithm of the density function.
*/


int unur_distr_cvec_set_mean( UNUR_DISTR *distribution, const double *mean );
/* 
   Set mean vector for multivariate @var{distribution}.
   @var{mean} must be a pointer to an array of size @code{dim}, where
   @code{dim} is the dimension returned by unur_distr_get_dim().
   A NULL pointer for @var{mean} is interpreted as the zero
   vector (0,@dots{},0).

   @strong{Important:} If the parameters of a distribution from the 
   UNU.RAN library of standard distributions
   (@pxref{Stddist,,Standard distributions})
   are changed, then neither its mode nor the normalization 
   constant are updated. Please use the respective calls
   unur_distr_cvec_upd_mode() and unur_distr_cvec_upd_pdfvol().
*/

const double *unur_distr_cvec_get_mean( const UNUR_DISTR *distribution );
/* 
   Get the mean vector of the @var{distribution}. The function returns a
   pointer to an array of size @code{dim}.
   If the mean vector is not marked as known the NULL pointer is
   returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 
   
   @emph{Important:} Do @strong{not} modify the array that holds the
   mean vector!
*/

int unur_distr_cvec_set_covar( UNUR_DISTR *distribution, const double *covar );
/* 
   Set covariance matrix for multivariate @var{distribution}.
   @var{covar} must be a pointer to an array of size
   @code{dim} x @code{dim}, where @code{dim} is the dimension returned
   by unur_distr_get_dim(). The rows of the matrix have to be stored
   consecutively in this array.

   @var{covar} must be a variance-covariance matrix of the
   @var{distribution}, i.e. it must be symmetric and positive definit and
   its diagonal entries (i.e. the variance of the components of the
   random vector) must be strictly positive.
   The Cholesky factor is computed (and stored) to verify the positive
   definiteness condition.
   Notice that the inverse of the given covariance matrix is
   automatically computed when it is requested by some routine.
   Notice that the computation of this inverse matrix is unstable in
   case of high correlations and/or high dimensions. Thus it might
   fail and methods that require this inverse cannot be used.
   As an alternative the inverse of the covariance matrix can be
   directly set by a unur_distr_cvec_set_covar_inv() call. 

   A NULL pointer for @var{covar} is interpreted as the
   identity matrix.

   @emph{Important:} This entry is abused in some methods which do not
   require the covariance matrix. It is then used to perform some
   transformation to obtain better performance.

   @emph{Important:} In case of an error (e.g. because @var{covar} is
   not a valid covariance matrix) an error code is returned.
   Moreover, the covariance matrix is not set and is marked as
   unknown. A previously set covariance matrix is then no longer
   available.

   @strong{Important:} If the parameters of a distribution from the 
   UNU.RAN library of standard distributions
   (@pxref{Stddist,,Standard distributions})
   are changed, then neither its mode nor the normalization 
   constant are updated. Please use the respective calls
   unur_distr_cvec_upd_mode() and unur_distr_cvec_upd_pdfvol().

   @emph{Remark:} UNU.RAN does not check whether the an eventually
   set covariance matrix and a rank-correlation matrix do not
   contradict each other.
*/

int unur_distr_cvec_set_covar_inv( UNUR_DISTR *distribution, const double *covar_inv );
/*
   Set inverse of the covariance matrix for multivariate @var{distribution}.
   @var{covar_inv} must be a pointer to an array of size
   @code{dim} x @code{dim}, where @code{dim} is the dimension returned
   by unur_distr_get_dim(). The rows of the matrix have to be stored
   consecutively in this array.

   @var{covar_inv} must be symmetric and positive definit. Only the
   symmetry of the matrix is checked. 

   A NULL pointer for @var{covar_inv} is interpreted as the identity matrix.

   @emph{Important:} In case of an error (because @var{covar_inv} is
   not symetric) an error code is returned.
   Moreover, the inverse of the covariance matrix is not set and is
   marked as unknown. A previously set inverse matrix is then no longer
   available.

   @emph{Remark:} UNU.RAN does not check whether the given matrix is
   positive definit.

   @emph{Remark:} UNU.RAN does not check whether the matrix
   @var{covar_inv} is the inverse of the eventually set covariance
   matrix.
*/


const double *unur_distr_cvec_get_covar( const UNUR_DISTR *distribution );
/* */

const double *unur_distr_cvec_get_cholesky( const UNUR_DISTR *distribution );
/* */

const double *unur_distr_cvec_get_covar_inv( UNUR_DISTR *distribution );
/*
   Get covariance matrix of @var{distribution}, its Cholesky factor,
   and its inverse, respectively. The function returns a
   pointer to an array of size @code{dim} x @code{dim}.
   The rows of the matrix are stored consecutively in this array.
   If the requested matrix is not marked as known the NULL
   pointer is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 

   @emph{Important:} Do @strong{not} modify the array that holds the
   covariance matrix!

   @emph{Remark:} The inverse of the covariance matrix is computed 
   if it is not already stored.
*/

int unur_distr_cvec_set_rankcorr( UNUR_DISTR *distribution, const double *rankcorr );
/* 
   Set rank-correlation matrix (Spearman's correlation) for
   multivariate @var{distribution}. 
   @var{rankcorr} must be a pointer to an array of size
   @code{dim} x @code{dim}, where @code{dim} is the dimension returned
   by unur_distr_get_dim(). The rows of the matrix have to be stored
   consecutively in this array.

   @var{rankcorr} must be a rank-correlation matrix of the
   @var{distribution}, i.e. it must be symmetric and positive definite
   and its diagonal entries must be equal to @code{1}.

   The Cholesky factor is computed (and stored) to verify the
   positive definiteness condition.

   A NULL pointer for @var{rankcorr} is interpreted as the identity matrix.

   @emph{Important:} In case of an error (e.g. because @var{rankcorr} is
   not a valid rank-correlation matrix) an error code is returned.
   Moreover, the rank-correlation matrix is not set and is marked as
   unknown. A previously set rank-correlation matrix is then no longer
   available.

   @emph{Remark:} UNU.RAN does not check whether the an eventually
   set covariance matrix and a rank-correlation matrix do not
   contradict each other.
*/

const double *unur_distr_cvec_get_rankcorr( const UNUR_DISTR *distribution );
/* */

const double *unur_distr_cvec_get_rk_cholesky( const UNUR_DISTR *distribution );
/*
   Get rank-correlation matrix and its cholesky factor, respectively,
   of @var{distribution}. The function
   returns a pointer to an array of size @code{dim} x @code{dim}.
   The rows of the matrix are stored consecutively in this array.
   If the requested matrix is not marked as known the NULL
   pointer is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_GET}. 

   @emph{Important:} Do @strong{not} modify the array that holds the
   rank-correlation matrix!
*/

int unur_distr_cvec_set_marginals( UNUR_DISTR *distribution, UNUR_DISTR *marginal );
/* 
   Sets marginal distributions of the given @var{distribution} to the
   same @var{marginal} distribution object. The @var{marginal}
   distribution must be an instance of a continuous univariate
   distribution object. Notice that the marginal distribution is
   copied into the @var{distribution} object.
*/

int unur_distr_cvec_set_marginal_array( UNUR_DISTR *distribution, UNUR_DISTR **marginals );
/* 
   Analogously to the above unur_distr_cvec_set_marginals() call.
   However, now an array @var{marginals} of the pointers to each of
   the marginal distributions must be given. It @strong{must} be an
   array of size @code{dim}, where @code{dim} is the dimension
   returned by unur_distr_get_dim(). 

   @emph{Notice}: Local copies for each of the entries are stored in
   the @var{distribution} object. If some of these entries are
   identical (i.e. contain the same pointer), then for each of these a
   new copy is made.
*/

int unur_distr_cvec_set_marginal_list( UNUR_DISTR *distribution, ... );
/* 
   Similar to the above unur_distr_cvec_set_marginal_array() call.
   However, now the pointers to the particular marginal distributions
   can be given as parameter and does not require an array of
   pointers. Additionally the given distribution objects are
   immediately destroyed. Thus calls like unur_distr_normal() can be 
   used as arguments. 
   (With unur_distr_cvec_set_marginal_array() the result of such call
   has to be stored in a pointer since it has to be freed afterwarts
   to avoid memory leaks!)

   The number of pointers to in the list of function arguments 
   @strong{must} be equal to the dimension of the @var{distribution},
   i.e. the dimension returned by unur_distr_get_dim(). 

   If one of the given pointer to marginal distributions is the NULL
   pointer then the marginal distributions of @var{distribution} are
   not set (or previous settings are not changed) and an error code is
   returned.

   @strong{Important:} All distribution objects given in the argument
   list are destroyed!
*/

const UNUR_DISTR *unur_distr_cvec_get_marginal( const UNUR_DISTR *distribution, int n );
/* 
   Get pointer to the @var{n}-th marginal distribution
   object from the given multivariate @var{distribution}. 
   If this does not exist, NULL is returned. 
   The marginal distributions are enumerated from @code{1}
   to @code{dim}, where @code{dim} is the dimension
   returned by unur_distr_get_dim(). 
*/

int unur_distr_cvec_set_pdfparams( UNUR_DISTR *distribution, const double *params, int n_params );
/* 
   Sets array of parameters for @var{distribution}. There is an upper limit
   for the number of parameters @code{n_params}. It is given by the
   macro @code{UNUR_DISTR_MAXPARAMS} in @file{unuran_config.h}. (It is set to
   5 by default but can be changed to any appropriate nonnegative number.)
   If @var{n_params} is negative or exceeds this limit no parameters
   are copied into the distribution object and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_NPARAMS}.

   For standard distributions from the UNU.RAN library the parameters 
   are checked. Moreover, the domain is updated automatically.
   If the given parameters are invalid for the standard distribution,
   then no parameters are set and an error code is returned.
   Notice that the given parameter list for such a distribution is
   handled in the same way as in the corresponding @command{new}
   calls, i.e. optional parameters for the PDF that are not present in
   the given list are (re-)set to their default values.

   @strong{Important:} If the parameters of a distribution from the 
   UNU.RAN library of standard distributions
   (@pxref{Stddist,,Standard distributions})
   are changed, then neither its mode nor the normalization 
   constant are updated. Please use the respective calls
   unur_distr_cvec_upd_mode() and unur_distr_cvec_upd_pdfvol().
*/

int unur_distr_cvec_get_pdfparams( const UNUR_DISTR *distribution, const double **params );
/* 
   Get number of parameters of the PDF and set pointer @var{params} to
   array of parameters. If no parameters are stored in the object, an
   error code is returned and @code{params} is set to NULL.
   
   @emph{Important:} Do @strong{not} change the entries in @var{params}!
*/


int unur_distr_cvec_set_pdfparams_vec( UNUR_DISTR *distribution, int par, const double *param_vec, int n_params );
/* 
   This function provides an interface for additional vector parameters for a
   multivariate @var{distribution} besides mean vector and covariance
   matrix which have their own calls.

   It sets the parameter with number @var{par}. 
   @var{par} indicates directly which of the parameters is set and
   must be a number between @code{0} and @code{UNUR_DISTR_MAXPARAMS}-1
   (the upper limit of possible parameters defined in
   @file{unuran_config.h}; it is set to 5 but can be changed to any 
   appropriate nonnegative number.)

   The entries of a this parameter are given by the array @var{param_vec}
   of size @var{n_params}. Notice that using this interface an
   An (@i{n} x @i{m})-matrix has to be stored in an array of length
   @var{n_params} = @i{n} times @i{m}; where the rows of the matrix
   are stored consecutively in this array.

   Due to great variety of possible parameters for a multivariate
   @var{distribution} there is no simpler interface.

   If @var{param_vec} is NULL then the corresponding entry is cleared.

   @strong{Important:} If the parameters of a distribution from the 
   UNU.RAN library of standard distributions
   (@pxref{Stddist,,Standard distributions})
   are changed, then neither its mode nor the normalization 
   constant are updated. Please use the respective calls
   unur_distr_cvec_upd_mode() and unur_distr_cvec_upd_pdfvol().

   If an error occurs no parameters are copied into the parameter
   object @code{unur_errno} is set to @code{UNUR_ERR_DISTR_DATA}.
*/

int unur_distr_cvec_get_pdfparams_vec( const UNUR_DISTR *distribution, int par, const double **param_vecs );
/* 
   Get parameter of the PDF with number @var{par}.
   The pointer to the parameter array is stored in @var{param_vecs}, its
   size is returned by the function.
   If the requested parameter is not set, then an error code is returned
   and @code{params} is set to NULL.

   @emph{Important:} Do @strong{not} change the entries in @var{param_vecs}!
*/

int unur_distr_cvec_set_domain_rect( UNUR_DISTR *distribution, const double *lowerleft, const double *upperright );
/* 
   Set rectangular domain for @var{distribution} with @var{lowerleft}
   and @var{upperright} vertices. Both must be pointer to an 
   array of the size returned by unur_distr_get_dim().
   A NULL pointer is interpreted as the zero vector (0,@dots{},0).
   For setting a coordinate of the boundary to @unurmath{\pm\infty}
   use @code{+/- UNUR_INFINITY}.
   The @var{lowerleft} vertex must be strictly smaller than
   @var{upperright} in each component. Otherwise no domain
   is set and @code{unur_errno} is set to @code{UNUR_ERR_DISTR_SET}.

   By default the domain of a distribution is unbounded. Thus one can
   use this call to truncate an existing distribution.

   @emph{Important:} Changing the domain of @var{distribution} 
   marks derived parameters like the mode or the center as unknown and
   must be set @emph{after} changing the domain. This is important for
   the already set (or default) value for the center does not 
   fall into the given domain.
   Notice that calls of the PDF and derived functions return @code{0.}
   when the parameter is not contained in the domain.
*/

int unur_distr_cvec_is_indomain( const double *x, const UNUR_DISTR *distribution );
/* 
   Check whether @var{x} falls into the domain of @var{distribution}.
*/

/*---------------------------------------------------------------------------*/

/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the
   essential parameters has been set or changed (and the parameter is
   required for the chosen method).
*/

int unur_distr_cvec_set_mode( UNUR_DISTR *distribution, const double *mode );
/* 
   Set mode of the @var{distribution}. @var{mode} must be a pointer to an
   array of the size returned by unur_distr_get_dim().
   A NULL pointer for @var{mode} is interpreted as the zero
   vector (0,@dots{},0).
*/

int unur_distr_cvec_upd_mode( UNUR_DISTR *distribution );
/* 
   Recompute the mode of the @var{distribution}. This call works
   properly for distribution objects from the UNU.RAN library of
   standard distributions when the corresponding function is
   available. If it failes @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
*/

const double *unur_distr_cvec_get_mode( UNUR_DISTR *distribution );
/* 
   Get mode of the @var{distribution}. The function returns a pointer to
   an array of the size returned by unur_distr_get_dim().
   If the mode is not marked as known the NULL pointer is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}. 
   (There is no difference between the case where no routine for
   computing the mode is available and the case where no mode exists
   for the @var{distribution} at all.)

   @emph{Important:} Do @strong{not} modify the array that holds the mode!
*/

int unur_distr_cvec_set_center( UNUR_DISTR *distribution, const double *center );
/* 
   Set center of the @var{distribution}. @var{center} must be a pointer to an
   array of the size returned by unur_distr_get_dim().
   A NULL pointer for @var{center} is interpreted as the zero
   vector (0,@dots{},0).

   The center is used by some methods to shift the distribution in
   order to decrease numerical round-off error.
   If not given explicitly a default is used.
   Moreover, it is used as starting point for several numerical search
   algorithm (e.g. for the mode). Then @var{center} must be a pointer
   where the call to the PDF returns a non-zero value.
   In particular @var{center} must contained in the domain of the distribution.

   Default: The mode, if given by a unur_distr_cvec_set_mode() call;
   else the mean, if given by a unur_distr_cvec_set_mean() call;
   otherwise the null vector (0,@dots{},0).
*/

const double *unur_distr_cvec_get_center( UNUR_DISTR *distribution );
/* 
   Get center of the @var{distribution}. The function returns a pointer to
   an array of the size returned by unur_distr_get_dim().
   It always returns some point as there always exists a default for
   the center, see unur_distr_cvec_set_center().

   @emph{Important:} Do @strong{not} modify the array that holds the center!
*/

int unur_distr_cvec_set_pdfvol( UNUR_DISTR *distribution, double volume );
/* 
   Set the volume below the PDF. If @var{vol} is non-positive, no
   volume is set and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_SET}. 
*/

int unur_distr_cvec_upd_pdfvol( UNUR_DISTR *distribution );
/*
   Recompute the volume below the PDF of the distribution. 
   It only works for distribution objects from the
   UNU.RAN library of standard distributions when the
   corresponding function is available. Otherwise @code{unur_errno} is
   set to @code{UNUR_ERR_DISTR_DATA}. 

   This call also sets the normalization constant such that the given
   PDF is the derivative of a given CDF, i.e. the volume is 1.
*/

double unur_distr_cvec_get_pdfvol( UNUR_DISTR *distribution );
/* 
   Get the volume below the PDF of the @var{distribution}. If this volume is
   not known,@* unur_distr_cont_upd_pdfarea() is called to compute
   it. If this is not successful @code{UNUR_INFINITY} is returned and
   @code{unur_errno} is set to @code{UNUR_ERR_DISTR_GET}.
*/

/* =END */

/*---------------------------------------------------------------------------*/

