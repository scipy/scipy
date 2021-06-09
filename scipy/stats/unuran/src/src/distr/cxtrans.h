/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: cxtrans.h                                                         *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         id   CXTRANS  (continuous distribution of transformed RV)         *
 *         type CONT     (continuous univariate distribution)                *
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
   =ExperimentalNODEX   CXTRANS   Continuous univariate distribution of transformed random variables

   =UP Distribution_objects [13]

   =DESCRIPTION
      Given a continuous univariate distribution object for a random
      variable @unurmath{X} with CDF @unurmath{F} and PDF @unurmath{f.}
      Then it can be transformed into a continuous univariate
      distribution object for @unurmath{Z=\phi(X)} for some strictly
      monotonically increasing continuous transformation
      @unurmath{\phi.} 
      Its respetive CDF and PDF are then given by

      @unurmathdisplay{G(z) = F(\phi^{-1}(z))}

      and 

      @unurmathdisplay{g(z) = f(\phi^{-1}(z))\cdot(\phi^{-1})'(z)}

      Currently, only power transformations, logarithmic and 
      exponential transformation are supported. Additionally,
      The random variable can be relocated and rescaled
      ("standardized"):

      @unurmathdisplay{
      \phi(X) = \log(Z) \\
      \phi(X) = sign(Z) |Z|^\alpha \\
      \phi(X) = \exp(Z) }

      where 

      @unurmathdisplay{
      Z = (X-\mu)/\sigma }

      @unurmath{\alpha} must be between 0 and @unurmath{\infty} where
      the logarithmic transformation is implemented as special case
      @unurmath{\alpha=0}, the exponential transformation as
      @unurmath{\alpha=\infty} and the standardization procedure as
      @unurmath{\alpha=1.}
      Moreover, @unurmath{\sigma} must be strictly positive.

      One of the most important applications of these transformations
      is to eliminate poles, i.e., where the PDF is not bounded from
      above. The implementation of a PDF must return
      @code{UNUR_INFINITY} at such points.
      By default the value of the PDF @i{g} of the transformed density
      is set to @code{0} for such points and its derivative dPDF to
      @code{UNUR_INFINITY}. These figures can be changed to more
      appropriate values. It must be again noted here that poles are
      (only) detected/indicated as points where the PDF @i{f} returns
      @code{UNUR_INFINITY}. 

      Notice that for @unurmath{\alpha > 1} these transformations are
      not defined at @unurmath{\mu.} It is difficult to check whether
      there exists a continuous continuation for the density of the
      transformed random variate at 0. Thus the above value for poles
      is used at 0.

      The distribution objects of transformed random variates are
      special cases of a continuous univariate distributions and thus
      they have most of these parameters (with the exception that
      functions cannot be changed). Additionally, 

      @itemize @minus
      @item there is a call to extract the underlying distribution,

      @item calls to handle the transformation.

      @end itemize

      Notice that by transformation @unurmath{\phi} the mode is
      usually not known. The domain of the transformed random variate
      is adjusted according to the transformation.

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous distributions of transformed
   random variables
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_cxtrans_new( const UNUR_DISTR *distribution );
/* 
   Create an object for the distribution of the transformed random
   variate. @var{distribution} must be a pointer to a univariate
   continuous distribution. 
   The resulting generator object is of the same type as of a
   unur_distr_cont_new() call.

   To select a transformation the unur_distr_cxtrans_set_alpha() call
   must be used. For the parameters for relocating and rescaling of the
   random variate unur_distr_cxtrans_set_rescale() has to be used.

   Without one of these additional calls the identity transformation
   is used, i.e. the random variate is not transformed at all.
*/

const UNUR_DISTR *unur_distr_cxtrans_get_distribution( const UNUR_DISTR *distribution );
/* 
   Get pointer to distribution object for underlying distribution.
*/



/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_cxtrans_set_alpha( UNUR_DISTR *distribution, double alpha );
/* 
   Change the parameter @var{alpha} for the power transformation. 
   One of the following transformations are possible:
   
   @multitable @columnfractions .25 .75 
   @headitem @var{alpha} @tab Transformation
   @item @code{UNUR_INFINITY}
   @tab  @unurmath{\phi(X)=\exp(Z)}
   @item @code{0}
   @tab  @unurmath{\phi(X)=\log(Z)}
   @item positive
   @tab  @unurmath{\phi(X) = sign(Z) |Z|^\alpha}
   @end multitable
   
   where @unurmath{Z = (X-\mu)/\sigma.}
   Negative values for @var{alpha} are not allowed.

   The relocating and rescaling parameters @unurmath{\mu} and
   @unurmath{\sigma} can be set by the respective call
   unur_distr_cxtrans_set_rescale().

   Default: @code{1} (i.e., identity).
*/

int unur_distr_cxtrans_set_rescale( UNUR_DISTR *distribution, double mu, double sigma );
/* 
   Change relocating and rescaling parameter. @var{sigma} must be
   strictly positive.

   Default: @var{mu} = @code{0.} and @var{sigma} = @code{1.} (i.e., identity)
 */

double unur_distr_cxtrans_get_alpha( const UNUR_DISTR *distribution );
/* */

double unur_distr_cxtrans_get_mu( const UNUR_DISTR *distribution );
/* */

double unur_distr_cxtrans_get_sigma( const UNUR_DISTR *distribution );
/* 
   Get parameters @unurmath{\alpha,} @unurmath{\mu,} and
   @unurmath{\sigma} for the power transformation, relocating and
   rescaling of the random variate.
*/

int unur_distr_cxtrans_set_logpdfpole( UNUR_DISTR *distribution, double logpdfpole, double dlogpdfpole );
/* 
   Set value for logarithm of PDF and its derivative that are used
   whenever a pole of the underlying distribution is detected (i.e.,
   when the PDF of the underlying distribution returns UNUR_INFINITY).

   Default: @code{-UNUR_INFINITY} and @code{UNUR_INFINITY}, respectively.
*/


/* ==DOC
   Additionally most of the set and get calls for continuous
   univariate distributions work. The most important exceptions are
   that the PDF and CDF cannot be changed and
   unur_distr_cont_upd_mode() uses in any way a (slow) numerical
   method that might fail.
*/


#define unur_distr_cxtrans_get_pdf(distr)   unur_distr_cont_get_pdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_cxtrans_get_pdf( UNUR_DISTR *distribution ); */

#define unur_distr_cxtrans_get_dpdf(distr)  unur_distr_cont_get_dpdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_cxtrans_get_dpdf( UNUR_DISTR *distribution ); */

#define unur_distr_cxtrans_get_cdf(distr)   unur_distr_cont_get_cdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_cxtrans_get_cdf( UNUR_DISTR *distribution ); */
/* 
   Get the respective pointer to the PDF, the derivative of the 
   PDF and the CDF of the distribution, respectively. The pointer is of type
   @code{double funct(double x, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
   See also unur_distr_cont_get_pdf().
   (Macro)
*/

#define unur_distr_cxtrans_eval_pdf(x,distr)  unur_distr_cont_eval_pdf((x),(distr))
/*  double unur_distr_cxtrans_eval_pdf( double x, UNUR_DISTR *distribution ); */

#define unur_distr_cxtrans_eval_dpdf(x,distr) unur_distr_cont_eval_dpdf((x),(distr))
/*  double unur_distr_cxtrans_eval_dpdf( double x, UNUR_DISTR *distribution ); */

#define unur_distr_cxtrans_eval_cdf(x,distr)  unur_distr_cont_eval_cdf((x),(distr))
/*  double unur_distr_cxtrans_eval_cdf( double x, UNUR_DISTR *distribution ); */
/* 
   Evaluate the PDF, derivative of the PDF. and the CDF,
   respectively, at @var{x}.
   Notice that @var{distribution} must not be the NULL pointer.
   If the corresponding function is not available for the distribution,
   @code{UNUR_INFINITY} is returned and @code{unur_errno} is set to
   @code{UNUR_ERR_DISTR_DATA}.
   See also unur_distr_cont_eval_pdf().
   (Macro)

   @emph{IMPORTANT:}
   In the case of a truncated standard distribution these calls always
   return the respective values of the @emph{untruncated} distribution!
*/

int unur_distr_cxtrans_set_domain( UNUR_DISTR *distribution, double left, double right );
/* 
   Set the left and right borders of the domain of the
   distribution. When the exponential transformation is used
   (i.e. @var{alpha} is set to UNUR_INFINITY), then the domain must be
   a subset of positive numbers.
   See unur_distr_cont_set_domain() for more details.

   @emph{Important:}
   The domain of the transformed random variate may be different
   from the domain of the underlying random variate.
*/

#define unur_distr_cxtrans_get_domain(distr,left,right)  unur_distr_cont_get_domain((distr),(left),(right))
/*  int unur_distr_cxtrans_get_domain( UNUR_DISTR *distribution, double *left, double *right ); */
/* 
   Get the left and right borders of the domain of the
   distribution. 
   See unur_distr_cont_get_domain() for details.
   (Macro)
*/


#define unur_distr_cxtrans_get_truncated(distr,left,right)  unur_distr_cont_get_truncated((distr),(left),(right))
/*  int unur_distr_cxtrans_get_truncated( UNUR_DISTR *distribution, double *left, double *right ); */
/* 
   Get the left and right borders of the (truncated) domain of the
   distribution.
   See unur_distr_cont_get_truncated() for details.
   (Macro)

   @emph{Important:}
   The domain of the transformed random variate may be different
   from the domain of the underlying random variate.
*/

/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the essential
   parameters has been set or changed (and the parameter is required
   for the chosen method).
*/

#define unur_distr_cxtrans_set_mode(distr,mode)   unur_distr_cont_set_mode((distr),(mode))
/*  int unur_distr_cxtrans_set_mode( UNUR_DISTR *distribution, double mode ); */
/* 
   Set mode of distribution. 
   See also unur_distr_cxtrans_set_mode().
   (Macro)
*/

#define unur_distr_cxtrans_upd_mode(distr)   unur_distr_cont_upd_mode((distr))
/*  double unur_distr_cxtrans_upd_mode( UNUR_DISTR *distribution ); */
/* 
   Recompute the mode of the distribution numerically. Notice that
   this routine is slow and might not work properly in every case.
   See also unur_distr_cont_upd_mode() for further details.
   (Macro)
*/

#define unur_distr_cxtrans_get_mode(distr)   unur_distr_cont_get_mode((distr))
/*  double unur_distr_cxtrans_get_mode( UNUR_DISTR *distribution ); */
/* 
   Get mode of distribution.
   See unur_distr_cont_get_mode() for details.
   (Macro)
*/


#define unur_distr_cxtrans_set_pdfarea(distr,area)   unur_distr_cont_set_pdfarea((distr),(area))
/*  int unur_distr_cxtrans_set_pdfarea( UNUR_DISTR *distribution, double area ); */
/* 
   Set the area below the PDF.
   See unur_distr_cont_set_pdfarea() for details.
   (Macro)
*/


#define unur_distr_cxtrans_upd_pdfarea(distr)   unur_distr_cont_upd_pdfarea((distr))
/*  double unur_distr_cxtrans_upd_pdfarea( UNUR_DISTR *distribution ); */
/*
   Recompute the area below the PDF of the distribution. 
   It only works for order statistics for distribution objects from
   the UNU.RAN library of standard distributions when the
   corresponding function is available.
   unur_distr_cont_upd_pdfarea() assumes that the PDF of the underlying
   distribution is normalized, i.e. it is the derivative of its CDF.
   Otherwise the computed area is wrong and there is @strong{no} warning
   about this failure.
   See unur_distr_cont_upd_pdfarea() for further details.
   (Macro)
*/

#define unur_distr_cxtrans_get_pdfarea(distr)   unur_distr_cont_get_pdfarea((distr))
/*  double unur_distr_cxtrans_get_pdfarea( UNUR_DISTR *distribution ); */
/* 
   Get the area below the PDF of the distribution.
   See unur_distr_cont_get_pdfarea() for details.
   (Macro)
*/

/* =END */

/*---------------------------------------------------------------------------*/
