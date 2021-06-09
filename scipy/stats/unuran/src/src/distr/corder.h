/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: corder.h                                                          *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for manipulating distribution objects of      *
 *         id   CORDER (continuous order statistics)                         *
 *         type CONT   (continuous univariate distribution)                  *
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
   =NODEX   CORDER   Continuous univariate order statistics

   =UP Distribution_objects [15]

   =DESCRIPTION
      These are special cases of a continuous univariate distributions
      and thus they have most of these parameters (with the exception
      that functions cannot be changed). Additionally,

      @itemize @minus
      @item there is a call to extract the underlying distribution,

      @item and a call to handle the @command{rank} of the order
      statistics.

      @end itemize

   =END
*/

/*---------------------------------------------------------------------------*/

/* 
   Routines for handling univariate continuous order statistics (CORDER).
*/

/* =ROUTINES */

UNUR_DISTR *unur_distr_corder_new( const UNUR_DISTR *distribution, int n, int k );
/* 
   Create an object for order statistics of sample size
   @var{n} and rank @var{k}.
   @var{distribution} must be a pointer to a univariate continuous
   distribution. 
   The resulting generator object is of the same type as of a
   unur_distr_cont_new() call.
   (However, it cannot be used to make an order statistics out of an
   order statistics.)

   To have a PDF for the order statistics, the given distribution
   object must contain a CDF and a PDF. Moreover, it is assumed that
   the given PDF is the derivative of the given CDF. Otherwise the
   area below the PDF of the order statistics is not computed correctly.

   @emph{Important:} There is no warning when the computed area below
   the PDF of the order statistics is wrong.
*/


const UNUR_DISTR *unur_distr_corder_get_distribution( const UNUR_DISTR *distribution );
/* 
   Get pointer to distribution object for underlying distribution.
*/



/* ==DOC
   @subsubheading Essential parameters
*/

int unur_distr_corder_set_rank( UNUR_DISTR *distribution, int n, int k );
/* 
   Change sample size @var{n} and rank @var{k} of order statistics.
   In case of invalid data, no parameters are changed.
   The area below the PDF can be set to that of the underlying
   distribution by a unur_distr_corder_upd_pdfarea() call.
*/

int unur_distr_corder_get_rank( const UNUR_DISTR *distribution, int *n, int *k );
/* 
   Get sample size @var{n} and rank @var{k} of order statistics.
   In case of error an error code is returned.
*/

/* ==DOC
   Additionally most of the set and get calls for continuous
   univariate distributions work. The most important exceptions are
   that the PDF and CDF cannot be changed and
   unur_distr_cont_upd_mode() uses in any way a (slow) numerical
   method that might fail.
*/


#define unur_distr_corder_get_pdf(distr)   unur_distr_cont_get_pdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_corder_get_pdf( UNUR_DISTR *distribution ); */

#define unur_distr_corder_get_dpdf(distr)  unur_distr_cont_get_dpdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_corder_get_dpdf( UNUR_DISTR *distribution ); */

#define unur_distr_corder_get_cdf(distr)   unur_distr_cont_get_cdf((distr))
/*  UNUR_FUNCT_CONT *unur_distr_corder_get_cdf( UNUR_DISTR *distribution ); */
/* 
   Get the respective pointer to the PDF, the derivative of the 
   PDF and the CDF of the distribution, respectively. The pointer is of type
   @code{double funct(double x, UNUR_DISTR *distr)}.
   If the corresponding function is not available for the distribution,
   the NULL pointer is returned.
   See also unur_distr_cont_get_pdf().
   (Macro)
*/

#define unur_distr_corder_eval_pdf(x,distr)  unur_distr_cont_eval_pdf((x),(distr))
/*  double unur_distr_corder_eval_pdf( double x, UNUR_DISTR *distribution ); */

#define unur_distr_corder_eval_dpdf(x,distr) unur_distr_cont_eval_dpdf((x),(distr))
/*  double unur_distr_corder_eval_dpdf( double x, UNUR_DISTR *distribution ); */

#define unur_distr_corder_eval_cdf(x,distr)  unur_distr_cont_eval_cdf((x),(distr))
/*  double unur_distr_corder_eval_cdf( double x, UNUR_DISTR *distribution ); */
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


#define unur_distr_corder_set_pdfparams(distr,params,n)  unur_distr_cont_set_pdfparams((distr),(params),(n))
/*  int unur_distr_corder_set_pdfparams(UNUR_DISTR *distribution,double *params,int n_params); */
/* 
   Set array of parameters for underlying distribution.
   See unur_distr_cont_set_pdfparams() for details.
   (Macro)
*/

#define unur_distr_corder_get_pdfparams(distr,params)  unur_distr_cont_get_pdfparams((distr),(params))
/*  int unur_distr_corder_get_pdfparams( UNUR_DISTR *distribution, double **params ); */
/* 
   Get number of parameters of the PDF of the underlying distribution
   and set pointer @var{params} to array of parameters. 
   See unur_distr_cont_get_pdfparams() for details.
   (Macro)
*/


#define unur_distr_corder_set_domain(distr,left,right)  unur_distr_cont_set_domain((distr),(left),(right))
/*  int unur_distr_corder_set_domain( UNUR_DISTR *distribution, double left, double right ); */
/* 
   Set the left and right borders of the domain of the
   distribution. 
   See unur_distr_cont_set_domain() for details.
   (Macro)
*/

#define unur_distr_corder_get_domain(distr,left,right)  unur_distr_cont_get_domain((distr),(left),(right))
/*  int unur_distr_corder_get_domain( UNUR_DISTR *distribution, double *left, double *right ); */
/* 
   Get the left and right borders of the domain of the
   distribution. 
   See unur_distr_cont_get_domain() for details.
   (Macro)
*/


#define unur_distr_corder_get_truncated(distr,left,right)  unur_distr_cont_get_truncated((distr),(left),(right))
/*  int unur_distr_corder_get_truncated( UNUR_DISTR *distribution, double *left, double *right ); */
/* 
   Get the left and right borders of the (truncated) domain of the
   distribution.
   See unur_distr_cont_get_truncated() for details.
   (Macro)
*/

/* ==DOC
   @subsubheading Derived parameters

   The following paramters @strong{must} be set whenever one of the essential
   parameters has been set or changed (and the parameter is required
   for the chosen method).
*/

#define unur_distr_corder_set_mode(distr,mode)   unur_distr_cont_set_mode((distr),(mode))
/*  int unur_distr_corder_set_mode( UNUR_DISTR *distribution, double mode ); */
/* 
   Set mode of distribution. 
   See also unur_distr_corder_set_mode().
   (Macro)
*/

#define unur_distr_corder_upd_mode(distr)   unur_distr_cont_upd_mode((distr))
/*  double unur_distr_corder_upd_mode( UNUR_DISTR *distribution ); */
/* 
   Recompute the mode of the distribution numerically. Notice that
   this routine is slow and might not work properly in every case.
   See also unur_distr_cont_upd_mode() for further details.
   (Macro)
*/

#define unur_distr_corder_get_mode(distr)   unur_distr_cont_get_mode((distr))
/*  double unur_distr_corder_get_mode( UNUR_DISTR *distribution ); */
/* 
   Get mode of distribution.
   See unur_distr_cont_get_mode() for details.
   (Macro)
*/


#define unur_distr_corder_set_pdfarea(distr,area)   unur_distr_cont_set_pdfarea((distr),(area))
/*  int unur_distr_corder_set_pdfarea( UNUR_DISTR *distribution, double area ); */
/* 
   Set the area below the PDF.
   See unur_distr_cont_set_pdfarea() for details.
   (Macro)
*/


#define unur_distr_corder_upd_pdfarea(distr)   unur_distr_cont_upd_pdfarea((distr))
/*  double unur_distr_corder_upd_pdfarea( UNUR_DISTR *distribution ); */
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

#define unur_distr_corder_get_pdfarea(distr)   unur_distr_cont_get_pdfarea((distr))
/*  double unur_distr_corder_get_pdfarea( UNUR_DISTR *distribution ); */
/* 
   Get the area below the PDF of the distribution.
   See unur_distr_cont_get_pdfarea() for details.
   (Macro)
*/

/* =END */

/*---------------------------------------------------------------------------*/
