/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: ninv.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method NINV                               *
 *         (Numerical INVersion of cumulative distribution function)         *
 *                                                                           *
 *   USAGE:                                                                  *
 *         only included in unuran.h                                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
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
   =METHOD NINV Numerical INVersion

   =UP  Methods_for_CONT

   =REQUIRED CDF

   =OPTIONAL PDF

   =SPEED Set-up: optional, Sampling: (very) slow

   =REINIT supported

   =DESCRIPTION
      NINV implementations of some methods for numerical inversion: 
      Newton's method, regula falsi (combined with interval
      bisectioning), and bisection method.
      Regula falsi and bisection method require only the CDF while
      Newton's method also requires the PDF.
      To speed up marginal generation times a table with suitable
      starting points can be created during the setup. 
      The performance of the algorithm can adjusted by the desired
      accuracy of the method.
      It is possible to use this method for generating from truncated
      distributions. The truncated domain can be changed for an
      existing generator object.

   =HOWTOUSE
      Method NINV generates random variates by numerical
      inversion and requires a continuous univariate distribution
      objects with given CDF. Three variants are available:

      @itemize @minus
      @item Regula falsi  [default]
      @item Newton's method
      @item Interval bisectioning
      @end itemize

      Newton's method additionally requires the PDF of the
      distribution and cannot be used otherwise (NINV automatically
      switches to regula falsi then).
      Default algorithm is regula falsi. It is slightly slower but
      numerically much more stable than Newton's algorithm.
      Interval bisectioning is the slowest method and should only be
      considered as a last resort when the other methods fails.

      It is possible to draw samples from truncated distributions. 
      The truncated domain can even be changed for an existing generator
      object by an unur_ninv_chg_truncated() call.
      
      Marginal generation times can be sped up by means of a table
      with suitable starting points which can be created during the
      setup. Using such a table can be switched on by means of a
      unur_ninv_set_table() call where the table size is given as a
      parameter. The table is still useful when the (truncated) domain
      is changed often, since it is computed for the 
      domain of the given distribution. (It is not possible to enlarge
      this domain.) If it is necessary to recalculate the table during
      sampling, the command unur_ninv_chg_table() can be used.
      As a rule of thumb using such a table is appropriate when the
      number of generated points exceeds the table size by a factor of
      100.
  
      The default number of iterations of NINV should be enough for all
      reasonable cases. Nevertheless, it is possible to adjust the maximal
      number of iterations with the commands
      unur_ninv_set_max_iter() and unur_ninv_chg_max_iter().
      In particular this might be necessary when the PDF has a pole or
      the distribution has extremely heavy tails.

      It is also possible to set/change the accuracy of the method
      (which also heavily influencies the generation time).
      We use two measures for the approximation error which can be
      used independently: x-error and u-error
      (@pxref{Inversion} for more details).
      It is possible to set the maximal tolerated error using
      with unur_ninv_set_x_resolution() and 
      with unur_ninv_set_u_resolution(), resp., and change it with the
      respective calls unur_ninv_chg_x_resolution() and
      unur_ninv_chg_x_resolution().
      The algorithm tries to satisfy @emph{both} accuracy goals (and 
      raises an error flag it this fails).
      One of these accuracy checks can be disabled by setting the
      accuracy goal to a negative value.
      
      NINV tries to use proper starting values for both the regula
      falsi and bisection method, and for Newton's method. Of course
      the user might have more knowledge about the properties of the
      target distribution and is able to share his wisdom with NINV
      using the respective commands 
      unur_ninv_set_start() and unur_ninv_chg_start().
      
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.
      The values given by the last unur_ninv_chg_truncated() call will be 
      then changed to the values of the domain of the underlying distribution
      object. It is important to note that for a distribution from the 
      UNU.RAN library of standard distributions
      (@pxref{Stddist,,Standard distributions})
      the normalization constant has to be updated using the 
      unur_distr_cont_upd_pdfarea() call whenever its parameters have been
      changed by means of a unur_distr_cont_set_pdfparams() call.

      It might happen that NINV aborts unur_sample_cont() without
      computing the correct value (because the maximal number
      iterations has been exceeded). Then the last approximate value
      for @i{x} is returned (with might be fairly false) and
      @code{unur_error} is set to @code{UNUR_ERR_GEN_SAMPLING}.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/*
  =ROUTINES
*/


UNUR_PAR *unur_ninv_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_ninv_set_useregula( UNUR_PAR *parameters );
/* 
   Switch to regula falsi combined with interval bisectioning.
   (This the default.)
*/

int unur_ninv_set_usenewton( UNUR_PAR *parameters );
/* 
   Switch to Newton's method.
   Notice that it is numerically less stable than regula falsi.
   It it is not possible to invert the CDF for a particular uniform random
   number @i{U} when calling unur_sample_cont(), @code{unur_error} is set
   to @code{UNUR_ERR_GEN_SAMPLING}.
   Thus it is recommended to check @code{unur_error} before
   using the result of the sampling routine.
*/

int unur_ninv_set_usebisect( UNUR_PAR *parameters );
/* 
   Switch to bisection method. This is a slow algorithm and should
   only be used as a last resort.
*/

int unur_ninv_set_max_iter( UNUR_PAR *parameters, int max_iter );
/* */

int unur_ninv_chg_max_iter(UNUR_GEN *generator, int max_iter);
/* 
   Set and change number of maximal iterations.
   Default is @code{100}.
*/


int unur_ninv_set_x_resolution( UNUR_PAR *parameters, double x_resolution);
/* */

int unur_ninv_chg_x_resolution(UNUR_GEN *generator, double x_resolution);
/* 
   Set and change the maximal tolerated relative x-error.
   If @var{x_resolution} is negative then checking of the x-error is
   disabled.

   Default is @code{1.e-8}.
*/

int unur_ninv_set_u_resolution( UNUR_PAR *parameters, double u_resolution);
/* */

int unur_ninv_chg_u_resolution(UNUR_GEN *generator, double u_resolution);
/*
   Set and change the maximal tolerated (abolute) u-error.
   If @var{u_resolution} is negative then checking of the u-error is
   disabled.

   Default is @code{-1} (disabled).
*/

int unur_ninv_set_start( UNUR_PAR *parameters, double left, double right);
/* 
   Set starting points.
   If not set, suitable values are chosen automatically.

   @multitable @columnfractions 0.2 0.2 0.6
   @item Newton:       @tab @var{left}:              @tab starting point
   @item Regula falsi: @tab @var{left}, @var{right}: @tab boundary of starting interval
   @end multitable

   If the starting points are not set then the follwing points are used by
   default:
   @multitable @columnfractions 0.2 0.2 0.6
   @item Newton:       @tab @var{left}:  @tab CDF(@var{left}) = 0.5
   @item Regula falsi: @tab @var{left}:  @tab CDF(@var{left}) = 0.1
   @item               @tab @var{right}: @tab CDF(@var{right}) = 0.9
   @end multitable

   If @var{left} == @var{right}, then  UNU.RAN always uses the default
   starting points!
*/
    
int unur_ninv_chg_start(UNUR_GEN *gen, double left, double right);
/* 
   Change the starting points for numerical inversion. 
   If left==right, then UNU.RAN uses the default starting points 
   (see unur_ninv_set_start()).
*/

int unur_ninv_set_table(UNUR_PAR *parameters, int no_of_points);
/* 
   Generates a table with @var{no_of_points} points containing
   suitable starting values for the iteration. The value of
   @var{no_of_points} must be at least 10 (otherwise it will be set
   to 10 automatically).

   The table points are chosen such that the CDF at these points
   form an equidistance sequence in the interval (0,1).

   If a table is used, then the starting points given by
   unur_ninv_set_start() are ignored.

   No table is used by default.
 */

int unur_ninv_chg_table(UNUR_GEN *gen, int no_of_points);
/*
   Recomputes a table as described in unur_ninv_set_table().
*/


int unur_ninv_chg_truncated(UNUR_GEN *gen, double left, double right);
/*
   Changes the borders of the domain of the (truncated) distribution. 

   Notice that the given truncated domain must be a subset of the
   domain of the given distribution. The generator always uses the
   intersection of the domain of the distribution and the truncated
   domain given by this call.
   Moreover the starting point(s) will not be changed.

   @emph{Important:} If the CDF is (almost) the same for @var{left} and 
   @var{right} and (almost) equal to @code{0} or @code{1}, then the truncated 
   domain is @emph{not} chanced and the call returns an error code.

   @emph{Notice:} If the parameters of the distribution has been changed by a 
   unur_distr_cont_set_pdfparams() call it is recommended to set the
   truncated domain again, since the former call might change the
   domain of the distribution but not update the values for the
   boundaries of the truncated distribution.
*/

double unur_ninv_eval_approxinvcdf( const UNUR_GEN *generator, double u );
/*
   Get approximate approximate value of inverse CDF at @var{u}.
   If @var{u} is out of the domain [0,1] then @code{unur_errno} is set
   to @code{UNUR_ERR_DOMAIN} and the respective bound of
   the domain of the distribution are returned (which is
   @code{-UNUR_INFINITY} or @code{UNUR_INFINITY} in the case of
   unbounded domains).

   @emph{Notice}: This function always evaluates the inverse CDF of
   the given distribution. A call to unur_ninv_chg_truncated() call
   has no effect.
*/

/* =END */
/*---------------------------------------------------------------------------*/
