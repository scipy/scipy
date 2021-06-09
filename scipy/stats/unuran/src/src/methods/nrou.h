/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: nrou.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method NROU                               *
 *         (Naive Ratio-Of-Uniforms method)                                  *
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
   =METHOD  NROU   Naive Ratio-Of-Uniforms method

   =UP  Methods_for_CONT

   =REQUIRED PDF

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =SPEED Set-up: slow or fast, Sampling: moderate

   =REINIT supported

   =REF  [HLD04: Sect.2.4 and Sect.6.4]

   =DESCRIPTION
      NROU is an implementation of the (generalized) ratio-of-uniforms
      method which uses (minimal) bounding rectangles, see
      @ref{Ratio-of-Uniforms}. It uses a positive control parameter
      @i{r} for adjusting the algorithm to the given distribution to
      improve performance and/or to make this method applicable.
      Larger values of @i{r} increase the class of distributions
      for which the method works at the expense of a higher rejection
      constant. For computational reasons @i{r=1} should be used if
      possible (this is the default).
      Moreover, this implementation uses the center @unurmath{\mu} of
      the distribution (see unur_distr_cont_get_center() for 
      details of its default values). 

      For the special case with @unurmath{r=1} the coordinates of the
      minimal bounding rectangles are given by 

      @unurmathdisplay{
      v^+ = \sup\limits_{x}         \sqrt{PDF(x)}, \\
      u^- = \inf\limits_{x} (x-\mu) \sqrt{PDF(x)}, \\
      u^+ = \sup\limits_{x} (x-\mu) \sqrt{PDF(x)}, }

      where @unurmath{\mu} is the center of the distribution.
      For other values of @i{r} we have

      @unurmathdisplay{
      v^+ = \sup\limits_{x}         (PDF(x))^{1/(r+1)}, \\
      u^- = \inf\limits_{x} (x-\mu) (PDF(x))^{r/(r+1)}, \\
      u^+ = \sup\limits_{x} (x-\mu) (PDF(x))^{r/(r+1)}. }

      These bounds can be given directly. Otherwise they are computed
      automatically by means of a (slow) numerical routine.
      Of course this routine can fail, especially when this rectangle
      is not bounded.

      It is important to note that the algorithm works with 
      @unurmath{PDF(x-\mu)} instead of @unurmath{PDF(x).}
      This is important as otherwise the acceptance region can become
      a very long and skinny ellipsoid along a diagonal of the (huge)
      bounding rectangle.

   =HOWTOUSE
      For using the NROU method UNU.RAN needs the PDF of the
      distribution. Additionally, the parameter @i{r} can be set via 
      a unur_vnrou_set_r() call. Notice that the acceptance
      probability decreases when @i{r} is increased. On the other
      hand is is more unlikely that the bounding rectangle does not
      exist if @i{r} is small.

      A bounding rectangle can be given by the
      unur_vnrou_set_u() and unur_vnrou_set_v() calls. 

      @emph{Important:} The bounding rectangle has to be
      provided for the function @unurmath{PDF(x-center)!} 
      Notice that @code{center} is the center of the given
      distribution, see unur_distr_cont_set_center().
      If in doubt or if this value is not optimal, it can be changed
      (overridden) by a unur_nrou_set_center() call.

      If the coordinates of the bounding rectangle are not provided by
      the user then the minimal bounding rectangle is computed
      automatically. 

      By means of unur_vnrou_set_verify() and unur_vnrou_chg_verify()
      one can run the sampling algorithm in a checking mode, i.e., in
      every cycle of the rejection loop it is checked whether the used
      rectangle indeed enclosed the acceptance region of the
      distribution. When in doubt (e.g., when it is not clear whether
      the numerical routine has worked correctly) this can be used to
      run a small Monte Carlo study.

      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.
      Notice, that derived parameters like the mode must also be (re-) set
      if the parameters or the domain has be changed.
      Notice, however, that then the values that has been set by 
      unur_vnrou_set_u() and unur_vnrou_set_v() calls are removed and 
      the coordinates of the bounding box are computed numerically.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_nrou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_nrou_set_u( UNUR_PAR *parameters, double umin, double umax );
/* 
   Sets left and right boundary of bounding rectangle.
   If no values are given, the boundary of the minimal bounding
   rectangle is computed numerically.
   
   @emph{Notice}: Computing the minimal bounding rectangle may fail
   under some circumstances. Moreover, for multimodal distributions
   the bounds might be too small as only local extrema are computed.
   Nevertheless, for @unurmath{T_c}-concave distributions with
   @unurmath{c=-1/2} it should work.

   @emph{Important:} The bounding rectangle that has to be
   provided is for the function @unurmath{PDF(x-center)!}

   Default: not set.
*/

int unur_nrou_set_v( UNUR_PAR *parameters, double vmax );
/* 
   Set upper boundary for bounding rectangle. If this value is not
   given then @unurmath{\sqrt{PDF(mode)}} is used instead.

   @emph{Notice}: When the mode is not given for the distribution
   object, then it will be computed numerically.

   Default: not set.
*/

int unur_nrou_set_r( UNUR_PAR *parameters, double r );
/*
   Sets the parameter @var{r} of the generalized ratio-of-uniforms 
   method.
   
   @emph{Notice}: This parameter must satisfy @var{r}>0.

   Default: @code{1}.
*/

int unur_nrou_set_center( UNUR_PAR *parameters, double center );
/* 
   Set the center (@unurmath{\mu}) of the PDF. 
   If not set the center of the given distribution object is used.

   Default: see unur_distr_cont_set_center().
*/

int unur_nrou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition @unurmath{PDF(x) \leq hat(x)} is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However, notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_nrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
