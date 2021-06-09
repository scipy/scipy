/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: vnrou.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method VNROU                              *
 *         (Vector Naive Ratio-Of-Uniforms method)                           *
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
   =METHOD  VNROU   Multivariate Naive Ratio-Of-Uniforms method

   =UP  Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =SPEED Set-up: fast or slow, Sampling: slow

   =REINIT supported

   =REF  [WGS91]

   =DESCRIPTION 
      VNROU is an implementation of the multivariate
      ratio-of-uniforms method which uses a (minimal) bounding
      hyper-rectangle, see also @ref{Ratio-of-Uniforms}.  It uses an
      additional parameter @i{r} that can be used for adjusting the
      algorithm to the given distribution to improve performance
      and/or to make this method applicable.  Larger values of
      @i{r} increase the class of distributions for which the
      method works at the expense of higher rejection
      constants. Moreover, this implementation uses the center 
      @unurmath{\mu} of the distribution (which is set to the mode or
      mean by default, see unur_distr_cvec_get_center() for details of
      its default values).

      The minimal bounding has then the coordinates

      @unurmathdisplay{
      v^+   = \sup\limits_{x}               (f(x))^{1/r\,d+1}, \\
      u^-_i = \inf\limits_{x_i} (x_i-\mu_i) (f(x))^{r/r\,d+1}, \\
      u^+_i = \sup\limits_{x_i} (x_i-\mu_i) (f(x))^{r/r\,d+1}, }

      where @unurmath{x_i} is the @i{i}-th coordinate of point @i{x};
      @unurmath{\mu_i} is the @i{i}-th coordinate of the center
      @unurmath{\mu.} 
      @unurmath{d} denotes the dimension of the distribution.
      These bounds can either be given directly, or are computed
      automatically by means of an numerical routine 
      by Hooke and Jeeves @unurbibref{HJa61} called direct search 
      (see @file{src/utils/hooke.c} for further references and
      details). Of course this algorithm can fail, especially when
      this rectangle is not bounded.

      It is important to note that the algorithm works with 
      @unurmath{PDF(x-center)} instead of 
      @unurmath{PDF(x),} i.e. the bounding rectangle has to be
      provided for @unurmath{PDF(x-center).}
      This is important as otherwise the acceptance region can become
      a very long and skinny ellipsoid along a diagonal of the (huge)
      bounding rectangle.

      VNROU is based on the rejection method (@pxref{Rejection}),
      and it is important to note that the acceptance probability
      decreases exponentially with dimension. Thus even for moderately
      many dimensions (e.g. 5) the number of repetitions to get one
      random vector can be prohibitively large and the algorithm seems
      to stay in an infinite loop.

   =HOWTOUSE
      For using the VNROU method UNU.RAN needs the PDF of the
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
      distribution, see unur_distr_cvec_set_center().
      If in doubt or if this value is not optimal, it can be changed
      (overridden) by a unur_distr_cvec_set_center() call.

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

      @strong{Important:}
      The rejection constant (i.e. the expected number of iterations
      for generationg one random vector) can be extremely high, in
      particular when the dimension is 4 or higher. 
      Then the algorithm will perform almost infinite loops.
      Thus it is recommended to read the volume below the hat function
      by means of the unur_vnrou_get_volumehat() call. The returned
      number divided by the volume below the PDF (which is 1 in case
      of a normalized PDF) gives the rejection constant.

      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.
      Notice, that the coordinates of a bounding rectangle given by
      unur_vnrou_set_u() and unur_vnrou_set_v() calls are used also
      when the generator is reused. These can be changed by means of
      unur_vnrou_chg_u() and unur_vnrou_chg_v() calls.
      (If no such coordinates have been given, then they are computed
      numerically during the reinitialization proceedure.) 

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_vnrou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_vnrou_set_u( UNUR_PAR *parameters, double *umin, double *umax );
/* 
   Sets left and right boundaries of bounding hyper-rectangle.
   If no values are given, the boundary of the minimal bounding
   hyper-rectangle is computed numerically.
  
   @strong{Important}: The boundaries are those of the density shifted
   by the center of the distribution, i.e., for the 
   function @unurmath{PDF(x-center)!}

   @emph{Notice}: Computing the minimal bounding rectangle may fail
   under some circumstances. Moreover, for multimodal distributions
   the bounds might be too small as only local extrema are computed.
   Nevertheless, for log-concave distributions it should work.

   Default: not set (i.e. computed automatically)
*/

int unur_vnrou_chg_u( UNUR_GEN *generator, double *umin, double *umax );
/* 
   Change left and right boundaries of bounding hyper-rectangle.
*/

int unur_vnrou_set_v( UNUR_PAR *parameters, double vmax );
/* 
   Set upper boundary for bounding hyper-rectangle. 
   If no values are given, the density at the mode is evaluated.
   If no mode is given for the distribution it is computed
   numerically (and might fail).
  
   Default: not set (i.e. computed automatically)
*/

int unur_vnrou_chg_v( UNUR_GEN *generator, double vmax );
/*
   Change upper boundary for bounding hyper-rectangle. 
*/

int unur_vnrou_set_r( UNUR_PAR *parameters, double r );
/* 
   Sets the parameter @var{r} of the generalized multivariate 
   ratio-of-uniforms method.

   @emph{Notice}: This parameter must satisfy @var{r}>0. 

   Default: @code{1}.
*/

int unur_vnrou_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_vnrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

double unur_vnrou_get_volumehat( const UNUR_GEN *generator );
/* 
   Get the volume of below the hat.
   For normalized densities, i.e. when the volume below PDF is 1, 
   this value equals the rejection constant for the vnrou method.

   In case of an error UNUR_INFINITY is returned.
*/


/* =END */
/*---------------------------------------------------------------------------*/
