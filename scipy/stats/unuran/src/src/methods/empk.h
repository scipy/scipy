/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: empk.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method EMPK                               *
 *         (EMPirical distribution with Kernel smoothing)                    *
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
   =METHOD  EMPK   EMPirical distribution with Kernel smoothing

   =UP  Methods_for_CEMP

   =REQUIRED observed sample

   =SPEED Set-up: slow (as sample is sorted),
          Sampling: fast (depends on kernel)

   =REINIT not implemented

   =REF  [HLa00] [HLD04: Sect.12.1.2]

   =DESCRIPTION
      EMPK generates random variates from an empirical distribution that is
      given by an observed sample. The idea is that simply choosing a random
      point from the sample and to return it with some added noise results
      in a method that has very nice properties, as it can be seen as sampling
      from a kernel density estimate. If the underlying distribution is
      continuous, especially the fine structur of the resulting empirical
      distribution is much better than using only resampling without noise.

      Clearly we have to decide about the density of the noise (called kernel)
      and about the standard deviation of the noise.
      The mathematical theory of kernel density estimation shows us that we
      are comparatively free in choosing the kernel. It also supplies us with
      a simple formula to compute the optimal standarddeviation of the noise,
      called bandwidth (or window width) of the kernel.

      The variance of the estimated density is slightly larger than
      that of the observed sample. However, this can be easily
      corrected if required.
      
      There is also a correction (mirroring technique) for
      distributions with non-negative support.

      A simple robust reference method is implemented to find a good
      standard deviation of the noise (i.e. the bandwidth of 
      kernel density estimation). For some cases 
      (e.g. densities with two or more sharp distinct peaks) there
      kernel density estimation can be adjusted by changing the
      smoothness factor and the so called beta factor.

   =HOWTOUSE
      EMPK uses empirical distributions. The main parameter is the
      choice if of kernel density. The most important kernels can be 
      set by unur_empk_set_kernel(). Additionally generators for other
      kernels can be used by using unur_empk_set_kernelgen() instead.
      Additionally variance correction and a correction for
      non-negative variates can be switched on.

      The two other parameters (smoothing factor and beta factor) are
      only useful for people knowing the theory of kernel density
      estimation. It is not necessary to change them if 
      the true underlying distribution is somehow comparable with a 
      bell-shaped curve, even skewed or with some not too sharp extra peaks.
      In all these cases the simple robust reference method implemented to 
      find a good standard deviation of the noise (i.e. the bandwidth of
      kernel density estimation) should give sensible results.
      However, it might be necessary to overwrite this automatic method
      to find the bandwidth eg. when resampling from data with
      two or more sharp distinct peaks. Then the distribution has nearly
      discrete components as well and our automatic method may
      easily choose too large a bandwidth which results in an
      empirical distribution which is oversmoothed (i.e. it has
      lower peaks than the original distribution). Then it
      is recommended to decrease the bandwidth using the 
      unur_empk_set_smoothing() call. A smoothing factor of @code{1}
      is the default. A smoothing factor of @code{0} leads to naive
      resampling of the data. Thus an appropriate value between these
      extremes should be choosen. We recommend to consult a reference
      on kernel smoothing when doing so; but it is not a simple problem
      to determine an optimal bandwidth for distributions with sharp peaks.

      In general, for most applications it is perfectly ok to use the
      default values offered. Unless you have some knowledge on
      density estimation we do not recommend to change anything. 
      There are two exceptions:

      @enumerate A
      @item
      In the case that the unknown underlying distribution is not continuous
      but discrete you should "turn off" the adding of the noise by setting:
      @smallexample
      unur_empk_set_smoothing(par, 0.)
      @end smallexample

      @item
      In the case that you are especially 
      interested in a fast sampling algorithm use the call
      @smallexample
      unur_empk_set_kernel(par, UNUR_DISTR_BOXCAR);
      @end smallexample
      to change the used noise distribution from the default Gaussian
      distribution to the uniform distribution.
      @end enumerate

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_empk_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_empk_set_kernel( UNUR_PAR *parameters, unsigned kernel);
/* 
   Select one of the supported kernel distributions. Currently the following 
   kernels are supported:

   @table @code
   @item UNUR_DISTR_GAUSSIAN
   Gaussian (normal) kernel
   @item UNUR_DISTR_EPANECHNIKOV
   Epanechnikov kernel
   @item UNUR_DISTR_BOXCAR
   Boxcar (uniform, rectangular) kernel
   @item UNUR_DISTR_STUDENT
   t3 kernel (Student's distribution with 3 degrees of freedom)
   @item UNUR_DISTR_LOGISTIC
   logistic kernel
   @end table

   For other kernels (including kernels with Student's distribution
   with other than 3 degrees of freedom) use the
   unur_empk_set_kernelgen() call.

   It is not possible to call unur_empk_set_kernel() twice.

   Default is the Gaussian kernel.
*/

int unur_empk_set_kernelgen( UNUR_PAR *parameters, const UNUR_GEN *kernelgen, double alpha, double kernelvar );
/* 
   Set generator for the kernel used for density estimation.

   @var{alpha} is used to compute the optimal bandwidth from the point of
   view of minimizing the mean integrated square error (MISE).
   It depends on the kernel K and is given by
   @smallexample
     alpha(K) = Var(K)^(-2/5)@{ \int K(t)^2 dt@}^(1/5)
   @end smallexample
   For standard kernels (see above) alpha is computed by the algorithm.

   @var{kernvar} is the variance of the used kernel. It is only required 
   for the variance corrected version of density estimation (which is 
   used by default); otherwise it is ignored.
   If @var{kernelvar} is nonpositive, variance correction is disabled.
   For standard kernels (see above) @var{kernvar} is computed by the
   algorithm.

   It is not possible to call unur_empk_set_kernelgen() after a standard kernel
   has been selected by a unur_empk_set_kernel() call.

   Notice that the uniform random number generator of the kernel
   generator is overwritten during the unur_init() call and at each
   unur_chg_urng() call with the uniform generator used for the empirical
   distribution.

   Default is the Gaussian kernel.
*/

int unur_empk_set_beta( UNUR_PAR *parameters, double beta );
/* 
   @var{beta} is used to compute the optimal bandwidth from the point
   of view of minimizing the mean integrated square error (MISE).
   @var{beta} depends on the (unknown) distribution of the sampled data
   points. 
   By default Gaussian distribution is assumed for the sample
   (@var{beta} = 1.3637439). There is no requirement to change
   @var{beta}.

   Default: @code{1.3637439}
*/

int unur_empk_set_smoothing( UNUR_PAR *parameters, double smoothing );
/* */

int unur_empk_chg_smoothing( UNUR_GEN *generator, double smoothing );
/* 
   Set and change the smoothing factor.
   The smoothing factor controlles how ``smooth'' the resulting density
   estimation will be. A smoothing factor equal to @code{0} results in naive
   resampling. A very large smoothing factor (together with the
   variance correction) results in a density which is approximately
   equal to the kernel.
   Default is 1 which results in a smoothing parameter minimising
   the MISE (mean integrated squared error) if the data are not too
   far away from normal. If a large smoothing factor is used, then
   variance correction must be switched on.

   Default: @code{1}
*/

int unur_empk_set_varcor( UNUR_PAR *parameters, int varcor );
/* */

int unur_empk_chg_varcor( UNUR_GEN *generator, int varcor );
/* 
   Switch variance correction in generator on/off.
   If @var{varcor} is TRUE then the variance of the used
   density estimation is the same as the sample variance. However this 
   increases the MISE of the estimation a little bit.

   Default is FALSE.
*/

int unur_empk_set_positive( UNUR_PAR *parameters, int positive );
/* 
   If @var{positive} is TRUE then only nonnegative random variates are
   generated. This is done by means of a mirroring technique.
   
   Default is FALSE.
*/

/* =END */

/*---------------------------------------------------------------------------*/

