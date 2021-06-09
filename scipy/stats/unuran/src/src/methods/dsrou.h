/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: dsrou.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method DSROU                              *
 *         (Discrete, Simple universal generator, Ratio-Of-Uniforms method)  *
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
   =METHOD  DSROU   Discrete Simple Ratio-Of-Uniforms method

   =UP  Methods_for_DISCR

   =REQUIRED T-concave PMF, mode, sum over PMF

   =SPEED Set-up: fast, Sampling: slow

   =REINIT supported

   =REF  [LJa01] [HLD04: Sect.10.3.2, Alg.10.6]

   =DESCRIPTION
      DSROU is based on the ratio-of-uniforms method
      (@pxref{Ratio-of-Uniforms}) but uses universal 
      inequalities for constructing a (universal) bounding rectangle.
      It works for all @i{T}-concave distributions with 
      @unurmath{T(x) = -1/\sqrt{x}}.

      The method requires the PMF, the (exact) location of the mode
      and the sum over the given PDF. The rejection constant is 4 for
      all @i{T}-concave distributions. Optionally the CDF at the mode
      can be given to increase the performance of the algorithm. Then
      the rejection constant is reduced to 2.
      
   =HOWTOUSE
      The method works for @i{T}-concave discrete distributions with
      given PMF. The sum over of the PMF or an upper bound of this sum
      must be known. 

      Optionally the CDF at the mode can be given to increase the
      performance using unur_dsrou_set_cdfatmode().
      However, this @strong{must not} be called if the sum over the
      PMF is replaced by an upper bound.
      
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.

      If any of mode, CDF at mode, or the sum over the PMF has been
      changed, then unur_reinit() must be executed. 
      (Otherwise the generator produces garbage).

      There exists a test mode that verifies whether the conditions
      for the method are satisfied or not while sampling. It can be
      switched on or off by calling unur_dsrou_set_verify() and
      unur_dsrou_chg_verify(), respectively.
      Notice however that sampling is (a little bit) slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_dsrou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_dsrou_set_cdfatmode( UNUR_PAR *parameters, double Fmode );
/* 
   Set CDF at mode. 
   When set, the performance of the algorithm is increased by factor 2.
   However, when the parameters of the distribution are changed
   unur_dsrou_chg_cdfatmode() has to be used to update this value.
   Notice that the algorithm detects a mode at the left boundary of
   the domain automatically and it is not necessary to use this call
   for a monotonically decreasing PMF. 

   Default: not set.
*/

int unur_dsrou_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_dsrou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PMF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/*...........................................................................*/

int unur_dsrou_chg_cdfatmode( UNUR_GEN *generator, double Fmode );
/* 
   Change CDF at mode of distribution.
   unur_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/
