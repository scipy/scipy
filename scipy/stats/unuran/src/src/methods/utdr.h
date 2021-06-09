/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: utdr.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method UTDR                               *
 *         (Universal Transformed Density Rejection; 3-point method)         *
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
   =METHOD  UTDR   Universal Transformed Density Rejection

   =UP  Methods_for_CONT

   =REQUIRED T-concave PDF, mode, approximate area

   =SPEED Set-up: moderate, Sampling: Moderate

   =REINIT supported

   =REF  [HWa95] [HLD04: Sect.4.5.4, Alg.4.4]

   =DESCRIPTION
      UTDR is based on the transformed density rejection and uses
      three almost optimal points for constructing hat and squeezes.
      It works for all @i{T}-concave distributions with 
      @unurmath{T(x) = -1/\sqrt(x).}
      
      It requires the PDF and the (exact) location of the mode.
      Notice that if no mode is given at all, a (slow) numerical mode
      finder will be used. 
      Moreover the approximate area below the given PDF is used.
      (If no area is given for the distribution the algorithm assumes that it
      is approximately 1.)
      The rejection constant is bounded from above by 4
      for all @i{T}-concave distributions.
      
   =HOWTOUSE
      UTDR works for any continuous univariate distribution object with
      given @i{T}-concave PDF (with @unurmath{T(x) = -1/\sqrt{x},)}
      mode and approximate area below PDF. 

      When the PDF does not change at the mode for varying parameters, then
      this value can be set with unur_utdr_set_pdfatmode() to avoid some 
      computations. Since this value will not be updated any more when the 
      parameters of the distribution are changed,
      the unur_utdr_chg_pdfatmode() call is necessary to do this manually.
      
      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.
      Notice, that derived parameters like the mode must also be (re-) set
      if the parameters or the domain has be changed.
      Moreover, if the PDF at the mode has been provided by a 
      unur_utdr_set_pdfatmode() call, additionally
      unur_utdr_chg_pdfatmode() must be used (otherwise this call is
      not necessary since then this figure is computed directly from
      the PDF). 

      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not. It can be switched on by calling 
      unur_utdr_set_verify() and unur_utdr_chg_verify(), respectively.
      Notice however that sampling is slower then.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_utdr_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_utdr_set_pdfatmode( UNUR_PAR *parameters, double fmode );
/* 
   Set pdf at mode. 
   When set, the PDF at the mode is never changed.          
   This is to avoid additional computations, when the PDF does not
   change when parameters of the distributions vary. 
   It is only useful when the PDF at the mode does not change with
   changing parameters for the distribution.

   Default: not set.
*/

int unur_utdr_set_cpfactor( UNUR_PAR *parameters, double cp_factor );
/* 
   Set factor for position of left and right construction point.
   The @var{cp_factor} is used to find almost optimal construction
   points for the hat function.
   There is no need to change this factor in almost all situations.

   Default is @code{0.664}.
*/

int unur_utdr_set_deltafactor( UNUR_PAR *parameters, double delta );
/* 
   Set factor for replacing tangents by secants.
   higher factors increase the rejection constant but reduces the risk of
   serious round-off errors.
   There is no need to change this factor it almost all situations.

   Default is @code{1.e-5}.
*/

int unur_utdr_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_utdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/*...........................................................................*/

int unur_utdr_chg_pdfatmode( UNUR_GEN *generator, double fmode );
/* 
   Change PDF at mode of distribution.
   unur_reinit() must be executed before sampling from the 
   generator again.
*/

/* =END */
/*---------------------------------------------------------------------------*/
