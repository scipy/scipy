/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: arou.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method AROU                               *
 *         (Adaptive Ratio-Of-Uniforms)                                      *
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
   =METHOD  AROU   Automatic Ratio-Of-Uniforms method

   =UP  Methods_for_CONT

   =REQUIRED T-concave PDF, dPDF

   =OPTIONAL mode

   =SPEED Set-up: slow, Sampling: fast

   =REINIT not implemented

   =REF  [LJa00]

   =DESCRIPTION
      AROU is a variant of the ratio-of-uniforms method that uses the
      fact that the transformed region is convex for many distributions.
      It works for all T-concave distributions with T(x) = -1/sqrt(x).
      
      It is possible to use this method for correlation induction by
      setting an auxiliary uniform random number generator via the
      unur_set_urng_aux() call. (Notice that this must be done after a
      possible unur_set_urng() call.)
      When an auxiliary generator is used then the number of used
      uniform random numbers that is used up for one generated random
      variate is constant and equal to 1.
      
      There exists a test mode that verifies whether the conditions for
      the method are satisfied or not while sampling. It can be
      switched on by calling unur_arou_set_verify() and
      unur_arou_chg_verify(), respectively.
      Notice however that sampling is (much) slower then.
      
      For densities with modes not close to 0 it is suggested to set
      either the mode or the center of the distribution by the
      unur_distr_cont_set_mode() or unur_distr_cont_set_center() call.
      The latter is the approximate location of the mode or the mean
      of the distribution. This location provides some information
      about the main part of the PDF and is used to avoid numerical
      problems.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_arou_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_arou_set_usedars( UNUR_PAR *parameters, int usedars );
/* 
   If @var{usedars} is set to TRUE, ``derandomized adaptive rejection
   sampling'' (DARS) is used in setup.
   Segments where the area between hat and squeeze is too
   large compared to the average area between hat and squeeze
   over all intervals are split.
   This procedure is repeated until the ratio between area below squeeze
   and area below hat exceeds the bound given by 
   unur_arou_set_max_sqhratio() call or the maximum number of segments is 
   reached. Moreover, it also aborts when no more segments can be
   found for splitting.

   Segments are split such that the angle of the segments are halved
   (corresponds to arc-mean rule of method TDR (@pxref{TDR})).

   Default is TRUE.
*/

int unur_arou_set_darsfactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for ``derandomized adaptive rejection sampling''.
   This factor is used to determine the segments that are ``too
   large'', that is, all segments where the area between squeeze and
   hat is larger than @var{factor} times the average area over all
   intervals between squeeze and hat.
   Notice that all segments are split when @var{factor} is set to
   @code{0.}, and that there is no splitting at all when @var{factor}
   is set to @code{UNUR_INFINITY}.

   Default is @code{0.99}. There is no need to change this parameter.
*/

int unur_arou_set_max_sqhratio( UNUR_PAR *parameters, double max_ratio );
/* 
   Set upper bound for the
   ratio (area inside squeeze) / (area inside envelope).
   It must be a number between 0 and 1.
   When the ratio exceeds the given number no further construction
   points are inserted via adaptive rejection sampling.
   Use @code{0} if no construction points should be added after the
   setup.
   Use @code{1} if adding new construction points should not be
   stopped until the maximum number of construction points is reached.

   Default is @code{0.99}.
*/

double unur_arou_get_sqhratio( const UNUR_GEN *generator );
/* 
   Get the current ratio (area inside squeeze) / (area inside envelope)
   for the generator. 
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_arou_get_hatarea( const UNUR_GEN *generator );
/* 
   Get the area below the hat for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_arou_get_squeezearea( const UNUR_GEN *generator );
/* 
   Get the area below the squeeze for the generator.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/


int unur_arou_set_max_segments( UNUR_PAR *parameters, int max_segs );
/* 
   Set maximum number of segements.
   No construction points are added @emph{after} the setup when the
   number of segments succeeds @var{max_segs}.

   Default is @code{100}.
*/


int unur_arou_set_cpoints( UNUR_PAR *parameters, int n_stp, const double *stp );
/* 
   Set construction points for enveloping polygon.
   If @var{stp} is NULL, then a heuristical rule of thumb is used to
   get @var{n_stp} construction points. 
   This is the default behavior when this routine is not called.
   The (default) number of construction points is @code{30}, then.
*/


int unur_arou_set_usecenter( UNUR_PAR *parameters, int usecenter );
/* 
   Use the center as construction point.
   Default is TRUE.
*/


int unur_arou_set_guidefactor( UNUR_PAR *parameters, double factor );
/* 
   Set factor for relative size of the guide table for indexed search
   (see also method DGT @ref{DGT}). It must be greater than or equal
   to @code{0}. 
   When set to @code{0}, then sequential search is used.

   Default is @code{2}.
*/


int unur_arou_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_arou_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/


int unur_arou_set_pedantic( UNUR_PAR *parameters, int pedantic );
/* 
   Sometimes it might happen that unur_init() has been executed
   successfully. But when additional construction points are added by
   adaptive rejection sampling, the algorithm detects that the
   PDF is not T-concave. 

   With @var{pedantic} being TRUE, the
   sampling routine is then exchanged by a routine that simply returns
   @code{UNUR_INFINITY}. Otherwise the new point is not added to the
   list of construction points. At least the hat function remains
   T-concave.

   Setting @var{pedantic} to FALSE allows sampling from a
   distribution which is ``almost'' T-concave and small errors are
   tolerated. However it might happen that the hat function cannot be
   improved significantly. When the hat function that has been
   constructed by the unur_init() call is extremely large then it
   might happen that the generation times are extremely high
   (even hours are possible in extremely rare cases).

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/



