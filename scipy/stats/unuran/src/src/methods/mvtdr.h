/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: mvtdr.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method MVTDR                              *
 *         (Multi-Variate Transformed Density Rejection)                     *
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
   =METHOD  MVTDR  Multi-Variate Transformed Density Rejection

   =UP  Methods_for_CVEC

   =REQUIRED log-concave (log)PDF, gradient of (log)PDF

   =OPTIONAL mode

   =SPEED Set-up: slow,
          Sampling: depends on dimension

   =REINIT not implemented

   =REF  [HLD04: Sect.11.3.4, Alg.11.15.] [LJa98]

   =DESCRIPTION
      MVTDR a multivariate version of the Transformed Density Rection
      (@pxref{TDR}) that works for log-concave densities.
      For this method the domain of the distribution is partitioned
      into cones with the mode (or the center) of the distribution as
      their (common) vertex. The hat function is then constructed as
      tangent planes of the transformed density in each of these
      cones. The respective construction points lie on the central
      lines in the cones through the vertex. The point is chosen such
      that the hat is minimal among all such points (see the given
      references for more details).

      The cones are created by starting with the orthants of the reals
      space. These are then iteratively split when the volume below
      the hat in such cones is too large. Thus an increasing number of
      cones results in a better fitting hat function.
      Notice however, that the required number of cones increases
      exponentially with the number of dimension.
      Moreover, due to the construction the rejection does not
      converge to 1 and remains strictly larger than 1.

      For distributions with bounded domains the cones are cut to
      pyramids that cover the domain.

   =HOWTOUSE
      Create a multivariate generator object that contains the PDF and
      its gradient. This object also should contain the mode of the
      distribution (or a point nearby should be provided as center of
      the distribution).  

      The method has three parameter to adjust the method for the given
      distribution:

      @table @code
      @item stepsmin
      Minimal number of iterations for splitting cones.
      Notice that we start with 2^dim initial cones and that we arrive
      at 2^(dim+stepsmin) cones after these splits. So this number
      must be set with care. It can be set by a 
      unur_mvtdr_set_stepsmin() call.

      @item boundsplitting
      Cones where the volume below the hat is relatively large
      (i.e. larger than the average volume over all cones times
      @code{boundsplitting} are further split.
      This parameter can set via a unur_mvtdr_set_boundsplitting() call.

      @item maxcones
      The maximum number of generated cones. When this number is
      reached, the initialization routine is stopped. Notice that the
      rejection constant can be still prohibitive large.
      This parameter can set via a unur_mvtdr_set_maxcones() call.

      @end table

      Setting of these parameter can be quite tricky. The default
      settings lead to hat functions where the volume below the hat is
      similar in each cone. However, there might be some problems with
      distributions with higher correlations, since then too few cones
      are created. Then it might be necessary to increase the values
      for @code{stepsmin} and @code{maxcones} and to set
      @code{boundsplitting} to @code{0}.

      The number of cones and the total volume below the hat can be
      controlled using the respective calls unur_mvtdr_get_ncones() and
      unur_mvtdr_get_hatvol(). Notice, that the rejection constant is
      bounded from below by some figure (larger than 1) that depends
      on the dimension.

      Unfortunately, the algorithm cannot detect the quality of the
      constructed hat. 

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_mvtdr_new( const UNUR_DISTR *distribution );
/* 
   Get parameters for generator.
*/

/*...........................................................................*/

int unur_mvtdr_set_stepsmin( UNUR_PAR *parameters, int stepsmin );
/* 
   Set minimum number of triangulation step for each starting cone.
   @var{stepsmin} must be nonnegative.

   Default: @code{5}.
*/

int unur_mvtdr_set_boundsplitting( UNUR_PAR *parameters, double boundsplitting );
/* 
   Set bound for splitting cones. All cones are split which have a
   volume below the hat that is greater than @var{bound_splitting} times
   the average over all volumes. However, the number given by the 
   unur_mvtdr_set_maxcones() is not exceeded.
   Notice that the later number is always reached 
   if @var{bound_splitting} is less than 1.

   Default: @code{1.5}
*/

int unur_mvtdr_set_maxcones( UNUR_PAR *parameters, int maxcones );
/* 
   Set maximum number of cones. 

   Notice that this number is always increased to 
   @unurmath{2^{dim+stepsmin}} where @i{dim} is the dimension of the
   distribution object and @i{stepsmin} the given mimimum number of
   triangulation steps.

   Notice: For higher dimensions and/or higher correlations between the
   coordinates of the random vector the required number of cones can
   be very high. A too small maximum number of cones can lead to 
   a very high rejection constant.

   Default: @code{10000}.
*/

int unur_mvtdr_get_ncones( const UNUR_GEN *generator );
/* 
   Get the number of cones used for the hat function of the 
   @var{generator}.
   (In case of an error @code{0} is returned.)
*/

double unur_mvtdr_get_hatvol( const UNUR_GEN *generator );
/* 
   Get the volume below the hat for the @var{generator}.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

int unur_mvtdr_set_verify( UNUR_PAR *parameters, int verify );
/* */

int unur_mvtdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.
   If the condition squeeze(@i{x}) <= PDF(@i{x}) <= hat(@i{x}) is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

/* =END */
/*---------------------------------------------------------------------------*/


