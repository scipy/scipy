/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: itdr.h                                                            *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method ITDR                               *
 *         (Inverse Transformed Density Rejection)                           *
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
   =METHOD  ITDR   Inverse Transformed Density Rejection

   =UP  Methods_for_CONT

   =REQUIRED monotone PDF, dPDF, pole

   =OPTIONAL splitting point between pole and tail region, c-values

   =SPEED Set-up: moderate, Sampling: moderate

   =REINIT supported

   =REF  [HLDa07]

   =DESCRIPTION
      ITDR is an acceptance/rejection method that works for monotone
      densities. It is especially designed for PDFs with a single
      pole. It uses different hat functions for the pole region and
      for the tail region. For the tail region @emph{Transformed Density
      Rejection} with a single construction point is used. 
      For the pole region a variant called @emph{Inverse Transformed
      Density Rejection} is used. The optimal splitting point between
      the two regions and the respective maximum local concavity and
      inverse local concavity (@pxref{Glossary}) that guarantee valid
      hat functions for each regions are estimated.
      This splitting point is set to the intersection point of local
      concavity and inverse local concavity.
      However, it is assumed that both, the local concavity and the
      inverse local concavity do not have a local minimum in the
      interior of the domain (which is the case for all standard
      distributions with a single pole).
      In other cases (or when the built-in search routines do not
      compute non-optimal values) one can provide the splitting point,
      and the @i{c}-values.

   =HOWTOUSE
      Method ITDR requires a distribution object with given PDF
      and its derivative and the location of the pole (or mode).
      The PDF must be monotone and may contain a pole.
      It must be set via the unur_distr_cont_set_pdf() and 
      unur_distr_cont_set_dpdf() calls. The PDF should return
      UNUR_INFINITY for the pole. Alternatively, one can also
      set the logarithm of the PDF and its derivative via the 
      unur_distr_cont_set_logpdf() and unur_distr_cont_set_dlogpdf()
      calls. This is in especially useful since then the setup and
      search routines are numerically more stable. Moreover, for many
      distributions computing the logarithm of the PDF is less
      expensive then computing the PDF directly.

      The pole of the distribution is given by a
      unur_distr_cont_set_mode() call. Notice that distributions with 
      ``heavy'' poles may have numerical problems caused by the
      resultion of the floating point numbers used by computers.
      While the minimal distance between two different floating point 
      numbers is about @code{1.e-320} near @code{0.} it increases
      to @code{1.e-16} near @code{1.} Thus any random variate
      generator implemented on a digital computer in fact draws samples
      from a discrete distribution that approximates the desired
      continuous distribution. For distributions with ``heavy'' poles
      not at 0 this approximation may be too crude and thus every
      goodness-of-fit test will fail.
      Besides this theoretic problem that cannot be resolved we 
      have to take into consideration that round-off errors occur more 
      frequently when we have PDFs with poles far away from
      @code{0.} Method ITDR tries to handles this situation as good as
      possible by moving the pole into @code{0.} 
      Thus do not use a wrapper for your PDF that hides this shift
      since the information about the resolution of the floating point
      numbers near the pole gets lost.

      Method ITDR uses different hats for the pole region and for the
      tail region. The splitting point between these two regions, the
      optimal @i{c}-value and design points for constructing the hats
      using Transformed Density Rejection are computed automatically.
      (The results of these computations can be read using the
      respective calls unur_itdr_get_xi(), unur_itdr_get_cp(), and 
      unur_itdr_get_ct() for the intersection point between local
      concavity and inverse local concavity, the @i{c}-value for the
      pole and the tail region.)
      However, one can also analyze the local concavity and inverse
      local concavity set the corresponding values using 
      unur_itdr_set_xi(), unur_itdr_set_cp(), and 
      unur_itdr_set_ct() calls.
      Notice, that @i{c}-values greater than -1/2 can be set to
      @code{-0.5}. Although this results in smaller acceptance
      probabities sampling from the hat distribution is much faster
      than for other values of @i{c}. Depending on the expenses of
      evaluating the PDF the resulting algorithm is usually faster.

      It is possible to change the parameters and the domain of the chosen 
      distribution and run unur_reinit() to reinitialize the generator object.
      However, the values given by unur_itdr_set_xi(), unur_itdr_set_cp(), 
      or unur_itdr_set_ct() calls are then ignored when unur_reinit() is
      called.
   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_itdr_new( const UNUR_DISTR *distribution );
/* 
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_itdr_set_xi( UNUR_PAR *parameters, double xi );
/* 
   Sets points where local concavity and inverse local concavity
   are (almost) equal. It is used to estimate the respective c-values
   for pole region and hat regions and to determine the splitting point @i{bx} 
   between pole and tail region. 
   If no such point is provided it will be computed automatically.

   Default: not set.
*/

int unur_itdr_set_cp( UNUR_PAR *parameters, double cp );
/* 
   Sets parameter @var{cp} for transformation T for inverse 
   density in pole region.
   It must be at most 0 and greater than -1. 
   A value of @code{-0.5} is treated separately and usually results in
   faster marginal generation time (at the expense of smaller
   acceptance probabilities.
   If no @var{cp}-value is given it is estimated automatically.

   Default: not set.
*/

int unur_itdr_set_ct( UNUR_PAR *parameters, double ct );
/* 
   Sets parameter @var{ct} for transformation T for 
   density in tail region.
   It must be at most 0. For densities with unbounded domain
   it must be greater than -1. 
   A value of @code{-0.5} is treated separately and usually results in
   faster marginal generation time (at the expense of smaller
   acceptance probabilities.
   If no @var{ct}-value is given it is estimated automatically.

   Default: not set.
*/

double unur_itdr_get_xi( UNUR_GEN *generator );
/* */

double unur_itdr_get_cp( UNUR_GEN *generator );
/* */

double unur_itdr_get_ct( UNUR_GEN *generator );
/* 
   Get intersection point @var{xi}, and c-values @var{cp} and @var{ct},
   respectively. 
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

double unur_itdr_get_area( UNUR_GEN *generator );
/* 
   Get area below hat.
   (In case of an error @code{UNUR_INFINITY} is returned.)
*/

int unur_itdr_set_verify( UNUR_PAR *parameters, int verify );
/* 
   Turn verifying of algorithm while sampling on/off.

   If the condition @unurmath{PDF(x) \leq hat(x)} is
   violated for some @i{x} then @code{unur_errno} is set to
   @code{UNUR_ERR_GEN_CONDITION}. However, notice that this might
   happen due to round-off errors for a few values of
   @i{x} (less than 1%).

   Default is FALSE.
*/

int unur_itdr_chg_verify( UNUR_GEN *generator, int verify );
/* 
   Change the verifying of algorithm while sampling on/off.
*/

/* =END */
/*---------------------------------------------------------------------------*/
