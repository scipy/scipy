/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: gibbs.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method GIBBS                              *
 *         (Markov Chain - GIBBS sampler)                                    *
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
   =METHOD  GIBBS   Markov Chain - GIBBS sampler

   =UP  MCMC_Methods_for_CVEC

   =REQUIRED T-concave logPDF, derivatives of logPDF

   =SPEED Set-up: fast, Sampling: moderate

   =REINIT not implemented

   =REF  [HLD04: Sect.14.1.2]

   =DESCRIPTION
      Method GIBBS implements a Gibbs sampler for a multivariate
      distribution with given joint density and its gradient.
      When running such a Markov chain all coordinates are updated
      cyclically using full conditional distributions. After each step
      the state of the chain is returned (i.e., a random point is
      returned whenever a single coordinate has been updated). 
      It is also possible to return only points after all coordinates
      have been updated by "thinning" the chain.
      Moreover, to reduce autocorrelation this thinning factor can be
      any integer. Notice, however, that the sampling time for a chain
      of given length is increased by the same factor, too.

      GIBBS also provides a variant of the Gibbs sampler where in
      each step a point from the full conditional distribution along
      some random direction is sampled. This direction is chosen
      uniformly from the sphere in each step.
      This method is also known as Hit-and-Run algorithm for
      non-uniform distributions.

      Our experiences shows that the original Gibbs sampler with
      sampling along coordinate axes is superior to random direction
      sampling as long as the correlations between the components of
      the random vector are not too high.

      For both variants transformed density rejection (see methods
      @pxref{TDR} and @pxref{ARS}) is used to
      sample from the full conditional distributions. In opposition to
      the univariate case, it is important that the factor @code{c} is
      as large as possible. I.e., for a log-concave density @code{c}
      must be set to @code{0.}, since otherwise numerical underflow
      might stop the algorithm.

      @emph{Important:} GIBBS does not generate independent random
      points. The starting point of the Gibbs chain must be in a
      "typical" region of the target distribution. If such a point is
      not known or would be too expensive, then the first part of the
      chain should be discarded (burn-in of the chain). 

   =HOWTOUSE
      For using the GIBBS method UNU.RAN needs the logarithm of the
      PDF of the multivariate joint distribution and its gradient or
      partial derivatives. 

      It provides two variants:
      @table @emph 
      @item coordinate direction sampling (Gibbs sampling) [default]
      The coordinates are updated cyclically.
      It requires the partial derivatives of the (logarithm of the)
      PDF of the target distribution, 
      see unur_distr_cvec_set_pdlogpdf(). 
      Otherwise, the gradient of the logPDF 
      (see unur_distr_cvec_set_dlogpdf()) 
      is used, which is more expensive. 

      This variant can be selected using
      unur_gibbs_set_variant_coordinate().

      @item random direction sampling (nonuniform Hit-and-Run algorithm)
      In each step is a direction is sampled uniformly from the sphere
      and the next point in the chain is sampled from the full
      conditional distribution along this direction.

      It requires the gradient of the logPDF and thus each step is
      more expensive than each step for coordinate direction sampling.

      This variant can be selected using
      unur_gibbs_set_variant_random_direction().
      @end table

      It is important that the @code{c} parameter for the TDR method
      is as large as possible. For logconcave distribution it must be
      set to @code{0}, since otherwise numerical underflow can cause
      the algorithm to stop.

      The starting point of the Gibbs chain must be "typical" for the
      target distribution. If such a point is not known or would be
      too expensive, then the first part of the chain should be
      discarded (burn-in of the chain). When using the 
      unur_gibbs_set_burnin() call this is done during the setup
      of the Gibbs sampler object.
      
      In case of a fatal error in the generator for conditional
      distributions the methods generates points that contain
      UNUR_INFINITY.

      @strong{Warning:} The algorithm requires that all full
      conditionals for the given distribution object are
      @i{T}-concave. However, this property is not checked. 
      If this property is not satisfied, then generation from the
      conditional distributions becomes (very) slow and might fail or
      (even worse) produces random vectors from an incorrect
      distribution.
      When using unur_gibbs_set_burnin() then the setup already
      might fail. Thus when in doubt whether GIBBS can be used for
      the targent distribution it is a good idea to use a burn-in for
      checking.  

      @emph{Remark:} It might happen (very rarely) that the chain
      becomes stuck due to numerical errors. (This is in particular the
      case when the given PDF does not fulfill the condition of this
      method.)
      When this happens during burn-in then the setup is aborted
      (i.e. it fails). Otherwise the chain restarts again from its
      starting point.

      @strong{Warning:} Be carefull with debugging flags. If it
      contains flag @code{0x01000000u} it produces a lot of output for
      each step in the algorithm.
      (This flag is switched of in the default debugging flags).

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_gibbs_new( const UNUR_DISTR *distribution );

/*...........................................................................*/

int unur_gibbs_set_variant_coordinate( UNUR_PAR *parameters );
/* 
   Coordinate Direction Sampling:
   Sampling along the coordinate directions (cyclic).

   This is the default.
*/

int unur_gibbs_set_variant_random_direction( UNUR_PAR *parameters );
/* 
   Random Direction Sampling:
   Sampling along the random directions.
*/

int unur_gibbs_set_c( UNUR_PAR *parameters, double c );
/* 
   Set parameter @var{c} for transformation @unurmath{T} of the
   transformed density rejection method.
   Currently only values between @code{0} and @code{-0.5} are
   allowed. If @code{c} is between @code{0} and @code{-0.5} it is set
   to @code{-0.5}. 

   For @var{c} @code{=0} (for logconcave densities) method ARS
   (@pxref{ARS}) is used which is very robust against badly
   normalized PDFs. For other values method TDR (@pxref{TDR}) is used.

   The value for @var{c} should be as large as possible to avoid 
   fatal numerical underflows. Thus for log-concave distributions
   @var{c} must be set to @code{0.}
 
   Default is @code{0}.
*/


int unur_gibbs_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
/* 
   Sets the starting point of the Gibbs sampler. @var{x0} must be 
   a "typical" point of the given distribution.
   If such a "typical" point is not known and a starting point is
   merely guessed, the first part of the Gibbs chain should be
   discarded (@emph{burn-in}), e.g.\ by mean of the
   unur_gibbs_set_burnin() call.

   Default is the result of unur_distr_cvec_get_center() for the 
   given distribution object.
*/

int unur_gibbs_set_thinning( UNUR_PAR *parameters, int thinning );
/*
   Sets the @var{thinning} parameter. When @var{thinning} is set to
   @i{k} then every @i{k}-th point from the iteration is returned by
   the sampling algorithm.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=1.

   Default: @code{1}.
*/

int unur_gibbs_set_burnin( UNUR_PAR *parameters, int burnin );
/*
   If a "typical" point for the target distribution is not known but
   merely guessed, the first part of the Gibbs chain should be
   discarded (@emph{burn-in}). This can be done during the
   initialization of the generator object. 
   The length of the burn-in can is then @var{burnin}.

   When method GIBBS is not applicable for the target distribution
   then the initialization already might fail during the burn-in.
   Thus this reduces the risk of running a generator that returns 
   UNUR_INFINITY cased by some fatal error during sampling.
   
   The thinning factor set by a unur_gibbs_set_thinning() call has
   no effect on the length of the burn-in, i.e., for the burn-in
   always a thinning factor @code{1} is used.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=0.

   Default: @code{0}.
*/

/*...........................................................................*/

const double *unur_gibbs_get_state( UNUR_GEN *generator );
/* */

int unur_gibbs_chg_state( UNUR_GEN *generator, const double *state );
/* 
   Get and change the current state of the Gibbs chain.
*/ 

int unur_gibbs_reset_state( UNUR_GEN *generator );
/* 
   Reset state of chain to starting point.

   @emph{Notice:} Currently this function does not reset
   the generators for conditional distributions. Thus it is not 
   possible to get the same Gibbs chain even when the underlying
   uniform random number generator is reset.
*/ 

/* =END */
/*---------------------------------------------------------------------------*/


