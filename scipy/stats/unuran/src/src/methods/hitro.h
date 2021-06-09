/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: hitro.h                                                           *
 *                                                                           *
 *   PURPOSE:                                                                *
 *         function prototypes for method HITRO                              *
 *         (Markov Chain - HIT-and-run sampler with Ratio-Of-uniforms)       *
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
   =METHOD  HITRO   Markov Chain - HIT-and-run sampler with Ratio-Of-uniforms

   =UP  MCMC_Methods_for_CVEC

   =REQUIRED PDF

   =OPTIONAL mode, center, bounding rectangle for acceptance region

   =REF [KLPa05]

   =SPEED Set-up: fast, Sampling: fast

   =REINIT not implemented

   =DESCRIPTION
      HITRO is an implementation of a hit-and-run sampler that runs on
      the acceptance region of the multivariate ratio-of-uniforms
      method, see @ref{Ratio-of-Uniforms}. 

      The Ratio-of-Uniforms transforms the region below the density
      into some region that we call "region of acceptance" in the
      following. The minimal bounding hyperrectangle of this region
      is given by 

      @unurmathdisplay{
      v^+   = \sup\limits_{x}               (f(x))^{1/r\,d+1}, \\
      u^-_i = \inf\limits_{x_i} (x_i-\mu_i) (f(x))^{r/r\,d+1}, \\
      u^+_i = \sup\limits_{x_i} (x_i-\mu_i) (f(x))^{r/r\,d+1}, }

      where @i{d} denotes the dimension of the distribution;
      @unurmath{x_i} is the @i{i}-th coordinate of point @i{x};
      @unurmath{\mu_i} is the @i{i}-th coordinate of the center
      @unurmath{\mu} of the distribution, i.e., a point in the 
      "main region" of the distribution.
      Using the center is important, since otherwise the acceptance
      region can become a very long and skinny ellipsoid along a
      diagonal of the (huge) bounding rectangle.

      For each step of the Hit-and-Run algorithm we have to choose
      some direction. This direction together with the current point
      of the chain determines a straight line. Then a point is sampled
      uniformly on intersection of this line and the region of
      acceptance. This is done by rejection from a uniform distribution
      on a line segment that covers it. 
      Depending of the chosen variant the endpoints of this covering
      line are computed either by means of a (not necessary minimal)
      bounding hyper-rectangle, or just the "covering plate" of the
      bounding hyper-rectangle.

      The required bounds of the hyper-rectable can be given directly
      by the user. Otherwise, these are computed automatically by
      means of a numerical routine by Hooke and Jeeves
      @unurbibref{HJa61} called direct search (see
      @file{src/utils/hooke.c} for further references and details). 
      However, this expensive computation can be avoided by determine
      these bounds "on the fly" by the following adaptive algorithm:
      Start with some (small) hyper-rectangle and enlarge it whenever
      the endpoints of the covering line segment are not contained in
      the acceptance region of the Ratio-of-Unfiorms method.
      This approach works reliable as long as the region of acceptance
      is convex. 

      The performance of the uniform sampling from the line segment is
      much improved if the covering line is adjusted (shortened)
      whenever a point is rejected (adaptive sampling). This technique
      reduces the expected number of iterations enormously.

      Method HITRO requires that the region of acceptance of the
      Ratio-of-Uniforms method is bounded. The shape of this region
      can be controlled by a parameter @i{r}. Higher values of @i{r}
      result in larger classes of distributions with bounded region
      of acceptance. (A distribution that has such a bounded region for
      some @i{r} also has a bounded region for every @i{r'} greater
      than @i{r}.) On the other hand the acceptance probability
      decreases with increasing @i{r}. Moreover, round-off errors are
      more likely and (for large values of @i{r}) might result in a
      chain with a stationary distribution different from the target
      distribution. 

      Method HITRO works optimal for distributions whose region of
      acceptance is convex. This is in particular the case for all
      log-concave distributions when we set @i{r} = @code{1}.
      For bounded but non-convex regions of acceptance convergence is
      yet not guarenteed by mathematical theory.


   =HOWTOUSE
      Method HITRO requires the PDF of the target distribution
      (derivatives are not necessary).

      The acceptance region of the Ratio-of-Uniforms transformation
      must be bounded. Its shape is controlled by parameter @i{r}. 
      By default this parameter is set to @code{1} as this guarentees
      a convex region of acceptance when the PDF of the given
      distribution is log-concave. It should only be set to a
      different (higher!) value using unur_vnrou_set_r() if otherwise
      @unurmath{x_i\,(f(x))^{r/r\,d+1}} were not 
      bounded for each coordinate.

      There are two variants of the HITRO sampler:
      @table @emph 
      @item coordinate direction sampling. [default]
      The coordinates are updated cyclically.
      This can be seen as a Gibbs sampler running on the acceptance
      region of the Ratio-of-Uniforms method.
      This variant can be selected using 
      unur_hitro_set_variant_coordinate().

      @item random direction sampling. 
      In each step is a direction is sampled uniformly from the
      sphere.

      This variant can be selected using
      unur_hitro_set_variant_random_direction().
      @end table

      Notice that each iteration of the coordinate direction sampler is
      cheaper than an iteration of the random direction sampler.

      Sampling uniformly from the line segment can be adjusted in
      several ways:

      @table @emph 
      @item Adaptive line sampling vs. simple rejection.
      When adaptive line sampling is switched on, the covering line is
      shortened whenever a point is rejected. However, when the region
      of acceptance is not convex the line segment from which we have
      to sample might not be connected. We found that the algorithm
      still works but at the time being there is no formal proof that
      the generated Markov chain has the required stationary
      distribution. 

      Adaptive line sampling can switch on/off by means of the
      unur_hitro_set_use_adaptiveline() call.

      @item Bounding hyper-rectangle vs. "covering plate".
      For computing the covering line we can use the bounding
      hyper-rectangle or just its upper bound.
      The latter saves computing time during the setup and 
      when computing the covering during at each iteration step
      at the expense of a longer covering line. When adaptive line
      sampling is used the total generation time for the entire chain
      is shorter when only the "covering plate" is used.

      @emph{Notice:} When coordinate sampling is used the entire
      bounding rectangle is used.

      Using the entire bounding hyper-rectangle can be switched on/off
      by means of the unur_hitro_set_use_boundingrectangle() call.

      @item Deterministic vs. adaptive bounding hyper-rectangle.
      A bounding rectangle can be given by the 
      unur_vnrou_set_u() and unur_vnrou_set_v() calls.
      Otherwise, the minimal bounding rectangle is computed
      automatically during the setup by means of a numerical
      algorithm. However, this is (very) slow especially in higher
      dimensions and it might happen that this algorithm (like
      any other numerical algorithm) does not return a correct result.

      Alternatively the bounding rectangle can be computed
      adaptively. In the latter case unur_vnrou_set_u() and
      unur_vnrou_set_v() can be used to provide a starting rectangle
      which must be sufficiently small. 
      Then both endpoints of the covering line segment are always
      check whether they are outside the acceptance region of the
      Ratio-of-Uniforms method. If they are not, then the line segment
      and the ("bounding") rectangle are enlarged using a factor that
      can be given using the unur_hitro_set_adaptive_multiplier() call.

      Notice, that running this method in the adaptive rectangle
      mode requires that the region of acceptance is convex when random
      directions are used, or the given PDF is unimodal when
      coordinate direction sampling is used.
      Moreover, it requires two additional calls to the PDF in each
      iteration step of the chain.

      Using addaptive bounding rectangles can be switched on/off
      by means of the unur_hitro_set_use_adaptiverectangle() call.

      @end table

      The algorithm takes of a bounded rectangular domain given by a
      unur_distr_cvec_set_domain_rect() call, i.e. the PDF is set to
      zero for every @i{x} outside the given domain.
      However, it is only the coordinate direction sampler where the
      boundary values are directly used to get the endpoins of the
      coverline line for the line sampling step.


      @emph{Important:} The bounding rectangle has to be
      provided for the function @unurmath{PDF(x-center)!}
      Notice that @code{center} is the center of the given
      distribution, see unur_distr_cvec_set_center().
      If in doubt or if this value is not optimal, it can be changed
      (overridden) by a unur_distr_cvec_set_center() call.

   =END
*/

/*---------------------------------------------------------------------------*/
/* Routines for user interface                                               */

/* =ROUTINES */

UNUR_PAR *unur_hitro_new( const UNUR_DISTR *distribution );
/*
   Get default parameters for generator.
*/

/*...........................................................................*/

int unur_hitro_set_variant_coordinate( UNUR_PAR *parameters );
/* 
   Coordinate Direction Sampling:
   Sampling along the coordinate directions (cyclic).

   @emph{Notice:} For this variant the entire bounding rectangle is
   always used independent of the
   unur_hitro_set_use_boundingrectangle() call.

   This is the default.
*/

int unur_hitro_set_variant_random_direction( UNUR_PAR *parameters );
/* 
   Random Direction Sampling:
   Sampling along the random directions.
*/

int unur_hitro_set_use_adaptiveline( UNUR_PAR *parameters, int adaptive );
/* 
   When @var{adaptive} is set to TRUE adaptive line sampling is
   applied, otherwise simple rejection is used.

   @emph{Notice:} When adaptive line sampling is switched off, 
   the entire bounding rectangle must be used since otherwise the
   sampling time can be arbitrarily slow.

   @emph{Warning:} When adaptive line sampling is switched off,
   sampling can be arbitrarily slow. In particular this happens
   when random direction sampling is used for distributions with
   rectangular domains. Then the algorithm can be trapped into 
   a vertex (or even edge).

   Default is TRUE.
*/

int unur_hitro_set_use_boundingrectangle( UNUR_PAR *parameters, int rectangle );
/* 
   When @var{rectangle} is set to TRUE the entire bounding rectangle is used
   for computing the covering line. Otherwise, only an upper bound for the 
   acceptance region is used.

   @emph{Notice:} When coordinate sampling is used the entire
   bounding rectangle has is always used and this call has no effect.

   Default: FALSE for random direction samplig, TRUE for coordinate
   direction sampling.
*/

int unur_hitro_set_use_adaptiverectangle( UNUR_PAR *parameters, int adaptive );
/* 
   When @var{adaptive} is set to FALSE the bounding rectangle is
   determined during the setup. Either, it is computed automatically by
   a (slow) numerical method, or it must be provided by 
   unur_vnrou_set_u() and unur_vnrou_set_v() calls.

   If @var{adaptive} is set to TRUE the bounding rectangle is computed
   adaptively. In this case the unur_vnrou_set_u() and
   unur_vnrou_set_v() calls can be used to provide a starting
   rectangle. This should be sufficiently small.
   If not given then we assume @unurmath{v_{max} = 1,}
   @unurmath{u_{min}=(-0.001,-0.001,\ldots,-0.001),} and
   @unurmath{u_{max}=(0.001,0.001,\ldots,0.001).}
   Adaptive enlargements of the bounding hyperrectangle can be
   controlled set setting an enlargement factor given 
   by a unur_hitro_set_adaptive_multiplier() call.
   
   Using adaptive computation of the bounding rectangle reduces the
   setup time significantly (when it is not given by the user) at the
   expense of two additional PDF evaluations during each iteration
   step.

   @emph{Important:} Using adaptive bounding rectangles requires that
   the region of acceptance is convex when random directions are used,
   or a unimodal PDF when coordinate direction sampling is used.

   Default: FALSE for random direction samplig, TRUE for coordinate
   direction sampling.
*/

/*...........................................................................*/

int unur_hitro_set_r( UNUR_PAR *parameters, double r );
/*
   Sets the parameter @var{r} of the generalized multivariate
   ratio-of-uniforms method.

   @emph{Notice}: This parameter must satisfy @var{r}>0.

   Default: @code{1}.
*/

int unur_hitro_set_v( UNUR_PAR *parameters, double vmax );
/*
   Set upper boundary for bounding hyper-rectangle.
   If not set not set the mode of the distribution is used.

   If adaptive bounding rectangles the value is used for the
   starting rectangle. If not given (and the mode of the distribution
   is not known) then @var{vmax}=@code{1e-3} is used.

   If deterministic bounding rectangles these values are the given
   values are used for the rectangle. If no value is given
   (and the mode of the distribution is not known), the upper
   bound of the minimal bounding hyper-rectangle is computed
   numerically (slow).

   Default: not set.
*/

int unur_hitro_set_u( UNUR_PAR *parameters, const double *umin, const double *umax );
/* 
   Sets left and right boundaries of bounding hyper-rectangle.

   If adaptive bounding rectangles these values are used for the
   starting rectangle. If not given then 
   @var{umin}=@code{@{-b,-b,@dots{},-b@}} and
   @var{umax}=@code{@{b,b,@dots{},b@}} with @code{b=1.e-3} is used.

   If deterministic bounding rectangles these values are the given
   values are used for the rectangle. If no values are given, the
   boundary of the minimal bounding hyper-rectangle is computed
   numerically (slow).

   @strong{Important}: The boundaries are those of the density shifted
   by the center of the distribution, i.e., for the
   function @unurmath{PDF(x-center)!}

   @emph{Notice}: Computing the minimal bounding rectangle may fail
   under some circumstances. Moreover, for multimodal distributions
   the bounds might be too small as only local extrema are computed.
   Nevertheless, for log-concave distributions it should work.

   Default: not set.
*/

int unur_hitro_set_adaptive_multiplier( UNUR_PAR *parameters, double factor );
/* 
   Adaptive enlargements of the bounding hyperrectangle can be
   controlled set setting the enlargement @var{factor}.
   This must be greater than 1. Values close to 1 result in small
   adaptive steps and thus reduce the risk of too large bounding
   rectangles. On the other hand many adaptive steps might be
   necessary. 

   @emph{Notice:} For practical reasons this call does not accept
   values for @var{factor} less than @code{1.0001}. If this value is
   UNUR_INFINITY this results in infinite loops.

   Default: @code{1.1}
*/ 

/*...........................................................................*/

int unur_hitro_set_startingpoint( UNUR_PAR *parameters, const double *x0 );
/* 
   Sets the starting point of the HITRO sampler in the original
   scale. @var{x0} must be a "typical" point of the given distribution.
   If such a "typical" point is not known and a starting point is
   merely guessed, the first part of the HITRO chain should be
   discarded (@emph{burn-in}), e.g.\ by mean of the
   unur_hitro_set_burnin() call.

   @emph{Important:} The PDF of the distribution must not vanish at
   the given point @var{x0}.

   Default is the result of unur_distr_cvec_get_center() for the 
   given distribution object.
*/

int unur_hitro_set_thinning( UNUR_PAR *parameters, int thinning );
/*
   Sets the @var{thinning} parameter. When @var{thinning} is set to
   @i{k} then every @i{k}-th point from the iteration is returned by
   the sampling algorithm.
   If thinning has to be set such that each coordinate is updated
   when using coordinate direction sampling, then @var{thinning}
   should be @code{dim+1} (or any multiple of it) where 
   @code{dim} is the dimension of the distribution object.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=1.

   Default: @code{1}.
*/

int unur_hitro_set_burnin( UNUR_PAR *parameters, int burnin );
/*
   If a "typical" point for the target distribution is not known but
   merely guessed, the first part of the HITRO chain should be
   discarded (@emph{burn-in}). This can be done during the
   initialization of the generator object. 
   The length of the burn-in can is then @var{burnin}.

   The thinning factor set by a unur_hitro_set_thinning() call has
   no effect on the length of the burn-in, i.e., for the burn-in
   always a thinning factor @code{1} is used.

   @emph{Notice}: This parameter must satisfy @var{thinning}>=0.

   Default: @code{0}.
*/

/*...........................................................................*/

const double *unur_hitro_get_state( UNUR_GEN *generator );
/* */

int unur_hitro_chg_state( UNUR_GEN *generator, const double *state );
/* 
   Get and change the current state of the HITRO chain.

   @emph{Notice:} The state variable contains the point in the
   @code{dim+1} dimensional point in the (tansformed) region of
   acceptance of the Ratio-of-Uniforms method. Its coordinate 
   are stored in the following order:
   @code{state[] = @{v, u1, u2, @dots{}, udim@}}.

   If the state can only be changed if the given @var{state} is inside
   this region.
*/ 

int unur_hitro_reset_state( UNUR_GEN *generator );
/* 
   Reset state of chain to starting point.

   @emph{Notice:} Currently this function does not reset
   the generators for conditional distributions. Thus it is not 
   possible to get the same HITRO chain even when the underlying
   uniform random number generator is reset.
*/ 

/* =END */

/*---------------------------------------------------------------------------*/

