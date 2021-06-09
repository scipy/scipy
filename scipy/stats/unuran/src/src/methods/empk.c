/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      empk.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generate from kernel estimation                              *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given observed sample.                                               *
 *      Produce a value x consistent with this sample                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to sample                                                    *
 *      alpha factor (is required if a kernel is given)                      *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      kernel                                                               *
 *      smoothing factor                                                     *
 *      beta factor                                                          *
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
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *                                                                           *
 *   [1] Devroye, L. and Gyorfi, L. (1985): Nonparametric Density            *
 *       Estimation: The $L_1$ View. New-York: John Wiley.                   *
 *                                                                           *
 *   [2] Devroye, L. (1986): Non-Uniform Random Variate Generation. (p. 765) *
 *       New York: Springer-Verlag.                                          *
 *                                                                           *
 *   [3] Hoermann, W. and Leydold, J. (2000): Automatic random variate       *
 *       generation for simulation input. Proceedings of the 2000 Winter     *
 *       Simulation Conference. (??? eds.)                                   *
 *                                                                           *
 *   [4] Silverman, B. (1986): Density Estimation for Statistics and         *
 *       Data Analysis. London: Chapman and Hall.                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   The method here is based on the fact that it is very easy to generate   *
 *   random variates that have as density the kernel density estimate of a   *
 *   given sample: It is enough to randomly select a point from the given    *
 *   sample and return it together with some additive noise (see [1], [2]    *
 *   [3] and [4]). Note that it is not necessary to compute the density      * 
 *   estimate itself! In fact it is much easier to generate variates from    *
 *   the density estimate than to compute the estimate itself!               *
 *                                                                           *
 *   For kernel density estimation we have to choose the kernel function     *
 *   (this is for us the density of the noise) and the smoothing parameter   *
 *   (i.e. the scale parameter of the kernel variate.) There is a lot        *
 *   of theory discussing the optimal smoothing parameter for density        *
 *   estimation. As we are mainly interested in a simple formula we          *
 *   use the formula explained in [4] p.46-47 that is minimizing the         *
 *   asymptotically MISE (mean integrated squared error) for a given         *
 *   kernel. We have the danger of oversmoothing the data if they are        *
 *   non-normal. Especially for multimodal distributions we get a better     *
 *   bandwidth if we replace the sample's stdev s in the formula by the      *
 *   minimum of s and the interquartilsrange R divided through 1.34.         *
 *                                                                           *
 *      bwidth_opt = alpha * beta *  min(s,R/1.34) * n^(-1/5)                *
 *                                                                           *
 *   If we assume normal data we have to take beta = 1.364, which            *
 *   is the default value of beta for the library.                           *
 *   To make it easy for the user to choose the bwidth as a fraction         *
 *   of the bwidth_opt explained above we introduced the parameter           *
 *   smoothing_factor and compute:                                           *
 *                                                                           *
 *      bwidth= bwidth_opt * smoothing_factor                                *
 *                                                                           *
 *   Alpha is a constant only depending on the kernel K(t) and can be        *
 *   computed as                                                             *
 *                                                                           *
 *      alpha(K) = Var(K)^(-2/5){ \int K(t)^2 dt}^(1/5)                      *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *   Algorithm EMPK:                                                         *
 *                                                                           *
 *   [Input]                                                                 *
 *   Sample X(i) for i=1,2,..,n                                              *
 *                                                                           *
 *   [Generate]                                                              *
 *   1: Generate a random variate I uniformly distributed on {1, 2,...,n}    *
 *   2: Generate a random variate W from the distribution with density K(t)  *
 *   3: Return X(I) + bwidth * W                                             *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *   As kernel we can take any unimodal density symmetric around 0 that has  *
 *   finite variance. Regarding the asymptotic MISE the so called            *
 *   Epanechnikov (or Bartlett) kernel (K(x) = 3/4(1-x^2) |x|<=1)            *
 *   is optimal, but several others have almost the same efficiency.         *
 *   From a statistical point of view it may be attractive to add            *
 *   Gaussian noise. If we are interested in the speed of the sampling       *
 *   procedure it is of course fastest to use the rectangular (or Boxcar)    *
 *   kernel as this means uniform noise.                                     *
 *                                                                           *
 *   The mean of the empirical distribution is equal to the sample mean.     *
 *   One disadvantage of this method lies in the fact that the variance of   *
 *   the empirical distribution generated by the kernel method is always     *
 *   higher than the sample standard-deviation s. As pointed out in [4]      *
 *   it is not difficult to correct this error. For the variance             *
 *   corrected version we have to compute the mean mux and the               *
 *   standardeviation s of the sample and we have to know the Variance of    *
 *   the kernel V(K). Then we can replace step (3) of the algorithm above    *
 *   by:                                                                     *
 *                                                                           *
 *   3': Return mux +(X(I) - mux + bwidth * W)/(1+bwidth^2 V(K)/s^2)^(1/2)   *
 *                                                                           *
 *   We can turn on or off variance correction with the function             *
 *                                                                           *
 *      unur_empk_set_varcor();                                              *
 *                                                                           *
 *   We also have problems with this method if we want to generate only      *
 *   positive random variates, as due to the added noise it can happen       *
 *   that the generated variate is negative, even if all values of the       *
 *   observed sample are positive. One simple possibility to get around      *
 *   this problem is the so called mirroring principle. If the generated     *
 *   variate X is negative we simply return -X.                              *
 *                                                                           *
 *   We can turn on or off the mirroring with the function                   *
 *                                                                           *
 *      unur_empk_set_positive();                                            *
 *                                                                           *
 *   It is not difficult to see that the mirroring principle is              *
 *   disturbing the variance correction, which can no longer guarantee       *
 *   that the variance of the generated distribution is equal to the         *
 *   sample variance.                                                        *
 *   The mirroring principle is also spoiling the nice property that the     *
 *   mean of the generated distribution is equal to the sample mean. Both    *
 *   problems are only of practical relevance if a high proportion of the    *
 *   sample is close to 0. There are no simple methods to get around these   *
 *   problems. It is a well-known fact that kernel density estimation        *
 *   performs poor around a point of discontinuity of the density! There     *
 *   are special kernels suggested in the literature to cope with this       *
 *   problem but they are not easily applicable to our algorithm.            *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <distributions/unur_stddistr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "cstd.h"
#include "empk.h"
#include "empk_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define EMPK_VARFLAG_VARCOR     0x001u   /* use variance correction          */
#define EMPK_VARFLAG_POSITIVE   0x002u   /* only positive values             */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define EMPK_DEBUG_PRINTDATA   0x00000100u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define EMPK_SET_KERNEL         0x010u    /* kernel density                  */
#define EMPK_SET_KERNGEN        0x020u    /* generator for kernel            */
#define EMPK_SET_KERNELVAR      0x001u    /* variance of kernel              */
#define EMPK_SET_ALPHA          0x002u    /* alpha factor                    */
#define EMPK_SET_BETA           0x004u    /* beta factor                     */
#define EMPK_SET_SMOOTHING      0x008u    /* smoothing factor                */

/*---------------------------------------------------------------------------*/

#define GENTYPE "EMPK"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empk_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empk_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empk_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_empk_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_empk_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

inline static int _unur_empk_comp_stddev( double *data, int n_data,
					  double *mean, double *stddev);
/*---------------------------------------------------------------------------*/
/* compute mean and standard deviation of data.                              */
/*---------------------------------------------------------------------------*/

inline static double _unur_empk_comp_iqrtrange( double *data, int n_data );
/*---------------------------------------------------------------------------*/
/* compute interquartile range.                                              */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_empk_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_empk_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cemp      /* data for distribution object      */

#define PAR       ((struct unur_empk_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_empk_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cemp /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

/*---------------------------------------------------------------------------*/

#define SQU(a) ((a)*(a))

/*---------------------------------------------------------------------------*/

/* compare two doubles (needed for sorting) */
inline static int 
compare_doubles (const void *a, const void *b)
{ 
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}

/*---------------------------------------------------------------------------*/

#define _unur_empk_getSAMPLE(gen)   (_unur_empk_sample)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_empk_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_par *par;

  /* check arguments */
  _unur_check_NULL( GENTYPE,distr,NULL );

  /* check distribution */
  if (distr->type != UNUR_DISTR_CEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CEMP,NULL);

  if (DISTR_IN.sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"observed sample"); return NULL; }
  if (DISTR_IN.n_sample < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"number of observed sample"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_empk_par) );
  COOKIE_SET(par,CK_EMPK_PAR);

  /* copy input */
  par->distr    = distr;          /* pointer to distribution object          */

  /* set default values */
  PAR->kernvar   = 1.;            /* variance of used kernel, only used if varcor == 1 */
  PAR->alpha     = 0.7763884;     /* alpha for Gaussian kernel, efficiency = 0.951 */
  PAR->beta      = 1.3637439;     /* optimal BETA if data follow normal distribution */
  PAR->smoothing = 1.;            /* determines how "smooth" the estimated density will be */

  /* kernel */
  PAR->kerngen   = NULL;          /* random variate generator for kernel (by user)  */
  PAR->kernel    = NULL;          /* random variate generator for kernel (internal) */

  par->method   = UNUR_METH_EMPK; /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init     = _unur_empk_init;

  return par;

} /* end of unur_empk_new() */

/*****************************************************************************/

int 
unur_empk_set_kernel( struct unur_par *par, unsigned kernel)
     /*----------------------------------------------------------------------*/
     /* set standard kernel and start kernel generator                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   kernel ... identifier for standard kernel                          */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  UNUR_DISTR *kerndist;
  double fpar[4];

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );

  /* It is not possible to call unur_empk_set_kernel() twice. */
  if (par->set & EMPK_SET_KERNEL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"Cannot overwrite kernel");
    return UNUR_ERR_PAR_SET;
  }

  /* check for standard distribution and set data */
  switch (kernel) {
  case UNUR_DISTR_EPANECHNIKOV:
    /* Epanechnikov (or Bartlett) kernel. efficiency 1.000
       minimizes asymptotic MISE.
       PDF: f(x)=3/4 (1-x^2) for -1<x<1
       (Beta distribution with p=2, q=2, a=-1, b=1) */
    fpar[0]     = 2.;
    fpar[1]     = 2.;
    fpar[2]     = -1.;
    fpar[3]     = 1.;
    kerndist    = unur_distr_beta( fpar, 4 );
    PAR->kernel  = unur_init( unur_arou_new ( kerndist ) );
    PAR->alpha   = 1.718771928;
    PAR->kernvar = 0.2;
    unur_distr_free( kerndist );
    break;

  case UNUR_DISTR_GAUSSIAN:
    /* Gaussian (normal) kernel. efficiency = 0.951 */
    kerndist    = unur_distr_normal( NULL, 0 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 0.7763884;
    PAR->kernvar = 1.;
    unur_distr_free( kerndist );
    break;

  case UNUR_DISTR_BOXCAR:
    /* Boxcar (uniform, rectangular) kernel. efficiency = 0.930 */
    fpar[0]     = -1.;
    fpar[1]     = 1.;
    kerndist    = unur_distr_uniform( fpar, 2 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 1.351;
    PAR->kernvar = 1./3.;
    unur_distr_free( kerndist );
    break;

  case UNUR_DISTR_STUDENT:
    /* t3 kernel (Student's distribution with 3 degrees of freedom).
       efficiency = 0.679 */
    fpar[0] = 3.;
    kerndist    = unur_distr_student( fpar, 1 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 0.48263;
    PAR->kernvar = 3.;
    unur_distr_free( kerndist );
    break;

  case UNUR_DISTR_LOGISTIC:
    /* logistic kernel. efficiency = efficiency = 0.887 */
    kerndist    = unur_distr_logistic( NULL, 0 );
    PAR->kernel  = unur_init( unur_cstd_new ( kerndist ) );
    PAR->alpha   = 0.434;
    PAR->kernvar = 3.289868133696; /* Pi^2/3 */
    unur_distr_free( kerndist );
    break;

  default:
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"Unknown kernel. make it manually");
    return UNUR_ERR_PAR_SET;
  }

  /* generation of kernel successful ? */
  if (PAR->kernel == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"Could not initialize kernel generator");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  /* changelog */
  par->set &= ~EMPK_SET_KERNGEN;   /* just for the case that ..._set_kerngen() has been called */
  par->set |= EMPK_SET_KERNEL | EMPK_SET_ALPHA | EMPK_SET_KERNELVAR;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_empk_set_kernel() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_kernelgen( struct unur_par *par, const struct unur_gen *kernelgen,
			 double alpha, double kernelvar )
     /*----------------------------------------------------------------------*/
     /* set generator for kernel distribution                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   kernelgen ... generator object that holds the kernel               */
     /*   alpha     ... parameter depending on kernel                        */
     /*   kernelvar ... variance of kernel                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, kernelgen, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );

  /* It is not possible to run this function after a 
     unur_empk_set_kernel() call. */
  if (par->set & EMPK_SET_KERNEL) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"Cannot overwrite kernel");
    return UNUR_ERR_PAR_SET;
  }
  /* check kernel distribution */
  if ( (kernelgen->method & UNUR_MASK_TYPE) != UNUR_METH_CONT ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return UNUR_ERR_DISTR_INVALID;
  }

  /* check new parameter for generator */
  if (alpha <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"alpha <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* set kernel distribution */
  PAR->kerngen = kernelgen;

  /* store alpha factor */
  PAR->alpha = alpha;

  /* changelog */
  par->set |= EMPK_SET_KERNGEN | EMPK_SET_ALPHA;

  /* set kernel variance */
  PAR->kernvar = kernelvar;
  
  if (kernelvar > 0.)
    par->set |= EMPK_SET_KERNELVAR;
  else
    /* variance correction disabled */
    par->set &= ~EMPK_SET_KERNELVAR;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_empk_set_kernelgen() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_beta( struct unur_par *par, double beta )
     /*----------------------------------------------------------------------*/
     /* set beta factor                                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   beta     ... parameter depending on sample                         */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );

  /* check new parameter for generator */
  if (beta <= 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"beta <= 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->beta = beta;

  /* changelog */
  par->set |= EMPK_SET_BETA;

  return UNUR_SUCCESS;

} /* end of unur_empk_set_beta() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_smoothing( struct unur_par *par, double smoothing )
     /*----------------------------------------------------------------------*/
     /* set smoothing factor                                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par       ... pointer to parameter for building generator object   */
     /*   smoothing ... smoothing factor                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );

  /* check new parameter for generator */
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->smoothing = smoothing;

  /* changelog */
  par->set |= EMPK_SET_SMOOTHING;

  return UNUR_SUCCESS;

} /* end of unur_empk_set_smoothing() */

/*---------------------------------------------------------------------------*/

int
unur_empk_chg_smoothing( struct unur_gen *gen, double smoothing )
     /*----------------------------------------------------------------------*/
     /* change smoothing factor                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*   smoothing ... smoothing factor                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, EMPK, UNUR_ERR_GEN_INVALID );
  
  /* no changelog required */

  /* check new parameter for generator */
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* recompute band width */
  GEN->bwidth = smoothing * GEN->bwidth_opt;

  /* recompute constant for variance corrected version */
  GEN->sconst = 1./sqrt(1. + GEN->kernvar * SQU( GEN->bwidth/GEN->stddev_observ ) );

  /* store smoothing factor */
  GEN->smoothing = smoothing;

  /* changelog */
  gen->set |= EMPK_SET_SMOOTHING;

  return UNUR_SUCCESS;

} /* end of unur_empk_chg_smoothing() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_varcor( struct unur_par *par, int varcor )
     /*----------------------------------------------------------------------*/
     /* turn variance correction on/off                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   varcor   ... 0 = no variance correction,  !0 = variance correction */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   variance correction is the default                                 */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );

  /* we do not know yet if the kernel variance known */

  /* we use a bit in variant */
  par->variant = (varcor) 
    ? (par->variant | EMPK_VARFLAG_VARCOR) 
    : (par->variant & (~EMPK_VARFLAG_VARCOR));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_empk_set_varcor() */

/*---------------------------------------------------------------------------*/

int
unur_empk_chg_varcor( struct unur_gen *gen, int varcor )
     /*----------------------------------------------------------------------*/
     /* turn variance correction on/off                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to generator object                           */
     /*   varcor   ... 0 = no variance correction,  !0 = variance correction */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, EMPK, UNUR_ERR_GEN_INVALID );
  
  /* kernel variance known ? */
  if (! (gen->set & EMPK_SET_KERNELVAR) ) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"variance correction disabled");
    return UNUR_ERR_PAR_SET;
  }

  /* no changelog required */

  /* we use a bit in variant */
  gen->variant = (varcor) 
    ? (gen->variant | EMPK_VARFLAG_VARCOR) 
    : (gen->variant & (~EMPK_VARFLAG_VARCOR));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_empk_chg_varcor() */

/*---------------------------------------------------------------------------*/

int
unur_empk_set_positive( struct unur_par *par, int positive )
     /*----------------------------------------------------------------------*/
     /* turn mirroring (produces positive random numbers only) on/off        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par      ... pointer to parameter for building generator object    */
     /*   positive ... 0 = no mirroring,  !0 = mirroring                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no mirroring is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, EMPK );

  /* we use a bit in variant */
  par->variant = (positive) 
    ? (par->variant | EMPK_VARFLAG_POSITIVE) 
    : (par->variant & (~EMPK_VARFLAG_POSITIVE));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_empk_set_positive() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_empk_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to paramters for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;
  double iqrtrange;       /* interquartile range of data */
  double sigma;           /* estimation (guess) for "real" standard deviation
			     of observed data. */

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_EMPK ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_EMPK_PAR,NULL);

  /* check kernel generator */
  if ( PAR->kerngen==NULL && PAR->kernel == NULL) {
    /* no kernel given. use default kernel */
    if ( unur_empk_set_kernel( par, UNUR_DISTR_GAUSSIAN )!=UNUR_SUCCESS ) {
      /* cannot make kernel generator (should not happen!) */
      _unur_par_free(par); return NULL; 
    }
  }

  /* create a new empty generator object */
  gen = _unur_empk_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* if variance correction is used, the variance of the kernel
     must be known and positive */
  if( (gen->variant & EMPK_VARFLAG_VARCOR) &&
      !( (gen->set & EMPK_SET_KERNELVAR) && GEN->kernvar > 0. )) {
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"variance correction disabled");
    gen->variant &= ~EMPK_SET_KERNELVAR;
  }

  /* set uniform random number generator */
  GEN->kerngen->urng = par->urng;

  /* copy debugging flags */
  GEN->kerngen->debug = par->debug;

  /* the observed data */

  /* sort entries */
  /** TODO: this sort can be removed after we have implemented a new 
      version of the function iqrtrange, that does not depend on sorting **/
  qsort( GEN->observ, (size_t)GEN->n_observ, sizeof(double), compare_doubles);

  /* compute mean and standard deviation of observed sample */
  _unur_empk_comp_stddev( GEN->observ, GEN->n_observ, &(GEN->mean_observ), &(GEN->stddev_observ) );

  /* compute interquartile range of the sample */
  iqrtrange = _unur_empk_comp_iqrtrange( GEN->observ, GEN->n_observ );

  /* get an estimation (guess) of the "real" standard deviation of 
     the observed data.
     For normal distributed data it is 1.34 times the interquartile range.
     If this is greater than the observed standard deviation we use that 
     instead of.  */
  sigma = iqrtrange / 1.34;
  if (GEN->stddev_observ < sigma) sigma = GEN->stddev_observ;

  /* compute band width (also called window width) */
  GEN->bwidth_opt = GEN->alpha * GEN->beta * sigma / exp(0.2 * log((double)GEN->n_observ));
  GEN->bwidth = GEN->smoothing * GEN->bwidth_opt;

  /* compute constant for variance corrected version */
  GEN->sconst = 1./sqrt(1. + GEN->kernvar * SQU( GEN->bwidth/GEN->stddev_observ ) );

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_empk_debug_init(gen);
#endif

  /* free parameters */
  _unur_par_free(par);

  return gen;

} /* end of _unur_empk_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_empk_create( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* allocate memory for generator                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to (empty) generator object with default settings          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_EMPK_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_empk_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_EMPK_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_empk_getSAMPLE(gen);
  gen->destroy = _unur_empk_free;
  gen->clone = _unur_empk_clone;

  /* copy observed data into generator object */
  GEN->observ   = DISTR.sample;          /* observations in distribution object */
  GEN->n_observ = DISTR.n_sample;        /* sample size */

  /* copy some parameters into generator object */
  GEN->smoothing = PAR->smoothing;    /* smoothing factor                      */
  GEN->alpha = PAR->alpha;            /* parameter for kernel width            */
  GEN->beta = PAR->beta;              /* parameter for kernel width            */

  /* copy kernel generator into generator object */
  if (PAR->kerngen)
    /* kernel provided by user */
    GEN->kerngen = _unur_gen_clone(PAR->kerngen);
  else
    /* kernel from UNURAN list of kernels */
    GEN->kerngen = PAR->kernel;

  /* variance of kernel */
  GEN->kernvar = PAR->kernvar;

  /* the kernel is an auxilliary generator for method EMPK, of course */
  gen->gen_aux = GEN->kerngen;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_empk_info;
#endif

/*   RETURN pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_empk_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_empk_clone( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* copy (clone) generator object                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of generator object                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
#define CLONE  ((struct unur_empk_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_EMPK_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy observed data into generator object */
  CLONE->observ = clone->distr->data.cemp.sample;   /* observations in distribution object */

  /* kernel generator is (also) stored as auxiliary generator */
  /* which has already been cloned by generic_clone.          */
  CLONE->kerngen = clone->gen_aux;

  return clone;

#undef CLONE
} /* end of _unur_empk_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_empk_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 

  /* check arguments */
  if( !gen ) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_EMPK ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_EMPK_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_empk_free() */

/*****************************************************************************/

double
_unur_empk_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{ 
  double U,K,X;
  int j;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_EMPK_GEN,INFINITY);

  /* select uniformly one of the observations */
  U = _unur_call_urng(gen->urng) * GEN->n_observ;
  j = (int) (U);

  /** TODO: recycle uniform random variate??
      maybe for Boxcar (Uniform) Kernel. **/
  /* U -= j;   u is now a "recycled" U(0,1) random variate, aber wie weiter verwenden, ? */

  /* sample from kernel distribution */
  K = unur_sample_cont( GEN->kerngen );
  
  if (gen->variant & EMPK_VARFLAG_VARCOR)
    /* use variance correction */
    X = GEN->mean_observ + (GEN->observ[j] - GEN->mean_observ + GEN->bwidth * K) * GEN->sconst;
  else
    /* no variance correction */
    X = GEN->observ[j] + GEN->bwidth * K;

  if (gen->variant & EMPK_VARFLAG_POSITIVE)
    /* use mirroring to avoid non-positive numbers */
    X = (X<0.) ? -X : X;

  return X;

} /* end of _unur_empk_sample() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_empk_comp_stddev( double *data, int n_data, double *mean, double *stddev)
     /*----------------------------------------------------------------------*/
     /* compute mean and standard deviation of data                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   data   ... pointer to array of data                                */
     /*   n_data ... number of data points                                   */
     /*   mean   ... pointer to store mean                                   */
     /*   stddev ... pointer to store stddev                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* REFERENCES:                                                          */
     /*   [1] Spicer C.C. (1972): Algorithm AS 52: Calculation of Power Sums */
     /*       of Deviations about the mean,                                  */
     /*       Applied Statistics 21(2), pp. 226-227.                         */
     /*----------------------------------------------------------------------*/
{
  double xsqu_sum;   /* sum of x[i]^2 */
  double dx;
  int n;

  if (n_data < 2)
    /* cannot compute standard deviation */
    return 0;

  /* initialize counters */
  *mean = 0.;
  xsqu_sum = 0.;

  /* compute sums */
  for (n=1; n <= n_data; n++) {
    dx = (data[n-1] - *mean) / n;
    xsqu_sum += n * (n - 1.) * dx * dx;
    *mean += dx;
  }

  /* compute standard deviation */
  *stddev = sqrt( xsqu_sum / (n_data - 1.));

  return UNUR_SUCCESS;
} /* end of _unur_empk_comp_stddev() */

/*---------------------------------------------------------------------------*/

/** TODO: implement a new version of interquartilsrang that does not depend
    on sorting. (Only important if we want to use really large samples) **/
double
_unur_empk_comp_iqrtrange( double *data, int n )
     /*----------------------------------------------------------------------*/
     /* compute interquartile range of sorted data (in ascending order)      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   data   ... pointer to array of data                                */
     /*   n      ... number of data points                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   interquartile range                                                */
     /*----------------------------------------------------------------------*/
{
  double lowerqrt,upperqrt;  /* lower and upper quartile */
  int j;

  /* data must be sorted in ascending order */

  j = n/2;

  if (j % 2) {
    lowerqrt = data[(j+1)/2-1];
    upperqrt = data[n-(j+1)/2];
  }
  else {
    lowerqrt = (data[j/2-1] + data[j/2+1-1])/2.;
    upperqrt = (data[n-j/2] + data[n-j/2-1])/2.;
  }
  
  return (upperqrt - lowerqrt);

} /* end of _unur_empk_comp_iqrange() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_empk_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_EMPK_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = EMPK (EMPirical distribution with Kernel smoothing)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cemp_debug( gen->distr, gen->genid, (gen->debug & EMPK_DEBUG_PRINTDATA));

  fprintf(LOG,"%s: sampling routine = _unur_empk_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: smoothing factor = %g",gen->genid, GEN->smoothing);
  _unur_print_if_default(gen,EMPK_SET_SMOOTHING); fprintf(LOG,"\n");
  if (gen->variant & EMPK_VARFLAG_POSITIVE)
    fprintf(LOG,"%s: positive random variable only; use mirroring \n",gen->genid);

  if (gen->variant & EMPK_VARFLAG_VARCOR)
    fprintf(LOG,"%s: use variance correction\n",gen->genid);
  else
    fprintf(LOG,"%s: no variance correction\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: Kernel:\n",gen->genid);

  fprintf(LOG,"%s:    type = %s  ",gen->genid,GEN->kerngen->distr->name);
  if (gen->set & EMPK_SET_KERNGEN)
    fprintf(LOG,"[kernel generator set]\n");
  else if (gen->set & EMPK_SET_KERNEL)
    fprintf(LOG,"[standard kernel]\n");
  else 
    fprintf(LOG,"[default kernel]\n");

  fprintf(LOG,"%s:    window width = %g\t(opt = %g)\n",gen->genid, GEN->bwidth, GEN->bwidth_opt);
  fprintf(LOG,"%s:    alpha = %g",gen->genid, GEN->alpha);
  _unur_print_if_default(gen,EMPK_SET_ALPHA); fprintf(LOG,"\n");
  if (gen->variant & EMPK_VARFLAG_VARCOR) {
    fprintf(LOG,"%s:    kernel variance = %g",gen->genid, GEN->kernvar);
    _unur_print_if_default(gen,EMPK_SET_KERNELVAR); fprintf(LOG,"\n");
    fprintf(LOG,"%s:    variance correction factor = %g\n",gen->genid, GEN->sconst);
  }

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: Data:\n",gen->genid);
  fprintf(LOG,"%s:    beta  = %g",gen->genid, GEN->beta);
  _unur_print_if_default(gen,EMPK_SET_BETA); fprintf(LOG,"\n");
  fprintf(LOG,"%s:    mean (data) = %g\n",gen->genid, GEN->mean_observ);
  fprintf(LOG,"%s:    stddev (data) = %g\n",gen->genid, GEN->stddev_observ);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_empk_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_empk_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = DATA  [length=%d]\n", GEN->n_observ);
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
  /* method */
  _unur_string_append(info,"method: EMPK (EMPirical distribution with Kernel smoothing)\n");

  _unur_string_append(info,"   kernel type = %s  (alpha=%g)  ", GEN->kerngen->distr->name, GEN->alpha);
  if (gen->set & EMPK_SET_KERNGEN)
    _unur_string_append(info,"[kernel generator set]\n");
  else if (gen->set & EMPK_SET_KERNEL)
    _unur_string_append(info,"[standard kernel]\n");
  else 
    _unur_string_append(info,"[default kernel]\n");

  _unur_string_append(info,"   window width = %g  (opt = %g)\n", GEN->bwidth, GEN->bwidth_opt);

  _unur_string_append(info,"   smoothing factor = %g\n", GEN->smoothing);

  if (gen->variant & EMPK_VARFLAG_POSITIVE)
    _unur_string_append(info,"   positive random variable only; use mirroring\n");

  if (gen->variant & EMPK_VARFLAG_VARCOR)
    _unur_string_append(info,"   variance correction factor = %g\n", GEN->sconst);
  else
    _unur_string_append(info,"   no variance correction\n");
  _unur_string_append(info,"\n");

  /* performance */
  /*   _unur_string_append(info,"performance characteristics:\n"); */
  /*   _unur_string_append(info,"\n"); */

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   smoothing = %g   %s\n", GEN->smoothing,
			(gen->set & EMPK_SET_SMOOTHING) ? "" : "[default]");

    if (gen->set & EMPK_SET_BETA)
      _unur_string_append(info,"   beta = %g\n", GEN->beta);

    if (gen->variant & EMPK_VARFLAG_VARCOR) 
      _unur_string_append(info,"   varcor = on\n");

    if (gen->variant & EMPK_VARFLAG_POSITIVE) 
      _unur_string_append(info,"   positive = on\n");

    _unur_string_append(info,"\n");

    /* Not displayed:
       int unur_empk_set_kernel( UNUR_PAR *parameters, unsigned kernel);
       int unur_empk_set_kernelgen( UNUR_PAR *parameters, const UNUR_GEN *kernelgen, double alpha, double kernelvar );
    */
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_empk_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
