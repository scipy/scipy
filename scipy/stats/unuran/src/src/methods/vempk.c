/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vempk.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    generate from kernel estimation                              *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given observed sample.                                               *
 *      Produce a vector x consistent with this sample                       *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to sample                                                    *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      smoothing factor                                                     *
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
 *      bwidth_opt = ????????????????????                                    *
 *                                                                           *
 *   Currently only Gaussian kernels with the same covariance matrix as the  *
 *   given sample are implemented.                                           *
 *                                                                           *
 *...........................................................................*
 *                                                                           *
 *   Algorithm VEMPK:                                                        *
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
 *   The mean of the empirical distribution is equal to the sample mean.     *
 *   One disadvantage of this method lies in the fact that the variance of   *
 *   the empirical distribution generated by the kernel method is always     *
 *   higher than the sample standard-deviation s. As pointed out in [4]      *
 *   it is not difficult to correct this error. For the variance             *
 *   corrected version we have to compute the mean mux and the               *
 *   covariance matrix of the sample and we have to know the covariance      *
 *   matrix of the kernel V(K). Then we can replace step (3) of the          *
 *   algorithm above by:                                                     *
 *                                                                           *
 *   3': Return ??????????????                                               *
 *                                                                           *
 *   We can turn on or off variance correction with the function             *
 *                                                                           *
 *      unur_vempk_set_varcor();                                             *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "mvstd.h"
#include "vempk.h"
#include "vempk_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define VEMPK_VARFLAG_VARCOR    0x001u   /* use variance correction          */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define VEMPK_DEBUG_PRINTDATA   0x00000100u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define VEMPK_SET_SMOOTHING      0x008u    /* smoothing factor                */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VEMPK"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vempk_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vempk_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vempk_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_vempk_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_vempk_sample_cvec( struct unur_gen *gen, double *result );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int compute_mean_covar( double *data, int n_data, int dim, double *xbar, double *S );
/*---------------------------------------------------------------------------*/
/* compute mean vector xbar and covariance matrix S of data.                 */
/*---------------------------------------------------------------------------*/

/*  inline static double _unur_empk_comp_iqrtrange( double *data, int n_data ); */
/*---------------------------------------------------------------------------*/
/* compute interquartile range.                                              */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_vempk_debug_init( const struct unur_par *par, const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_vempk_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvemp      /* data for distribution object      */

#define PAR       ((struct unur_vempk_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_vempk_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvemp /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/

#define _unur_vempk_getSAMPLE(gen)  ( _unur_vempk_sample_cvec )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_vempk_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CVEMP) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEMP,NULL);

  if (DISTR_IN.sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"observed sample"); return NULL; }
  if (DISTR_IN.n_sample < 2) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"size of observed sample"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_vempk_par) );
  COOKIE_SET(par,CK_VEMPK_PAR);

  /* copy input */
  par->distr    = distr;         /* pointer to distribution object           */

  /* set default values */
  PAR->smoothing = 1.;            /* determines how "smooth" the estimated density will be */

  par->method   = UNUR_METH_VEMPK; /* method                                 */
  par->variant  = 0u;              /* default variant                        */

  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init     = _unur_vempk_init;

  return par;

} /* end of unur_vempk_new() */

/*****************************************************************************/

int
unur_vempk_set_smoothing( struct unur_par *par, double smoothing )
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
  _unur_check_par_object( par, VEMPK );

  /* check new parameter for generator */
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->smoothing = smoothing;

  /* changelog */
  par->set |= VEMPK_SET_SMOOTHING;

  return UNUR_SUCCESS;

} /* end of unur_vempk_set_smoothing() */

/*---------------------------------------------------------------------------*/

int
unur_vempk_chg_smoothing( struct unur_gen *gen, double smoothing )
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
  _unur_check_gen_object( gen, VEMPK, UNUR_ERR_GEN_INVALID );
  
  /* no changelog required */

  /* check new parameter for generator */
  if (smoothing < 0.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"smoothing factor < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store smoothing factor */
  GEN->smoothing = smoothing;

  /* recompute band width */
  GEN->hact = GEN->hopt * GEN->smoothing;

  /* recompute constant for variance corrected version */
  GEN->corfac = 1./sqrt( 1. + GEN->hact * GEN->hact);

  /* changelog */
  gen->set |= VEMPK_SET_SMOOTHING;

  return UNUR_SUCCESS;

} /* end of unur_vempk_chg_smoothing() */

/*---------------------------------------------------------------------------*/

int
unur_vempk_set_varcor( struct unur_par *par, int varcor )
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
  _unur_check_par_object( par, VEMPK );

  /* we use a bit in variant */
  par->variant = (varcor) 
    ? (par->variant | VEMPK_VARFLAG_VARCOR) 
    : (par->variant & (~VEMPK_VARFLAG_VARCOR));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vempk_set_varcor() */

/*---------------------------------------------------------------------------*/

int
unur_vempk_chg_varcor( struct unur_gen *gen, int varcor )
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
  _unur_check_gen_object( gen, VEMPK, UNUR_ERR_GEN_INVALID );
  
  /* no changelog required */

  /* we use a bit in variant */
  gen->variant = (varcor) 
    ? (gen->variant | VEMPK_VARFLAG_VARCOR) 
    : (gen->variant & (~VEMPK_VARFLAG_VARCOR));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_vempk_chg_varcor() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_vempk_init( struct unur_par *par )
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
  double *S;                  /* covariance matrix of sample */
  UNUR_DISTR *kernel_distr;   /* kernel distribution */

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_VEMPK ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VEMPK_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_vempk_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* compute mean vector and covariance matrix of sample */
  GEN->xbar = _unur_xmalloc( GEN->dim * sizeof(double) );
  S  = _unur_xmalloc( GEN->dim * GEN->dim * sizeof(double) );
  compute_mean_covar( DISTR.sample, DISTR.n_sample, GEN->dim, GEN->xbar, S );
  
  /* make a distribution object for multinormal distribution */
  kernel_distr = unur_distr_multinormal( GEN->dim, NULL, S );

  /* create kernel generator */
  GEN->kerngen = unur_init( unur_mvstd_new( kernel_distr ) );
  if (GEN->kerngen==NULL) {
    _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    _unur_par_free (par); free (S); _unur_vempk_free(gen);
    return NULL;
  }

  /* set uniform random number generator */
  GEN->kerngen->urng = par->urng;

  /* copy debugging flags */
  GEN->kerngen->debug = par->debug;

  /* the kernel is an auxiliary generator for method VEMPK, of course */
  gen->gen_aux = GEN->kerngen;
  
  /* compute bandwith:
     Silverman (1986), p.86, formular 4.14 for normal kernel */
  GEN->hopt = (exp((1./(GEN->dim+4.))*log(4./(GEN->dim+2.)))*
	      exp(-1./(GEN->dim+4.) *log((double)GEN->n_observ)));

  /* hopt is optimal for the multi-normal distribution
     as it is oversmoothing especially multimodal distributions */
  
  /* compute used band width */
  GEN->hact = GEN->hopt * GEN->smoothing;
  
  /* compute constant for variance corrected version */
  GEN->corfac = 1./sqrt(1. + GEN->hact * GEN->hact);
  
#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_vempk_debug_init(par,gen);
#endif
  
  /* free parameters */
  _unur_par_free(par);

  /* we do not need mu and S any more */
  free(S);
  unur_distr_free(kernel_distr);

  return gen;

} /* end of _unur_vempk_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_vempk_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VEMPK_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_vempk_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_VEMPK_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim; 

  /* copy observed data into generator object */
  GEN->observ   = DISTR.sample;          /* observations in distribution object */
  GEN->n_observ = DISTR.n_sample;        /* sample size */

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_vempk_getSAMPLE(gen);
  gen->destroy = _unur_vempk_free;
  gen->clone = _unur_vempk_clone;

  /* copy some parameters into generator object */
  GEN->smoothing = PAR->smoothing;    /* smoothing factor                      */

  /* initialize pointer */
  GEN->kerngen = NULL;               /* generator for kernel distribution     */
  GEN->xbar = NULL;                  /* mean vector of sample                 */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_vempk_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_vempk_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vempk_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_vempk_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VEMPK_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy additional data for generator object */
  CLONE->observ = clone->distr->data.cvemp.sample;   /* observations in distribution object */

  if (GEN->xbar) {
    CLONE->xbar = _unur_xmalloc( GEN->dim * sizeof(double) );
    memcpy( CLONE->xbar, GEN->xbar, GEN->dim * sizeof(double) );
  }

  /* kernel generator is (also) stored as auxiliary generator */
  /* which has already been cloned by generic_clone.          */
  CLONE->kerngen = clone->gen_aux;

  return clone;

#undef CLONE
} /* end of _unur_vempk_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_vempk_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_VEMPK ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VEMPK_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->xbar)   free( GEN->xbar );

  _unur_generic_free(gen);

} /* end of _unur_vempk_free() */

/*****************************************************************************/

int
_unur_vempk_sample_cvec( struct unur_gen *gen, double *result )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   result ... random vector (result)                                  */
     /*----------------------------------------------------------------------*/
{ 
  #define idx(a,b) (a*GEN->dim+b)
  double U;
  int j,k;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);  COOKIE_CHECK(gen,CK_VEMPK_GEN,UNUR_ERR_COOKIE);

  /* select uniformly one of the observations */
  U = _unur_call_urng(gen->urng) * GEN->n_observ;
  j = (int) (U);

  /* sample from kernel distribution */
  unur_sample_vec( GEN->kerngen, result );

  if (gen->variant & VEMPK_VARFLAG_VARCOR)
    /* use variance correction */
    for (k=0; k<GEN->dim; k++) 
      result[k] = GEN->xbar[k] + (GEN->observ[idx(j,k)] - GEN->xbar[k] + result[k]*GEN->hact) * GEN->corfac;
  
  else
    /* no variance correction */
    for (k=0; k<GEN->dim; k++) 
      result[k] = GEN->hact * result[k] + GEN->observ[idx(j,k)];
  
  return UNUR_SUCCESS;
#undef idx
} /* end of _unur_vempk_sample() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
compute_mean_covar( double *data, int n_data, int dim, 
		    double *xbar, double *S ) 
     /*----------------------------------------------------------------------*/
     /* compute mean vector xbar and covariance matrix S of data             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   data   ... pointer to array of data                                */
     /*   n_data ... number of data points                                   */
     /*   xbar   ... pointer to store mean vector                            */
     /*   S      ... pointer to store covariance matrix                      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) (a*dim+b)

  int i,j,k;
  double *x;
  
  /* allocate working array */
  x = malloc(dim*sizeof(double));
  
  /* initialization */
  
  for(j=0; j<dim; j++) {
    xbar[j] = 0.;
    for(k=0; k<dim; k++)
      S[idx(j,k)] = 0.;
  }
  
  /* the mean vector xbar of all observations is computed */
  for (i=0; i<n_data; i++)
    for(j=0; j<dim; j++)
      xbar[j] += data[idx(i,j)];
  for(j=0; j<dim; j++)
    xbar[j] /= n_data;
  
  /* the variance covariance matrix S of all observations is computed */
  for (i=0; i<n_data; i++) {
    for(j=0; j<dim; j++)
      x[j] = data[idx(i,j)] - xbar[j];
    for(j=0; j<dim; j++)
      for(k=0; k<=j; k++) 
	S[idx(j,k)] += x[j] * x[k];
  }

  for (j=dim-1; j>=0; j--)
    for (k=0; k<=j; k++) {
      S[idx(j,k)] /= (n_data-1);
      if (k!=j)
	S[idx(k,j)] = S[idx(j,k)];   /* covariance matrix must be symmetric */
    }
  
  /* free working array */
  free(x);

  /* o.k. */
  return UNUR_SUCCESS;

} /* compute_mean_covar() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_vempk_debug_init( const struct unur_par *par, const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_VEMPK_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VEMPK_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = VEMPK ((Vector) EMPirical distribution with Kernel smoothing)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);
  _unur_distr_cvemp_debug( gen->distr, gen->genid, (gen->debug & VEMPK_DEBUG_PRINTDATA));

  fprintf(LOG,"%s:\tmean vector =\n",gen->genid);
  fprintf(LOG,"%s:\t   ( %g",gen->genid,GEN->xbar[0]);
  for (i=1; i<GEN->dim; i++) 
    fprintf(LOG,", %g",GEN->xbar[i]);
  fprintf(LOG,")\n%s:\n",gen->genid);

  fprintf(LOG,"%s:\tcovariance matrix = [see %s]\n",gen->genid, GEN->kerngen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: sampling routine = _unur_vempk_sample_cvec()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: smoothing factor = %g",gen->genid, PAR->smoothing);
  _unur_print_if_default(par,VEMPK_SET_SMOOTHING); fprintf(LOG,"\n");
  fprintf(LOG,"%s:\n",gen->genid);


  fprintf(LOG,"%s: bandwith hopt = %g\n",gen->genid, GEN->hopt);
  fprintf(LOG,"%s: (used)   hact = %g\n",gen->genid, GEN->hact);
  fprintf(LOG,"%s:\n",gen->genid);

  if (gen->variant & VEMPK_VARFLAG_VARCOR) {
    fprintf(LOG,"%s: use variance correction\n",gen->genid);
    fprintf(LOG,"%s:\tcorrection factor = %g\n",gen->genid, GEN->corfac);
  }
  else
    fprintf(LOG,"%s: no variance correction\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_vempk_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_vempk_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = DATA  [length=%d]\n", GEN->n_observ);
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
  /* method */
  _unur_string_append(info,"method: VEMPK (EMPirical distribution with Kernel smoothing)\n");

  _unur_string_append(info,"   kernel type = multinormal\n");
  _unur_string_append(info,"   smoothing factor = %g\n", GEN->smoothing);
  _unur_string_append(info,"   bandwith = %g\n", GEN->hact);

  if (gen->variant & VEMPK_VARFLAG_VARCOR)
    _unur_string_append(info,"   variance correction factor = %g\n", GEN->corfac);
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
			(gen->set & VEMPK_SET_SMOOTHING) ? "" : "[default]");
    if (gen->variant & VEMPK_VARFLAG_VARCOR)
      _unur_string_append(info,"   varcor = on\n");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_vempk_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
