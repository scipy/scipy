/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dsrou.c                                                      *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    discrete simple universal method (ratio-of-uniforms method)  *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given PMF and mode of a T_{-1/2}-concave distribution                *
 *      produce a value X consistent with its PMF                            *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the mass function                                         *
 *      mode of distribution                                                 *
 *      sum over PMF                                                         *
 *                                                                           *
 *   OPTIONAL:                                                               *
 *      CDF at mode                                                          *
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
 *   [1] Leydold J. (2001): A simple universal generator for continuous and  *
 *       discrete univariate T-concave distributions,                        *
 *       ACM Trans. Math. Software 27(1), pp. 66--82.                        *
 *                                                                           *
 *   [2] Kinderman, A.J. and Monahan, F.J. (1977): Computer generation of    *
 *       random variables using the ratio of uniform deviates,               *
 *       ACM Trans. Math. Software 3(3), pp. 257--260.                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * The ratio-of-uniforms method introduced in [2] is a flexible method that  *
 * is based on the following theorem:                                        *
 *                                                                           *
 * THEOREM:                                                                  *
 *    Let X be a random variable with density function f(x) = g(x) / G,      *
 *    where g(x) is a positive integrable function with support (x_0,x_1)    *
 *    not necessarily finite and G = integral g(x) dx.                       *
 *    If (V,U) is uniformly distributed in                                   *
 *       A = {(v,u): 0 < u <= sqrt(g(v/u)), x_0 < v/u < x_1},                *
 *    then X = V/U has probability density function f(x).                    *
 *                                                                           *
 * Generating point (V,U) uniformly distributed in A is done by rejection    *
 * from an enveloping region, usually from the minimal bounding rectangle.   *
 *                                                                           *
 * For discrete random variates the continuous PDF                           *
 *    PDF(x) = PMF(floor(x))                                                 *
 * is used.                                                                  *
 *                                                                           *
 * The implemented algorithm uses the fact, that for many distributions,     *
 * the polygon having the "spikes" of A as its verticeds is convex.          *
 * Then we can find the follow bounding rectangles:                          *
 * (For simplicity we have assumed that the sum over the PMF is 1)           *
 *                                                                           *
 * Define                                                                    *
 *    R = {(v,u):  -1/u_l <= v <= 0, 0 <= u <= u_l} \cup                     *
 *        {(v,u):  0 <= v <= 1/u_r, 0 <= u <= u_r}                           *
 *    Q = {(v,u):  v_l <= v <= 0, 0 <= u <= u_l} \cup                        *
 *        {(v,u):  0 <= v <= v_r, 0 <= u <= u_r}                             *
 * where                                                                     *
 *    u_l = sqrt(PMF(mode-1)), u_r = sqrt(PMF(mode)),                        *
 *    v_l = -F(mode-1)/u_l, v_r = (1-F(mode-1))/u_r                          *
 * Then                                                                      *
 *    A subset R subset Q                                                    *
 *                                                                           *
 * Thus we can use R to generate whenever the CDF F(mode-1) at the mode      *
 * is known, and Q otherwise.                                                *
 * Notice, that the rection constant is 2 in the first case and 4 and the    *
 * latter.                                                                   *
 *                                                                           *
 * Distributions with a convex set A are characterized by the following      *
 * theorem that shows a connection to transformed density rejection TDR.     *
 *                                                                           *
 * THEOREM:                                                                  *
 *    A is convex if and only if g is T-concave with transformation          *
 *    T(x) = -1/sqrt(x), i.e., -1/sqrt(g(x)) is a concave function.          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dsrou.h"
#include "dsrou_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define DSROU_VARFLAG_VERIFY   0x002u  /* run verify mode                    */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DSROU_DEBUG_REINIT    0x00000010u  /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DSROU_SET_CDFMODE     0x001u   /* CDF at mode is known               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "DSROU"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dsrou_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dsrou_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dsrou_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dsrou_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dsrou_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_dsrou_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dsrou_sample( struct unur_gen *gen );
static int _unur_dsrou_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_dsrou_rectangle( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute universal bounding rectangle.                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_dsrou_debug_init( const struct unur_gen *gen, int is_reinit );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_dsrou_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       ((struct unur_dsrou_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dsrou_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */     

#define PMF(x)    _unur_discr_PMF((x),(gen->distr))   /* call to PMF         */

/*---------------------------------------------------------------------------*/

#define _unur_dsrou_getSAMPLE(gen) \
   ( ((gen)->variant & DSROU_VARFLAG_VERIFY) \
     ? _unur_dsrou_sample_check : _unur_dsrou_sample )

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

#define SQRT2     (M_SQRT2)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_dsrou_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_DISCR) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);

  if (DISTR_IN.pmf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PMF"); 
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_dsrou_par) );
  COOKIE_SET(par,CK_DSROU_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  PAR->Fmode     = -1.;                /* CDF at mode (unknown yet)           */

  par->method   = UNUR_METH_DSROU;    /* method and default variant          */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dsrou_init;

  return par;

} /* end of unur_dsrou_new() */

/*****************************************************************************/

int 
unur_dsrou_set_cdfatmode( struct unur_par *par, double Fmode )
     /*----------------------------------------------------------------------*/
     /* set value of cdf at mode                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   Fmode  ... CDF at mode                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DSROU );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
   }
 
  /* store date */
  PAR->Fmode = Fmode;

  /* changelog */
  par->set |= DSROU_SET_CDFMODE;

  return UNUR_SUCCESS;
 
} /* end of unur_dsrou_set_cdfatmode() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_set_verify( struct unur_par *par, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DSROU );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | DSROU_VARFLAG_VERIFY) : (par->variant & (~DSROU_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_ERR_NULL;
} /* end of unur_dsrou_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_dsrou_chg_verify( struct unur_gen *gen, int verify )
     /*----------------------------------------------------------------------*/
     /* turn verifying of algorithm while sampling on/off                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   verify ... 0 = no verifying,  !0 = verifying                       */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no verifying is the default                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );
 
  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_discr_error) 
    return UNUR_FAILURE;

  if (verify)
    /* turn verify mode on */
    gen->variant |= DSROU_VARFLAG_VERIFY;
  else
    /* turn verify mode off */
    gen->variant &= ~DSROU_VARFLAG_VERIFY;

  SAMPLE = _unur_dsrou_getSAMPLE(gen); 

  /* o.k. */
  return UNUR_SUCCESS;
 
} /* end of unur_dsrou_chg_verify() */

/*****************************************************************************/

int
unur_dsrou_chg_cdfatmode( struct unur_gen *gen, double Fmode )
     /*----------------------------------------------------------------------*/
     /* change value of cdf at mode                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen    ... pointer to generator object                             */
     /*   Fmode  ... CDF at mode                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, DSROU, UNUR_ERR_GEN_INVALID );

  /* check new parameter for generator */
  if (Fmode < 0. || Fmode > 1.) {
    _unur_warning(gen->genid,UNUR_ERR_PAR_SET,"CDF(mode)");
    return UNUR_ERR_PAR_SET;
  }
  
  /* copy parameters */
  GEN->Fmode = Fmode;

  /* changelog */
  gen->set |= DSROU_SET_CDFMODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_dsrou_chg_cdfatmode() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_dsrou_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params ... pointer to paramters for building generator object      */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_DSROU ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DSROU_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_dsrou_create(par);

  /* free parameters */
  _unur_par_free(par);

  if (!gen) return NULL;

  /* check parameters */
  if (_unur_dsrou_check_par(gen) != UNUR_SUCCESS) {
    _unur_dsrou_free(gen); return NULL;
  }

  /* compute universal bounding rectangle */
  if ( _unur_dsrou_rectangle(gen)!=UNUR_SUCCESS ) {
    _unur_dsrou_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_dsrou_debug_init(gen, FALSE);
#endif

  return gen;

} /* end of _unur_dsrou_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dsrou_reinit( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* re-initialize (existing) generator.                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* check parameters */
  if ( (result = _unur_dsrou_check_par(gen)) != UNUR_SUCCESS)
    return result;

  /* compute universal bounding rectangle */
  if ( (result = _unur_dsrou_rectangle(gen)) != UNUR_SUCCESS)
    return result;

  /* (re)set sampling routine */
  SAMPLE = _unur_dsrou_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & DSROU_DEBUG_REINIT) _unur_dsrou_debug_init(gen,TRUE);
#endif

  return UNUR_SUCCESS;
} /* end of _unur_dsrou_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dsrou_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DSROU_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_dsrou_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DSROU_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_dsrou_getSAMPLE(gen);
  gen->destroy = _unur_dsrou_free;
  gen->clone = _unur_dsrou_clone;
  gen->reinit = _unur_dsrou_reinit;

  /* copy some parameters into generator object */
  GEN->Fmode = PAR->Fmode;            /* CDF at mode                           */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_dsrou_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dsrou_create() */

/*---------------------------------------------------------------------------*/

int
_unur_dsrou_check_par( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* check parameters of given distribution and method                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check for required data: mode */
  if (!(gen->distr->set & UNUR_DISTR_SET_MODE)) {
    _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode: try finding it (numerically)"); 
    if (unur_distr_discr_upd_mode(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mode"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* check for required data: sum over PMF */
  if (!(gen->distr->set & UNUR_DISTR_SET_PMFSUM))
    if (unur_distr_discr_upd_pmfsum(gen->distr)!=UNUR_SUCCESS) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"sum over PMF");
      return UNUR_ERR_DISTR_REQUIRED;
    }

  /* mode must be in domain */
  if ( (DISTR.mode < DISTR.BD_LEFT) ||
       (DISTR.mode > DISTR.BD_RIGHT) ) {
    /* there is something wrong.
       assume: user has change domain without changing mode.
       but then, she probably has not updated area and is to large */
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,"area and/or CDF at mode");
    DISTR.mode = _unur_max(DISTR.mode,DISTR.BD_LEFT);
    DISTR.mode = _unur_min(DISTR.mode,DISTR.BD_RIGHT);
  }

  return UNUR_SUCCESS;
} /* end of _unur_dsrou_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dsrou_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_dsrou_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DSROU_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_dsrou_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_dsrou_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_DSROU ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DSROU_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_dsrou_free() */

/*****************************************************************************/

int
_unur_dsrou_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   int (sample from random variate)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,X;
  int I;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DSROU_GEN,INT_MAX);

  while (1) {
    /* generate point uniformly in union of rectangles */
    V = GEN->al + _unur_call_urng(gen->urng) * (GEN->ar - GEN->al);
    V /= (V<0.) ? GEN->ul : GEN->ur;    /* if ul==0. then al==0. and thus V>=0. */

    while ( _unur_iszero(U = _unur_call_urng(gen->urng)));
    U *= (V<0.) ? GEN->ul : GEN->ur;

    /* ratio */
    X = floor(V/U) + DISTR.mode;

    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;

    /* convert to int */
    I = (int) X;

    /* accept or reject */
    if (U*U <= PMF(I))
      return I;
  }
} /* end of _unur_dsrou_sample() */

/*---------------------------------------------------------------------------*/

int
_unur_dsrou_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify that method can be used             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   int (sample from random variate)                                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{ 
  double U,V,pI,VI,X;
  double um2, vl, vr;
  int I;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DSROU_GEN,INT_MAX);

  while (1) {
    /* generate point uniformly in union of rectangles */
    V = GEN->al + _unur_call_urng(gen->urng) * (GEN->ar - GEN->al);
    V /= (V<0.) ? GEN->ul : GEN->ur;

    while ( _unur_iszero(U = _unur_call_urng(gen->urng)));
    U *= (V<0.) ? GEN->ul : GEN->ur;

    /* ratios */
    X = floor(V/U) + DISTR.mode;

    /* inside domain ? */
    if ( (X < DISTR.BD_LEFT) || (X > DISTR.BD_RIGHT) )
      continue;

    /* convert to int */
    I = (int) X;

    /* values of PMF and v-coordinate of point */
    pI = PMF(I);
    VI = V/U * sqrt(pI);

    /* values of boundary of rectangle          */
    /* (avoid roundoff error with FP registers) */
    um2 = (2.+4.*DBL_EPSILON) * ((V<0) ? GEN->ul*GEN->ul : GEN->ur*GEN->ur);
    vl = (GEN->ul>0.) ? (1.+UNUR_EPSILON) * GEN->al/GEN->ul : 0.;
    vr = (1.+UNUR_EPSILON) * GEN->ar/GEN->ur;

    /* check hat */
    if ( pI > um2 || VI < vl || VI > vr ) {
      /*        printf("pI = %g < %g     VI = %g < %g < %g\n", */
      /*  	     pI, ur2, vl, VI, vr); */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PMF(x) > hat(x)");
    }

    /* accept or reject */
    if (U*U <= pI)
      return I;
  }
} /* end of _unur_dsrou_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_dsrou_rectangle( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding rectangle                                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double pm, pbm;               /* PMF at mode and mode-1                    */

  /* check arguments */
  CHECK_NULL( gen, UNUR_ERR_NULL );
  COOKIE_CHECK( gen,CK_DSROU_GEN, UNUR_ERR_COOKIE );

  /* compute PMF at mode and mode-1 */
  pm = PMF(DISTR.mode);
  pbm = (DISTR.mode-1 < DISTR.BD_LEFT) ? 0. : PMF(DISTR.mode-1);

  /* pm and pbm must be positive */
  if (pm <= 0. || pbm < 0.) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"PMF(mode) <= 0.");
    return UNUR_ERR_GEN_DATA;
  }

  /* heights of rectangles */
  GEN->ul = sqrt(pbm);
  GEN->ur = sqrt(pm);

  /* areas of rectangle */
  if (_unur_iszero(GEN->ul)) {
    /* PMF monotonically decreasing */
    GEN->al = 0.;
    GEN->ar = DISTR.sum;
  }
  else if (gen->set & DSROU_SET_CDFMODE) {
    /* CDF at mode known */
    GEN->al = -(GEN->Fmode * DISTR.sum)+pm;
    GEN->ar = DISTR.sum + GEN->al;
  }
  else {
    GEN->al = -(DISTR.sum - pm);
    GEN->ar = DISTR.sum;
  }    

  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of _unur_dsrou_rectangle() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_dsrou_debug_init( const struct unur_gen *gen, int is_reinit )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen       ... pointer to generator object                          */
     /*   is_reinit ... if TRUE the generator has been reinitialized         */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DSROU_GEN,RETURN_VOID);

  LOG = unur_get_stream();
  
  fprintf(LOG,"%s:\n",gen->genid);
  if (!is_reinit) {
    fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
    fprintf(LOG,"%s: method  = dsrou (discrete simple universal ratio-of-uniforms)\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: reinit!\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);
 
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );
  
  fprintf(LOG,"%s: sampling routine = _unur_dsrou_sample",gen->genid);
  if (gen->variant & DSROU_VARFLAG_VERIFY)
    fprintf(LOG,"_check");
  fprintf(LOG,"()\n%s:\n",gen->genid);
 
  if (gen->set & DSROU_SET_CDFMODE)
    fprintf(LOG,"%s: CDF(mode) = %g\n",gen->genid,GEN->Fmode);
  else
    fprintf(LOG,"%s: CDF(mode) unknown\n",gen->genid);

  fprintf(LOG,"%s: no (universal) squeeze\n",gen->genid);
  fprintf(LOG,"%s: no mirror principle\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: Rectangles:\n",gen->genid);
  if (GEN->ul > 0.)
    fprintf(LOG,"%s:    left upper point  = (%g,%g) \tarea = %g   (%5.2f%%)\n",
	    gen->genid,GEN->al/GEN->ul,GEN->ul,fabs(GEN->al),100.*fabs(GEN->al)/(-GEN->al+GEN->ar));
  else
    fprintf(LOG,"%s:    left upper point  = (0,0) \tarea = 0   (0.00%%)\n",gen->genid);

  fprintf(LOG,"%s:    right upper point = (%g,%g) \tarea = %g   (%5.2f%%)\n",
	  gen->genid,GEN->ar/GEN->ur,GEN->ur,GEN->ar,100.*GEN->ar/(-GEN->al+GEN->ar));

  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_dsrou_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_dsrou_info( struct unur_gen *gen, int help )
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
  struct unur_distr *distr = gen->distr;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PMF\n");
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   mode      = %d   %s\n", DISTR.mode,
                      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   sum(PMF)  = %g\n", DISTR.sum);
  if (gen->set & DSROU_SET_CDFMODE)
    _unur_string_append(info,"   F(mode)   = %g\n", GEN->Fmode);
  else
    _unur_string_append(info,"   F(mode)   = [unknown]\n");
  _unur_string_append(info,"\n");

  if (help) {
    if ( distr->set & UNUR_DISTR_SET_MODE_APPROX ) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may provide the \"mode\"");
    _unur_string_append(info,"\n");
  }

  /* method */
  _unur_string_append(info,"method: DSROU (Discrete Simple Ratio-Of-Uniforms)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");

  _unur_string_append(info,"   enveloping rectangle = (%g,%g) x (%g,%g)  [left]\n",
		      (GEN->ul > 0.)?GEN->al/GEN->ul:0., 0.,
		      0., (GEN->ul > 0.)?GEN->ul:0.);
  _unur_string_append(info,"                          (%g,%g) x (%g,%g)  [right]\n",
		      0.,GEN->ar/GEN->ur,  0., GEN->ur);

  _unur_string_append(info,"   area(hat) = %g + %g = %g\n",
		      fabs(GEN->al), GEN->ar, -GEN->al+GEN->ar);

  _unur_string_append(info,"   rejection constant = %g\n",
		      2. * (-GEN->al+GEN->ar) / DISTR.sum);
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    if (gen->set & DSROU_SET_CDFMODE)
      _unur_string_append(info,"   cdfatmode = %g\n", GEN->Fmode);
    else
      _unur_string_append(info,"   cdfatmode = [not set]\n");

    if (gen->variant & DSROU_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");

    _unur_string_append(info,"\n");
  }

  /* Hints */
  if (help) {
    if ( !(gen->set & DSROU_SET_CDFMODE))
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"cdfatmode\" to reduce the rejection constant.");
    _unur_string_append(info,"\n");
  }

} /* end of _unur_dsrou_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
