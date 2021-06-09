/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      cstd.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generators for standard distribution                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold             *
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
 * CSTD is a wrapper for special generator for Continuous univariate         *
 * STandarD distributions. It only works for distributions in the UNURAN     *
 * library of distributions or those continuous distributions that have      *
 * the inverse CDF implemented. Otherwise it refuses to work.                *
 *                                                                           *
 * It calls the initialization routine provided by the distribution object.  *
 * This routine has to do all setup steps for the special generator.         *
 * If no such routine is given (i.e. distr->init==NULL) or when it does not  *
 * implement the inversion method then the inverse CDF stored in the         *
 * distribution object (if available) is used.                               *
 * If neither is available, then unur_cstd_new() does not work and the       *
 * NULL pointer is returned instead of the pointer to a parameter object.    *
 *                                                                           *
 * Notice that using a truncated distribution (this can be constructed by    *
 * changing the default domain of a distribution by means of an              *
 * unur_distr_cont_set_domain() call) is only allowed if the inversion       *
 * method is used. Otherwise no parameter object is returned by the          *
 * unur_cstd_new() call.                                                     *
 *                                                                           *
 * Variants (different algorithms for the same distribution) are possible    *
 * and can be selected by unsigned integers using the                        *
 * unur_cstd_set_variant() call.                                             *
 * For possible variants see the generator files for each distribution in    *
 * the distributions directory. However the following are common to all      *
 * distributions:                                                            *
 *                                                                           *
 *    UNUR_STDGEN_DEFAULT   ... the default generator                        *
 *    UNUR_STDGEN_INVERSION ... the inversion method (if available)          *
 *    UNUR_STDGEN_FAST      ... the fasted available special generator       *
 *                                                                           *
 * unur_cstd_set_variant() return 0 if a variant is not implemented, and 1   *
 * otherwise. In the first case the selected variant is not changed.         *
 *                                                                           *
 * The domain of a (truncated) distribution can be changed without building  *
 * a new generator object by means of the unur_cstd_chg_truncated() call.    *
 * Notice that this only works when the inversion method is used. Otherwise  *
 * nothing happens to the domain and an error message is produced.           *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cont.h>
#include <distributions/unur_stddistr.h>
#include <distributions/unur_distributions_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "cstd.h"
#include "cstd_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define CSTD_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */
#define CSTD_DEBUG_CHG       0x00001000u   /* print changed parameters       */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define CSTD_SET_VARIANT          0x01u

/*---------------------------------------------------------------------------*/

#define GENTYPE "CSTD"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_cstd_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_cstd_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_cstd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_cstd_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_cstd_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* There are no sampling routines, since every distribution has its own.     */
/* Sampling routines are defined in ../distributions/ for each distributions.*/
/* double _unur_cstd_sample( UNUR_GEN *gen ); does not exist!                */
/*---------------------------------------------------------------------------*/

static double _unur_cstd_sample_inv( struct unur_gen *gen ); 
/*---------------------------------------------------------------------------*/
/* Generic inversion method.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_cstd_inversion_init( struct unur_par *par, struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Initialize special generator for inversion method.                        */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_chg_pdfparams( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print new (changed) parameters of distribution                            */
/*---------------------------------------------------------------------------*/

static void _unur_cstd_debug_chg_truncated( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print new (changed) domain of (truncated) distribution                    */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_cstd_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_cstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_cstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define CDF(x)    _unur_cont_CDF((x),(gen->distr))    /* call to CDF         */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_cstd_new( const struct unur_distr *distr )
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
  _unur_check_NULL(GENTYPE,distr,NULL);

  /* check distribution */
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.init == NULL && DISTR_IN.invcdf == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"init() for special generators or inverse CDF");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_cstd_par) );
  COOKIE_SET(par,CK_CSTD_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  par->method   = UNUR_METH_CSTD;   /* method                                */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* indicate default parameters           */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for initializing generator */
  par->init = _unur_cstd_init;

  return par;

} /* end of unur_cstd_new() */

/*****************************************************************************/

int 
unur_cstd_set_variant( struct unur_par *par, unsigned variant )
     /*----------------------------------------------------------------------*/
     /* set variant of method                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par     ... pointer to parameter for building generator object     */
     /*   variant ... indicator for variant                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  unsigned old_variant;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, par->distr, UNUR_ERR_NULL );
  _unur_check_par_object( par, CSTD );

  /* store date */
  old_variant = par->variant;
  par->variant = variant;

  /* check variant. run special init routine only in test mode */
  if ( (par->DISTR_IN.init == NULL || par->DISTR_IN.init(par,NULL)!=UNUR_SUCCESS) &&
       _unur_cstd_inversion_init(par,NULL)!=UNUR_SUCCESS ) {
    /* variant not valid */
    _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
    par->variant = old_variant;
    return UNUR_ERR_PAR_VARIANT;
  }

  /* changelog */
  par->set |= CSTD_SET_VARIANT;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end if unur_cstd_set_variant() */

/*---------------------------------------------------------------------------*/

int 
unur_cstd_chg_truncated( struct unur_gen *gen, double left, double right )
     /*----------------------------------------------------------------------*/
     /* change the left and right borders of the domain of the               */
     /* (truncated) distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen   ... pointer to generator object                              */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   the new boundary points may be +/- INFINITY                        */
     /*----------------------------------------------------------------------*/
{
  double Umin, Umax;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, CSTD, UNUR_ERR_GEN_INVALID );

  /* domain can only be changed for inversion method! */
  if ( ! GEN->is_inversion ) { 
    /* this is not the inversion method */
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"truncated domain for non inversion method");
    return UNUR_ERR_GEN_DATA;
  }

  /* CDF required ! */
  if (DISTR.cdf == NULL) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_DATA,"truncated domain, CDF required");
    return UNUR_ERR_GEN_DATA;
  }

  /* check new parameter for generator */
  /* (the truncated domain must be a subset of the domain) */
  if (left < DISTR.domain[0]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    left = DISTR.domain[0];
  }
  if (right > DISTR.domain[1]) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"truncated domain too large");
    right = DISTR.domain[1];
  }

  if (left >= right) {
    _unur_warning(NULL,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* compute umin and umax */
  Umin = (left > -INFINITY) ? CDF(left) : 0.;
  Umax = (right < INFINITY) ? CDF(right) : 1.;

  /* check result */
  if (Umin > Umax) {
    /* this is a serios error that should not happen */
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  if (_unur_FP_equal(Umin,Umax)) {
    /* CDF values very close */
    _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values very close");
    if (_unur_iszero(Umin) || _unur_FP_same(Umax,1.)) {
      /* this is very bad */
      _unur_warning(gen->genid,UNUR_ERR_DISTR_SET,"CDF values at boundary points too close");
      return UNUR_ERR_DISTR_SET;
    }
  }


  /* copy new boundaries into generator object */
  DISTR.trunc[0] = left;
  DISTR.trunc[1] = right;
  GEN->Umin = Umin;
  GEN->Umax = Umax;

  /* changelog */
  gen->distr->set |= UNUR_DISTR_SET_TRUNCATED;

  /* indicate that we have a truncated distribution.
     (do not have the standard domain any more) */
  gen->distr->set &= ~UNUR_DISTR_SET_STDDOMAIN;

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & CSTD_DEBUG_CHG) 
    _unur_cstd_debug_chg_truncated( gen );
#endif
  
  /* o.k. */
  return UNUR_SUCCESS;
  
} /* end of unur_cstd_chg_truncated() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_cstd_init( struct unur_par *par )
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

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_CSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_cstd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* run special init routine for generator */
  GEN->is_inversion = FALSE;   /* reset flag for inversion method */
  if ( (DISTR.init == NULL || DISTR.init(NULL,gen)!=UNUR_SUCCESS) &&
       _unur_cstd_inversion_init(NULL,gen)!=UNUR_SUCCESS ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    _unur_cstd_free(gen); return NULL; 
  }

  /* check parameters */
  if (_unur_cstd_check_par(gen) != UNUR_SUCCESS) {
    _unur_cstd_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_cstd_debug_init(gen);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_cstd_init() */

/*---------------------------------------------------------------------------*/

int
_unur_cstd_reinit( struct unur_gen *gen )
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
  int rcode;

  /* run special init routine for generator */
  GEN->is_inversion = FALSE;   /* reset flag for inversion method */
  if ( (DISTR.init == NULL || DISTR.init(NULL,gen)!=UNUR_SUCCESS) &&
       _unur_cstd_inversion_init(NULL,gen)!=UNUR_SUCCESS ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return UNUR_ERR_GEN_DATA;
  }

  /* check parameters */
  if ( (rcode = _unur_cstd_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug & CSTD_DEBUG_REINIT) 
      _unur_cstd_debug_chg_pdfparams( gen );
#endif

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_cstd_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_cstd_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_CSTD_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_cstd_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_CSTD_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = NULL;      /* will be set in _unur_cstd_init() */
  gen->destroy = _unur_cstd_free;
  gen->clone = _unur_cstd_clone;
  gen->reinit = _unur_cstd_reinit;

  /* defaults */
  GEN->gen_param = NULL;        /* parameters for the generator              */
  GEN->n_gen_param = 0;         /* (computed in special GEN->init()          */
  GEN->is_inversion = FALSE;    /* method not based on inversion             */
  GEN->sample_routine_name = NULL ;  /* name of sampling routine             */

  /* copy some parameters into generator object */
  GEN->Umin = 0.;               /* cdf at left boundary of domain            */
  GEN->Umax = 1.;               /* cdf at right boundary of domain           */

  /* GEN->is_inversion is set in _unur_cstd_inversion_init() */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_cstd_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_cstd_create() */

/*---------------------------------------------------------------------------*/

int
_unur_cstd_check_par( struct unur_gen *gen )
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
  /* domain valid for special generator ?? */
  if (!(gen->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    /* domain has been modified */
    gen->distr->set &= UNUR_DISTR_SET_TRUNCATED;
    DISTR.trunc[0] = DISTR.domain[0];
    DISTR.trunc[1] = DISTR.domain[1];

    if ( ! GEN->is_inversion ) {
      /* this is not the inversion method */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed for non inversion method");
      return UNUR_ERR_GEN_DATA;
    }

    if (DISTR.cdf == NULL) {
      /* using a truncated distribution requires a CDF */
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"domain changed, CDF required");
      return UNUR_ERR_GEN_DATA;
    }

    /* compute Umin and Umax */
    GEN->Umin = (DISTR.trunc[0] > -INFINITY) ? CDF(DISTR.trunc[0]) : 0.;
    GEN->Umax = (DISTR.trunc[1] < INFINITY)  ? CDF(DISTR.trunc[1]) : 1.;
  }

  return UNUR_SUCCESS;
} /* end of _unur_cstd_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_cstd_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_cstd_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_CSTD_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy parameters for special generators */
  if (GEN->gen_param) {
    CLONE->gen_param = _unur_xmalloc( GEN->n_gen_param * sizeof(double) );
    memcpy( CLONE->gen_param, GEN->gen_param, GEN->n_gen_param * sizeof(double) );
  }

  return clone;

#undef CLONE
} /* end of _unur_cstd_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_cstd_free( struct unur_gen *gen )
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

  /* magic cookies */
  COOKIE_CHECK(gen,CK_CSTD_GEN,RETURN_VOID);

  /* check input */
  if ( gen->method != UNUR_METH_CSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->gen_param)  free(GEN->gen_param);

  _unur_generic_free(gen);

} /* end of _unur_cstd_free() */

/*****************************************************************************/

/** 
    double _unur_cstd_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

/*****************************************************************************/

double
_unur_cstd_sample_inv( struct unur_gen *gen ) 
     /*----------------------------------------------------------------------*/
     /* generic inversion method.                                            */
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
  double U;

  if (!DISTR.invcdf) return INFINITY;

  /* sample from uniform random number generator */
  while (_unur_iszero(U = GEN->Umin + _unur_call_urng(gen->urng) * (GEN->Umax-GEN->Umin)));

  /* compute inverse CDF */
  return (DISTR.invcdf(U,gen->distr));

} /* _unur_cstd_sample_inv() */

/*---------------------------------------------------------------------------*/

double
unur_cstd_eval_invcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate inverse CDF at u.                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   u   ... argument for inverse CDF (0<=u<=1, no validation!)         */
     /*                                                                      */
     /* return:                                                              */
     /*   double (inverse CDF)                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double x;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( gen->method != UNUR_METH_CSTD ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY;
  }
  COOKIE_CHECK(gen,CK_CSTD_GEN,INFINITY);

  if (!DISTR.invcdf) {
    /* no inverse CDF available */
    _unur_error(gen->genid,UNUR_ERR_NO_QUANTILE,"inversion CDF required");
    return UNUR_INFINITY;
  } 

  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.trunc[0];
    if (u>=1.) return DISTR.trunc[1];
    return u;  /* = NaN */
  }
  
  /* rescale given u */
  u = GEN->Umin + u * (GEN->Umax - GEN->Umin);

  /* compute inverse CDF */
  x = DISTR.invcdf(u,gen->distr);

  /* validate range */
  if (x<DISTR.trunc[0]) x = DISTR.trunc[0];
  if (x>DISTR.trunc[1]) x = DISTR.trunc[1];

  return x;

} /* end of unur_cstd_eval_invcdf() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_cstd_inversion_init( struct unur_par *par, struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for inversion method                    */
     /* when inverse CDF is available for the distribution.                  */
     /* if gen == NULL then only check existance of variant.                 */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* one of par and gen must not be the NULL pointer */
  switch ((par) ? par->variant : gen->variant) {

  case 0:  /* DEFAULT */
  case UNUR_STDGEN_INVERSION:   /* inversion method */
    if (gen) {
      if (DISTR.invcdf) {
	GEN->is_inversion = TRUE;
	_unur_cstd_set_sampling_routine(gen,_unur_cstd_sample_inv);
	return UNUR_SUCCESS;
      }
    }
    else {
      if ((par->distr->data.cont).invcdf) {
	return UNUR_SUCCESS;
      }
    }

  default: /* no such generator */
    if (gen) _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return UNUR_FAILURE;
  }
  
} /* end of _unur_cstd_inversion_init() */

/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_cstd_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_CSTD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_cont_debug( gen->distr, gen->genid );

  /* sampling routine */
  fprintf(LOG,"%s: sampling routine = ",gen->genid);
  if (GEN->sample_routine_name)
    fprintf(LOG,"%s()",GEN->sample_routine_name);
  else
    fprintf(LOG,"(Unknown)");
  if (GEN->is_inversion)
    fprintf(LOG,"   (Inversion)");
  fprintf(LOG,"\n%s:\n",gen->genid);

  if (!(gen->distr->set & UNUR_DISTR_SET_STDDOMAIN)) {
    fprintf(LOG,"%s: domain has been changed. U in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);
    fprintf(LOG,"%s:\n",gen->genid);
  }

} /* end of _unur_cstd_debug_init() */

/*---------------------------------------------------------------------------*/

void 
_unur_cstd_debug_chg_pdfparams( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) parameters of distribution                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_CSTD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: parameters of distribution changed:\n",gen->genid);
  for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\tparam[%d] = %g\n",gen->genid,i,DISTR.params[i]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED)
    fprintf(LOG,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);

} /* end of _unur_cstd_debug_chg_pdfparams() */

/*---------------------------------------------------------------------------*/

void 
_unur_cstd_debug_chg_truncated( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print new (changed) domain of (truncated) distribution               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_CSTD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: domain of truncated distribution changed:\n",gen->genid);
  fprintf(LOG,"%s:\tdomain = (%g, %g)\n",gen->genid, DISTR.trunc[0], DISTR.trunc[1]);
  fprintf(LOG,"%s:\tU in (%g,%g)\n",gen->genid,GEN->Umin,GEN->Umax);

} /* end of _unur_cstd_debug_chg_truncated() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_cstd_info( struct unur_gen *gen, int help )
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
  int samplesize = 10000;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

  /* method */
  _unur_string_append(info,"method: CSTD (special generator for Continuous STandarD distribution)\n");
  _unur_string_append(info,"   variant = %d  %s\n", gen->variant,
		      (GEN->is_inversion)?"[implements inversion method]" : "");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#urn] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   variant = %d  %s\n", gen->variant,
			(gen->set & CSTD_SET_VARIANT) ? "" : "[default]");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_cstd_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
