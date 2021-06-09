/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      hrd.c                                                        *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    Hazard Rate Decreasing                                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Rejection with dynamic thinning.                                     *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to the hazard rate                                           *
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
 *   Evrim Ozgul (2002): The generation of random variates with a given      *
 *      hazard rate, M.Sc. thesis, Department of Industrial Engineering,     *
 *      Bogazici University, Istanbul.                                       *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Rejection from majorizing hazard rate using dynamic thinning.           *
 *   This means that rejection from a constant hazard rate is used.          *
 *   When a point Xi is generated but rejected a new constant majorizing     *
 *   hazard rate with value HR(Xi) starting at Xi is used.                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "hrd.h"
#include "hrd_struct.h"

#ifdef UNUR_ENABLE_INFO
#  include <tests/unuran_tests.h>
#endif

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define HRD_VARFLAG_VERIFY     0x01u    /* flag for verifying mode           */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define HRD_DEBUG_REINIT    0x00000010u   /* print parameters after reinit  */
#define HRD_DEBUG_SAMPLE       0x01000000u    /* trace sampling
						 (only if verify mode is on) */
                                                
/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "HRD"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hrd_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_hrd_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hrd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_hrd_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hrd_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_hrd_free( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_hrd_sample( struct unur_gen *gen );
static double _unur_hrd_sample_check( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_hrd_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_hrd_debug_sample( const struct unur_gen *gen, double x, int i );
/*---------------------------------------------------------------------------*/
/* trace sampling.                                                           */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_hrd_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cont      /* data for distribution object      */

#define PAR       ((struct unur_hrd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_hrd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define HR(x)     _unur_cont_HR((x),(gen->distr))   /* call to hazard rate   */

/*---------------------------------------------------------------------------*/

#define _unur_hrd_getSAMPLE(gen) \
   ( ((gen)->variant & HRD_VARFLAG_VERIFY) \
     ? _unur_hrd_sample_check : _unur_hrd_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_hrd_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CONT,NULL);

  if (DISTR_IN.hr == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"HR"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_hrd_par) );
  COOKIE_SET(par,CK_HRD_PAR);

  /* copy input */
  par->distr   = distr;         /* pointer to distribution object            */

  /* set default values */

  par->method   = UNUR_METH_HRD;  /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                      /* inidicate default parameters   */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_hrd_init;

  return par;

} /* end of unur_hrd_new() */

/*****************************************************************************/

int
unur_hrd_set_verify( struct unur_par *par, int verify )
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
  _unur_check_par_object( par, HRD );

  /* we use a bit in variant */
  par->variant = (verify) ? (par->variant | HRD_VARFLAG_VERIFY) : (par->variant & (~HRD_VARFLAG_VERIFY));

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hrd_set_verify() */

/*---------------------------------------------------------------------------*/

int
unur_hrd_chg_verify( struct unur_gen *gen, int verify )
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
  _unur_check_gen_object( gen, HRD, UNUR_ERR_GEN_INVALID );

  /* we must not change this switch when sampling has been disabled by
     using a pointer to the error producing routine                          */
  if (SAMPLE == _unur_sample_cont_error) 
    return UNUR_FAILURE;

  /* we use a bit in variant */
  gen->variant = (verify) 
    ? (gen->variant | HRD_VARFLAG_VERIFY) 
    : (gen->variant & (~HRD_VARFLAG_VERIFY));

  /* sampling routine */
  SAMPLE = _unur_hrd_getSAMPLE(gen);

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_hrd_chg_verify() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_hrd_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   params  pointer to paramters for building generator object         */
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
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_HRD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HRD_PAR,NULL);

  /* create a new empty generator object */    
  gen = _unur_hrd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if (_unur_hrd_check_par(gen) != UNUR_SUCCESS) {
    _unur_hrd_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_hrd_debug_init(gen);
#endif
  
  return gen;

} /* end of _unur_hrd_init() */

/*---------------------------------------------------------------------------*/

int
_unur_hrd_reinit( struct unur_gen *gen )
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

  /* check parameters */
  if ( (rcode = _unur_hrd_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* (re)set sampling routine */
  SAMPLE = _unur_hrd_getSAMPLE(gen);

  /* nothing to do */
  return UNUR_SUCCESS;
} /* end of _unur_hrd_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hrd_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HRD_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_hrd_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_HRD_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_hrd_getSAMPLE(gen);
  gen->destroy = _unur_hrd_free;
  gen->clone = _unur_hrd_clone;
  gen->reinit = _unur_hrd_reinit;

  /* default values */

  /* initialize variables */
  GEN->left_border = 0.;             /* left border of domain                 */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_hrd_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_hrd_create() */

/*---------------------------------------------------------------------------*/

int
_unur_hrd_check_par( struct unur_gen *gen )
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
  /* set left border and check domain */
  if (DISTR.domain[0] < 0.)       DISTR.domain[0] = 0.;
  if (DISTR.domain[1] < INFINITY) DISTR.domain[1] = INFINITY;
  GEN->left_border = DISTR.domain[0];

  /* compute upper bound for hazard rate (at left border) */
  GEN->upper_bound = HR(GEN->left_border);
  if (GEN->upper_bound <= 0. || _unur_FP_is_infinity(GEN->upper_bound)) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"no valid upper bound for HR at left boundary");
    return UNUR_ERR_GEN_CONDITION;
  }

  return UNUR_SUCCESS;
} /* end of _unur_hrd_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hrd_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_hrd_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HRD_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_hrd_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_hrd_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_HRD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HRD_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_hrd_free() */

/*****************************************************************************/

double
_unur_hrd_sample( struct unur_gen *gen )
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
  double U,V,E,X,hrx;
  double lambda;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HRD_GEN,INFINITY);

  /* parameter for majorizing hazard rate */
  lambda = GEN->upper_bound;

  /* starting point */
  X = GEN->left_border;

  for(;;) {
    /* sample from U(0,1) */
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );

    /* sample from exponential distribution with scale parameter lambda */
    E = -log(U) / lambda;

    /* Remark: by this construction E is monotically increasing with
       the uniform random number generated by the _unur_call_urng() call */

    /* next step */
    X += E;

    /* hazard rate at generated point */
    hrx = HR(X);

    /* reject or accept */
    V =  lambda * _unur_call_urng(gen->urng);
    if( V <= hrx ) 
      return X;      /* accept */

    /* else: reject */

    /* update majorizing hazard rate */
    if (hrx > 0.)
      lambda = hrx;   
    else {
      /* the given function is not a hazard rate */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not valid");
      return INFINITY;
    } 
  }

} /* end of _unur_hrd_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_hrd_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (verify mode)                                  */
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
  double U,V,E,X,hrx;
  double lambda;
  int i;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HRD_GEN,INFINITY);

  /* parameter for majorizing hazard rate */
  lambda = GEN->upper_bound;

  /* starting point */
  X = GEN->left_border;

  for(i=1;;i++) {
    /* sample from U(0,1) */
    while ( _unur_iszero(U = 1.-_unur_call_urng(gen->urng)) );

    /* sample from exponential distribution with scale parameter lambda */
    E = -log(U) / lambda;

    /* Remark: by this construction E is monotically increasing with
       the uniform random number generated by the _unur_call_urng() call */

    /* next step */
    X += E;

    /* hazard rate at generated point */
    hrx = HR(X);

    /* verify upper bound */
    if ( (1.+UNUR_EPSILON) * lambda < hrx )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not decreasing");

    /* reject or accept */
    V =  lambda * _unur_call_urng(gen->urng);
    if( V <= hrx ) {
#ifdef UNUR_ENABLE_LOGGING
      /* write info into LOG file */
      if (gen->debug & HRD_DEBUG_SAMPLE)
	_unur_hrd_debug_sample( gen, X, i );
#endif
      return X;
    }

    /* else: reject */

    /* update majorizing hazard rate */
    if (hrx > 0.)
      lambda = hrx;   
    else {
      /* the given function is not a hazard rate */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"HR not valid");
      return INFINITY;
    } 
  }

} /* end of _unur_hrd_sample_check() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_hrd_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HRD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = HRD (Hazard Rate Decreasing)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_hrd_sample",gen->genid);
  if (gen->variant & HRD_VARFLAG_VERIFY)
    fprintf(LOG,"_check()\n");
  else
    fprintf(LOG,"()\n");
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_hrd_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_hrd_debug_sample( const struct unur_gen *gen, double x, int i )
     /*----------------------------------------------------------------------*/
     /* write info about generated point into LOG file                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   x   ... generated point                                            */
     /*   i   ... number of iterations                                       */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HRD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: X = %g\t #iterations = %d\n",gen->genid,x,i);

} /* end of _unur_hrd_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_hrd_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   functions = HR\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: HRD (Hazard Rate Decreasing)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E[#iterations] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
    if (help) {
      _unur_string_append(info,"parameters:\n");
      if (gen->variant & HRD_VARFLAG_VERIFY)
	_unur_string_append(info,"   verify = on\n");
      _unur_string_append(info,"\n");
    }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_hrd_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
