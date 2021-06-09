/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvstd.c                                                      *
 *                                                                           *
 *   TYPE:      multivariate continuous random variate                       *
 *   METHOD:    generators for standard distribution                         *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold             *
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
 * MVSTD is a wrapper for special generator for MultiVariate continuous      *
 * STandarD distributions. It only works for distributions in the UNURAN     *
 * library of distributions and will refuse to work otherwise.               *
 * (In detail it rejects a distribution if its id is equal to DISTR_GENERIC, *
 * the id inserted my the unur_distr_cvec_new() call.)                       *
 *                                                                           *
 * It calls the initialzation routine provided by the distribution object.   *
 * This routine has to do all setup steps for the special generator.         *
 * If no such routine is given, i.e. distr->init==NULL, then                 *
 * unur_mvstd_new() does not work and the NULL pointer is returned instead   *
 * of the pointer to a parameter object.                                     *
 *                                                                           *
 * Notice that using a truncated distribution (this can be constructed by    *
 * changing the default domain of a distribution by means of an              *
 * unur_distr_cvec_set_domain_rect() call) is not possible.                  *
 * method is used. Otherwise no parameter object is returned by the          *
 * unur_cstd_new() call.                                                     *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distributions/unur_stddistr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "mvstd.h"
#include "mvstd_struct.h"

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

#define MVSTD_DEBUG_REINIT   0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MVSTD"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvstd_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_mvstd_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvstd_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_mvstd_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mvstd_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_mvstd_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* There are sampling routines, since every distribution has its own.        */
/* Sampling routines are defined in ../distributions/ for each distributions.*/
/* double _unur_mvstd_sample( UNUR_GEN *gen ); does not exist!               */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mvstd_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_mvstd_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_mvstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mvstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mvstd_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (distr->id == UNUR_DISTR_GENERIC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"standard distribution");
    return NULL;
  }
  if (DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"init() for special generators");
    return NULL;
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mvstd_par) );
  COOKIE_SET(par,CK_MVSTD_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  par->method   = UNUR_METH_MVSTD;  /* method                                */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* indicate default parameters           */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for initializing generator */
  par->init = _unur_mvstd_init;

  return par;

} /* end of unur_mvstd_new() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_mvstd_init( struct unur_par *par )
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
  /* check for required data: initializing routine for special generator */
  if (par->DISTR_IN.init == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_NULL,"");
    return NULL;
  }

  /* check input */
  if ( par->method != UNUR_METH_MVSTD ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_MVSTD_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mvstd_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  
  /* run special init routine for generator */
  if ( DISTR.init(gen)!=UNUR_SUCCESS ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"variant for special generator");
    _unur_mvstd_free(gen); return NULL; 
  }

  /* check parameters */
  if (_unur_mvstd_check_par(gen) != UNUR_SUCCESS) {
    _unur_mvstd_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_mvstd_debug_init(gen);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_mvstd_init() */

/*---------------------------------------------------------------------------*/

int
_unur_mvstd_reinit( struct unur_gen *gen )
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
  /* run special init routine for generator */
  if ( DISTR.init(gen)!=UNUR_SUCCESS ) {
    /* init failed --> could not find a sampling routine */
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"parameters");
    return UNUR_ERR_GEN_DATA;
  }

  /* check parameters */
  return _unur_mvstd_check_par(gen);
} /* end of _unur_mvstd_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mvstd_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MVSTD_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mvstd_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MVSTD_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = NULL;      /* will be set in _unur_mvstd_init() */
  gen->destroy = _unur_mvstd_free;
  gen->clone = _unur_mvstd_clone;
  gen->reinit = _unur_mvstd_reinit;

  /* defaults */
  GEN->sample_routine_name = NULL ;  /* name of sampling routine             */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_mvstd_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_mvstd_create() */

/*---------------------------------------------------------------------------*/

int
_unur_mvstd_check_par( struct unur_gen *gen )
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
  /* we do not allow truncated multinormal distributions */
  if (gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"truncated domain");
    return UNUR_ERR_GEN_CONDITION;
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_mvstd_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mvstd_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_mvstd_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MVSTD_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_mvstd_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_mvstd_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_MVSTD_GEN,RETURN_VOID);

  /* check input */
  if ( gen->method != UNUR_METH_MVSTD ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_mvstd_free() */

/*****************************************************************************/

/** 
    double _unur_mvstd_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are defined in ../distributions/ for each distributions.
**/

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
_unur_mvstd_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MVSTD_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = generator for standard distribution\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_cvec_debug( gen->distr, gen->genid );

  /* sampling routine */
  fprintf(LOG,"%s: sampling routine = ",gen->genid);
  if (GEN->sample_routine_name)
    fprintf(LOG,"%s()\n",GEN->sample_routine_name);
  else
    fprintf(LOG,"(Unknown)\n");
  fprintf(LOG,"%s:\n",gen->genid);

  fflush(LOG);

} /* end of _unur_mvstd_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_mvstd_info( struct unur_gen *gen, int help )
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
  int dim = gen->distr->dim;
  int samplesize = 10000;
  double E_urn;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",dim);
  _unur_distr_cvec_info_domain(gen);
  _unur_string_append(info,"\n\n");
  
  /*   if (help) { */
  /*   _unur_string_append(info,"\n"); */
  /*   } */

  /* method */
  _unur_string_append(info,"method: MVSTD (special generator for MultiVariate continuous STandarD distribution)\n");
  /*   _unur_string_append(info,"   variant = %d\n", gen->variant); */
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");

  E_urn = unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize);
  _unur_string_append(info,"   E [#urn] = %.2f x %d = %.2f  [approx.]\n",
		      E_urn / dim, dim, E_urn);
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    /*     _unur_string_append(info,"   variant = %d  %s\n", gen->variant, */
    /* 			(gen->set & MVSTD_SET_VARIANT) ? "" : "[default]"); */
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_mvstd_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
