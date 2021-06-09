/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dext.c                                                       *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    wrapper for external generator                               *
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
 * DEXT is a wrapper that allows to use external generators within the       *
 * framework of UNURAN.                                                      *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "dext.h"
#include "dext_struct.h"

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

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "DEXT"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dext_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dext_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dext_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dext_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_dext_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* There are sampling routines, since every distribution has its own.        */
/* Sampling routines are defined in ../distributions/ for each distributions.*/
/* double _unur_dext_sample( UNUR_GEN *gen ); does not exist!                */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dext_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_dext_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr     /* data for distribution object      */

#define PAR       ((struct unur_dext_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dext_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_dext_new( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object (can be NULL!)            */
     /*                                                                      */
     /* return:                                                              */
     /*   default parameters (pointer to structure)                          */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_par *par;

  /* check distribution */
  if (distr != NULL) {
    if (distr->type != UNUR_DISTR_DISCR) {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
    COOKIE_CHECK(distr,CK_DISTR_DISCR,NULL);
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_dext_par) );
  COOKIE_SET(par,CK_DEXT_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  PAR->init     = NULL;   /* pointer to initialization routine for external generator */
  PAR->sample   = NULL;   /* pointer to sampling routine for external generator */
  
  par->method   = UNUR_METH_DEXT;   /* method                                */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* inidicate default parameters          */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for initializing generator */
  par->init = _unur_dext_init;

  return par;

} /* end of unur_dext_new() */

/*****************************************************************************/

int
unur_dext_set_init( struct unur_par *par, int (*init)(struct unur_gen *gen) )
     /*----------------------------------------------------------------------*/
     /* set initialization routine for external generator                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par   ... pointer to parameter for building generator object       */
     /*   init  ... pointer to initialization routine for external generator */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, DEXT );

  /* store date */
  PAR->init = init;

  return UNUR_SUCCESS;

} /* end if unur_dext_set_init() */

/*---------------------------------------------------------------------------*/

int
unur_dext_set_sample( struct unur_par *par, int (*sample)(struct unur_gen *gen) )
     /*----------------------------------------------------------------------*/
     /* set sampling routine for external generator                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   sample ... pointer to sampling routine for external generator      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, sample, UNUR_ERR_NULL );

  /* check input */
  _unur_check_par_object( par, DEXT );

  /* store date */
  PAR->sample = sample;

  return UNUR_SUCCESS;

} /* end if unur_dext_set_sample() */

/*---------------------------------------------------------------------------*/

void *
unur_dext_get_params( struct unur_gen *gen, size_t size )
     /*----------------------------------------------------------------------*/
     /* Get pointer to memory block for storing parameters of external       */
     /* generator. The memory block is (re-) allocated if necessary.         */
     /* If size is set 0, then only the pointer stored in the generator      */
     /* object is returned.                                                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   size ... size if the memory block                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to memory block that contains parameters                   */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  _unur_check_NULL( GENTYPE, gen, NULL );
  COOKIE_CHECK(gen, CK_DEXT_GEN, NULL);

  /* Is it really necessary to check the method? No! */
  /*   if ( gen->method != UNUR_METH_DEXT ) {  */
  /*     _unur_error((gen)->genid,UNUR_ERR_GEN_INVALID,""); */
  /*     return NULL;  */
  /*   } */

  if (size && size != GEN->size_param) {
    /* allocate memory block */
    GEN->param = _unur_xrealloc(GEN->param, size);
    GEN->size_param = size;
  }

  return GEN->param;
} /* end of unur_dext_get_params() */


/*---------------------------------------------------------------------------*/

double  *
unur_dext_get_distrparams( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Get pointer to array of parameters of underlying distribution        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to double array                                            */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  CHECK_NULL(gen, NULL);
  COOKIE_CHECK(gen, CK_DEXT_GEN, NULL);

  return DISTR.params;
} /* end of unur_dext_get_distrparams() */

/*---------------------------------------------------------------------------*/

int
unur_dext_get_ndistrparams( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Get size of array of parameters of underlying distribution           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*                                                                      */
     /* return:                                                              */
     /*   length of double array                                             */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check input */
  CHECK_NULL(gen, 0);
  COOKIE_CHECK(gen, CK_DEXT_GEN, 0);

  return DISTR.n_params;
} /* end of unur_dext_get_ndistrparams() */

/*---------------------------------------------------------------------------*/


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_dext_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_DEXT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL;
  }
  COOKIE_CHECK(par,CK_DEXT_PAR,NULL);

  /* we need a sampling routine */
  if (PAR->sample == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_CONDITION,"sampling routine missing");
    return NULL;
  }

  /* create a new empty generator object */
  gen = _unur_dext_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;
  
  /* run special init routine for external generator */
  if (GEN->init != NULL) {
    if (GEN->init(gen) != UNUR_SUCCESS) {
      /* init failed --> could not find a sampling routine */
      _unur_error(GENTYPE,UNUR_FAILURE,"init for external generator failed");
      _unur_dext_free(gen); return NULL;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_dext_debug_init(gen);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_dext_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dext_reinit( struct unur_gen *gen )
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
  /* run special init routine for external generator */
  if (GEN->init != NULL) {
    if (GEN->init(gen) != UNUR_SUCCESS) {
      /* init failed --> could not find a sampling routine */
      _unur_error(GENTYPE,UNUR_FAILURE,"init for external generator failed");
      return UNUR_FAILURE;
    }
  }

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_dext_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dext_create( struct unur_par *par )
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
  struct unur_distr *distr = NULL;

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DEXT_PAR,NULL);

  /* need a distribution object */
  if (par->distr == NULL)
    par->distr = distr = unur_distr_discr_new();

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_dext_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DEXT_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = PAR->sample;      /* will be set in _unur_dext_init() */
  gen->destroy = _unur_dext_free;
  gen->clone = _unur_dext_clone;
  gen->reinit = _unur_dext_reinit;

  /* copy data */
  GEN->init = PAR->init;
  GEN->sample = PAR->sample;

  /* defaults */
  GEN->param    = NULL;   /* parameters for the generator                    */
  GEN->size_param  = 0;   /* size of parameter object                        */

  /* clean up */
  if (distr) _unur_distr_free(distr);

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_dext_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_dext_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dext_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_dext_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DEXT_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy parameters for special generators */
  if (GEN->param) {
    CLONE->param = _unur_xmalloc( GEN->size_param );
    memcpy( CLONE->param, GEN->param, GEN->size_param );
  }

  return clone;

#undef CLONE
} /* end of _unur_dext_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_dext_free( struct unur_gen *gen )
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
  COOKIE_CHECK(gen,CK_DEXT_GEN,RETURN_VOID);

  /* check input */
  if ( gen->method != UNUR_METH_DEXT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return;
  }

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->param)  free(GEN->param);

  _unur_generic_free(gen);

} /* end of _unur_dext_free() */

/*****************************************************************************/

/**
    double _unur_dext_sample( struct unur_gen *gen ) {}
    Does not exists !!!
    Sampling routines are set by user via unur_dext_set_sample().
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
_unur_dext_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DEXT_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = wrapper for external generator\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  /* distribution */
  _unur_distr_discr_debug( gen->distr, gen->genid, FALSE );

} /* end of _unur_dext_info_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_dext_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
  /* method */
  _unur_string_append(info,"method: DEXT (wrapper for Discrete EXTernal generators)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#urn] = %.2f  [approx.]\n",
		      unur_test_count_urn(gen,samplesize,0,NULL)/((double)samplesize));
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
    /* Not displayed:
       int unur_dext_set_init( UNUR_PAR *parameters, int (*init)(UNUR_GEN *gen) );
       int unur_dext_set_sample( UNUR_PAR *parameters, int (*sample)(UNUR_GEN *gen) );
    */
  }


  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_dext_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
