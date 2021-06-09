/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      empl.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generate from empirical CDF with linear interpolation        *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given observed sample.                                               *
 *      Produce a value x consistent with this sample                        *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to sample                                                    *
 *                                                                           *
 *   OPTIONAL:                                                               *
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
 *   [1] Law, A.M. and Kelton, D. (2000): Simulation Modeling and Analysis.  *
 *       New York: McGraw-Hill.                                              *
 *                                                                           *
 *   [2] Bratley, P., Fox, B.L., and Schrage E.L. (1987): A Guide to         *
 *       Simulation. New York: Springer-Verlag.                              *
 *                                                                           *
 *   [3] Hoermann, W. and Leydold, J. (2000): Automatic random variate       *
 *       generation for simulation input. Proceedings of the 2000 Winter     *
 *       Simulation Conference. (??? eds.)                                   *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "empl.h"
#include "empl_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define EMPL_DEBUG_PRINTDATA   0x00000100u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "EMPL"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empl_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empl_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_empl_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_empl_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_empl_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_empl_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_empl_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cemp      /* data for distribution object      */

#define PAR       ((struct unur_empl_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_empl_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cemp /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

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

#define _unur_empl_getSAMPLE(gen)   (_unur_empl_sample)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_empl_new( const struct unur_distr *distr )
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
  par = _unur_par_new( sizeof(struct unur_empl_par) );
  COOKIE_SET(par,CK_EMPL_PAR);

  /* copy input */
  par->distr    = distr;          /* pointer to distribution object          */

  /* set default values */

  par->method   = UNUR_METH_EMPL; /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init     = _unur_empl_init;

  return par;

} /* end of unur_empl_new() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_empl_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_EMPL ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_EMPL_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_empl_create(par);
  _unur_par_free(par);
  if (!gen) return NULL; 

  /* the observed data */

  /* sort entries */
  qsort( GEN->observ, (size_t)GEN->n_observ, sizeof(double), compare_doubles);

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_empl_debug_init(gen);
#endif

  return gen;

} /* end of _unur_empl_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_empl_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_EMPL_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_empl_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_EMPL_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_empl_getSAMPLE(gen);
  gen->destroy = _unur_empl_free;
  gen->clone = _unur_empl_clone;

  /* copy observed data into generator object */
  GEN->observ   = DISTR.sample;          /* observations in distribution object */
  GEN->n_observ = DISTR.n_sample;        /* sample size */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_empl_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_empl_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_empl_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_empl_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_EMPL_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy observed data into generator object */
  CLONE->observ = clone->distr->data.cemp.sample;   /* observations in distribution object */

  return clone;

#undef CLONE
} /* end of _unur_empl_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_empl_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_EMPL ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_EMPL_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_empl_free() */

/*****************************************************************************/

double
_unur_empl_sample( struct unur_gen *gen )
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
  double U,X;
  int J;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_EMPL_GEN,INFINITY);

  /* select uniformly an interval */
  U = _unur_call_urng(gen->urng) * (GEN->n_observ-1);
  J = (int) (U);

  /* linear CDF between obeservation points */
  X = GEN->observ[J] + (U-J)*(GEN->observ[J+1] - GEN->observ[J]);

  return X;

} /* end of _unur_empl_sample() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_empl_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_EMPL_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = EMPL (EMPirical distribution with Linear interpolation)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cemp_debug( gen->distr, gen->genid, (gen->debug & EMPL_DEBUG_PRINTDATA));

  fprintf(LOG,"%s: sampling routine = _unur_empl_sample()\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_empl_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_empl_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"method: EMPL (EMPirical distribution with Linear interpolation)\n");
  _unur_string_append(info,"\n");

  /* performance */
  /*   _unur_string_append(info,"performance characteristics:\n"); */
  /*   _unur_string_append(info,"\n"); */

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters: none\n");
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_empl_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
