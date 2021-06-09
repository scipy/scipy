/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      hist.c                                                       *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    generate from histogram                                      *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given histogram of observed sample.                                  *
 *      Produce a value x consistent with this histogram (use inversion)     *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointer to histogram                                                 *
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

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen_source.h"
#include "hist.h"
#include "hist_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants: none                                                            */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define HIST_DEBUG_PRINTHIST   0x00000100u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "HIST"         /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hist_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hist_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_hist_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_hist_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_hist_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_hist_create_tables( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create (allocate) tables                                                  */
/*---------------------------------------------------------------------------*/

static int _unur_hist_make_guidetable( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create table for indexed search                                           */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/
static void _unur_hist_debug_init( const struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_hist_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cemp      /* data for distribution object      */

#define PAR       ((struct unur_hist_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_hist_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cemp /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/
/* constants                                                                 */

/*---------------------------------------------------------------------------*/

#define _unur_hist_getSAMPLE(gen)   (_unur_hist_sample)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_hist_new( const struct unur_distr *distr )
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

  if (DISTR_IN.hist_prob == NULL || !(distr->set & UNUR_DISTR_SET_DOMAIN)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"histogram"); return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_hist_par) );
  COOKIE_SET(par,CK_HIST_PAR);

  /* copy input */
  par->distr    = distr;          /* pointer to distribution object          */

  /* set default values */
  par->method   = UNUR_METH_HIST; /* method                                  */
  par->variant  = 0u;             /* default variant                         */

  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init     = _unur_hist_init;

  return par;

} /* end of unur_hist_new() */

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_hist_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_HIST ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_HIST_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_hist_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* compute guide table */
  if ( (_unur_hist_create_tables(gen) != UNUR_SUCCESS) ||
       (_unur_hist_make_guidetable(gen) != UNUR_SUCCESS) ) {
    _unur_hist_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_hist_debug_init(gen);
#endif

  return gen;

} /* end of _unur_hist_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hist_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_HIST_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_hist_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_HIST_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_hist_getSAMPLE(gen);
  gen->destroy = _unur_hist_free;
  gen->clone = _unur_hist_clone;

  /* make sure that the domain coincides with bin data      */
  if (DISTR.hist_bins) {
    DISTR.hmin = DISTR.hist_bins[0];
    DISTR.hmax = DISTR.hist_bins[DISTR.n_hist];
  }

  /* copy observed data into generator object */
  GEN->n_hist = DISTR.n_hist;      /* size of histogram     */
  GEN->prob   = DISTR.hist_prob;   /* probabilities of bins */
  GEN->hmin   = DISTR.hmin;        /* lower ...             */
  GEN->hmax   = DISTR.hmax;        /* ... and upper bound   */
  GEN->hwidth = (DISTR.hmax - DISTR.hmin) / DISTR.n_hist;
  GEN->bins   = (DISTR.hist_bins) ? DISTR.hist_bins : NULL;

  /* set all pointers to NULL */
  GEN->sum = 0.;
  GEN->cumpv = NULL;
  GEN->guide_table = NULL;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_hist_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_hist_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_hist_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_hist_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_HIST_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy histrogram into generator object */
  CLONE->prob = clone->distr->data.cemp.hist_prob;   /* probabilities of bins */
  CLONE->bins = clone->distr->data.cemp.hist_bins;   /* location of bins */

  /* copy data for distribution */
  CLONE->cumpv = _unur_xmalloc( GEN->n_hist * sizeof(double) );
  memcpy( CLONE->cumpv, GEN->cumpv, GEN->n_hist * sizeof(double) );
  CLONE->guide_table = _unur_xmalloc( GEN->n_hist * sizeof(int) );
  memcpy( CLONE->guide_table, GEN->guide_table, GEN->n_hist * sizeof(int) );

  return clone;

#undef CLONE
} /* end of _unur_hist_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_hist_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_HIST ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_HIST_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free two auxiliary tables */
  if (GEN->guide_table) free(GEN->guide_table);
  if (GEN->cumpv)       free(GEN->cumpv);

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_hist_free() */

/*****************************************************************************/

double
_unur_hist_sample( struct unur_gen *gen )
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
  double U;
  int J;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_HIST_GEN,INFINITY);

  /* sample from U(0,1) */
  U = _unur_call_urng(gen->urng);

  /* look up in guide table ... */
  J = GEN->guide_table[(int)(U * GEN->n_hist)];
  /* ... and search */
  U *= GEN->sum;
  while (GEN->cumpv[J] < U) J++;

  /* reuse of uniform random number: U ~ U(0,1) */
  U = (U - (J ? GEN->cumpv[J-1] : 0.)) / GEN->prob[J];

  /* sample uniformly from bin */
  if (GEN->bins) 
    return (U * GEN->bins[J+1] + (1.-U) * GEN->bins[J]);
  else
    return (GEN->hmin + (U+J)*GEN->hwidth);

} /* end of _unur_hist_sample() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_hist_create_tables( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create (allocate) tables                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  /* allocation for cummulated probabilities */
  GEN->cumpv = _unur_xrealloc( GEN->cumpv, GEN->n_hist * sizeof(double) );

  /* allocate memory for the guide table */
  GEN->guide_table = _unur_xrealloc( GEN->guide_table, GEN->n_hist * sizeof(int) );

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_hist_create_tables() */

/*---------------------------------------------------------------------------*/

int
_unur_hist_make_guidetable( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* create guide table                                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{ 
  double *pv;                   /* pointer to probability vector */
  int n_pv;                     /* length of probability vector */
  double pvh;                   /* aux variable for computing cumulated sums */
  double gstep;                 /* step size when computing guide table */
  int i,j;
  
  /* probability vector */
  pv = GEN->prob;
  n_pv = GEN->n_hist;

  /* computation of cumulated probabilities */
  for( i=0, pvh=0.; i<n_pv; i++ ) {
    GEN->cumpv[i] = ( pvh += pv[i] );
    /* ... and check probability vector */
    if (pv[i] < 0.) {
      _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"probability < 0");
      return UNUR_ERR_GEN_DATA;
    }
  }
  GEN->sum = GEN->cumpv[n_pv-1];

  /* computation of guide-table */
  gstep = GEN->sum / GEN->n_hist;
  pvh = 0.;
  for( j=0, i=0; j<GEN->n_hist;j++ ) {
    while (GEN->cumpv[i] < pvh) 
      i++;
    if (i >= n_pv) {
      _unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
      break;
    }
    GEN->guide_table[j] = i;
    pvh += gstep;
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->n_hist; j++ )
    GEN->guide_table[j] = n_pv - 1;
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_hist_make_urntable() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_hist_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_HIST_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = HIST (HISTogram of empirical distribution)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cemp_debug( gen->distr, gen->genid, (gen->debug & HIST_DEBUG_PRINTHIST));

  fprintf(LOG,"%s: sampling routine = _unur_hist_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_hist_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_hist_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   functions = DATA  [histogram of size=%d]\n", DISTR.n_hist);
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
  /* method */
  _unur_string_append(info,"method: HIST (HISTogramm of empirical distribution)\n");
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

} /* end of _unur_hist_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
