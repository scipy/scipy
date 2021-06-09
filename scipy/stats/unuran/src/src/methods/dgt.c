/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      dgt.c                                                        *
 *                                                                           *
 *   TYPE:      discrete univariate random variate                           *
 *   METHOD:    guide table (indexed search)                                 *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given N discrete events with different probabilities P[k]            *
 *      produce a value k consistent with its probability.                   *
 *                                                                           *
 *   REQUIRED:  pointer to probability vector                                *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold             *
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
 *   [1] Chen, H. C. and Asau, Y. (1974): On generating random variates      *
 *       from an empirical distribution, AIIE Trans. 6, pp. 163-166          *
 *                                                                           *
 *   SUMMARY:                                                                *
 *   [2] Devroye, L. (1986): Non-Uniform Random Variate Generation, New-York *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   This method a varient of the inversion method, i.e.                     *
 *                                                                           *
 *   (1) Generate a random number U ~ U(0,1).                                *
 *   (2) Find largest integer I such that F(I) = P(X<=I) <= U.               *
 *                                                                           *
 *   Step (2) is the crucial step. Using sequential search requires O(N)     *
 *   comparisons. Indexed search however uses a guide table to jump to some  *
 *   I' <= I near I. In a preprossing step (0,1) is partitioned into N       *
 *   equal intervals and P(X<=I) <= (k-1)/N is solved for all k=1,...,N,     *
 *   and the solutions are stored in a table (guide table). Setup this       *
 *   table can be done in O(N). [2] has shown that the expected number of    *
 *   of comparisons is at most 2.                                            *
 *                                                                           *
 *   In the current implementation we use a variant that allows a different  *
 *   size for the guide table. Bigger guide tables reduce the expected       *
 *   number of comparisions but needs more memory and need more setup time.  *
 *   On the other hand, the guide table can be made arbitrarily small to     *
 *   save memory and setup time. Indeed, for size = 1 we have sequential     *
 *   search again that requires no preprocessing.                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   VARIANTS:                                                               *
 *                                                                           *
 *   We have three variants for the setup procedure:                         *
 *                                                                           *
 *   1 ... compute (k-1)/N for each k=1,...,N                                *
 *   2 ... compute (k-1)/N by summing up 1/N                                 *
 *   0 ... use 1 for large N and 2 for small N (default)                     *
 *                                                                           *
 *   variant 2 is faster but is more sensitive to roundoff errors when       *
 *   N is large.                                                             *
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
#include "dgt.h"
#include "dgt_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define DGT_VARFLAG_DIV     0x01u     /* compute guide table by division n/k */
#define DGT_VARFLAG_ADD     0x02u     /* compute guide table by adding       */

#define DGT_VAR_THRESHOLD   1000      /* above this value: use variant 1, else 2 */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define DGT_DEBUG_REINIT       0x00000010u  /* print parameters after reinit */
#define DGT_DEBUG_PRINTVECTOR  0x00000100u
#define DGT_DEBUG_TABLE        0x00000200u

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define DGT_SET_GUIDEFACTOR    0x010u
#define DGT_SET_VARIANT        0x020u

/*---------------------------------------------------------------------------*/

#define GENTYPE "DGT"         /* type of generator                           */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dgt_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dgt_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dgt_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_dgt_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_dgt_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_dgt_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_dgt_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_dgt_create_tables( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create (allocate) tables                                                  */
/*---------------------------------------------------------------------------*/

static int _unur_dgt_make_guidetable( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* create table for indexed search                                           */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_dgt_debug_init( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_dgt_debug_table( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print data for guide table.                                               */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_dgt_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.discr      /* data for distribution object      */

#define PAR       ((struct unur_dgt_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_dgt_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.discr /* data for distribution in generator object */

#define SAMPLE    gen->sample.discr     /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/

#define _unur_dgt_getSAMPLE(gen)  (_unur_dgt_sample)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_dgt_new( const struct unur_distr *distr )
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

  if (DISTR_IN.pv == NULL) {
    /* There is no PV try to compute it.                         */
    if ( DISTR_IN.pmf
	 && ( (((unsigned)DISTR_IN.domain[1] - (unsigned)DISTR_IN.domain[0]) < UNUR_MAX_AUTO_PV)
	      || ( (distr->set & UNUR_DISTR_SET_PMFSUM) && DISTR_IN.domain[0] > INT_MIN ) ) ) {
      /* However this requires a PMF and either a bounded domain   */
      /* or the sum over the PMF.                                  */
      _unur_warning(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV. Try to compute it.");
    }
    else {
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); return NULL;
    }
  }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_dgt_par) );
  COOKIE_SET(par,CK_DGT_PAR);

  /* copy input */
  par->distr       = distr;          /* pointer to distribution object       */

  /* set default values */
  PAR->guide_factor = 1.;            /* use same size for guide table        */

  par->method      = UNUR_METH_DGT;  /* method                               */
  par->variant     = 0u;             /* default variant                      */
  par->set         = 0u;             /* inidicate default parameters         */    
  par->urng        = unur_get_default_urng(); /* use default urng            */
  par->urng_aux    = NULL;                    /* no auxilliary URNG required */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_dgt_init;

  return par;

} /* end of unur_dgt_new() */

/*---------------------------------------------------------------------------*/

int
unur_dgt_set_variant( struct unur_par *par, unsigned variant )
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
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DGT );

  /* check new parameter for generator */
  if (variant != DGT_VARFLAG_ADD && variant != DGT_VARFLAG_DIV) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_VARIANT,"");
    return UNUR_ERR_PAR_VARIANT;
  }

  /* changelog */
  par->set |= DGT_SET_VARIANT;

  par->variant = variant;

  return UNUR_SUCCESS;
} /* end of unur_dgt_set_variant() */

/*---------------------------------------------------------------------------*/

int
unur_dgt_set_guidefactor( struct unur_par *par, double factor )
     /*----------------------------------------------------------------------*/
     /* set factor for relative size of guide table                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   factor ... relative size of table                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, DGT );

  /* check new parameter for generator */
  if (factor < 0) {
    _unur_warning(GENTYPE,UNUR_ERR_PAR_SET,"relative table size < 0");
    return UNUR_ERR_PAR_SET;
  }

  /* store date */
  PAR->guide_factor = factor;

  /* changelog */
  par->set |= DGT_SET_GUIDEFACTOR;

  return UNUR_SUCCESS;

} /* end of unur_dgt_set_guidefactor() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_dgt_init( struct unur_par *par )
     /*----------------------------------------------------------------------*/
     /* initialize new generator                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par ... pointer to parameter for building generator object         */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to generator object                                        */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_gen *gen;         /* pointer to generator object */
  
  /* check arguments */
  CHECK_NULL(par,NULL);
  
  /* check input */
  if ( par->method != UNUR_METH_DGT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_DGT_PAR,NULL);
  
  /* create a new empty generator object */
  gen = _unur_dgt_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check parameters */
  if ( _unur_dgt_check_par(gen) != UNUR_SUCCESS ) {
    _unur_dgt_free(gen); return NULL;
  }

  /* compute guide table */
  if ( (_unur_dgt_create_tables(gen) != UNUR_SUCCESS) ||
       (_unur_dgt_make_guidetable(gen) != UNUR_SUCCESS) ) {
    _unur_dgt_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_dgt_debug_init(gen);
#endif

  return gen;
} /* end of _unur_dgt_init() */

/*---------------------------------------------------------------------------*/

int
_unur_dgt_reinit( struct unur_gen *gen )
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
  if ( (rcode = _unur_dgt_check_par(gen)) != UNUR_SUCCESS)
    return rcode;

  /* compute table */
  if ( ((rcode = _unur_dgt_create_tables(gen)) != UNUR_SUCCESS) ||
       ((rcode = _unur_dgt_make_guidetable(gen)) != UNUR_SUCCESS) ) {
    return rcode;
  }

  /* (re)set sampling routine */
  SAMPLE = _unur_dgt_getSAMPLE(gen);

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & DGT_DEBUG_REINIT) _unur_dgt_debug_init(gen);
#endif

  return UNUR_SUCCESS;
} /* end of _unur_dgt_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dgt_create( struct unur_par *par )
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
  struct unur_gen *gen;       /* pointer to generator object */

  /* check arguments */
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_DGT_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_dgt_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_DGT_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_dgt_getSAMPLE(gen);
  gen->destroy = _unur_dgt_free;
  gen->clone = _unur_dgt_clone;
  gen->reinit = _unur_dgt_reinit;

  /* copy some parameters into generator object */
  GEN->guide_factor = PAR->guide_factor;

  /* set all pointers to NULL */
  GEN->cumpv = NULL;
  GEN->guide_table = NULL;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_dgt_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_dgt_create() */

/*---------------------------------------------------------------------------*/

int
_unur_dgt_check_par( struct unur_gen *gen )
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
  /* we need a PV */
  if (DISTR.pv == NULL) {
    /* try to compute PV */
    if (unur_distr_discr_make_pv( gen->distr ) <= 0) {
      /* not successful */
      _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"PV"); 
      return UNUR_ERR_DISTR_REQUIRED;
    }
  }

  /* default variant? */
  if (gen->variant == 0)   /* default variant */
    gen->variant = (DISTR.n_pv > DGT_VAR_THRESHOLD) 
      ? DGT_VARFLAG_DIV : DGT_VARFLAG_ADD;

  return UNUR_SUCCESS;
} /* end of _unur_dgt_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_dgt_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_dgt_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_DGT_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* copy data for distribution */
  CLONE->cumpv = _unur_xmalloc( DISTR.n_pv * sizeof(double) );
  memcpy( CLONE->cumpv, GEN->cumpv, DISTR.n_pv * sizeof(double) );
  CLONE->guide_table = _unur_xmalloc( GEN->guide_size * sizeof(int) );
  memcpy( CLONE->guide_table, GEN->guide_table, GEN->guide_size * sizeof(int) );

  return clone;

#undef CLONE
} /* end of _unur_dgt_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_dgt_free( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* deallocate generator object                                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{ 

  /* check arguments */
  if (!gen) /* nothing to do */
    return;

  /* check input */
  if ( gen->method != UNUR_METH_DGT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_DGT_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free two auxiliary tables */
  if (GEN->guide_table) free(GEN->guide_table);
  if (GEN->cumpv)       free(GEN->cumpv);

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_dgt_free() */

/*****************************************************************************/

int
_unur_dgt_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   integer (sample from random variate)                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{ 
  int j;
  double u;

  /* check arguments */
  CHECK_NULL(gen,INT_MAX);  COOKIE_CHECK(gen,CK_DGT_GEN,INT_MAX);

  /* sample from U(0,1) */
  u = _unur_call_urng(gen->urng);

  /* look up in guide table ... */
  j = GEN->guide_table[(int)(u * GEN->guide_size)];
  /* ... and search */
  u *= GEN->sum;
  while (GEN->cumpv[j] < u) j++;

  return (j + DISTR.domain[0]);

} /* end of _unur_dgt_sample() */


/*---------------------------------------------------------------------------*/

int
unur_dgt_eval_invcdf_recycle( const struct unur_gen *gen, double u, double *recycle )
     /*----------------------------------------------------------------------*/
     /* evaluate inverse CDF at u.                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   u       ... argument for inverse CDF (0<=u<=1)                     */
     /*   recycle ... if not NULL then store recycled 'u'                    */
     /*                                                                      */
     /* return:                                                              */
     /*   integer (inverse CDF)                                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{
  int j;

  /* set default */
  if (recycle) *recycle = 0.;

  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INT_MAX );
  if ( gen->method != UNUR_METH_DGT ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INT_MAX;
  }
  COOKIE_CHECK(gen,CK_DGT_GEN,INT_MAX);

  /* check range of u */
  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return INT_MAX;  /* u == NaN */
  }

  /* look up in guide table ... */
  j = GEN->guide_table[(int)(u * GEN->guide_size)];
  /* ... and search */
  u *= GEN->sum;
  while (GEN->cumpv[j] < u) j++;

  if (recycle) {
    *recycle = 1. - (GEN->cumpv[j] - u) / DISTR.pv[j];
  }

  j+=DISTR.domain[0];

  /* validate range */
  if (j<DISTR.domain[0]) j = DISTR.domain[0];
  if (j>DISTR.domain[1]) j = DISTR.domain[1];

  return j;

} /* end of unur_dgt_eval_invcdf_recycle() */

/*---------------------------------------------------------------------------*/

int
unur_dgt_eval_invcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate inverse CDF at u.                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   u       ... argument for inverse CDF (0<=u<=1)                     */
     /*                                                                      */
     /* return:                                                              */
     /*   integer (inverse CDF)                                              */
     /*                                                                      */
     /* error:                                                               */
     /*   return INT_MAX                                                     */
     /*----------------------------------------------------------------------*/
{
  return unur_dgt_eval_invcdf_recycle(gen,u,NULL);
} /* end of unur_dgt_eval_invcdf() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_dgt_create_tables( struct unur_gen *gen )
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
  /* size of guide table */
  GEN->guide_size = (int)( DISTR.n_pv * GEN->guide_factor);
  if (GEN->guide_size <= 0)
    /* do not use a guide table whenever params->guide_factor is 0 or less */
    GEN->guide_size = 1;

  /* allocation for cummulated probabilities */
  GEN->cumpv = _unur_xrealloc( GEN->cumpv, DISTR.n_pv * sizeof(double) );

  /* allocate memory for the guide table */
  GEN->guide_table = _unur_xrealloc( GEN->guide_table, GEN->guide_size * sizeof(int) );

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_dgt_create_tables() */

/*---------------------------------------------------------------------------*/

int
_unur_dgt_make_guidetable( struct unur_gen *gen )
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
  pv = DISTR.pv;
  n_pv = DISTR.n_pv;

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
  
  if (gen->variant == DGT_VARFLAG_DIV) {
    GEN->guide_table[0] = 0;
    for( j=1, i=0; j<GEN->guide_size ;j++ ) {
      while( GEN->cumpv[i]/GEN->sum < ((double)j)/GEN->guide_size ) 
	i++;
      if (i >= n_pv) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN->guide_table[j]=i;
    }
  }

  else { /* gen->variant == DGT_VARFLAG_ADD */
    gstep = GEN->sum / GEN->guide_size;
    pvh = 0.;
    for( j=0, i=0; j<GEN->guide_size ;j++ ) {
      while (GEN->cumpv[i] < pvh) 
	i++;
      if (i >= n_pv) {
	_unur_warning(gen->genid,UNUR_ERR_ROUNDOFF,"guide table");
	break;
      }
      GEN->guide_table[j] = i;
      pvh += gstep;
    }
  }

  /* if there has been an round off error, we have to complete the guide table */
  for( ; j<GEN->guide_size ;j++ )
    GEN->guide_table[j] = n_pv - 1;
  /* o.k. */
  return UNUR_SUCCESS;

} /* end of _unur_dgt_make_guidetable() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_dgt_debug_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DGT_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = discrete univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = indexed search (guide table)\n",gen->genid);

  fprintf(LOG,"%s: variant = %u ",gen->genid,gen->variant);
  _unur_print_if_default(gen,DGT_SET_VARIANT);
  fprintf(LOG,"\n%s:\n",gen->genid);

  _unur_distr_discr_debug( gen->distr,gen->genid,(gen->debug & DGT_DEBUG_PRINTVECTOR));

  fprintf(LOG,"%s: sampling routine = _unur_dgt_sample()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s: length of probability vector = %d\n",gen->genid,DISTR.n_pv);
  fprintf(LOG,"%s: length of guide table = %d   (rel. = %g%%",
	  gen->genid,GEN->guide_size,100.*GEN->guide_factor);
  _unur_print_if_default(gen,DGT_SET_GUIDEFACTOR);
  if (GEN->guide_size == 1) 
    fprintf(LOG,") \t (-->sequential search");
  fprintf(LOG,")\n%s:\n",gen->genid);

  fprintf(LOG,"%s: sum over PMF (as computed) = %#-20.16g\n",gen->genid,GEN->sum);

  if (gen->debug & DGT_DEBUG_TABLE)
    _unur_dgt_debug_table(gen);

} /* end of _unur_dgt_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_dgt_debug_table( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write guide table into LOG file                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{   
  FILE *LOG;
  int i,j,m;
  int n_asts;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_DGT_GEN,RETURN_VOID);

  LOG = unur_get_stream();
  
  fprintf(LOG,"%s: guide table:\n", gen->genid); 
  fprintf(LOG,"%s:\n", gen->genid);
  n_asts = 0;
  for (i=0; i<GEN->guide_size; i++){
    fprintf(LOG,"%s: [%5d] -> %5d ", gen->genid, i, GEN->guide_table[i]);
    /* print row of asterisks */
    if (i == GEN->guide_size-1)
      j = GEN->guide_size - GEN->guide_table[i];
    else
      j = GEN->guide_table[i+1] - GEN->guide_table[i] + 1;
    for (m=0; m<j && m<10; m++ ) {
      fprintf(LOG," *");
      ++n_asts;
    }
    /* too many asterisks print */
    if (m<j){
      n_asts += j-m;
      fprintf(LOG," ... %d", j);
    }
    fprintf(LOG,"\n");
  }

  /* print expected number of comparisons */
  fprintf(LOG,"%s:\n", gen->genid);
  fprintf(LOG,"%s: expected number of comparisons = %g\n",gen->genid,
          ((double)n_asts)/GEN->guide_size);
  fprintf(LOG,"%s:\n", gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);

} /*  end of _unur_dgt_debug_table() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_dgt_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   functions = PV  [length=%d%s]\n",
		      DISTR.domain[1]-DISTR.domain[0]+1,
		      (DISTR.pmf==NULL) ? "" : ", created from PMF");
  _unur_string_append(info,"   domain    = (%d, %d)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: DGT (Guide Table)\n");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   E [#look-ups] = %g\n", 1+1./GEN->guide_factor);
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   guidefactor = %g  %s\n", GEN->guide_factor,
			(gen->set & DGT_SET_GUIDEFACTOR) ? "" : "[default]");
    if (gen->set & DGT_SET_VARIANT)
      _unur_string_append(info,"   variant = %d\n", gen->variant);
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_dgt_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
