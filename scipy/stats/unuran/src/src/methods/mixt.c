/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mixt.c                                                        *
 *                                                                           *
 *   TYPE:      (continuous) univariate random variate                       *
 *   METHOD:    mixture of distributions (meta method)                       *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given an array of components and array of probabilities              *
 *      produce a value X consistent with the mixture distribution           *
 *                                                                           *
 *   REQUIRED:                                                               *
 *      pointers to generators of components                                 *
 *      array of probabilities                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2010 Wolfgang Hoermann and Josef Leydold                  *
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
#include <distr/cont.h>
#include <distr/discr.h>
#include <urng/urng.h>
#include <utils/unur_fp_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "dgt.h"
#include "dgt_struct.h"
#include "mixt.h"
#include "mixt_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants:                                                                 */

#define MIXT_VARFLAG_INVERSION   0x004u    /* use inversion method (if possible) */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define MIXT_SET_USEINVERSION     0x001u    /* use inverion method            */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MIXT"          /* type of generator                         */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static int _unur_mixt_check_par( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Check parameters of given distribution and method                         */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_mixt_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static double _unur_mixt_sample( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static double _unur_mixt_sample_inv( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* sample from generator by inversion                                        */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mixt_indexgen( const double *prob, int n_prob );
/*---------------------------------------------------------------------------*/
/* create generator for index.                                               */
/*---------------------------------------------------------------------------*/

static int _unur_mixt_get_boundary( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* compute boundary of mixture and check for overlapping domains.            */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mixt_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_mixt_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_mixt_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mixt_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cont /* data for distribution in generator object */

#define SAMPLE    gen->sample.cont      /* pointer to sampling routine       */

#define BD_LEFT   domain[0]             /* left boundary of domain of distribution */
#define BD_RIGHT  domain[1]             /* right boundary of domain of distribution */

/* shortcuts for auxiliary generators */
#define INDEX     gen_aux

#define PROB      gen_aux->distr->data.discr.pv
#define COMP      gen_aux_list
#define N_COMP    n_gen_aux_list

/*---------------------------------------------------------------------------*/

#define _unur_mixt_getSAMPLE(gen) \
   ( ((gen)->variant & MIXT_VARFLAG_INVERSION) \
     ? _unur_mixt_sample_inv : _unur_mixt_sample )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mixt_new( int n, const double *prob, struct unur_gen **comp )
     /*----------------------------------------------------------------------*/
     /* get default parameters                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   n    ... number of components                                      */
     /*   prob ... probabilities for components                              */
     /*   comp ... array of pointers to components (generators)              */
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
  _unur_check_NULL( GENTYPE, prob, NULL );
  _unur_check_NULL( GENTYPE, comp, NULL );
  if (n<1) { _unur_error(GENTYPE,UNUR_ERR_DISTR_DOMAIN,"n < 1"); return NULL; }
  /* checking type of generator objects in 'comp' is delayed to init */

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mixt_par) );
  COOKIE_SET(par,CK_MIXT_PAR);

  /* copy input */
  par->distr    = NULL;      /* pointer to distribution object               */

  /* copy data */
  PAR->n_comp   = n;         /* number of components                         */
  PAR->prob     = prob;      /* probabilities for components                 */
  PAR->comp     = comp;      /* array of pointers to components (generators) */

  /* set default values */
  par->method   = UNUR_METH_MIXT;   /* method and default variant            */
  par->variant  = 0u;               /* default variant                       */
  par->set      = 0u;               /* inidicate default parameters          */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_mixt_init;

  return par;

} /* end of unur_mixt_new() */

/*****************************************************************************/

int
unur_mixt_set_useinversion( struct unur_par *par, int useinversion )
     /*----------------------------------------------------------------------*/
     /* set flag for using inversion method (default: off)                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   par    ... pointer to parameter for building generator object      */
     /*   useinversion ... !0 = try inversion method,  0 = no inversion      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   no squeeze is the default                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MIXT );

  /* we use a bit in 'variant' */
  par->variant = (useinversion)
    ? (par->variant | MIXT_VARFLAG_INVERSION)
    : (par->variant & (~MIXT_VARFLAG_INVERSION));

  /* changelog */
  par->set |= MIXT_SET_USEINVERSION;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_mixt_set_useinversion() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_mixt_init( struct unur_par *par )
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
  int i;

  /* check arguments */
  CHECK_NULL(par,NULL);

  /* check input */
  if ( par->method != UNUR_METH_MIXT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MIXT_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mixt_create(par);
  if (!gen) { _unur_par_free(par); return NULL; }

  /* probabilities */
  gen->INDEX = _unur_mixt_indexgen(PAR->prob,PAR->n_comp);

  /* components */
  gen->N_COMP = PAR->n_comp;    /* number of components                         */
  gen->COMP = _unur_xmalloc( gen->N_COMP * sizeof(struct unur_gen *));
  for (i=0; i<gen->N_COMP; i++)
    gen->COMP[i] = unur_gen_clone(PAR->comp[i]);

  /* free parameters */
  _unur_par_free(par);

  /* check parameters */
  if (_unur_mixt_check_par(gen) != UNUR_SUCCESS) {
    _unur_mixt_free(gen); return NULL;
  }

  /* compute boundary of mixture and check for overlapping domains */
  if ( _unur_mixt_get_boundary(gen) != UNUR_SUCCESS ) {
    _unur_mixt_free(gen); return NULL;
  }

  /* set name of distribution */
  unur_distr_set_name(gen->distr, "(mixture)");

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file */
    if (gen->debug) _unur_mixt_debug_init(gen);
#endif

  return gen;
} /* end of _unur_mixt_init() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mixt_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MIXT_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mixt_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MIXT_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* object 'par' does not contain a distribution object */
  /* so we create one.                                   */
  gen->distr = unur_distr_cont_new();

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mixt_getSAMPLE(gen);
  gen->destroy = _unur_mixt_free;
  gen->clone = _unur_mixt_clone;
  gen->reinit = NULL;    /* reinit not implemented ! */

  /* copy some parameters into generator object */
  GEN->is_inversion = (gen->variant & MIXT_VARFLAG_INVERSION) ? TRUE : FALSE;

  /* initialize parameters: none */

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_mixt_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_mixt_create() */

/*---------------------------------------------------------------------------*/

int
_unur_mixt_check_par( struct unur_gen *gen )
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
  int i;
  int type;

  /* check probabilities */
  if (gen->INDEX == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"invalid probabilities");
    return UNUR_ERR_GEN_DATA;
  }

  /* check generator objects */
  for (i=0; i<gen->N_COMP; i++) {

    /* generators must not be NULL */
    if (gen->COMP[i] == NULL) {
      _unur_error(gen->genid,UNUR_ERR_NULL,"component is NULL");
      return UNUR_ERR_NULL;
    }

    /* all generators must sample from univariate distributions */
    type = gen->COMP[i]->method & UNUR_MASK_TYPE;
    if ( type != UNUR_METH_DISCR && 
	 type != UNUR_METH_CONT  &&
	 type != UNUR_METH_CEMP  ) {
      _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"component not univariate");
      return UNUR_ERR_GEN_INVALID;
    }

    /* we only can use inversion method if all generators use inversion method */
    if (GEN->is_inversion && (! _unur_gen_is_inversion (gen->COMP[i]))) {
      _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"component does not implement inversion");
      return UNUR_ERR_GEN_INVALID;
    }
  }

  return UNUR_SUCCESS;
} /* end of _unur_mixt_check_par() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mixt_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_mixt_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MIXT_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  return clone;

#undef CLONE
} /* end of _unur_mixt_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_mixt_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_MIXT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MIXT_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  _unur_generic_free(gen);

} /* end of _unur_mixt_free() */

/*****************************************************************************/

double
_unur_mixt_sample( struct unur_gen *gen )
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
  struct unur_gen *comp;
  int J;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_MIXT_GEN,INFINITY);

  /* sample index */
  J = unur_sample_discr(gen->INDEX);
  
  /* get component */
  comp = gen->COMP[J];

  /* sample from selected component */
  switch(comp->method & UNUR_MASK_TYPE) {
  case UNUR_METH_DISCR:
    return ((double) comp->sample.discr(comp));
  case UNUR_METH_CONT:
  case UNUR_METH_CEMP:
  default:
    return (comp->sample.cont(comp));
  }

} /* end of _unur_mixt_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_mixt_sample_inv( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator by inversion                                   */
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
  double U, recycle;
  int J;

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_MIXT_GEN,INFINITY);

  /* sample index */
  U = _unur_call_urng(gen->urng);
  J =unur_dgt_eval_invcdf_recycle( gen->INDEX, U, &recycle );

  /* the resolution of recycle is less than that of U. */
  /* the values 0. and 1. may be result in INFINITY.   */
  /* thus we make a small perturbation in this case.   */
  if (_unur_iszero(recycle)) recycle = DBL_MIN;
  if (_unur_isone(recycle))  recycle = 1. - DBL_EPSILON;

  /* sample from component */
  return unur_quantile(gen->COMP[J], recycle);

} /* end of _unur_mixt_sample_inv() */

/*---------------------------------------------------------------------------*/

double
unur_mixt_eval_invcdf( const struct unur_gen *gen, double u )
     /*----------------------------------------------------------------------*/
     /* evaluate inverse CDF at u.                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   u   ... argument for inverse CDF                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   double (inverse CDF)                                               */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*----------------------------------------------------------------------*/
{
  double recycle;
  int J;
  
  /* check arguments */
  _unur_check_NULL( GENTYPE, gen, INFINITY );
  if ( ! (gen->method == UNUR_METH_MIXT && GEN->is_inversion) ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return INFINITY;
  }
  COOKIE_CHECK(gen,CK_MIXT_GEN,INFINITY);

  if ( ! (u>0. && u<1.)) {
    if ( ! (u>=0. && u<=1.)) {
      _unur_warning(gen->genid,UNUR_ERR_DOMAIN,"U not in [0,1]");
    }
    if (u<=0.) return DISTR.domain[0];
    if (u>=1.) return DISTR.domain[1];
    return u;  /* = NaN */
  }

  /* get index */
  J =unur_dgt_eval_invcdf_recycle( gen->INDEX, u, &recycle );

  /* the resolution of recycle is less than that of U. */
  /* the values 0. and 1. may be result in INFINITY.   */
  /* thus we make a small perturbation in this case.   */
  if (_unur_iszero(recycle)) recycle = DBL_MIN;
  if (_unur_isone(recycle))  recycle = 1. - DBL_EPSILON;

  /* get from component */
  return unur_quantile(gen->COMP[J], recycle);

} /* end of unur_mixt_eval_invcdf() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

struct unur_gen *
_unur_mixt_indexgen( const double *prob, int n_prob )
/*---------------------------------------------------------------------------*/
/* create generator for index using method DGT.                              */
/*                                                                           */
/* parameters:                                                               */
/*   prob   ... probability vector                                           */
/*   n_prob ... length of probability vector                                 */
/*                                                                           */
/* return:                                                                   */
/*   pointer to generator object                                             */
/*                                                                           */
/* error:                                                                    */
/*   return NULL                                                             */
/*---------------------------------------------------------------------------*/
{
  struct unur_distr *distr;
  struct unur_par *par;
  struct unur_gen *igen;

  /* create generator */
  distr = unur_distr_discr_new();
  unur_distr_discr_set_pv(distr, prob, n_prob);
  par = unur_dgt_new(distr);
  igen = unur_init(par);

  /* clear working space */
  unur_distr_free(distr);
  
  return igen;

} /* end of _unur_mixt_indexgen() */

/*---------------------------------------------------------------------------*/

int
_unur_mixt_get_boundary( struct unur_gen *gen )
/*---------------------------------------------------------------------------*/
/* compute boundary of mixture and check for overlapping domains.            */
/*                                                                           */
/* parameters:                                                               */
/*   gen ... pointer to generator object                                     */
/*                                                                           */
/* return:                                                                   */
/*   UNUR_SUCCESS ... when everything is o.k.                                */
/*   UNUR_NULL    ... when a component is NULL (which should not happen)     */
/*   UNUR_ERR_GEN_INVALID .. if domains of components overlap or not ordered */
/*---------------------------------------------------------------------------*/
{
  int i;
  int overlap = FALSE;
  double comp_left, comp_right;
  double bd_left, bd_right;
  struct unur_gen *comp;

  /* initialize boundary values for mixture */
  bd_left = INFINITY;
  bd_right = -INFINITY;

  for (i=0; i<gen->N_COMP; i++) {
    comp = gen->COMP[i];

    /* comp has already been check. but we want to be on the safe side. */
    CHECK_NULL(comp,UNUR_ERR_NULL);

    /* domain of component [i] */
    switch (comp->method & UNUR_MASK_TYPE) {
    case UNUR_METH_CONT:
      comp_left  = comp->distr->data.cont.BD_LEFT;
      comp_right = comp->distr->data.cont.BD_RIGHT;
      break;

    case UNUR_METH_DISCR:
      comp_left  = (double) (comp->distr->data.discr.BD_LEFT);
      comp_right = (double) (comp->distr->data.discr.BD_RIGHT);
      break;
      
    default:
      /* cannot estimate boundary */
      comp_left = -INFINITY;
      comp_right = INFINITY;
    }

    /* check for overlapping domains */
    if ( _unur_FP_less(comp_left,bd_right) )
      overlap = TRUE;

    /* update domain of mixture */
    bd_left = _unur_min(bd_left, comp_left);
    bd_right = _unur_max(bd_right, comp_right);
  }

  /* overlap or unordered domains? */
  if (GEN->is_inversion && overlap) {
    _unur_error(gen->genid,UNUR_ERR_GEN_INVALID,"domains of components overlap or are unsorted");
    return UNUR_ERR_GEN_INVALID;
  }

  /* store boundary */
  unur_distr_cont_set_domain(gen->distr, bd_left, bd_right);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_mixt_get_boundary() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_mixt_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  struct unur_gen *comp;
  int i;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MIXT_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = MIXT (MIXTure of distributions -- meta method)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cont_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: sampling routine = _unur_mixt_sample",gen->genid);
  if (GEN->is_inversion) fprintf(LOG,"_inv");
  fprintf(LOG,"()\n%s:\n",gen->genid);

  fprintf(LOG,"%s: use inversion = %s",gen->genid,
	  (GEN->is_inversion) ? "on" : "off");
  _unur_print_if_default(gen,MIXT_SET_USEINVERSION);
  fprintf(LOG,"\n%s:\n",gen->genid);

  /* probabilities */
  fprintf(LOG,"%s: probabilities (%d) = \n",gen->genid, gen->N_COMP);
  fprintf(LOG,"%s:   %g",gen->genid, (gen->PROB)[0]);
  for (i=1; i<gen->N_COMP; i++)
    fprintf(LOG,", %g", (gen->PROB)[i]);
  fprintf(LOG,"\n%s:\n",gen->genid);

  /* components */
  fprintf(LOG,"%s: components (%d):\n",gen->genid, gen->N_COMP);
  for (i=0; i<gen->N_COMP; i++) {
    comp = gen->COMP[i];
    fprintf(LOG,"%s:   [%d]: %s\n",gen->genid, i, comp->genid);
    fprintf(LOG,"%s:\t type = ",gen->genid); 
    switch (comp->distr->type) {
    case UNUR_DISTR_CONT:
    case UNUR_DISTR_CEMP:
      fprintf(LOG,"continuous\n");
      break;
    case UNUR_DISTR_DISCR:
      fprintf(LOG,"discrete\n");
      break;
    default:
      fprintf(LOG,"[unknown]\n");
    }
    fprintf(LOG,"%s:\t name = %s\n",gen->genid, comp->distr->name);
  }
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_mixt_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_mixt_info( struct unur_gen *gen, int help )
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
  struct unur_gen *comp;
  int i;
  double sum;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   domain    = (%g, %g)\n", DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   # components = %d\n", gen->N_COMP);

  if (help) {
    sum = ((struct unur_dgt_gen*)gen->INDEX->datap)->sum;
    _unur_string_append(info,"   probabilities = (%g", gen->PROB[0] / sum);
    for (i=1; i<gen->N_COMP; i++)
      _unur_string_append(info,", %g", gen->PROB[i] / sum);
    _unur_string_append(info,")\n");
    
    _unur_string_append(info,"   components = \n");
    for (i=0; i<gen->N_COMP; i++) {
      comp = gen->COMP[i];
      _unur_string_append(info,"\t[%d] %s - ",i, comp->genid);
      switch (comp->distr->type) {
      case UNUR_DISTR_CONT:
      case UNUR_DISTR_CEMP:
	_unur_string_append(info,"continuous");
	break;
      case UNUR_DISTR_DISCR:
	_unur_string_append(info,"discrete");
	break;
      default:
	_unur_string_append(info,"[unknown]");
      }
      _unur_string_append(info,": %s\n",comp->distr->name);
    }
  }
  _unur_string_append(info,"\n");

  /* method */
  _unur_string_append(info,"method: MIXT (MIXTure of distributions -- meta method)\n");
  _unur_string_append(info,"   select component = method DGT\n");
  _unur_string_append(info,"   inversion method = %s\n",
		      (GEN->is_inversion) ? "TRUE" : "FALSE");
  _unur_string_append(info,"\n");

  /* performance */
  _unur_string_append(info,"performance characteristics: depends on components\n");
  _unur_string_append(info,"\n");

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   useinversion = ");
    if (gen->variant & MIXT_VARFLAG_INVERSION)
      _unur_string_append(info,"on\n");
    else
      _unur_string_append(info,"off  [default]\n");
  }

} /* end of _unur_mixt_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
