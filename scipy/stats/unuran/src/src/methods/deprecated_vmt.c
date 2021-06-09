/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      deprecated_vmt.h                                             *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    generated random vector with independent components with     *
 *              given marginal distribution and use linear transformation    *
 *              of vector.                                                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      multivariate distribution with given mean vector and                 *
 *      covariance matrix.                                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  THIS METHOD AND THE CORRESPONDING ROUTINES SHOULD NOT BE USED ANY MORE!  *
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
 *   [1] Devroye, L. (1986): Non-Uniform Random Variate Generation, New-York *
 *                                                                           *
 *   [2] Hoermann, W., J. Leydold, and G. Derflinger (2004):                 *
 *       Automatic Nonuniform Random Variate Generation, Springer, Berlin.   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  VMT generates random vectors for distributions with given mean           *
 *  vector mu and covariance matrix Sigma. It produces random vectors        *
 *  of the form X = L Y + mu, where L is the Cholesky factor of Sigma,       *
 *  i.e. L L^t = Sigma, and Y has independent components of the same         *
 *  distribution with mean 0 and standard deviation 1.                       *
 *                                                                           *
 *  See [2], Sect.11.1.6, Alg.11.3.                                          *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "auto.h"
#include "deprecated_vmt.h"
#include "deprecated_vmt_struct.h"

/*---------------------------------------------------------------------------*/
#ifdef USE_DEPRECATED_CODE
/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define VMT_DEBUG_REINIT     0x00000010u   /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "VMT"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vmt_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_vmt_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vmt_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_vmt_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_vmt_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_vmt_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_vmt_make_marginal_gen( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* make generators for marginal distributions                                */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_vmt_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_vmt_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_vmt_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

/*---------------------------------------------------------------------------*/

#define _unur_vmt_getSAMPLE(gen) ( _unur_vmt_sample_cvec )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_vmt_new( const struct unur_distr *distr )
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
  if (distr->type != UNUR_DISTR_CVEC) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (!(distr->set & UNUR_DISTR_SET_MEAN)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"mean"); return NULL; }
  if (!(distr->set & UNUR_DISTR_SET_COVAR)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"covariance matrix");
    return NULL; }
  if (!(distr->set & UNUR_DISTR_SET_STDMARGINAL)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"standardized marginals");
    return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_vmt_par) );
  COOKIE_SET(par,CK_VMT_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_VMT ;     /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_vmt_init;

  return par;

} /* end of unur_vmt_new() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_vmt_init( struct unur_par *par )
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
  _unur_check_NULL( GENTYPE,par,NULL );

  /* check input */
  if ( par->method != UNUR_METH_VMT ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_VMT_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_vmt_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* initialize generators for marginal distribution */
  if (_unur_vmt_make_marginal_gen(gen) != UNUR_SUCCESS) {
    _unur_vmt_free(gen); return NULL;
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_vmt_debug_init(gen);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_vmt_init() */

/*---------------------------------------------------------------------------*/

int
_unur_vmt_reinit( struct unur_gen *gen )
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
  /* (re)set sampling routine */
  SAMPLE = _unur_vmt_getSAMPLE(gen);

  /* free list of (old) marginal generators */
  if (gen->gen_aux_list)
    _unur_gen_list_free( gen->gen_aux_list, gen->n_gen_aux_list );

  /* initialize generators for marginal distribution */
  return _unur_vmt_make_marginal_gen(gen);
} /* end of _unur_vmt_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vmt_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_VMT_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_vmt_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_VMT_GEN);

  /* dimension of distribution */
  GEN->dim = gen->distr->dim; 

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_vmt_getSAMPLE(gen);
  gen->destroy = _unur_vmt_free;
  gen->clone = _unur_vmt_clone;
  gen->reinit = _unur_vmt_reinit;

  /* cholesky factor of covariance matrix */
  GEN->cholesky =  DISTR.cholesky;  

  /* initialize pointer */
  GEN->marginalgen_list = NULL;

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_vmt_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_vmt_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_vmt_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_VMT_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* cholesky factor of covariance matrix */
  CLONE->cholesky =  clone->distr->data.cvec.cholesky;

  /* marginal generators are (also) stored as auxiliary generator */
  /* which has already been cloned by generic_clone.              */
  CLONE->marginalgen_list = clone->gen_aux_list;

  return clone;

#undef CLONE
} /* end of _unur_vmt_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_vmt_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_VMT ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_VMT_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  _unur_generic_free(gen);

} /* end of _unur_vmt_free() */

/*****************************************************************************/

int
_unur_vmt_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) (a*GEN->dim+b)
  int j,k;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_VMT_GEN,UNUR_ERR_COOKIE);

  while (1) {

    /* generate random vector with independent components */
    for (j=0; j<GEN->dim; j++)
      vec[j] = unur_sample_cont(GEN->marginalgen_list[j]);
    
    /* 
       transform to desired covariance structure: 
       X = L.Y + mu 
       where
       L  ... cholesky factor of the covariance matrix
       Y  ... vector with indenpent components (generated above)
       mu ... mean vector
       (notice that L is a lower triangular matrix)
    */
    for (k=GEN->dim-1; k>=0; k--) {
      vec[k] *= GEN->cholesky[idx(k,k)];
      for (j=k-1; j>=0; j--)
	vec[k] += vec[j] * GEN->cholesky[idx(k,j)];
      vec[k] += DISTR.mean[k];
    }

    if ( !(gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED) || 
	 _unur_distr_cvec_is_indomain( vec, gen->distr) )
      return UNUR_SUCCESS;
  }

#undef idx
} /* end of _unur_vmt_sample_cvec() */

/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_vmt_make_marginal_gen( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* make generators for marginal distributions                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  if (_unur_distr_cvec_marginals_are_equal(DISTR.stdmarginals, GEN->dim)) {
    /* we can use the same generator object for all marginal distributions */
    struct unur_gen *marginalgen = unur_init( unur_auto_new( DISTR.stdmarginals[0] ) );
    if (marginalgen)
      gen->gen_aux_list = _unur_gen_list_set(marginalgen,GEN->dim);
      gen->n_gen_aux_list = GEN->dim;
  }

  else {
    int i,j;
    int failed = FALSE;
    struct unur_gen **marginalgens = _unur_xmalloc( GEN->dim * sizeof(struct unur_gen*) );
    for (i=0; i<GEN->dim; i++) {
      marginalgens[i] = unur_init( unur_auto_new( DISTR.stdmarginals[i] ) );
      if (marginalgens[i]==NULL) {
	failed=TRUE; break; 
      }
    }
    if (failed) {
      for (j=0; j<i; j++) _unur_free(marginalgens[j]);
      free (marginalgens);
    }
    else {
      gen->gen_aux_list = marginalgens;
      gen->n_gen_aux_list = GEN->dim;
    }
  }

  /* the marginal generator is an auxiliary generator for method VMT, of course */
  GEN->marginalgen_list = gen->gen_aux_list;
  
  /* verify initialization of marginal generators */
  if (GEN->marginalgen_list == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GENERIC,"init of marginal generators failed");
    return UNUR_FAILURE;
  }

  return UNUR_SUCCESS;
} /* end of _unur_vmt_make_marginal_gen() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_vmt_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  int i;
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_VMT_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = VMT (Vector Matrix Transformation)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

  fprintf(LOG,"%s: generators for standardized marginal distributions = \n",gen->genid);
  fprintf(LOG,"%s:\t",gen->genid);
  for (i=0; i<GEN->dim; i++)
    fprintf(LOG,"[%s] ", GEN->marginalgen_list[i]->genid);
  fprintf(LOG,"\n%s:\n",gen->genid);

  fprintf(LOG,"%s: sampling routine = _unur_vmt_sample_cvec()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: INIT completed **********************\n",gen->genid);

} /* end of _unur_vmt_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/
#endif   /* USE_DEPRECATED_CODE */
/*---------------------------------------------------------------------------*/
