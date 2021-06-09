/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mcorr.c                                                      *
 *                                                                           *
 *   TYPE:      random matrix                                                *
 *   METHOD:    Matrix -- COORelation matrix                                 *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      random correlation matrix.                                           *
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
 *   [1] Devroye, L. (1986): Non-Uniform Random Variate Generation,          *
 *       New-York, Sect.6.1, p.605.                                          *
 *                                                                           *
 *   [2] Marsaglia, G. and I. Olkin (1984):                                  *
 *       Generating Correlation Matrices.                                    *
 *       SIAM J. Sci. Stat. Comput 5, 470-475.                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * Methods MCORR generates random correlation matrices.                      *
 * It implements two algorithms:                                             *
 *                                                                           *
 * (1) HH: Generate a random matrix H where all rows are independent and     *
 *     uniformly distributed on the sphere and returns HH'; see ref. [1].    *
 *                                                                           *
 * (2) eigen: The algorithm by Marsaglia and Olkin [2] generates a           *
 *     correlation at random with a given set of eigenvalues.                *
 *                                                                           *
 * MCORR uses Algorithm (2) when the set of eigenvalues is given; and        *
 * Algorithm (1) otherwise.
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/matr.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include <utils/unur_fp_source.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "arou.h"
#include "mcorr.h"
#include "mcorr_struct.h"

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

#define MCORR_DEBUG_REINIT    0x00000010u  /* print parameters after reinit  */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

#define MCORR_SET_EIGENVALUES  0x001u   /* set eigenvalues of corr-matrix    */

/*---------------------------------------------------------------------------*/

#define GENTYPE "MCORR"        /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcorr_init( struct unur_par *par );
static int _unur_mcorr_init_HH( struct unur_gen *gen );
static int _unur_mcorr_init_eigen( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_mcorr_reinit( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Reinitialize generator.                                                   */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcorr_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_mcorr_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_mcorr_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/
static int _unur_mcorr_sample_matr_HH( struct unur_gen *gen, double *mat );
/* Algorithm (1) */
static int _unur_mcorr_sample_matr_eigen( struct unur_gen *gen, double *mat );
/* Algorithm (2) */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_mcorr_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_mcorr_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.matr      /* data for distribution object      */

#define PAR       ((struct unur_mcorr_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mcorr_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.matr /* data for distribution in generator object */

#define SAMPLE    gen->sample.matr      /* pointer to sampling routine       */

/*---------------------------------------------------------------------------*/
#define NORMAL    gen->gen_aux        /* pointer to normal variate generator */
/*---------------------------------------------------------------------------*/

#define _unur_mcorr_getSAMPLE(gen) \
   ( ((gen)->set & MCORR_SET_EIGENVALUES) \
     ? _unur_mcorr_sample_matr_eigen : _unur_mcorr_sample_matr_HH )

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_mcorr_new( const struct unur_distr *distr )
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
  if ( !(distr->type == UNUR_DISTR_MATR &&
	 distr->id == UNUR_DISTR_MCORRELATION) ) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_MATR,NULL);

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_mcorr_par) );
  COOKIE_SET(par,CK_MCORR_PAR);

  /* copy input */
  par->distr    = distr;      /* pointer to distribution object              */

  /* set default values */
  par->method   = UNUR_METH_MCORR;    /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* number of rows and columns (dimension of distribution). */
  /* do not confuse with distr->dim which is the size of     */
  /* the array that stores the matrix.                       */
  PAR->dim = distr->data.matr.n_rows;

  PAR->eigenvalues = NULL; /* (optional) eigenvalues of correlation matrix */

  /* routine for starting generator */
  par->init = _unur_mcorr_init;

  return par;

} /* end of unur_mcorr_new() */

/*---------------------------------------------------------------------------*/

int
unur_mcorr_set_eigenvalues( UNUR_PAR *par, const double *eigenvalues )
     /*----------------------------------------------------------------------*/
     /* sets the (optional) eigenvalues of the correlation matrix            */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check arguments */
  _unur_check_NULL( GENTYPE, par, UNUR_ERR_NULL );
  _unur_check_par_object( par, MCORR );
  _unur_check_NULL( GENTYPE, eigenvalues, UNUR_ERR_NULL );

  /* check for eigenvalues */
  for (i=0; i<PAR->dim; i++)
    if (eigenvalues[i] <= 0.) {
      _unur_error(GENTYPE, UNUR_ERR_PAR_SET,"eigenvalue <= 0");
      return UNUR_ERR_PAR_SET;
    }

  /* store date */
  PAR->eigenvalues = eigenvalues;

  /* changelog */
  par->set |= MCORR_SET_EIGENVALUES;

  return UNUR_SUCCESS;
} /* unur_mcorr_set_eigenvalues() */


/*---------------------------------------------------------------------------*/

int
unur_mcorr_chg_eigenvalues( UNUR_GEN *gen, const double *eigenvalues )
     /*----------------------------------------------------------------------*/
     /* change the (optional) eigenvalues of the correlation matrix          */
     /*----------------------------------------------------------------------*/
{
  int i;

  /* check input */
  _unur_check_NULL( GENTYPE, gen, UNUR_ERR_NULL );
  _unur_check_gen_object( gen, MCORR, UNUR_ERR_GEN_INVALID );
  _unur_check_NULL( GENTYPE, eigenvalues, UNUR_ERR_NULL );

  /* check for eigenvalues */
  for (i=0; i<GEN->dim; i++)
    if (eigenvalues[i] <= 0.) {
      _unur_error(GENTYPE, UNUR_ERR_PAR_SET,"eigenvalue <= 0");
      return UNUR_ERR_PAR_SET;
    }

  /* store date */
  if (GEN->eigenvalues == NULL)
    GEN->eigenvalues = _unur_xmalloc(GEN->dim * sizeof(double));
  memcpy(GEN->eigenvalues, eigenvalues, GEN->dim * sizeof(double));

  /* changelog */
  gen->set |= MCORR_SET_EIGENVALUES;

  return UNUR_SUCCESS;
} /* unur_mcorr_chg_eigenvalues() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_mcorr_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_MCORR ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_MCORR_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_mcorr_create(par);
  _unur_par_free(par);
  if (!gen) return NULL; 

  /* run special initialize routines */
  if (gen->set && MCORR_SET_EIGENVALUES) {
    if (_unur_mcorr_init_eigen(gen) != UNUR_SUCCESS) {
      _unur_mcorr_free(gen); return NULL;
    }
  }
  else {
    if (_unur_mcorr_init_HH(gen) != UNUR_SUCCESS) {
      _unur_mcorr_free(gen); return NULL;
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_mcorr_debug_init(gen);
#endif

  /* o.k. */
  return gen;

} /* end of _unur_mcorr_init() */

/*---------------------------------------------------------------------------*/

int
_unur_mcorr_init_HH( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize generator for Algorithm (1)                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* we need a generator for standard normal distributons */

  if (NORMAL==NULL) {
    struct unur_distr *normaldistr = unur_distr_normal(NULL,0);
    struct unur_par   *normalpar = unur_arou_new( normaldistr );

    unur_arou_set_usedars( normalpar, TRUE );
    NORMAL = unur_init( normalpar );
    _unur_distr_free( normaldistr );
    if (NORMAL == NULL) {
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"Cannot create aux Gaussian generator");
      return UNUR_FAILURE;
    }
    /* need same uniform random number generator and debugging flags */
    NORMAL->urng = gen->urng;
    NORMAL->debug = gen->debug;
  }

  return UNUR_SUCCESS;

} /* end of _unur_mcorr_init_HH() */

/*---------------------------------------------------------------------------*/

int
_unur_mcorr_init_eigen( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize generator for Algorithm (2)                               */
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
  double sum_eigenvalues = 0.;

  /* allocate working array */
  GEN->M = _unur_xrealloc(GEN->M, (5*GEN->dim + 2*GEN->dim*GEN->dim)*sizeof(double));

  /* we have to normalize the eigenvalues:    */
  /*   sum(eigenvalues) == dim                */

  /* check for eigenvalues */
  for (i=0; i<GEN->dim; i++) {
    if (GEN->eigenvalues[i] <= 0.) {
      _unur_error(GENTYPE, UNUR_ERR_SHOULD_NOT_HAPPEN,"eigenvalue <= 0");
      return UNUR_FAILURE;
    }
    sum_eigenvalues += GEN->eigenvalues[i];
  }
  
  /* scaling values */
  if (!_unur_FP_equal(sum_eigenvalues, (double) GEN->dim))
    _unur_warning(GENTYPE, UNUR_ERR_GENERIC,"scaling sum(eigenvalues) -> dim");
  for (i=0; i<GEN->dim; i++)
    GEN->eigenvalues[i] *= GEN->dim / sum_eigenvalues;

  return UNUR_SUCCESS;

} /* end of _unur_mcorr_init_eigen() */

/*---------------------------------------------------------------------------*/

int
_unur_mcorr_reinit( struct unur_gen *gen )
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
  SAMPLE = _unur_mcorr_getSAMPLE(gen);

  /* run special initialize routines */
  if (gen->set && MCORR_SET_EIGENVALUES)
    return _unur_mcorr_init_eigen(gen);
  else
    return _unur_mcorr_init_HH(gen);

} /* end of _unur_mcorr_reinit() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mcorr_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_MCORR_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_mcorr_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_MCORR_GEN);

  /* number of rows and columns (dimension of distribution). */
  /* do not confuse with distr->dim which is the size of     */
  /* the array that stores the matrix.                       */
  GEN->dim = DISTR.n_rows;

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_mcorr_getSAMPLE(gen);
  gen->destroy = _unur_mcorr_free;
  gen->clone = _unur_mcorr_clone;
  gen->reinit = _unur_mcorr_reinit;

  /* initialize pointers */
  GEN->M = NULL;
  GEN->H = NULL;
  GEN->eigenvalues = NULL;

  /* copy optional eigenvalues of the correlation matrix */
  if (gen->set && MCORR_SET_EIGENVALUES) {
    GEN->eigenvalues = _unur_xmalloc(GEN->dim * sizeof(double));
    memcpy(GEN->eigenvalues, PAR->eigenvalues, GEN->dim * sizeof(double));
  }

  /* allocate working array */
  if (gen->set && MCORR_SET_EIGENVALUES) {
    GEN->M = _unur_xmalloc((5*GEN->dim + 2*GEN->dim*GEN->dim) * sizeof(double));
  }
  else {
    GEN->H = _unur_xmalloc(GEN->dim * GEN->dim * sizeof(double));
  }

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_mcorr_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;

} /* end of _unur_mcorr_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_mcorr_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_mcorr_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_MCORR_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate new working array */
  if (GEN->M) 
    CLONE->M = _unur_xmalloc((5*GEN->dim + 2*GEN->dim*GEN->dim) * sizeof(double));
  if (GEN->H)
    CLONE->H = _unur_xmalloc(GEN->dim * GEN->dim * sizeof(double));

  /* copy optional eigenvalues */
  if (GEN->eigenvalues) {
    CLONE->eigenvalues = _unur_xmalloc(GEN->dim * sizeof(double));
    memcpy(CLONE->eigenvalues, GEN->eigenvalues, GEN->dim * sizeof(double));
  }
  
  return clone;

#undef CLONE
} /* end of _unur_mcorr_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_mcorr_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_MCORR ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  /* free memory */
  if (GEN->eigenvalues) free(GEN->eigenvalues);
  if (GEN->H)           free(GEN->H);
  if (GEN->M)           free(GEN->M);

  _unur_generic_free(gen);

} /* end of _unur_mcorr_free() */

/*****************************************************************************/

int
_unur_mcorr_sample_matr_HH( struct unur_gen *gen, double *mat )
     /*----------------------------------------------------------------------*/
     /* sample from generator - Algorithm (1)                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   mat ... random matrix (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*(GEN->dim)+(b))
  int i,j,k;
  double sum, norm, x;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MCORR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(mat,UNUR_ERR_NULL);

  /* generate rows vectors of matrix H uniformly distributed in the unit sphere */
  /** TODO: sum != 0 and all columns must be independent **/

  for (i=0; i<GEN->dim; i++) {
    sum=0.;
    for (j=0; j<GEN->dim; j++) {
      x = _unur_sample_cont(NORMAL);
      GEN->H[idx(i,j)] = x;
      sum += x * x;
    }
    norm = sqrt(sum);
    for (j=0; j<GEN->dim; j++) GEN->H[idx(i,j)] /= norm;
  }

  /* Compute HH' */
  for (i=0; i<GEN->dim; i++)
    for (j=0; j<GEN->dim; j++) {
      if (j<i)
	mat[idx(i,j)] = mat[idx(j,i)];
      else if(j==i)
	mat[idx(i,j)] = 1.;
      else {
	sum=0.;
	for (k=0; k<GEN->dim; k++)
	  sum += GEN->H[idx(i,k)]*GEN->H[idx(j,k)];
	mat[idx(i,j)] = sum;
      }
    }

  return UNUR_SUCCESS;
#undef idx
} /* end of _unur_mcorr_sample_matr_HH() */

/*--------------------------------------------------------------------------*/

int
_unur_mcorr_sample_matr_eigen( struct unur_gen *gen, double *mat )
     /*----------------------------------------------------------------------*/
     /* sample from generator - Algorithm (2)                                */
     /* (eigenvalues given, Marsaglia-Olkin method)                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   mat ... random matrix (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))
  int i,j,k, dim;
  double *E, *P;
  double *x, *y, *z, *w, *r; /* misc vectors used in the marsaglia-olkin method */
  double a, b, c, e, e2;
  int s; /* random sign +-1 */

  /* check parameters */
  CHECK_NULL(gen, UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MCORR_GEN,UNUR_ERR_COOKIE);
  CHECK_NULL(mat, UNUR_ERR_NULL);

  dim = GEN->dim; 
  
  if (dim<1) {
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"dimension < 1");
    return UNUR_ERR_SHOULD_NOT_HAPPEN;
  }

  /* initialization steps */
  /* setting working arrays */
  x = GEN->M + (0*dim);
  y = GEN->M + (1*dim);
  z = GEN->M + (2*dim);
  w = GEN->M + (3*dim);
  r = GEN->M + (4*dim);
  E = GEN->M + (5*dim);
  P = GEN->M + (5*dim+dim*dim);

  /* initially E is an identity matrix */
  for (i=0; i<dim; i++)
    for (j=0; j<dim; j++)
      E[idx(i,j)] = (i==j) ? 1 : 0;

  for (k=0; k<dim-1; k++) {
    /* w is a random vector */
    for (i=0; i<dim; i++) w[i] = _unur_call_urng(gen->urng);
    /* x = E*w */
    for (i=0; i<dim; i++) {
      x[i]=0;
      for (j=0; j<dim; j++) {
        x[i] += E[idx(i,j)]*w[j];
      }
    }

    /* a=sum{(1-lambda_i)*x_i*x_i} */
    a=0;
    for (i=0; i<dim; i++)
      a += (1-GEN->eigenvalues[i])*x[i]*x[i];

    /* check if all eigenvalues are ~1 */
    if (fabs(a)<DBL_EPSILON) {
      /* return identity matrix */
      for (i=0; i<dim; i++) {
	for (j=0; j<dim; j++) {
	  mat[idx(i,j)] = (i==j) ? 1: 0;
	}}
      _unur_warning(gen->genid, UNUR_ERR_GEN_CONDITION,"all eigenvalues are ~1 -> identity matrix");
      
      return UNUR_ERR_GEN_CONDITION;
    }

    do {

      /* z is a random vector */
      for (i=0; i<dim; i++) z[i] = _unur_call_urng(gen->urng);

      /* y = E*z */
      for (i=0; i<dim; i++) {
        y[i]=0;
        for (j=0; j<dim; j++) {
          y[i] += E[idx(i,j)]*z[j];
        }
      }

      /* b=sum{(1-lambda_i)*x_i*y_i} */
      /* c=sum{(1-lambda_i)*y_i*y_i} */
      b=0; c=0;
      for (i=0; i<dim; i++) {
        b += (1-GEN->eigenvalues[i])*x[i]*y[i];
        c += (1-GEN->eigenvalues[i])*y[i]*y[i];
      }

      /* e^2 = b^2 - a*c */
      e2 = b*b - a*c;

    } while (e2<0);

    e=sqrt(e2);


    /* random sign */
    s = ( _unur_call_urng(gen->urng) >.5) ? 1: -1 ;

    /* r=x*(b+s*e)/a - y */
    for (i=0; i<dim; i++) r[i] = x[i]*(b+s*e)/a - y[i];

    /* another random sign */
    s = ( _unur_call_urng(gen->urng) >.5) ? 1: -1 ;

    /* pk=s*r/norm(r) */
    _unur_vector_normalize(dim, r);
    for (i=0; i<dim; i++) P[idx(k,i)] = s * r[i];

    /* E = E - r r^T */
    for (i=0; i<dim; i++) {
      for (j=0; j<dim; j++) {
        E[idx(i,j)] -= r[i]*r[j];
      }
    }

  } /* next k */

  /* w is a random vector */
  for (i=0; i<dim; i++) w[i] = _unur_call_urng(gen->urng);

  /* x = E*w */
  for (i=0; i<dim; i++) {
    x[i]=0;
    for (j=0; j<dim; j++) {
      x[i] += E[idx(i,j)]*w[j];
    }
  }

  _unur_vector_normalize(dim, x);

  /* last row of the orthogonal matrix P */
  for (i=0; i<dim; i++) {
    P[idx(dim-1,i)] = x[i];
  }

  /* mat = P L P^T, where L diagonal containing the eigenvalues */
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      mat[idx(i,j)] = 0;
      for (k=0; k<dim; k++) {
        mat[idx(i,j)] += P[idx(i,k)] * GEN->eigenvalues[k] * P[idx(j,k)];
      }
    }
  }

  /* symmetrization (necessary due to rounding-errors) */
  for (i=0; i<dim; i++) {
    for (j=(i+1); j<dim; j++) {
      mat[idx(i,j)] = (mat[idx(i,j)]+mat[idx(j,i)])/2.;
      mat[idx(j,i)] = mat[idx(i,j)];
    }
  }  

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_mcorr_sample_eigen() */

/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

static void
_unur_mcorr_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_MCORR_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = random matrix\n",gen->genid);
  fprintf(LOG,"%s: method  = MCORR (Matrix - CORRELATION matrix)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_matr_debug( gen->distr, gen->genid );

  /* eigenvalues */
  if (gen->set && MCORR_SET_EIGENVALUES)
    _unur_matrix_print_vector( GEN->dim, GEN->eigenvalues, "eigenvalues =", LOG, gen->genid, "\t   ");

  if (gen->set && MCORR_SET_EIGENVALUES)
    fprintf(LOG,"%s: sampling routine = _unur_mcorr_sample_matr_eigen()\n",gen->genid);
  else
    fprintf(LOG,"%s: sampling routine = _unur_mcorr_sample_matr_HH()\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_mcorr_debug_init() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_mcorr_info( struct unur_gen *gen, int help )
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
  _unur_string_append(info,"   dimension = %d x %d   (= %d)\n", 
		      DISTR.n_rows, DISTR.n_cols, gen->distr->dim);

  if (gen->set && MCORR_SET_EIGENVALUES) {
    _unur_string_append(info,"   eigenvalues = ");
    _unur_distr_info_vector( gen, GEN->eigenvalues, GEN->dim);
    _unur_string_append(info,"\n");
  }
  _unur_string_append(info,"\n");

  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */
  
  /* method */
  _unur_string_append(info,"method: MCORR (Random CORRelation matrix)\n");
  if (gen->set && MCORR_SET_EIGENVALUES)
    _unur_string_append(info,"   generate correlation matrix with given eigenvalues\n");
  _unur_string_append(info,"\n");

  /* performance */
  /*   _unur_string_append(info,"performance characteristics:\n"); */
  /*   _unur_string_append(info,"\n"); */

  /* parameters */
  if (help) {
    _unur_string_append(info,"parameters: \n");
    if (gen->set && MCORR_SET_EIGENVALUES) {
      _unur_string_append(info,"   eigenvalues = ");
      _unur_distr_info_vector( gen, GEN->eigenvalues, GEN->dim);
      _unur_string_append(info,"\n");
    }
    _unur_string_append(info,"\n");
  }

  /* Hints */
  /*   if (help) { */
  /*     _unur_string_append(info,"\n"); */
  /*   } */

} /* end of _unur_mcorr_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
