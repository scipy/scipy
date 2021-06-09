/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      norta.c                                                      *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    generates random vector with given covariance matrix and     *
 *              given marginal distribution.                                 *
 *                                                                           *
 *   DESCRIPTION:                                                            *
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
 *   [1] Hoermann, W., J. Leydold, and G. Derflinger (2004):                 *
 *       Automatic Nonuniform Random Variate Generation, Springer, Berlin.   *
 *                                                                           *
 *   [2] Gosh, S. (????): Eigenvector Correction method                      *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   NORTA (NORmal to anything) is a model to get random vectors with        *
 *   given marginal distributions and rank correlation. (Notice that this    *
 *   does not uniquely define the multivariate distribution.)                *
 *                                                                           *
 *   In the NORTA model multinormal random variates with the given           *
 *   rank (Spearman's) correlations are generated. For this task the rank    *
 *   correlations are transformed into the corresponding                     *
 *   (Pearson) correlation matrix and the VMT method is used to sample from  *
 *   the resulting multinormal distribution (by means of the Cholesky        *
 *   decomposition of the covariance matrix).                                *
 *   In a second step the (standard normal distributed) marginal variates    *
 *   are transformed by means of the CDF of the normal distribution to get   *
 *   uniform marginals. The resulting random vectors then have uniform       *
 *   marginals and the desired rank correlation between its components.      *
 *   Such a random vector is called 'copula'.                                *                   
 *                                                                           *
 *   By means of the inverse CDF the uniform marginals are transformed into  *
 *   the target marginal distributions. This transformation does not change  *
 *   the rank correlation.                                                   *
 *                                                                           *
 *   It can happen that the desired rank correlation matrix is not feasible, *
 *   i.e., it cannot occur as rank correlation matrix of a multinormal       *
 *   distribution. For this case we use Gosh's "eigenvector correction       *
 *   method" [2]. In this case a spectral decomposition of the corresponding *
 *   (Pearson) correlation matrix is carried out. If there are some          *
 *   non-positive eigenvalues (i.e. this matrix is not a true correlation    *
 *   matrix) then these are set to a small value and a new (positive         *
 *   definite) correlation matrix is created and used.                       *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <distr/cont.h>
#include <distributions/unur_distributions.h>
#include <urng/urng.h>
#include <utils/matrix_source.h>
#include "unur_methods_source.h"
#include "x_gen.h"
#include "x_gen_source.h"
#include "cstd.h"
#include "hinv.h"
#include "ninv.h"
#include "pinv.h"
#include "mvstd.h"
#include "norta.h"
#include "norta_struct.h"

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* smallest eigenvalue allowed for correlation matrix                        */
#define UNUR_NORTA_MIN_EIGENVALUE  (1.e-10)

/*---------------------------------------------------------------------------*/
/* Variants                                                                  */

/*---------------------------------------------------------------------------*/
/* Debugging flags                                                           */
/*    bit  01    ... pameters and structure of generator (do not use here)   */
/*    bits 02-12 ... setup                                                   */
/*    bits 13-24 ... adaptive steps                                          */
/*    bits 25-32 ... trace sampling                                          */

#define NORTA_DEBUG_SIGMA_Y     0x00000010u   /* print sigma_y for normal    */

/*---------------------------------------------------------------------------*/
/* Flags for logging set calls                                               */

/*---------------------------------------------------------------------------*/

#define GENTYPE "NORTA"          /* type of generator                          */

/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_norta_init( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* Initialize new generator.                                                 */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_norta_create( struct unur_par *par );
/*---------------------------------------------------------------------------*/
/* create new (almost empty) generator object.                               */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_norta_clone( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* copy (clone) generator object.                                            */
/*---------------------------------------------------------------------------*/

static void _unur_norta_free( struct unur_gen *gen);
/*---------------------------------------------------------------------------*/
/* destroy generator object.                                                 */
/*---------------------------------------------------------------------------*/

static int _unur_norta_sample_cvec( struct unur_gen *gen, double *vec );
/*---------------------------------------------------------------------------*/
/* sample from generator                                                     */
/*---------------------------------------------------------------------------*/

static int _unur_norta_nortu_setup( struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* Compute parameter for NORTU (normal to uniform).                          */
/*---------------------------------------------------------------------------*/

static int _unur_norta_make_correlationmatrix( int dim, double *M);
/*---------------------------------------------------------------------------*/
/* make correlation matrix by transforming symmetric positive definit matrix */
/*---------------------------------------------------------------------------*/

static struct unur_gen *_unur_norta_make_marginalgen( const struct unur_gen *gen,
						      const struct unur_distr *marginal );
/*---------------------------------------------------------------------------*/
/* make generator for marginal distribution                                  */
/*---------------------------------------------------------------------------*/

#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/
/* the following functions print debugging information on output stream,     */
/* i.e., into the LOG file if not specified otherwise.                       */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_init( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print after generator has been initialized has completed.                 */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_sigma_y( const struct unur_gen *gen, 
				       const double *sigma_y, 
				       const char *comment );
/*---------------------------------------------------------------------------*/
/* print sigma_y of corresponding normal distribution.                       */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_eigensystem( const struct unur_gen *gen,
					   const double *eigenvalues,
					   const double *eigenvectors );
/*---------------------------------------------------------------------------*/
/* print eigensystem of sigma_y.                                             */
/*---------------------------------------------------------------------------*/

static void _unur_norta_debug_nmgenerator( const struct unur_gen *gen );
/*---------------------------------------------------------------------------*/
/* print genid of multinormal generator.                                     */
/*---------------------------------------------------------------------------*/
#endif

#ifdef UNUR_ENABLE_INFO
static void _unur_norta_info( struct unur_gen *gen, int help );
/*---------------------------------------------------------------------------*/
/* create info string.                                                       */
/*---------------------------------------------------------------------------*/
#endif

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define DISTR_IN  distr->data.cvec      /* data for distribution object      */

#define PAR       ((struct unur_norta_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_norta_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

#define SAMPLE    gen->sample.cvec      /* pointer to sampling routine       */     

#define MNORMAL   gen->gen_aux          /* pointer to multinormal generator  */

/*---------------------------------------------------------------------------*/

#define _unur_norta_getSAMPLE(gen)   (_unur_norta_sample_cvec)

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**  Public: User Interface (API)                                           **/
/*****************************************************************************/

struct unur_par *
unur_norta_new( const struct unur_distr *distr )
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
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(distr,CK_DISTR_CVEC,NULL);

  if (!(distr->set & UNUR_DISTR_SET_RANKCORR)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"rank correlation matrix");
    return NULL; }

  if (!(distr->set & UNUR_DISTR_SET_MARGINAL)) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_REQUIRED,"marginals");
    return NULL; }

  /* allocate structure */
  par = _unur_par_new( sizeof(struct unur_norta_par) );
  COOKIE_SET(par,CK_NORTA_PAR);

  /* copy input */
  par->distr    = distr;            /* pointer to distribution object        */

  /* set default values */
  par->method   = UNUR_METH_NORTA ;   /* method                              */
  par->variant  = 0u;                 /* default variant                     */
  par->set      = 0u;                 /* inidicate default parameters        */    
  par->urng     = unur_get_default_urng(); /* use default urng               */
  par->urng_aux = NULL;                    /* no auxilliary URNG required    */

  par->debug    = _unur_default_debugflag; /* set default debugging flags    */

  /* routine for starting generator */
  par->init = _unur_norta_init;

  return par;

} /* end of unur_norta_new() */


/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

struct unur_gen *
_unur_norta_init( struct unur_par *par )
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
  if ( par->method != UNUR_METH_NORTA ) {
    _unur_error(GENTYPE,UNUR_ERR_PAR_INVALID,"");
    return NULL; }
  COOKIE_CHECK(par,CK_NORTA_PAR,NULL);

  /* create a new empty generator object */
  gen = _unur_norta_create(par);
  _unur_par_free(par);
  if (!gen) return NULL;

  /* check for truncated domain */
  if ( gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED ) {
    if ( DISTR.domainrect == NULL ) {
      /* cannot handle non-rectangular domain */
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"cannot handle non-rectangular domain");
      _unur_norta_free(gen); return NULL;
    }
    else { /* rectangular domain */
      if (_unur_distr_cvec_marginals_are_equal(DISTR.marginals, GEN->dim)) {
	/* we have to handle these equal marginal distributions independently */
	if ( _unur_distr_cvec_duplicate_firstmarginal(gen->distr) != UNUR_SUCCESS ) {
	  _unur_norta_free(gen); return NULL;
	}
      }
    }
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug) _unur_norta_debug_init(gen);
#endif

  /* compute parameters for NORTU (normal to uniform) */
  if (_unur_norta_nortu_setup(gen) != UNUR_SUCCESS) {
    _unur_norta_free(gen); return NULL;
  }

  /* distribution object for standard normal distribution */
  GEN->normaldistr = unur_distr_normal(NULL,0);

  if (gen->distr->id != UNUR_DISTR_COPULA) {
    /* we have to initialize generator for marginal distributions */

    if (_unur_distr_cvec_marginals_are_equal(DISTR.marginals, GEN->dim)) {
      /* we can use the same generator object for all marginal distribuitons */
      struct unur_gen *marginalgen = _unur_norta_make_marginalgen( gen, DISTR.marginals[0] );
      if (marginalgen)
	GEN->marginalgen_list = _unur_gen_list_set(marginalgen,GEN->dim);
    }

    else {
      /* we need different generator objects for all marginal distribuitons */
      int i,j;
      int failed = FALSE;
      struct unur_gen **marginalgens = _unur_xmalloc( GEN->dim * sizeof(struct unur_gen*) );

      /* check for truncated domain */
      if ( gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED ) {
	/* set domains for each marginal distribution */
	for (i=0; i<GEN->dim && !failed; i++) {
	  if ( (unur_distr_cont_set_domain(DISTR.marginals[i],
					   DISTR.domainrect[2*i], DISTR.domainrect[2*i+1]))
	       != UNUR_SUCCESS) {
	    failed = TRUE; break;
	  }
	}
      }
      
      /* create generator for each marginal distribution */
      for (i=0; i<GEN->dim && !failed; i++) {
	marginalgens[i] = _unur_norta_make_marginalgen( gen, DISTR.marginals[i] );
	if (marginalgens[i]==NULL) {
	  failed=TRUE; break;
	}
      }

      /* check creating of generators */
      if (failed) {
	for (j=0; j<i; j++) _unur_free(marginalgens[j]);
	free (marginalgens);
      }
      else
	GEN->marginalgen_list = marginalgens;
    }

    /* verify initialization of marginal generators */
    if (GEN->marginalgen_list == NULL) {
      _unur_error(gen->genid,UNUR_ERR_GENERIC,"init of marginal generators failed");
      _unur_norta_free(gen);
      return NULL;
    }
  }

  /* o.k. */
  return gen;

} /* end of _unur_norta_init() */

/*---------------------------------------------------------------------------*/

static struct unur_gen *
_unur_norta_create( struct unur_par *par )
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
  CHECK_NULL(par,NULL);  COOKIE_CHECK(par,CK_NORTA_PAR,NULL);

  /* create new generic generator object */
  gen = _unur_generic_create( par, sizeof(struct unur_norta_gen) );

  /* magic cookies */
  COOKIE_SET(gen,CK_NORTA_GEN);

  /* set generator identifier */
  gen->genid = _unur_set_genid(GENTYPE);

  /* routines for sampling and destroying generator */
  SAMPLE = _unur_norta_getSAMPLE(gen);
  gen->destroy = _unur_norta_free;
  gen->clone = _unur_norta_clone;

  /* dimension of distribution */
  GEN->dim = gen->distr->dim;

  /* allocate array for auxiliary copula */
  GEN->copula = _unur_xmalloc(sizeof(double)*GEN->dim);

  /* initialize pointer */
  MNORMAL = NULL;
  GEN->normaldistr = NULL;
  GEN->marginalgen_list = NULL;

#ifdef UNUR_ENABLE_INFO
  /* set function for creating info string */
  gen->info = _unur_norta_info;
#endif

  /* return pointer to (almost empty) generator object */
  return gen;
  
} /* end of _unur_norta_create() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_norta_clone( const struct unur_gen *gen )
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
#define CLONE  ((struct unur_norta_gen*)clone->datap)

  struct unur_gen *clone;

  /* check arguments */
  CHECK_NULL(gen,NULL);  COOKIE_CHECK(gen,CK_NORTA_GEN,NULL);

  /* create generic clone */
  clone = _unur_generic_clone( gen, GENTYPE );

  /* allocate array for auxiliary copula */
  CLONE->copula = _unur_xmalloc(sizeof(double)*GEN->dim);

  /* clone marginal distribution */
  CLONE->normaldistr = _unur_distr_clone(GEN->normaldistr);

  /* marginal generators */
  if (GEN->marginalgen_list)
    CLONE->marginalgen_list = _unur_gen_list_clone( GEN->marginalgen_list, GEN->dim );

  return clone;

#undef CLONE
} /* end of _unur_norta_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_norta_free( struct unur_gen *gen )
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
  if ( gen->method != UNUR_METH_NORTA ) {
    _unur_warning(gen->genid,UNUR_ERR_GEN_INVALID,"");
    return; }
  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  /* free auxiliary arrays */
  if (GEN->copula) free (GEN->copula);

  /* free normal distribution object */
  if (GEN->normaldistr) _unur_distr_free (GEN->normaldistr);

  /* free marginal generators */
  if (GEN->marginalgen_list)
    _unur_gen_list_free( GEN->marginalgen_list, GEN->dim);

  /* we cannot use this generator object any more */
  SAMPLE = NULL;   /* make sure to show up a programming error */

  _unur_generic_free(gen);

} /* end of _unur_norta_free() */

/*****************************************************************************/

int
_unur_norta_sample_cvec( struct unur_gen *gen, double *vec )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*GEN->dim+(b))
  int j;
  double *u;

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_NORTA_GEN,UNUR_ERR_COOKIE);

  /* pointer to auxiliary array of uniforms */
  u = GEN->copula;

  /* sample from multinormal distribution */
  _unur_sample_vec(MNORMAL,u);

  /* make copula */
  for (j=0; j<GEN->dim; j++)
    vec[j] = unur_distr_cont_eval_cdf( u[j], GEN->normaldistr );

  if (gen->distr->id == UNUR_DISTR_COPULA)
    /* we want to have a normal copula --> just return data */
    return UNUR_SUCCESS;
  
  /* else non-uniform marginals */
  for (j=0; j<GEN->dim; j++) {
    /* call marginal generator (using this U-value) */
    vec[j] = unur_quantile(GEN->marginalgen_list[j], vec[j]);
  }

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_norta_sample_cvec() */


/*****************************************************************************/
/**  Auxilliary Routines                                                    **/
/*****************************************************************************/

int
_unur_norta_nortu_setup( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* Compute parameter for NORTU (Normal to uniform)                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS if computed successfully                              */
     /*   error code   otherwise                                             */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int dim = GEN->dim;    /* dimension of distribution */
  double *sigma_y;      /* correlation matrix for corresponding normal distr. */
  double *eigenvalues;  /* eigenvalues of sigma_y */
  double *eigenvectors; /* eigenvectors of sigma_y */
  int eigenvalues_positive; /* boolean indicating whether all eigenvalues are 
			       strictly positive */
  struct unur_distr *mn_distr; /* multinormal distribution */ 
  struct unur_gen   *mn_gen;   /* generator for multinormal distribution */ 
  int i,j;

  /* setup correlation matrix for corresponding normal distribution */
  sigma_y = _unur_xmalloc(dim * dim * sizeof(double));
  for(i=0; i<dim; i++) {
    /* left lower part: make matrix symmetric */
    for(j=0; j<i; j++)
      sigma_y[idx(i,j)] = sigma_y[idx(j,i)];
    /*   diagonal */
    sigma_y[idx(i,i)] = 1.;
    /* right upper part */
    for(j=i+1; j<dim; j++)
      sigma_y[idx(i,j)] = 2.*sin(DISTR.rankcorr[idx(i,j)]*(M_PI/6.));  
  }

#ifdef UNUR_ENABLE_LOGGING
  /* write info into LOG file */
  if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
    _unur_norta_debug_sigma_y( gen, sigma_y, "NORTU setup:" );
#endif

  /* compute eigenvectors and eigenvalues of sigma_y */
  eigenvalues = _unur_xmalloc(dim * sizeof(double));
  eigenvectors = _unur_xmalloc(dim * dim * sizeof(double));
  if (_unur_matrix_eigensystem(dim, sigma_y, eigenvalues, eigenvectors) != UNUR_SUCCESS) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"cannot compute eigenvalues for given sigma_y");
    free(sigma_y); free(eigenvalues); free(eigenvectors);
    return UNUR_ERR_GEN_DATA;
  }

#ifdef UNUR_ENABLE_LOGGING
  if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
    _unur_norta_debug_eigensystem( gen, eigenvalues, eigenvectors );
#endif

  /* check if all eigenvalues are positive */
  /* otherwise set to small values close to 0 */
  eigenvalues_positive = TRUE;
  for(i=0; i<dim; i++)
    if(eigenvalues[i] < UNUR_NORTA_MIN_EIGENVALUE) {
      eigenvalues[i] = UNUR_NORTA_MIN_EIGENVALUE;
      eigenvalues_positive = FALSE;
    }

  /* make corrected correlation matrix */

  if (!eigenvalues_positive) {
    _unur_matrix_transform_diagonal(dim,eigenvectors,eigenvalues,sigma_y);
    _unur_norta_make_correlationmatrix(dim,sigma_y);
    _unur_warning(GENTYPE,UNUR_ERR_GEN_DATA,
		  "sigma_y not positive definite -> corrected matrix");
#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
      _unur_norta_debug_sigma_y( gen, sigma_y, "\tEigenvalues < 0 --> correction required" );
#endif
  }

  /* clear working arrays */
  free(eigenvalues);
  free(eigenvectors);

  /* make generator for multinormal distribution */
  mn_distr = unur_distr_multinormal(dim, NULL, sigma_y);
  mn_gen = NULL;
  if (mn_distr) {
    mn_gen = unur_init(unur_mvstd_new(mn_distr));
    _unur_distr_free(mn_distr);
  }
  if (mn_gen == NULL) {
    _unur_error(GENTYPE,UNUR_ERR_GEN_DATA,"(corrected) sigma_y not positive definit");
    free(sigma_y);
    return UNUR_ERR_GEN_DATA;
  }
  MNORMAL = mn_gen;
  /* need same uniform random number generator as NORTA generator */
  MNORMAL->urng = gen->urng;
  /* copy debugging flags */
  MNORMAL->debug = gen->debug;

#ifdef UNUR_ENABLE_LOGGING
    if (gen->debug & NORTA_DEBUG_SIGMA_Y) 
      _unur_norta_debug_nmgenerator( gen );
#endif

  /* clear working arrays */
  free(sigma_y);

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_norta_nortu_setup() */

/*---------------------------------------------------------------------------*/

int
_unur_norta_make_correlationmatrix( int dim, double *M)
     /*----------------------------------------------------------------------*/
     /* make correlation matrix by transforming symmetric positive definit   */
     /* matrix.                                                              */
     /*                                                                      */
     /* There is no checking whether M fulfills the conditions!              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   dim ... dimension of matrix M                                      */
     /*   M   ... symmetric square matrix                                    */
     /*----------------------------------------------------------------------*/
{
#define idx(a,b) ((a)*dim+(b))

  int i,j;

  /* diagonal is used to store the roots of the diagonal elements */
  for (i=0; i<dim; i++)
    M[idx(i,i)] = sqrt(M[idx(i,i)]);

  for (i=0; i<dim; i++)
    for (j=i; j<dim; j++)
      if(i!=j) {
	M[idx(i,j)] /= M[idx(i,i)] * M[idx(j,j)];
	M[idx(j,i)] = M[idx(i,j)];
      }

  /* the diagonal elements are set to 1. */
  for (i=0; i<dim; i++) 
    M[idx(i,i)] = 1.;

  return UNUR_SUCCESS;
#undef idx
} /* end of _unur_norta_make_correlationmatrix() */

/*---------------------------------------------------------------------------*/

struct unur_gen *
_unur_norta_make_marginalgen( const struct unur_gen *gen,
			      const struct unur_distr *marginal )
     /*----------------------------------------------------------------------*/
     /* make generator for marginal distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen      ... pointer to multivariate generator object              */
     /*   marginal ... distribution object for marginal distribution         */
     /*                                                                      */
     /* return:                                                              */
     /*   generator for marginal distribution                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_gen *marginalgen;
  struct unur_par *par;

  /* check arguments */
  CHECK_NULL(gen,NULL);      COOKIE_CHECK(gen,CK_NORTA_GEN,NULL);
  CHECK_NULL(marginal,NULL);

  /* check distribution */
  if (marginal->type != UNUR_DISTR_CONT) {
    _unur_error(GENTYPE,UNUR_ERR_DISTR_INVALID,""); return NULL; }
  COOKIE_CHECK(marginal,CK_DISTR_CONT,NULL);

  /* try PINV, CSTD+Inversion, HINV, and NINV */
  do {
    /* PINV (Polynomial interpolation based INVersion of CDF) */
    par = unur_pinv_new( marginal );
    if ( (marginalgen = _unur_init(par)) != NULL )
      break;

    /* CSTD + inversion */
    par = unur_cstd_new( marginal );
    if (unur_cstd_set_variant( par, UNUR_STDGEN_INVERSION)==UNUR_SUCCESS) {
      marginalgen = _unur_init(par);
      break;
    }
    else {
      _unur_par_free(par);
    }

    /* HINV (Hermite interpolation based INVersion of CDF) */
    par = unur_hinv_new( marginal );
    if ( (marginalgen = _unur_init(par)) != NULL )
      break;

    /* NINV (Numerical inversion with regula falsi */
    par = unur_ninv_new( marginal );
    unur_ninv_set_table( par, 100 );
    if ( (marginalgen = _unur_init(par)) != NULL ) 
      break;

    /* no inversion method avaiblable */
    _unur_error(gen->genid,UNUR_ERR_DISTR_REQUIRED,"data for (numerical) inversion of marginal missing");
    return NULL;
    
  } while (1);

  /* set debugging flags */
  marginalgen->debug = gen->debug;

  /* return generator object */
  return marginalgen;

} /* end of _unur_norta_make_marginalgen() */


/*****************************************************************************/
/**  Debugging utilities                                                    **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_init( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* write info about generator into LOG file                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*----------------------------------------------------------------------*/
{
/*   int i; */
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s:\n",gen->genid);
  fprintf(LOG,"%s: type    = continuous multivariate random variates\n",gen->genid);
  fprintf(LOG,"%s: method  = NORTA (Vector Matrix Transformation)\n",gen->genid);
  fprintf(LOG,"%s:\n",gen->genid);

  _unur_distr_cvec_debug( gen->distr, gen->genid );

} /* end of _unur_norta_debug_init() */

/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_sigma_y( const struct unur_gen *gen, 
			   const double *sigma_y, 
			   const char *comment )
     /*----------------------------------------------------------------------*/
     /* print sigma_y of corresponding normal distribution.                  */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen     ... pointer to generator object                            */
     /*   sigma_y ... pointer to correlation matrix                          */
     /*   comment ... additional string printed                              */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: %s\n",gen->genid,comment);
  fprintf(LOG,"%s:\n",gen->genid);
  _unur_matrix_print_matrix( GEN->dim, sigma_y, "\tsigma_y =", 
			     LOG, gen->genid, "\t   ");

} /* end of _unur_norta_debug_sigma_y() */

/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_eigensystem( const struct unur_gen *gen,
			       const double *eigenvalues,
			       const double *eigenvectors )
     /*----------------------------------------------------------------------*/
     /* print eigensystem of sigma_y.                                        */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*   eigenvalues  ... eigenvalues                                       */
     /*   eigenvectors ... eigenvalues                                       */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  _unur_matrix_print_vector( GEN->dim, eigenvalues, 
			     "\teigenvalues of sigma_y =", 
			     LOG, gen->genid, "\t   ");
  _unur_matrix_print_matrix( GEN->dim, eigenvectors, 
			     "\teigenvectors of sigma_y [rows] =", 
			     LOG, gen->genid, "\t   ");

} /* end of _unur_norta_debug_eigensystem() */

/*---------------------------------------------------------------------------*/

void
_unur_norta_debug_nmgenerator( const struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* print genid of multinormal generator.                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen          ... pointer to generator object                       */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;

  /* check arguments */
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_NORTA_GEN,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: generator for multinormal auxiliary distribution = %s\n", gen->genid,
	  MNORMAL->genid );
  fprintf(LOG,"%s:\n",gen->genid);

} /* end of _unur_norta_debug_nmgenerator() */

/*---------------------------------------------------------------------------*/
#endif   /* UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_norta_info( struct unur_gen *gen, int help )
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
  struct unur_distr *distr = gen->distr;
  int i;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  
  /* distribution */
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = MARGINAL distributions\n");

  _unur_string_append(info,"   marginals =");
  for (i=0; i<distr->dim; i++)
    _unur_string_append(info," %s", distr->data.cvec.marginals[i]->name);
  _unur_string_append(info,"\n\n");
  
  /*   if (help) { */
  /*   _unur_string_append(info,"\n"); */
  /*   } */

  /* method */
  _unur_string_append(info,"method: NORTA (NORmal To Anything)\n");
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

} /* end of _unur_norta_info() */

/*---------------------------------------------------------------------------*/
#endif   /* end UNUR_ENABLE_INFO */
/*---------------------------------------------------------------------------*/
