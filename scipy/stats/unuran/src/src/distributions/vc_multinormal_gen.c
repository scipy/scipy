/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      vc_multinormal_gen.c                                         *
 *                                                                           *
 *   Special generators for MultiNormal distribution                         *
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
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr_source.h>
#include <methods/x_gen.h>
#include <methods/cstd.h>
#include <methods/mvstd.h>
#include <methods/mvstd_struct.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"

/*---------------------------------------------------------------------------*/
/* init routines for special generators                                      */

static int _unur_stdgen_init_multinormal_cholesky( struct unur_gen *gen );

/*---------------------------------------------------------------------------*/
/* abbreviations */

#define PAR       ((struct unur_mvstd_par*)par->datap) /* data for parameter object */
#define GEN       ((struct unur_mvstd_gen*)gen->datap) /* data for generator object */
#define DISTR     gen->distr->data.cvec /* data for distribution in generator object */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Inititialize                                                           **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_stdgen_multinormal_init( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* initialize special generator for multinormal distribution            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check parameters */
  /* the generator does not work correctly when the domain has been truncated*/
  if ( gen->distr->set & UNUR_DISTR_SET_DOMAINBOUNDED ) {
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"truncated domain not allowed");
    return UNUR_FAILURE;
  }

  /* set sampling routine */
  gen->sample.cvec = _unur_stdgen_sample_multinormal_cholesky;

  /* set routine name */
  GEN->sample_routine_name = "_unur_stdgen_sample_multinormal_cholesky";

  /* run special init routine */
  return _unur_stdgen_init_multinormal_cholesky(gen);
} /* end of _unur_stdgen_multinormal_init() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/**  Special generators                                                     **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

/*****************************************************************************
 *                                                                           *
 * MultiNormal Distribution: Cholesky decomposition                          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FUNCTION:   - samples a random number from the                            *
 *               MultiNormal distribution N(mu,sigma).                       *
 *                                                                           *
 *****************************************************************************
 * UNURAN (c) 2007  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien *
 *****************************************************************************/

#define NORMAL  gen->gen_aux        /* pointer to normal variate generator   */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_init_multinormal_cholesky( struct unur_gen *gen )
{
  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MVSTD_GEN,UNUR_ERR_COOKIE);

  /* we do not allow truncated multinormal distributions.     */
  /* [ currently checked in routine _unur_mvstd_check_par() ] */

  /* make a normal variate generator */
  if (NORMAL==NULL) {
    /* make generator object */
    struct unur_distr *distr = unur_distr_normal(NULL,0);
    NORMAL = unur_init( unur_cstd_new( distr ) );
    _unur_check_NULL( gen->genid, NORMAL, UNUR_ERR_NULL );
    /* use same uniform random number generator */
    NORMAL->urng = gen->urng;
    /* copy debugging flags */
    NORMAL->debug = gen->debug;
    /* we do not need the distribution object any more */
    _unur_distr_free( distr );
  }
  /* else we are in the re-init mode 
     --> there is no necessity to make the generator object again */

  return UNUR_SUCCESS;
} /* end of _unur_stdgen_init_multinormal_cholesky() */

/*---------------------------------------------------------------------------*/

int
_unur_stdgen_sample_multinormal_cholesky( struct unur_gen *gen, double *X )
{
#define idx(a,b) ((a)*dim+(b))
  int j,k;

  int dim = gen->distr->dim;     /* dimension of distribution */
  double *L = DISTR.cholesky;    /* cholesky factor of covariance matrix  */
  double *mean = DISTR.mean;     /* mean (mode) of distribution */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MVSTD_GEN,UNUR_ERR_COOKIE);

  /* generate random vector with independent components */
  for (j=0; j<dim; j++)
    X[j] = unur_sample_cont(NORMAL);
    
  /*
    transform to desired covariance structure:
    X = L.Y + mu
    where
    L  ... cholesky factor of the covariance matrix
    Y  ... vector with indenpent components (generated above)
    mu ... mean vector
    (notice that L is a lower triangular matrix)
  */
  for (k=dim-1; k>=0; k--) {
    X[k] *= L[idx(k,k)];
    for (j=k-1; j>=0; j--)
      X[k] += X[j] * L[idx(k,j)];
    X[k] += mean[k];
  }

  return UNUR_SUCCESS;

#undef idx
} /* end of _unur_stdgen_sample_multinormal_cholesky() */

/*---------------------------------------------------------------------------*/
#undef NORMAL
/*---------------------------------------------------------------------------*/
