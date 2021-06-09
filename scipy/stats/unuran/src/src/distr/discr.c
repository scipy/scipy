/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: discr.c                                                           *
 *                                                                           *
 *   manipulate univariate discrete distribution objects                     *
 *                                                                           *
 *   return:                                                                 *
 *     UNUR_SUCCESS ... on success                                           *
 *     error code   ... on error                                             *
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
#include <distributions/unur_stddistr.h>
#include <parser/functparser_source.h>
#include "distr_source.h"
#include "distr.h"
#include "discr.h"

/*---------------------------------------------------------------------------*/

#define DISTR distr->data.discr

/*---------------------------------------------------------------------------*/
/* Constants */

/* maximum size of domain for which the pmfsum is computed automatically */
#define MAX_PMF_DOMAIN_FOR_UPD_PMFSUM   (1000)

/*---------------------------------------------------------------------------*/

static double _unur_distr_discr_eval_pmf_tree( int k, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for PMF.                                           */
/*---------------------------------------------------------------------------*/

static double _unur_distr_discr_eval_cdf_tree( int k, const struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* evaluate function tree for CDF.                                           */
/*---------------------------------------------------------------------------*/

static void _unur_distr_discr_free( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* destroy distribution object.                                              */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/

static int _unur_distr_discr_find_mode( struct unur_distr *distr );
/*---------------------------------------------------------------------------*/
/* find mode of unimodal probability vector numerically by bisection.        */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/**                                                                         **/
/** univariate discrete distributions                                       **/
/**                                                                         **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

struct unur_distr *
unur_distr_discr_new( void )
     /*----------------------------------------------------------------------*/
     /* create a new (empty) distribution object                             */
     /* type: univariate discete                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   none                                                               */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to distribution object                                     */
     /*                                                                      */
     /* error:                                                               */
     /*   return NULL                                                        */
     /*----------------------------------------------------------------------*/
{
  register struct unur_distr *distr;
  register int i;

  /* get empty distribution object */
  distr = _unur_distr_generic_new();
  if (!distr) return NULL;

  /* set magic cookie */
  COOKIE_SET(distr,CK_DISTR_DISCR);

  /* set type of distribution */
  distr->type = UNUR_DISTR_DISCR;

  /* set id to generic distribution */
  distr->id = UNUR_DISTR_GENERIC;

  /* dimension of random vector */
  distr->dim = 1;   /* univariant */

  /* destructor */
  distr->destroy = _unur_distr_discr_free;

  /* clone */
  distr->clone = _unur_distr_discr_clone;

  /* set defaults                                                            */

  /* finite probability vector */
  DISTR.pv        = NULL;          /* probability vector (PV)                */
  DISTR.n_pv      = 0;             /* length of PV                           */

  /* probability mass function */
  DISTR.pmf       = NULL;          /* pointer to PMF                         */
  DISTR.cdf       = NULL;          /* pointer to CDF                         */
  DISTR.invcdf    = NULL;          /* pointer to inverse CDF                 */

  DISTR.init      = NULL;          /* pointer to special init routine        */

  DISTR.set_params= NULL;          /* funct for setting parameters and domain*/

  DISTR.n_params  = 0;             /* number of parameters of the pmf        */
  /* initialize parameters of the PMF                                        */
  for (i=0; i<UNUR_DISTR_MAXPARAMS; i++)
    DISTR.params[i] = 0.;

  DISTR.norm_constant = 1.;        /* (log of) normalization constant for PMF
				      (initialized to avoid accidently floating
				      point exception                        */

  DISTR.trunc[0] = DISTR.domain[0] = 0;         /* left boundary of domain   */
  DISTR.trunc[1] = DISTR.domain[1] = INT_MAX;   /* right boundary of domain  */

  DISTR.mode     = 0;              /* location of mode                       */
  DISTR.upd_mode = _unur_distr_discr_find_mode;  /* funct for computing mode */

  DISTR.sum     = 1.;              /* sum over PMF                           */
  DISTR.upd_sum = NULL;            /* funct for computing sum                */

  DISTR.pmftree    = NULL;         /* pointer to function tree for PMF       */
  DISTR.cdftree    = NULL;         /* pointer to function tree for CDF       */

  /* return pointer to object */
  return distr;

} /* end of unur_distr_discr_new() */

/*---------------------------------------------------------------------------*/

struct unur_distr *
_unur_distr_discr_clone( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* copy (clone) distribution object                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to source distribution object                    */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to clone of distribution object                            */
     /*----------------------------------------------------------------------*/
{
#define CLONE clone->data.discr

  struct unur_distr *clone;

  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  /* allocate memory */
  clone = _unur_xmalloc( sizeof(struct unur_distr) );
  
  /* copy distribution object into clone */
  memcpy( clone, distr, sizeof( struct unur_distr ) );

  /* copy function trees into generator object (when there is one) */
  CLONE.pmftree  = (DISTR.pmftree) ? _unur_fstr_dup_tree(DISTR.pmftree) : NULL;
  CLONE.cdftree  = (DISTR.cdftree) ? _unur_fstr_dup_tree(DISTR.cdftree) : NULL;

  /* copy probability vector into generator object (when there is one) */
  if (DISTR.pv) {
    CLONE.pv = _unur_xmalloc( DISTR.n_pv * sizeof(double) );
    memcpy( CLONE.pv, DISTR.pv, DISTR.n_pv * sizeof(double) );
  }

  /* copy user name for distribution */
  if (distr->name_str) {
    size_t len = strlen(distr->name_str) + 1;
    clone->name_str = _unur_xmalloc(len);
    memcpy( clone->name_str, distr->name_str, len );
    clone->name = clone->name_str;
  }

  return clone;

#undef CLONE
} /* end of _unur_distr_discr_clone() */

/*---------------------------------------------------------------------------*/

void
_unur_distr_discr_free( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* free distribution object                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  if( distr == NULL ) /* nothing to do */
    return;
  _unur_check_distr_object( distr, DISCR, RETURN_VOID );

  if (DISTR.pmftree)  _unur_fstr_free(DISTR.pmftree);
  if (DISTR.cdftree)  _unur_fstr_free(DISTR.cdftree);

  if (DISTR.pv) free( DISTR.pv );

  /* user name for distribution */
  if (distr->name_str) free(distr->name_str);

  COOKIE_CLEAR(distr);
  free( distr );

} /* end of unur_distr_discr_free() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pv( struct unur_distr *distr, const double *pv, int n_pv )
     /*----------------------------------------------------------------------*/
     /* set probability vector for distribution                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr   ... pointer to distribution object                         */
     /*   pv      ... pointer to PV                                          */
     /*   n_pv    ... length of PV                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* it is not possible to set a PV when a PMF is given. */
  if (DISTR.pmf != NULL || DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PMF/CDF given, cannot set PV");
    return UNUR_ERR_DISTR_SET;
  }

  /* check new parameter for distribution */
  if (n_pv < 0) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"length of PV");
    return UNUR_ERR_DISTR_SET;
  }

  /* n_pv must not be too large */
  if ( (DISTR.domain[0] > 0) && ((unsigned)DISTR.domain[0] + (unsigned)n_pv > INT_MAX) ) {
    /* n_pv too large, causes overflow */
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"length of PV too large, overflow");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.domain[1] = DISTR.domain[0] + n_pv - 1;

  /* we do not check non-negativity of p.v.
     (it is cheaper to do it when unur_init() is called */

  /* allocate memory for probability vector */
  DISTR.pv = _unur_xrealloc( DISTR.pv, n_pv * sizeof(double) );
  if (!DISTR.pv) return UNUR_ERR_MALLOC;

  /* copy probability vector */
  memcpy( DISTR.pv, pv, n_pv * sizeof(double) );
  DISTR.n_pv = n_pv;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_pv() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_make_pv( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* compute probability vector                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr   ... pointer to distribution object                         */
     /*                                                                      */
     /* return:                                                              */
     /*   length of probability vector                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  double *pv;          /* pointer to probability vector */
  int n_pv;            /* length of PV */
  double cdf, cdf_old; /* cumulated sum of PV */
  double thresh_cdf;   /* threshold for truncating PV */
  int valid;           /* whether cumputed PV is valid */
  int i;

  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  /* PMF or CDF required */
  if ( DISTR.pmf == NULL && DISTR.cdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_GET,"PMF or CDF");
    return 0;
  }

  /* if there exists a PV, it has to be removed */
  if (DISTR.pv != NULL) {
    free(DISTR.pv); DISTR.n_pv = 0;
  }

  /* compute PV */

  if ((unsigned)DISTR.domain[1] - (unsigned)DISTR.domain[0] < UNUR_MAX_AUTO_PV ) {

    /* first case: bounded domain */
    n_pv = DISTR.domain[1] - DISTR.domain[0] + 1;
    pv = _unur_xmalloc( n_pv * sizeof(double) );
    if (DISTR.pmf) {
      for (i=0; i<n_pv; i++)
	pv[i] = _unur_discr_PMF(DISTR.domain[0]+i,distr);
    }
    else if (DISTR.cdf) {
      cdf_old = 0.;
      for (i=0; i<n_pv; i++) {
	cdf = _unur_discr_CDF(DISTR.domain[0]+i,distr);
	pv[i] = cdf - cdf_old;
	cdf_old = cdf;
      }
    }
    valid = TRUE;
  }

  else {
    /* second case: domain too big but sum over PMF given       */
    /* we chop off the trailing part of the distribution        */
#define MALLOC_SIZE 1000 /* allocate 1000 doubles at once       */

    int n_alloc;         /* number of doubles allocated         */
    int max_alloc;       /* maximal number of allocated doubles */
    int size_alloc;      /* size of allocated blocks            */

    /* get maximal size of PV */
    if ( (DISTR.domain[0] <= 0) || (INT_MAX - DISTR.domain[0] >= UNUR_MAX_AUTO_PV - 1) ) {
      /* we can have a PV of length UNUR_MAX_AUTO_PV */
      size_alloc = MALLOC_SIZE;
      max_alloc = UNUR_MAX_AUTO_PV;
    }
    else { /* length of PV must be shorter than UNUR_MAX_AUTO_PV */
      size_alloc = max_alloc = INT_MAX - DISTR.domain[0];
    }

    /* init counter */
    n_pv = 0;
    pv = NULL;
    valid = FALSE;  /* created PV is empty yet and not valid */
    cdf = 0.;       /* cumulated sum of PV                   */
    cdf_old = 0.;   /* cumulated sum of PV in last iteration */
    /* threshold for truncating PV */
    thresh_cdf = (distr->set & UNUR_DISTR_SET_PMFSUM) ? (1.-1.e-8)*DISTR.sum : INFINITY;

    /* compute PV */
    for (n_alloc = size_alloc; n_alloc <= max_alloc; n_alloc += size_alloc) {
      pv = _unur_xrealloc( pv, n_alloc * sizeof(double) );

      if (DISTR.pmf) {
	for (i=0; i<size_alloc; i++) {
	  cdf += pv[n_pv] = _unur_discr_PMF(DISTR.domain[0]+n_pv,distr);
	  n_pv++;
	  if (cdf > thresh_cdf) { valid = TRUE; break; }
	}
      }
      else if (DISTR.cdf) {
	for (i=0; i<size_alloc; i++) {
	  cdf = _unur_discr_CDF(DISTR.domain[0]+n_pv,distr);
	  pv[n_pv] = cdf - cdf_old;
	  cdf_old = cdf;
	  n_pv++;
	  if (cdf > thresh_cdf) { valid = TRUE; break; }
	}
      }	  
      if (cdf > thresh_cdf) break;
    }

    if (distr->set & UNUR_DISTR_SET_PMFSUM) {
      /* make a warning if computed PV might not be valid */
      if (valid != TRUE)
	/* not successful */
	_unur_warning(distr->name,UNUR_ERR_DISTR_GET,"PV truncated");
    }
    else { /* PMFSUM not known */
      /* assume we have the important part of distribution */
      valid = TRUE;
      DISTR.sum = cdf;
      distr->set |= UNUR_DISTR_SET_PMFSUM;
    }
    
#undef MALLOC_SIZE
  }

  /* store vector */
  DISTR.pv = pv;
  DISTR.n_pv = n_pv;
  DISTR.domain[1] = DISTR.domain[0] + n_pv - 1;

  /* o.k. */
  return (valid) ? n_pv : -n_pv;
} /* end of unur_distr_discr_make_pv() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_get_pv( const struct unur_distr *distr, const double **pv )
     /*----------------------------------------------------------------------*/
     /* get length of probability vector and set pointer to probability      */
     /* vector                                                               */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   pv       ... pointer to probability vector                         */
     /*                                                                      */
     /* return:                                                              */
     /*   length of probability vector                                       */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  *pv = (DISTR.pv) ? DISTR.pv : NULL;
  return DISTR.n_pv;

} /* end of unur_distr_discr_get_pv() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_pv( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* returns the value of the probability vector at k or, if there is no  */
     /* probability vector defined, evaluates  the pmf                       */
     /*                                                                      */
     /* parampeters:                                                         */
     /*  k     ... argument for probability vector of pmf                    */
     /*  distr ... pointer to distribution object                            */
     /*                                                                      */
     /* return:                                                              */
     /*   pv[k] or pmf(k)                                                    */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  if (DISTR.pv != NULL) {
    /* use probability vector */
    if (k < DISTR.domain[0] || k > DISTR.domain[1])
      return 0.;
    else
      return (DISTR.pv[k-DISTR.domain[0]]);
  }

  if (DISTR.pmf != NULL) {
    /* use PMF */
    double px = _unur_discr_PMF(k,distr);
    if (_unur_isnan(px)) {
      _unur_warning(distr->name,UNUR_ERR_DISTR_DATA,"PMF returns NaN");
      return 0.;
    }
    else
      return px;
  }

  /* else: data missing */
  _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
  return INFINITY;

} /* end of unur_distr_discr_eval_pv() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmf( struct unur_distr *distr, UNUR_FUNCT_DISCR *pmf )
     /*----------------------------------------------------------------------*/
     /* set PMF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   pmf   ... pointer to PMF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr,  UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, pmf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* it is not possible to set both a PMF and a PV */
  if (DISTR.pv != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"delete exisiting PV");
    free(DISTR.pv); DISTR.n_pv = 0;
  }

  /* we do not allow overwriting a PMF */
  if (DISTR.pmf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PMF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  DISTR.pmf = pmf;
  return UNUR_SUCCESS;

} /* end of unur_distr_discr_set_pmf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_cdf( struct unur_distr *distr, UNUR_FUNCT_DISCR *cdf )
     /*----------------------------------------------------------------------*/
     /* set CDF of distribution                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   cdf   ... pointer to CDF                                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, cdf, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  
  /* it is not possible to set both a CDF and a PV */
  if (DISTR.pv != NULL) {
    _unur_warning(distr->name,UNUR_ERR_DISTR_SET,"delete exisiting PV");
    free(DISTR.pv); DISTR.n_pv = 0;
  }

  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, sum, etc. might be wrong now! */

  DISTR.cdf = cdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_invcdf( struct unur_distr *distr, UNUR_IFUNCT_DISCR *invcdf )
     /*----------------------------------------------------------------------*/
     /* set inverse CDF of distribution                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   invcdf ... pointer to inverse CDF                                  */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_NULL( distr->name, invcdf,UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  
  /* we do not allow overwriting an inverse cdf */
  if (DISTR.invcdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of inverse CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_INVALID;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  DISTR.invcdf = invcdf;
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_invcdf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_pmf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to PMF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to PMF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.pmf;
} /* end of unur_distr_discr_get_pmf() */

/*---------------------------------------------------------------------------*/

UNUR_FUNCT_DISCR *
unur_distr_discr_get_cdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to CDF of distribution                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to CDF                                                     */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.cdf;
} /* end of unur_distr_discr_get_cdf() */

/*---------------------------------------------------------------------------*/

UNUR_IFUNCT_DISCR *
unur_distr_discr_get_invcdf( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get pointer to inverse CDF of distribution                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to inverse CDF                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );

  return DISTR.invcdf;
} /* end of unur_distr_discr_get_invcdf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_pmf( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate PMF of distribution at k                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for pmf                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   pmf(k)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  if (DISTR.pmf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_discr_PMF(k,distr);
} /* end of unur_distr_discr_eval_pmf() */

/*---------------------------------------------------------------------------*/

double
unur_distr_discr_eval_cdf( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate CDF of distribution at k                                    */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF(k)                                                             */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  if (DISTR.cdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INFINITY;
  }

  return _unur_discr_CDF(k,distr);
} /* end of unur_distr_discr_eval_cdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_eval_invcdf( double u, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate inverse CDF of distribution at u                            */
     /*                                                                      */
     /* parameters:                                                          */
     /*   u     ... argument for inverse CDF                                 */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   invcdf(u)                                                          */
     /*   INT_MAX    ... on error                                            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INT_MAX );
  _unur_check_distr_object( distr, DISCR, INT_MAX );

  if (DISTR.invcdf == NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return INT_MAX;
  }

  if (u<=0.)
    return DISTR.domain[0];
  if (u>=1.)
    return DISTR.domain[1];
  else
    return _unur_discr_invCDF(u,distr);

} /* end of unur_distr_discr_eval_invcdf() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfstr( struct unur_distr *distr, const char *pmfstr )
     /*----------------------------------------------------------------------*/
     /* set PMF of distribution via a string interface                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   pmfstr ... string that describes function term of PMF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, pmfstr, UNUR_ERR_NULL );

  /* it is not possible to set a PMF when a PV is given. */
  if (DISTR.pv != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"PV given, cannot set PMF");
    return UNUR_ERR_DISTR_SET;
  }

  /* we do not allow overwriting a PMF */
  if (DISTR.pmf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of PMF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_DATA;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse PMF string */
  if ( (DISTR.pmftree = _unur_fstr2tree(pmfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }
  DISTR.pmf  = _unur_distr_discr_eval_pmf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_pmfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_cdfstr( struct unur_distr *distr, const char *cdfstr )
     /*----------------------------------------------------------------------*/
     /* set CDF of distribution via a string interface                       */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*   cdfstr ... string that describes function term of CDF              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  _unur_check_NULL( NULL, cdfstr, UNUR_ERR_NULL );

  /* we do not allow overwriting a CDF */
  if (DISTR.cdf != NULL) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Overwriting of CDF not allowed");
    return UNUR_ERR_DISTR_SET;
  }

  /* for derived distributions (e.g. order statistics) not possible */
  if (distr->base) return UNUR_ERR_DISTR_DATA;

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* parse string */
  if ( (DISTR.cdftree = _unur_fstr2tree(cdfstr)) == NULL ) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"Syntax error in function string");
    return UNUR_ERR_DISTR_SET;
  }

  /* set evaluation function */
  DISTR.cdf  = _unur_distr_discr_eval_cdf_tree;

  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_cdfstr() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_discr_eval_pmf_tree( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for PMF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for PMF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   PMF at k                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  return ((DISTR.pmftree) ? _unur_fstr_eval_tree(DISTR.pmftree,(double)k) : 0.);
} /* end of _unur_distr_discr_eval_pmf_tree() */

/*---------------------------------------------------------------------------*/

double
_unur_distr_discr_eval_cdf_tree( int k, const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* evaluate function tree for CDF.                                      */
     /*                                                                      */
     /* parameters:                                                          */
     /*   k     ... argument for CDF                                         */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   CDF at k                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  return ((DISTR.cdftree) ? _unur_fstr_eval_tree(DISTR.cdftree,(double)k) : 0.);
} /* end of _unur_distr_discr_eval_cdf_tree() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_discr_get_pmfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get PMF string that is given via the string interface                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  _unur_check_NULL( NULL, DISTR.pmftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.pmftree,"x","PMF",TRUE);
} /* end of unur_distr_discr_get_pmfstr() */

/*---------------------------------------------------------------------------*/

char *
unur_distr_discr_get_cdfstr( const struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get CDF string that is given via the string interface                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr  ... pointer to distribution object                          */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to resulting string.                                       */
     /*                                                                      */
     /* comment:                                                             */
     /*   This string should be freed when it is not used any more.          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, NULL );
  _unur_check_distr_object( distr, DISCR, NULL );
  _unur_check_NULL( NULL, DISTR.cdftree, NULL );

  /* make and return string */
  return _unur_fstr_tree2string(DISTR.cdftree,"x","CDF",TRUE);
} /* end of unur_distr_discr_get_cdfstr() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfparams( struct unur_distr *distr, const double *params, int n_params )
     /*----------------------------------------------------------------------*/
     /* set array of parameters for distribution                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... list of arguments                                     */
     /*   n_params ... number of arguments                                   */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );
  if (n_params>0) _unur_check_NULL(distr->name, params, UNUR_ERR_NULL);

  /* first check number of new parameter for the distribution */
  if (n_params < 0 || n_params > UNUR_DISTR_MAXPARAMS ) {
    _unur_error(NULL,UNUR_ERR_DISTR_NPARAMS,"");
    return UNUR_ERR_DISTR_NPARAMS;
  }

  /* changelog */
  distr->set &= ~UNUR_DISTR_SET_MASK_DERIVED;
  /* derived parameters like mode, area, etc. might be wrong now! */

  /* even if the set routine fails, the derived parameters are
     marked as unknown. but this is o.k. since in this case something
     has been wrong. */

  /* use special routine for setting parameters
     (if there is one) */

  if (DISTR.set_params)
    return (DISTR.set_params(distr,params,n_params));

  /* otherwise simply copy parameters */

  DISTR.n_params = n_params;
  if (n_params) memcpy( DISTR.params, params, n_params*sizeof(double) );

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_pmfparams( const struct unur_distr *distr, const double **params )
     /*----------------------------------------------------------------------*/
     /* get number of pmf parameters and sets pointer to array params[] of   */
     /* parameters                                                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr    ... pointer to distribution object                        */
     /*   params   ... pointer to list of arguments                          */
     /*                                                                      */
     /* return:                                                              */
     /*   number of pmf parameters                                           */
     /*                                                                      */
     /* error:                                                               */
     /*   return 0                                                           */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, 0 );
  _unur_check_distr_object( distr, DISCR, 0 );

  *params = (DISTR.n_params) ? DISTR.params : NULL;
  return DISTR.n_params;

} /* end of unur_distr_discr_get_pmfparams() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_domain( struct unur_distr *distr, int left, int right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity            */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (left >= right) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"domain, left >= right");
    return UNUR_ERR_DISTR_SET;
  }

  /* store data */
  DISTR.trunc[0] = DISTR.domain[0] = left;
  DISTR.trunc[1] = DISTR.domain[1] = (DISTR.pv == NULL) ? right : left+DISTR.n_pv-1;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_DOMAIN;

  /* if distr is an object for a standard distribution, this   */
  /* not the original domain of it. (not a "standard domain")  */
  /* However, since we have changed the domain, we assume      */
  /* that this is not a truncated distribution.                */
  /* At last we have to mark all derived parameters as unknown */
  distr->set &= ~(UNUR_DISTR_SET_STDDOMAIN |
		  UNUR_DISTR_SET_TRUNCATED | 
		  UNUR_DISTR_SET_MASK_DERIVED );

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_discr_set_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_domain( const struct unur_distr *distr, int *left, int *right )
     /*----------------------------------------------------------------------*/
     /* set the left and right borders of the domain of the distribution     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   left  ... left boundary point                                      */
     /*   right ... right boundary point                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* comment:                                                             */
     /*   INT_MIN and INT_MAX are interpreted as (minus) infinity            */
     /*   if no boundaries have been set [INT_MIN, INT_MAX] is returned.     */
     /*----------------------------------------------------------------------*/
{
  /* in case of error the boundaries are set to +/- INFINITY */
  *left = INT_MIN;
  *right = INT_MAX;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* o.k. */
  *left  = DISTR.domain[0];
  *right = DISTR.domain[1];

  return UNUR_SUCCESS;
} /* end of unur_distr_discr_get_domain() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_mode( struct unur_distr *distr, int mode )
     /*----------------------------------------------------------------------*/
     /* set mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   mode  ... mode of PMF                                              */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  DISTR.mode = mode;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_MODE;

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of unur_distr_discr_set_mode() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_upd_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute mode of distribution (if possible)                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  if (DISTR.upd_mode == NULL) {
    /* no function to compute mode available */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

  /* compute mode */
  if ((DISTR.upd_mode)(distr)==UNUR_SUCCESS) {
    /* changelog */
    distr->set |= UNUR_DISTR_SET_MODE;
    return UNUR_SUCCESS;
  }
  else {
    /* computing of mode failed */
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"");
    return UNUR_ERR_DISTR_DATA;
  }

} /* end of unur_distr_discr_upd_mode() */
  
/*---------------------------------------------------------------------------*/

int
unur_distr_discr_get_mode( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get mode of distribution                                             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   mode of distribution                                               */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INT_MAX );
  _unur_check_distr_object( distr, DISCR, INT_MAX );

  /* mode known ? */
  if ( !(distr->set & UNUR_DISTR_SET_MODE) ) {
    /* try to compute mode */
    if (DISTR.upd_mode == NULL) {
      /* no function to compute mode available */
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
      return INT_MAX;
    }
    else {
      /* compute mode */
      if (unur_distr_discr_upd_mode(distr)!=UNUR_SUCCESS) {
	/* finding mode not successfully */
	_unur_error(distr->name,UNUR_ERR_DISTR_GET,"mode");
	return INT_MAX;
      }
    }
  }

  return DISTR.mode;

} /* end of unur_distr_discr_get_mode() */

/*---------------------------------------------------------------------------*/

int
unur_distr_discr_set_pmfsum( struct unur_distr *distr, double sum )
     /*----------------------------------------------------------------------*/
     /* set sum over PMF                                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   sum   ... sum over PMF                                             */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* check new parameter for distribution */
  if (sum <= 0.) {
    _unur_error(distr->name,UNUR_ERR_DISTR_SET,"pmf sum <= 0");
    return UNUR_ERR_DISTR_SET;
  }

  DISTR.sum = sum;

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PMFSUM;

  /* o.k. */
  return UNUR_SUCCESS;

} /* end of unur_distr_discr_set_pmfsum() */

/*---------------------------------------------------------------------------*/

int 
unur_distr_discr_upd_pmfsum( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* (re-) compute sum over PMF of distribution (if possible)             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
  double sum = 0.;
  int k, left, right, length;

  /* check arguments */
  _unur_check_NULL( NULL, distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_SET );

  /* changelog */
  distr->set |= UNUR_DISTR_SET_PMFSUM;

  if (DISTR.upd_sum != NULL) {
    /* try function given by distribution */
    if ((DISTR.upd_sum)(distr)==UNUR_SUCCESS)
      return UNUR_SUCCESS;
  }

  /* no function to compute sum available */
  left  = DISTR.domain[0];
  right = DISTR.domain[1];
  length = right - left;
  /* remark: length < 0 if right-left overflows */


  if (DISTR.cdf != NULL) {
    /* use CDF */
    if (left > INT_MIN) left -= 1;
    DISTR.sum = _unur_discr_CDF(right,distr) - _unur_discr_CDF(left,distr);
    return UNUR_SUCCESS;
  }

  if (DISTR.pv != NULL) {
    for (k = 0; k<= length; k++)
      /* use probability vector */
      sum += DISTR.pv[k];
    DISTR.sum = sum;
    return UNUR_SUCCESS;
  }

  if (DISTR.pmf != NULL && length > 0 && length <= MAX_PMF_DOMAIN_FOR_UPD_PMFSUM) {
    /* use PMF */
    for (k = left; k<= right; k++)
      sum += _unur_discr_PMF(k,distr);
    DISTR.sum = sum;
    return UNUR_SUCCESS;
  }

  /* error: data missing */
  distr->set &= ~UNUR_DISTR_SET_PMFSUM;
  _unur_error(distr->name,UNUR_ERR_DISTR_DATA,"Cannot compute sum");
  return UNUR_ERR_DISTR_DATA;

} /* end of unur_distr_discr_upd_pmfsum() */
  
/*---------------------------------------------------------------------------*/

double
unur_distr_discr_get_pmfsum( struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* get sum over PMF of distribution                                     */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*                                                                      */
     /* return:                                                              */
     /*   sum over PMF of distribution                                       */
     /*----------------------------------------------------------------------*/
{
  /* check arguments */
  _unur_check_NULL( NULL, distr, INFINITY );
  _unur_check_distr_object( distr, DISCR, INFINITY );

  /* sum known ? */
  if ( !(distr->set & UNUR_DISTR_SET_PMFSUM) ) {
    /* try to compute sum */
    if ( unur_distr_discr_upd_pmfsum(distr) != UNUR_SUCCESS ) {
      _unur_error(distr->name,UNUR_ERR_DISTR_GET,"sum");
      return INFINITY;
    }
  }

  return DISTR.sum;

} /* end of unur_distr_discr_get_pmfsum() */

/*****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_LOGGING
/*---------------------------------------------------------------------------*/

void
_unur_distr_discr_debug( const struct unur_distr *distr, const char *genid, unsigned printvector )
     /*----------------------------------------------------------------------*/
     /* write info about distribution into LOG file                          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   distr ... pointer to distribution object                           */
     /*   genid ... pointer to generator id                                  */
     /*   printvector ... whether the probability vector (if given)          */
     /*----------------------------------------------------------------------*/
{
  FILE *LOG;
  int i;

  /* check arguments */
  CHECK_NULL(distr,RETURN_VOID);
  COOKIE_CHECK(distr,CK_DISTR_DISCR,RETURN_VOID);

  LOG = unur_get_stream();

  fprintf(LOG,"%s: distribution:\n",genid);
  fprintf(LOG,"%s:\ttype = discrete univariate distribution\n",genid);
  fprintf(LOG,"%s:\tname = %s\n",genid,distr->name);

  if ( DISTR.pmf ) {
    /* have probability mass function */
    fprintf(LOG,"%s:\tPMF with %d argument(s)\n",genid,DISTR.n_params);
    for( i=0; i<DISTR.n_params; i++ )
      fprintf(LOG,"%s:\t\tparam[%d] = %g\n",genid,i,DISTR.params[i]);
  }

  if (DISTR.n_pv>0) {
    /* have probability vector */
    fprintf(LOG,"%s:\tprobability vector of length %d",genid,DISTR.n_pv);
    if (printvector) {
      for (i=0; i<DISTR.n_pv; i++) {
	if (i%10 == 0)
	  fprintf(LOG,"\n%s:\t",genid);
	fprintf(LOG,"  %.5f",DISTR.pv[i]);
      }
    }
    fprintf(LOG,"\n%s:\n",genid);
  }

  /* domain */
  if ( DISTR.pmf ) {
    /* have probability mass function */
    fprintf(LOG,"%s:\tdomain for pmf = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[1]);
    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
    fprintf(LOG,"\n%s:\n",genid);
  }

  if (DISTR.n_pv>0) {
    /* have probability vector */
    fprintf(LOG,"%s:\tdomain for pv = (%d, %d)",genid,DISTR.domain[0],DISTR.domain[0]-1+DISTR.n_pv);
    _unur_print_if_default(distr,UNUR_DISTR_SET_DOMAIN);
    fprintf(LOG,"\n%s:\n",genid);
  }

  if (distr->set & UNUR_DISTR_SET_MODE)
    fprintf(LOG,"%s:\tmode = %d\n",genid,DISTR.mode);
  else
    fprintf(LOG,"%s:\tmode unknown\n",genid);
  
  if (distr->set & UNUR_DISTR_SET_PMFSUM)
    fprintf(LOG,"%s:\tsum over PMF = %g\n",genid,DISTR.sum);
  else
    fprintf(LOG,"%s:\tsum over PMF unknown\n",genid);
  fprintf(LOG,"%s:\n",genid);
  
} /* end of _unur_distr_discr_debug() */

/*---------------------------------------------------------------------------*/
#endif    /* end UNUR_ENABLE_LOGGING */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int 
_unur_distr_discr_find_mode(struct unur_distr *distr )
     /*----------------------------------------------------------------------*/
     /* find mode of a probability vector by bisection.                      */
     /*                                                                      */
     /* Assumptions: PMF is unimodal, no "inflexion" ("saddle") points.      */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*----------------------------------------------------------------------*/
{
#define N_TRIALS  (100)

  int x[3];                       /* bracket for mode x[0] < x[1] < x[2]     */
  double fx[3];                   /* PMF values at bracket points            */
  int xnew;                       /* new point in iteration                  */
  double fxnew;                   /* PMF value at new point                  */
  int step;                       /* auxiliary variable for searching point  */
  int this, other;                /* side of bracket under consideration     */
  int cutthis;                    /* which side of bracket should be cut     */

  const double r = (sqrt(5.)-1.)/2.;     /* sectio aurea                     */

  /* check arguments */
  CHECK_NULL( distr, UNUR_ERR_NULL );
  _unur_check_distr_object( distr, DISCR, UNUR_ERR_DISTR_INVALID );

  /* left and right boundary point of bracket */
  x[0] = DISTR.domain[0];
  x[2] = DISTR.domain[1];
  fx[0] = unur_distr_discr_eval_pv(x[0], distr);
  fx[2] = unur_distr_discr_eval_pv(x[2], distr);

  /* we assume for our algorithm that there are at least three points */
  if (x[2] <= x[0] + 1) {
    /* we have one or two points only */
    DISTR.mode = (fx[0] <= fx[2]) ? x[2] : x[0];
    distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
    return UNUR_SUCCESS;
  }

  /* --- middle point of bracket ------------------------------------------- */
  /* we have to find a middle point where the PMF is non-zero */

  /* first trial: mean of boundary points */
  x[1]  = (x[0]/2) + (x[2]/2);
  if (x[1]<=x[0]) x[1]++;
  if (x[1]>=x[2]) x[1]--;
  fx[1] = unur_distr_discr_eval_pv(x[1], distr); 

  /* second trial: start search from left boundary */
  if ( !(fx[1]>0.)) {
    xnew = (DISTR.domain[0]!=INT_MIN) ? DISTR.domain[0] : 0;
    for (step = 1; step < N_TRIALS; step++) {
      xnew += step;
      if (xnew >= DISTR.domain[1]) break;
      if ((fxnew = unur_distr_discr_eval_pv(xnew,distr)) > 0.) {
	x[1] = xnew; fx[1] = fxnew; break;
      }
    }
  }

  /* third trial: start search from 0 */
  if ( !(fx[1]>0.) && DISTR.domain[0]!=0) {
    xnew = 0;
    for (step = 1; step < N_TRIALS; step++) {
      xnew += step;
      if (xnew >= DISTR.domain[1]) break;
      if ((fxnew = unur_distr_discr_eval_pv(xnew,distr)) > 0.) {
	x[1] = xnew; fx[1] = fxnew; break;
      }
    }
  }

  /* forth trial: start search from right boundary */
  if ( !(fx[1]>0.) && DISTR.domain[1]!=INT_MAX) {
    xnew = DISTR.domain[1];
    for (step = 1; step < N_TRIALS; step++) {
      xnew -= step;
      if (xnew <= DISTR.domain[0]) break;
      if ((fxnew = unur_distr_discr_eval_pv(xnew,distr)) > 0.) {
	x[1] = xnew; fx[1] = fxnew; break;
      }
    }
  }

  if ( !(fx[1]>0.)) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA,
		"find_mode(): no positive entry in PV found");
    return UNUR_ERR_DISTR_DATA;
  }
  if (fx[1]<fx[0] && fx[1]<fx[2]) {
    _unur_error(distr->name,UNUR_ERR_DISTR_DATA, "find_mode(): PV not unimodal");
    return UNUR_ERR_DISTR_DATA;
  }

  /* --- middle point of bracket found ------------------------------------- */

  /* --- interate until maximum is found ----------------------------------- */

  while (1) {

    /* fprintf(stderr,"x = %d, %d, %d\n",x[0],x[1],x[2]); */
    /* fprintf(stderr,"fx = %g, %g, %g\n",fx[0],fx[1],fx[2]); */

    /* terminate */
    if (x[0]+1 >= x[1] && x[1] >= x[2]-1) {
      DISTR.mode = (fx[0]>fx[2]) ? x[0] : x[2];
      if (fx[1]>DISTR.mode) DISTR.mode = x[1];
      distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
      return UNUR_SUCCESS;
    } 

    /* new point */
    xnew  = (int) (r*x[0] + (1.-r)*x[2]);
    if (xnew == x[0])  ++xnew;
    if (xnew == x[2])  --xnew;
    if (xnew == x[1])  xnew += (x[1]-1==x[0]) ? 1 : -1;

    /* side of bracket */
    if (xnew < x[1]) {
      this = 0; other = 2; } /* l.h.s. of bracket */
    else {
      this = 2; other = 0; } /* r.h.s. of bracket */

    /* value at new point */
    fxnew = unur_distr_discr_eval_pv(xnew,distr);
    if ( fxnew < fx[0] && fxnew < fx[2] ) {
      _unur_error(distr->name,UNUR_ERR_DISTR_DATA, "find_mode(): PV not unimodal");
      return UNUR_ERR_DISTR_DATA;
    }
    
    do {

      if (!_unur_FP_same(fxnew,fx[1])) {
	cutthis = (fxnew > fx[1]) ? FALSE : TRUE;
	break;
      }

      /* else: fxnew == fx[1] */

      if (fx[this]  > fx[1]) { cutthis = FALSE; break; }
      if (fx[other] > fx[1]) { cutthis = TRUE;  break; }

      /* else: fx[0] < fxnew && fx[1] == fxnew && fx[2] < fxnew */

      for (step = 1; step < N_TRIALS && xnew >= x[0]  && xnew <= x[2]; step++) {
	xnew += (this==0) ? -1 : 1;
	fxnew = unur_distr_discr_eval_pv(xnew,distr);
	if (_unur_FP_less(fxnew,fx[1])) {
	  DISTR.mode = x[1];
	  distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ; 
	  return UNUR_SUCCESS;
	}
      }

      _unur_error(distr->name,UNUR_ERR_DISTR_DATA, "find_mode(): PV not unimodal");
      return UNUR_ERR_DISTR_DATA;
      
    } while (0);

    if (cutthis) {
      x[this] = xnew; fx[this] = fxnew;
    }
    else {
      x[other] = x[1]; fx[other] = fx[1];
      x[1] = xnew; fx[1] = fxnew;
    }

  }   /* --- end while(1) --- */

  /* changelog */
  /*   distr->set |= UNUR_DISTR_SET_MODE | UNUR_DISTR_SET_MODE_APPROX ;  */

  /* o.k. */
  /*   return UNUR_SUCCESS; */

} /* end of _unur_distr_discr_find_mode() */

/*---------------------------------------------------------------------------*/
#undef DISTR
/*---------------------------------------------------------------------------*/
