/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mrou_rectangle.c                                             *
 *                                                                           *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      The bounding rectangle for the multivariate RoU-methods is computed  *
 *      numerically.                                                         *
 *      Only power transformations with parameter r are considered.          *
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
 *   Notation follows [2].                                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   REFERENCES:                                                             *
 *   [1] Wakefield J.C., Gelfand A.E., Smith A.F.M.                          *
 *       Efficient generation of random variates via the ratio-of-uniforms   *
 *       method.                                                             *
 *       Statistics and Computing (1991) 1, pp (129-133)                     *
 *                                                                           *
 *   [2] Hoermann, W., Leydold J., and Derflinger, G. (2004):                *
 *       Automatic non-uniform random variate generation, Springer, Berlin.  *
 *       Section 2.4, Algorithm 2.9 (RoU), p.35                              *
 *                                                                           *
 *   [3] Hooke, R. and Jeeves, T.A. (1961):                                  *
 *       Direct Search Solution of Numerical and Statistical Problems.       *
 *       Journal of the ACM, Vol. 8, April 1961, pp. 212-229.                *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/cvec.h>
#include <utils/fmax_source.h>
#include <utils/hooke_source.h>
#include <utils/matrix_source.h>
#include <utils/unur_fp_source.h>
#include <utils/mrou_rectangle_struct.h>
#include <utils/mrou_rectangle_source.h>

/*---------------------------------------------------------------------------*/
/* Constants                                                                 */

/* Convergence parameters for the hooke optimization algorithm */

#define MROU_HOOKE_RHO     (0.5)
#define MROU_HOOKE_EPSILON (1.e-7)
#define MROU_HOOKE_MAXITER (1000L)

/* Scaling factor for the computed minimum bounding rectangle.               */
/* The computed rectangle  (0, vmax)x(umin[d], umax[d]) is scaled by this    */
/* factor, i.e. :                                                            */
/* vmax = vmax * ( 1+ MROU_RECT_SCALING)                                     */
/* umin[i] = umin[i] - (umax[i]-umin[i])*MROU_RECT_SCALING/2.                */
/* umax[i] = umax[i] + (umax[i]-umin[i])*MROU_RECT_SCALING/2.                */
#define MROU_RECT_SCALING (1.e-4)


static double _unur_mrou_rectangle_aux_vmax(double *x, void *p );
static double _unur_mrou_rectangle_aux_umin(double *x, void *p );
static double _unur_mrou_rectangle_aux_umax(double *x, void *p );
/*---------------------------------------------------------------------------*/
/* Auxiliary functions used in the computation of the bounding rectangle     */
/*---------------------------------------------------------------------------*/

#define PDF(x)    _unur_cvec_PDF((x),(distr))    /* call to PDF              */

/*---------------------------------------------------------------------------*/

double
_unur_mrou_rectangle_aux_vmax(double *x, void *p )
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  struct MROU_RECTANGLE *rr;
  rr = p; /* typecast from void* to unur_rou_rectangle* */

  return -pow( _unur_cvec_PDF((x),(rr->distr)) ,
	       1./(1.+ rr->r * rr->dim) );
}

/*---------------------------------------------------------------------------*/

double
_unur_mrou_rectangle_aux_umin(double *x, void *p)
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  struct MROU_RECTANGLE *rr;
  rr = p; /* typecast from void* to unur_rou_rectangle* */

  return ( (x[rr->aux_dim] - rr->center[rr->aux_dim])
	   * pow( _unur_cvec_PDF((x),(rr->distr)),
		  rr->r / (1.+ rr->r * rr->dim) ) );
}

/*---------------------------------------------------------------------------*/

double
_unur_mrou_rectangle_aux_umax(double *x, void *p)
     /*----------------------------------------------------------------------*/
     /* Auxiliary function used in the computation of the bounding rectangle */
     /*----------------------------------------------------------------------*/
{
  return (- _unur_mrou_rectangle_aux_umin(x,p)) ;
}

/*---------------------------------------------------------------------------*/

struct MROU_RECTANGLE *
_unur_mrou_rectangle_new( void )
     /*----------------------------------------------------------------------*/
     /* create empty MROU rectangle object.                                  */
     /*                                                                      */
     /* parameters: none                                                     */
     /*                                                                      */
     /* return:                                                              */
     /*   pointer to allocated and initialized MROU_RECTANGLE object         */
     /*----------------------------------------------------------------------*/
{
  struct MROU_RECTANGLE *rr;

  rr = _unur_xmalloc(sizeof(struct MROU_RECTANGLE ));

  rr->distr  = NULL;
  rr->dim    = 0;
  rr->umin   = NULL;
  rr->umax   = NULL;
  rr->r      = 1;
  rr->bounding_rectangle = 1;
  rr->center = NULL;
  rr->genid  = "";
  
  return rr;
} /* end of _unur_mrou_rectangle_new() */

/*---------------------------------------------------------------------------*/

int
_unur_mrou_rectangle_compute( struct MROU_RECTANGLE *rr )
     /*----------------------------------------------------------------------*/
     /* compute universal bounding hyper-rectangle                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   rr ... mRoU object                                                 */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS ... on success                                        */
     /*   error code   ... on error                                          */
     /*                                                                      */
     /* remark:                                                              */
     /*   This routine computes two parameters of the bounding rectangle:    */
     /*      vmax  .... the upper bound of the v-coordinate                  */
     /*      umin, umax ... left lower and right upper vertex of rectangle   */
     /*                     in u-hyperplane                                  */
     /*                                                                      */
     /*   vmax is computed by means of the mode of the given distribution.   */
     /*   If this is not available, a method by Hooke and Jeeves [3] is      */
     /*   used.                                                              */
     /*                                                                      */
     /*   umin and umax is computed if flag rr->bounding_rectangle is        */
     /*   nonzero. Again the algorithm by Hooke and Jeeves [3] is used.      */
     /*   In this case rr->umin and rr->umax must contain points to          */
     /*   to double arrays each of size (at least) (rr->dim)+1.              */
     /*   Otherwise, these two pointers are not needed and can be NULL.      */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_funct_vgeneric faux; /* function to be minimized/maximized    */
  double *xstart, *xend, *xumin, *xumax; /* coordinate arrays used in maximum/minimum calculations */
  int d, dim;             /* index used in dimension loops (0 <= d < dim) */
  int hooke_iters_vmax;   /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umin;   /* actual number of min/max iterations = return value of hooke()*/
  int hooke_iters_umax;   /* actual number of min/max iterations = return value of hooke()*/
  double scaled_epsilon;  /* to be used in the hooke algorithm */
  int flag_finite = TRUE; /* sanity flag FALSE when some of the rectangle (u,v)-values are not finite */

  /* dimension of the distribution */
  dim = rr->dim;

  /* allocate memory for the coordinate vectors */
  xstart = _unur_xmalloc(dim * sizeof(double));
  xend   = _unur_xmalloc(dim * sizeof(double));
  xumin  = _unur_xmalloc(dim * sizeof(double));
  xumax  = _unur_xmalloc(dim * sizeof(double));

  /* --- compute vmax --- */

  if ( (rr->distr->set & UNUR_DISTR_SET_MODE) && (rr->distr->data.cvec.mode != NULL)) { 
    /* position of mode is known ... vmax = f(mode)^(1/r*dim+1)) */
    faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_vmax;
    faux.params = rr;

    rr->vmax = -faux.f(rr->distr->data.cvec.mode, faux.params);
  }
  else {
    /* calculation of vmax */
    faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_vmax;
    faux.params = rr;
    
    /* starting point */
    memcpy(xstart, rr->center, dim * sizeof(double));
    
    hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend,
				    MROU_HOOKE_RHO, MROU_HOOKE_EPSILON, MROU_HOOKE_MAXITER);
    
    rr->vmax = -faux.f(xend, faux.params);
    
    if (hooke_iters_vmax >= MROU_HOOKE_MAXITER) {
      scaled_epsilon = MROU_HOOKE_EPSILON * rr->vmax;
      if (scaled_epsilon>MROU_HOOKE_EPSILON) scaled_epsilon=MROU_HOOKE_EPSILON;
      
      /* recalculating extremum with scaled_epsilon and new starting point */
      memcpy(xstart, xend, dim * sizeof(double));
      hooke_iters_vmax = _unur_hooke( faux, dim, xstart, xend,
                                      MROU_HOOKE_RHO, scaled_epsilon , MROU_HOOKE_MAXITER);
      rr->vmax = -faux.f(xend, faux.params);
      if (hooke_iters_vmax >= MROU_HOOKE_MAXITER) {
        _unur_warning(rr->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (vmax)");
      }
    }

    /* additional scaling of boundary rectangle */
    rr->vmax = rr->vmax * ( 1+ MROU_RECT_SCALING);
  }

  /* check for finite results */
  flag_finite = _unur_isfinite(rr->vmax);

  /* --- compute umin and umax --- */
  
  if (rr->bounding_rectangle) {

    /* check pointers to avoid segfault */
    if (rr->umin == NULL || rr->umax == NULL) {
      free(xstart); free(xend); free(xumin); free(xumax);
      _unur_error(rr->genid,UNUR_ERR_NULL,"");
      return UNUR_ERR_NULL;
    }

    /* calculation of umin and umax */
    for (d=0; d<dim; d++) {
      
      /* setting coordinate dimension to be used by the auxiliary functions */
      rr->aux_dim  = d;
      
      /* starting point at center */
      memcpy(xstart, rr->center, dim * sizeof(double));
      
      /*-----------------------------------------------------------------------------*/
      /* calculation for umin */
      
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umin;
      faux.params = rr;
      
      hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend,
				      MROU_HOOKE_RHO, MROU_HOOKE_EPSILON, MROU_HOOKE_MAXITER);
      rr->umin[d] = faux.f(xend, faux.params);
      
      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumin, xend, dim * sizeof(double));
      
      /*-----------------------------------------------------------------------------*/
      /* and now, an analogue calculation for umax */
      
      faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umax;
      faux.params = rr;
      
      hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend,
				      MROU_HOOKE_RHO, MROU_HOOKE_EPSILON, MROU_HOOKE_MAXITER);
      rr->umax[d] = -faux.f(xend, faux.params);
      
      /* storing actual endpoint in case we need a recalculation */
      memcpy(xumax, xend, dim * sizeof(double));
      
      /*-----------------------------------------------------------------------------*/
      /* checking if we need to recalculate umin */
      if (hooke_iters_umin >= MROU_HOOKE_MAXITER) {      
	scaled_epsilon = MROU_HOOKE_EPSILON * (rr->umax[d]-rr->umin[d]);
	if (scaled_epsilon>MROU_HOOKE_EPSILON) scaled_epsilon=MROU_HOOKE_EPSILON;
	
	/* recalculating extremum with scaled_epsilon and new starting point */
	faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umin;
	faux.params = rr;
	
	memcpy(xstart, xumin, dim * sizeof(double));
	hooke_iters_umin = _unur_hooke( faux, dim, xstart, xend,
					MROU_HOOKE_RHO, scaled_epsilon , MROU_HOOKE_MAXITER);
	rr->umin[d] = faux.f(xend, faux.params);
	if (hooke_iters_umin >= MROU_HOOKE_MAXITER) {
	  _unur_warning(rr->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umin)");
	}
      }
      
      /* checking if we need to recalculate umax */
      if (hooke_iters_umax >= MROU_HOOKE_MAXITER) {
	scaled_epsilon = MROU_HOOKE_EPSILON * (rr->umax[d]-rr->umin[d]);
	if (scaled_epsilon>MROU_HOOKE_EPSILON) scaled_epsilon=MROU_HOOKE_EPSILON;
	
	/* recalculating extremum with scaled_epsilon and new starting point */
	faux.f = (UNUR_FUNCT_VGENERIC*) _unur_mrou_rectangle_aux_umax;
	faux.params = rr;
	
	memcpy(xstart, xumax, dim * sizeof(double));
	hooke_iters_umax = _unur_hooke( faux, dim, xstart, xend,
					MROU_HOOKE_RHO, scaled_epsilon , MROU_HOOKE_MAXITER);
	rr->umin[d] = faux.f(xend, faux.params);
	if (hooke_iters_umax >= MROU_HOOKE_MAXITER) {
	  _unur_warning(rr->genid , UNUR_ERR_GENERIC, "Bounding rect uncertain (umax)");
	}
      }
      
      /*-----------------------------------------------------------------------------*/
      /* additional scaling of boundary rectangle */
      rr->umin[d] = rr->umin[d] - (rr->umax[d]-rr->umin[d])*MROU_RECT_SCALING/2.;
      rr->umax[d] = rr->umax[d] + (rr->umax[d]-rr->umin[d])*MROU_RECT_SCALING/2.;
      
      /* check for finite results */
      flag_finite = flag_finite && _unur_isfinite(rr->umin[d]) && _unur_isfinite(rr->umax[d]);
    }
  }
  
  /* free working arrays */
  free(xstart); free(xend); free(xumin); free(xumax);

  if (rr->vmax <= 0.) {
    /* vmax must be strictly positive! */
    _unur_error("RoU",UNUR_ERR_DISTR_DATA,"cannot find bounding rectangle");
    return UNUR_ERR_DISTR_DATA;
  }
  
  /* return status of computation */
  return (flag_finite ? UNUR_SUCCESS : UNUR_ERR_INF);

} /* end of _unur_mrou_rectangle() */


#undef PDF
#undef MROU_HOOKE_RHO
#undef MROU_HOOKE_EPSILON
#undef MROU_HOOKE_MAXITER
#undef MROU_RECT_SCALING

/*---------------------------------------------------------------------------*/

