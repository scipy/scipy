/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      mvtdr_sample.c                                               *
 *                                                                           *
 *   TYPE:      continuous multivariate random variate                       *
 *   METHOD:    multivariate transformed density rejection                   *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Given (logarithm of the) PDF of a log-concave distribution;          *
 *      produce a value x consistent with its density.                       *
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

/*****************************************************************************/
/**  Private                                                                **/
/*****************************************************************************/

/*****************************************************************************/
/**  Sampling                                                               **/
/*****************************************************************************/

/*---------------------------------------------------------------------------*/

int
_unur_mvtdr_sample_cvec( struct unur_gen *gen, double *rpoint )
     /*----------------------------------------------------------------------*/
     /* sample from generator                                                */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   vec ... random vector (result)                                     */
     /*----------------------------------------------------------------------*/
{
  CONE *c;       /* cone for generating point */
  double gx;     /* distance of random point */
  double U;      /* uniformly distributed random number */
  double f, h;   /* value of density and hat at random point */
  int i,j;

  double *S = GEN->S;  /* working array for storing point on simples */

  /* check arguments */
  CHECK_NULL(gen,UNUR_ERR_NULL);
  COOKIE_CHECK(gen,CK_MVTDR_GEN,UNUR_ERR_COOKIE);

  /* loop until random point is accepted */
  while( 1 ) { 

    /*.......................................................................*/
    /** find a cone **/

    U = _unur_call_urng(gen->urng);      /* sample from uniform distribution */

    /* look up in guide table and search for cone */
    c = (GEN->guide)[(int) (U * GEN->guide_size)]; 
    U *= GEN->Htot;
    while (c->next!=NULL && c->Hsum < U) 
      c = c->next;

    /*.......................................................................*/
    /** get random point and distance of hyper plane **/

    /* get x value for marginal distribution of hat --> hyperplane */
    if (GEN->has_domain)
      /* update domain of gamma distribution */
      /* remark: this is rather expensive    */
      unur_tdr_chg_truncated(GEN_GAMMA, 0., c->beta * c->height );
    gx = unur_sample_cont(GEN_GAMMA) / (c->beta);
      
    /* nonnegative uniform random numbers with sum u_i = 1 */
    _unur_mvtdr_simplex_sample(gen, S);
      
    /* move point into center */
    for( i=0; i<GEN->dim; i++ ) rpoint[i] = GEN->center[i];
      
    /* calculate random point on chosen hyper-plane */
    for( j=0; j<GEN->dim; j++ ) {
      double x = gx * S[j] / c->gv[j];
      for( i=0; i<GEN->dim; i++ )
	rpoint[i] += x * (c->v[j])->coord[i];
    }

    /*.......................................................................*/
    /** accept or reject **/

    f = PDF(rpoint);                        /* density */
    h = T_inv( c->alpha - c->beta * gx );   /* hat */

    /* verify hat function */
    if ( (gen->variant & MVTDR_VARFLAG_VERIFY) &&
	 ((1.+UNUR_EPSILON) * h < f ) )
      _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF(x) > hat(x)");

    /* accept point */
    if( _unur_call_urng(gen->urng) * h <= f )
      return UNUR_SUCCESS;
  }

} /* end of _unur_mvtdr_sample_cvec() */

/*-----------------------------------------------------------------*/

int
_unur_mvtdr_simplex_sample( const struct unur_gen *gen, double *U )
     /*----------------------------------------------------------------------*/
     /* sample point uniformly on standard simplex                           */
     /* point in standard simplex 0 <= u[0] <= ... <= u[dim-2] <= 1          */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   U   ... array for storing result                                   */
     /*----------------------------------------------------------------------*/
{
  int dim = GEN->dim;

  /*................................................................*/
  if (dim == 2) {
    U[0] = _unur_call_urng(gen->urng);
    U[1] = 1. - U[0];
    return UNUR_SUCCESS;
  }
  /*................................................................*/
  if (dim == 3) {
    U[0] = _unur_call_urng(gen->urng);
    U[1] = _unur_call_urng(gen->urng);
    if( U[0] > U[1] ) {
      U[2] = U[0]; U[0] = U[1]; U[1] = U[2];
    }
    U[2] = 1. - U[1];
    U[1] = U[1] - U[0];
    return UNUR_SUCCESS;
  }
  /*................................................................*/
  if (dim >3) {
    int i,j;
    double U_aux;

    /* generate dim-1 numbers */
    for( i=0; i<dim-1; i++ )
      U[i] = _unur_call_urng(gen->urng);

    /* sort numbers (insertion sort) */
    /** TODO!! replace by exp random variates!! **/
    for( i=1; i<dim-1; i++ ) {
      U_aux = U[i];
      for( j=i; j>0 && U[j-1] > U_aux; j-- )
	U[j] = U[j-1];
      U[j] = U_aux;
    }
    
    /* transform to affine coordinates */
    U[dim-1] = 1.;
    for( i=dim-1; i>0; i-- )
      U[i] -= U[i-1];

    return UNUR_SUCCESS;
  }
  /*................................................................*/
  /* else --> make error message!! */
  _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return UNUR_FAILURE;

} /* end of _unur_mvtdr_simplex_sample() */

/*---------------------------------------------------------------------------*/
