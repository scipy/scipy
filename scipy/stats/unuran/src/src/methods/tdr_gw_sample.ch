/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_gw_sample.c                                              *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Sampling routines for variant of Gilks and Wild.                     *
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
/**  Sampling routines                                                      **/
/*****************************************************************************/

double
_unur_tdr_gw_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (original variant by Gilks & Wild)             */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*                                                                      */
     /* return:                                                              */
     /*   double (sample from random variate)                                */
     /*                                                                      */
     /* error:                                                               */
     /*   return INFINITY                                                    */
     /*                                                                      */
     /*======================================================================*/
     /* comment:                                                             */
     /*   x   ... random point                                               */
     /*   x0  ... left construction point in interval                        */
     /*   x1  ... right construction point in interval                       */
     /*   f   ... PDF                                                        */
     /*   Tf  ... transformed PDF                                            */
     /*   dTf ... derivative of transformed PDF                              */
     /*   sq  ... slope of squeeze in interval                               */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   if (Tf)'(x0) == 0:                                                 */
     /*   X = x0 + U / f(x0)                                                 */
     /*   U ~ U(0,area below hat)                                            */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   log(x):                                                            */
     /*   squeeze(x) = f(x0) * exp(sq * (x-x0))                              */
     /*                                                                      */
     /*   left hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                   */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = f(x1) * exp( (Tf)'(x1) *  (x-x1) )                  */
     /*   generation:                                                        */
     /*      X = x1 + 1/(Tf)'(x1) * \log( (Tf)'(x1)/f(x1) * U + 1 )          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   squeeze(x) = 1 / (Tf(x0) + sq * (x-x0))^2                          */
     /*                                                                      */
     /*   left hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                  */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(0,area below left hat)                                    */
     /*                                                                      */
     /*   right hat(x) = 1 / (Tf(x1) + (Tf)'(x1) * (x-x1))^2                 */
     /*   generation:                                                        */
     /*      X = x1 + (Tf(x1)^2 * U) / (1 - Tf(x1) * (Tf)'(x1) * U)          */
     /*      U ~ U(- area below right hat,0)                                 */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv, *pt;
  double U, V;                 /* uniform random number                      */
  double X;                    /* generated point                            */
  double fx, sqx, hx;          /* values of density, squeeze, and hat at X   */
  double Tsqx, Thx;            /* values of transformed squeeze and hat at X */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN->Umin + _unur_call_urng(urng) * (GEN->Umax - GEN->Umin);

    /* look up in guide table and search for segment */
    iv =  GEN->guide[(int) (U * GEN->guide_size)];
    U *= GEN->Atotal;
    while (iv->Acum < U) {
      iv = iv->next;
    }

    /* reuse of uniform random number */
    U -= iv->Acum;    /* result: U in (-A_hat, 0) */

    /* left or right side of hat */
    if (-U < iv->Ahatr) { /* right */
      pt = iv->next;
      /* u unchanged */
    }
    else {                /* left */
      pt = iv;
      U += iv->Ahat;
    }

    /* we have three different types of transformations */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      /* random variate */
      if (_unur_iszero(pt->dTfx))
	X = pt->x + U / pt->fx;
      else
	{
	  double t = pt->dTfx * U / pt->fx;
	  if (fabs(t) > 1.e-6)
	    X = pt->x + log(t + 1.) * U / (pt->fx * t);
	  /* x = pt->x + log(t + 1.) / pt->dTfx; is cheaper but numerical unstable */
	  else if (fabs(t) > 1.e-8)
	    /* use Taylor series */
	    X = pt->x + U / pt->fx * (1 - t/2. + t*t/3.);
	  else
	    X = pt->x + U / pt->fx * (1 - t/2.);
	}

      /* accept or reject */
      hx = pt->fx * exp(pt->dTfx*(X - pt->x));    /* value of hat at x */   
      V = _unur_call_urng(urng) * hx;  /* a random point between 0 and hat at x */
      
      /* below mininum of density in interval ? */
      if (V <= iv->fx && V <= iv->next->fx)
	return X;

      /* below squeeze ? */
      sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(X - iv->x)) : 0.;     /* value of squeeze at x */
      if (V <= sqx)
	return X;

      break;

    case TDR_VAR_T_SQRT:
      /* random variate */
      if (_unur_iszero(pt->dTfx))
	X = pt->x + U /pt->fx;
      else {
	/* it would be less expensive to use:
	   X = pt->x + pt->Tfx/pt->dTfx * (1. - 1./(1. + pt->dTfx * pt->Tfx * U) )
	   however, this is unstable for small pt->dTfx */
	X = pt->x + (pt->Tfx*pt->Tfx*U) / (1.-pt->Tfx*pt->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }

      /* accept or reject */
      Thx = pt->Tfx + pt->dTfx * (X - pt->x);     /* transformed hat at x */ 
      hx = 1./(Thx*Thx);
      V = _unur_call_urng(urng) * hx;  /* a random point between 0 and hat at x */

      /* below mininum of density in interval ? */
      if (V <= iv->fx && V <= iv->next->fx)
	return X;

      /* below squeeze ? */
      Tsqx = (iv->Asqueeze > 0.) ? (iv->Tfx + iv->sq * (X - iv->x)) : -INFINITY; /* transformed squeeze at x */ 
      sqx = (iv->Asqueeze > 0.) ? 1./(Tsqx*Tsqx) : 0.;
      if (V <= sqx)
	return X;
      
      break;

    case TDR_VAR_T_POW:
      /** TODO **/

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return INFINITY;

    } /* end switch */

    /* value of PDF at x */
    fx = PDF(X);
    /** TODO: for the changing parameter case it would be better
        not to added the new construction point when we accept
        the generated point! **/

    /* being above squeeze is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_gw_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    if (V <= fx)
      /* between PDF and squeeze */
      return X;

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }

} /* end of _unur_tdr_gw_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_gw_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results                             */
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
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv, *pt;
  double U, V;                 /* uniform random number                      */
  double X;                    /* generated point                            */
  double fx, sqx, hx;          /* value of density, squeeze, and hat at X    */
  int error = 0;               /* indicates error                            */

  /* check arguments */
  CHECK_NULL(gen,INFINITY);  COOKIE_CHECK(gen,CK_TDR_GEN,INFINITY);

  if (GEN->iv == NULL) {
    _unur_error(gen->genid,UNUR_ERR_GEN_DATA,"empty generator object");
    return INFINITY;
  } 

  /* main URNG */
  urng = gen->urng;

  while (1) {

    /* sample from U( Umin, Umax ) */
    U = GEN->Umin + _unur_call_urng(urng) * (GEN->Umax - GEN->Umin);

    /* evaluate inverse of hat CDF */
    X = _unur_tdr_gw_eval_invcdfhat( gen, U, &hx, &fx, &sqx, &iv, &pt );

    /* check result */
    if (X < DISTR.BD_LEFT || X > DISTR.BD_RIGHT) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (_unur_FP_greater(fx,hx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
      error = 1;
    }
    if (_unur_FP_less(fx,sqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
      error = 1;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file (in case error) */
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_gw_debug_sample( gen, iv, pt, X, fx, hx, sqx ); 
#endif

    /* accept or reject */
    V = _unur_call_urng(urng) * hx;  /* a random point between 0 and hat at X */

    /* below mininum of density in interval ? */
    if (V <= iv->fx && V <= iv->next->fx)
      return X;

    /* below squeeze ? */
    if (V <= sqx)
      return X;

    /* being above squeeze is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_gw_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    if (V <= fx)
      /* between PDF and squeeze */
      return X;

    /* reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }

} /* end of _unur_tdr_gw_sample_check() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_gw_eval_invcdfhat( const struct unur_gen *gen, double U, 
			     double *hx, double *fx, double *sqx,
			     struct unur_tdr_interval **ivl,
			     struct unur_tdr_interval **cpt )
     /*----------------------------------------------------------------------*/
     /* evaluate the inverse of the hat CDF at u                             */
     /* (original variant by Gilks & Wild)                                   */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen ... pointer to generator object                                */
     /*   U   ... argument for inverse CDF (0<=U<=1, no validation!)         */
     /*   hx  ... pointer for storing hat at sampled X                       */
     /*   fx  ... pointer for storing squeeze at sampled X                   */
     /*   sqx ... pointer for storing density at sampled X                   */
     /*   ivl ... pointer to interval for hat                                */
     /*   cpt ... pointer to interval where hat construction point is stored */
     /*                                                                      */
     /* return:                                                              */
     /*   inverse hat CDF.                                                   */
     /*                                                                      */
     /*   values of hat, density, squeeze (computation is suppressed if      */
     /*   corresponding pointer is NULL).                                    */
     /*                                                                      */
     /*   ivl and cpt can be NULL. then no pointers are stored.              */
     /*----------------------------------------------------------------------*/
{ 
  struct unur_tdr_interval *iv, *pt;        /* pointer to intervals          */
  double X;                                 /* inverse of hat CDF at U       */  
  double Tsqx, Thx;                         /* transformed hat and squeeze   */
  double t;                                 /* aux variable                  */

  /** -1- Compute inverse of hat CDF at U **/ 

  /* look up in guide table and search for interval */
  iv =  GEN->guide[(int) (U * GEN->guide_size)];
  U *= GEN->Atotal;
  while (iv->Acum < U) {
    iv = iv->next;
  }

  /* transform u such that u in (-A_hat, 0) */
  /* (reuse of uniform random number)       */
  U -= iv->Acum;

  /* left or right hand side of hat */
  if (-U < iv->Ahatr) { /* r.h.s. */
    pt = iv->next;
    /* U unchanged */
  }
  else {                /* l.h.s. */
    pt = iv;
    U += iv->Ahat;
  }
    
  /* we have three different types of transformations */
  switch (gen->variant & TDR_VARMASK_T) {
    
  case TDR_VAR_T_LOG:
    /* inverse of CDF (random variate) */
    if (_unur_iszero(pt->dTfx))
      X = pt->x + U / pt->fx;
    else {
      t = pt->dTfx * U / pt->fx;
      if (fabs(t) > 1.e-6)
	X = pt->x + log(t + 1.) * U / (pt->fx * t);
      /* x = pt->x + log(t + 1.) / pt->dTfx; is cheaper but numerical unstable */
      else if (fabs(t) > 1.e-8)
	/* use Taylor series */
	X = pt->x + U / pt->fx * (1 - t/2. + t*t/3.);
      else
	X = pt->x + U / pt->fx * (1 - t/2.);
    }
    break;

  case TDR_VAR_T_SQRT:
    /* inverse of CDF (random variate) */
    if (_unur_iszero(pt->dTfx))
      X = pt->x + U / (pt->fx);
    else {
      /* it would be less expensive to use:
	 X = pt->x + pt->Tfx/(pt->dTfx) * (1. - 1./(1. + pt->dTfx * pt->Tfx * U) )
	 however, this is unstable for small pt->dTfx */
      X = pt->x + (pt->Tfx*(pt->Tfx)*U) / (1.-(pt->Tfx)*(pt->dTfx)*U);  
      /* It cannot happen, that the denominator becomes 0 ! */
    }
    break;

  case TDR_VAR_T_POW:
    /** TODO **/
    
  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    X = INFINITY;
    
  } /* end switch */

  /** -2- Evaluate hat at X **/

  if (hx != NULL) {
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      *hx = pt->fx * exp(pt->dTfx*(X - pt->x));
      break;
    case TDR_VAR_T_SQRT:
      Thx = pt->Tfx + pt->dTfx * (X - pt->x);     /* transformed hat at x */ 
      *hx = 1./(Thx*Thx);
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      *hx = INFINITY;
    }
  }

  /** -3- Evaluate density at X **/
  if (fx != NULL) {
    *fx = PDF(X);
  }

  /** -4- Evaluate squeeze at X **/
  if (sqx != NULL) {
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      *sqx = (iv->Asqueeze > 0.) ? iv->fx * exp(iv->sq*(X - iv->x)) : 0.;
      break;
    case TDR_VAR_T_SQRT:
      if (iv->Asqueeze > 0.) {
	Tsqx = iv->Tfx + iv->sq * (X - iv->x);
	*sqx = 1./(Tsqx*Tsqx);
      }
      else
	*sqx = 0.;
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
    default:  /* this should not happen */
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      *sqx = 0.;
    }
  }


  /* store interval pointer */
  if (ivl) *ivl = iv;
  if (cpt) *cpt = pt;

  /* o.k. */
  return X;

} /* end of _unur_tdr_gw_eval_invcdfhat() */

/*---------------------------------------------------------------------------*/

int
_unur_tdr_gw_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv,
			  double x, double fx )
     /*----------------------------------------------------------------------*/
     /* improve hat function by splitting interval                           */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen        ... pointer to generator object                         */
     /*   iv         ... pointer to interval that has to be split            */
     /*   x          ... splitting point                                     */
     /*   fx         ... value of PDF at splitting point                     */
     /*                                                                      */
     /* return:                                                              */
     /*   UNUR_SUCCESS    ... improving hat successful                       */
     /*   others          ... error: PDF not monotone in interval            */
     /*----------------------------------------------------------------------*/
{
  int result;

  /* is there any reason to improve hat ? */
  if (! (GEN->max_ratio * GEN->Atotal > GEN->Asqueeze) ) {
    /* no more construction points (avoid calling this function any more) */
    GEN->max_ivs = GEN->n_ivs;
    return UNUR_SUCCESS;
  }

  /* add construction point */
  result = _unur_tdr_gw_interval_split(gen, iv, x, fx);
  if (result!=UNUR_SUCCESS && result!=UNUR_ERR_SILENT && result!=UNUR_ERR_INF) {
    /* condition for PDF is violated! */
    _unur_error(gen->genid,UNUR_ERR_GEN_CONDITION,"");
    if (gen->variant & TDR_VARFLAG_PEDANTIC || result == UNUR_ERR_ROUNDOFF) {
      /* replace sampling routine by dummy routine that just returns INFINITY */
      SAMPLE = _unur_sample_cont_error;
      return UNUR_ERR_GEN_CONDITION;
    }
  }

  /* splitting successful --> update guide table */
  /** TODO: it is not necessary to update the guide table every time. 
      But then (1) some additional bookkeeping is required and
      (2) the guide table method requires a acc./rej. step. **/
  _unur_tdr_make_guide_table(gen);

  /* o.k. */
  return UNUR_SUCCESS;
} /* end of _unur_tdr_gw_improve_hat() */

/*****************************************************************************/
