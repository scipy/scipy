/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_ps_sample.c                                              *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Sampling routines for variant with proportional squeezes.            *
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
_unur_tdr_ps_sample( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator (proportional squeeze)                         */
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
     /*                                                                      */
     /*   hat(x) = f(x0) * exp( (Tf)'(x0) *  (x-x0) )                        */
     /*   generation:                                                        */
     /*      X = x0 + 1/(Tf)'(x0) * \log( (Tf)'(x0)/f(x0) * U + 1 )          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*                                                                      */
     /*----------------------------------------------------------------------*/
     /*   T(x) = -1/sqrt(x):                                                 */
     /*                                                                      */
     /*   hat(x) = 1 / (Tf(x0) + (Tf)'(x0) * (x-x0))^2                       */
     /*   generation:                                                        */
     /*      X = x0 + (Tf(x0)^2 * U) / (1 - Tf(x0) * (Tf)'(x0) * U)          */
     /*      U ~ U(-area below left hat, area below left hat)                */
     /*----------------------------------------------------------------------*/
{ 
  UNUR_URNG *urng;             /* pointer to uniform random number generator */
  struct unur_tdr_interval *iv;
  double U, V;                 /* uniform random number                      */
  double X;                    /* generated point                            */
  double fx;                   /* value of density at X                      */
  double Thx;                  /* value of transformed hat at X              */

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
    U -= iv->Acum - iv->Ahatr;    /* result: U in (-A_hatl, A_hatr) */

    /* generate from hat distribution */
    switch (gen->variant & TDR_VARMASK_T) {

    case TDR_VAR_T_LOG:
      if (_unur_iszero(iv->dTfx))
	X = iv->x + U / iv->fx;
      else {
	double t = iv->dTfx * U / iv->fx;
	if (fabs(t) > 1.e-6)
	  /* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	  X = iv->x + log(t + 1.) * U / (iv->fx * t);
	else if (fabs(t) > 1.e-8)
	  /* use Taylor series */
	  X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
	else
	  X = iv->x + U / iv->fx * (1 - t/2.);
      }
      break;

    case TDR_VAR_T_SQRT:
      if (_unur_iszero(iv->dTfx))
	X = iv->x + U / iv->fx;
      else {
	/* it would be less expensive to use:
	   X = iv->x + iv->Tfx/iv->dTfx * (1. - 1./(1. + iv->dTfx * iv->Tfx * U) )
	   however, this is unstable for small iv->dTfx */
	X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
	/* It cannot happen, that the denominator becomes 0 ! */
      }
      break;

    case TDR_VAR_T_POW:
      /** TODO **/

    default:
      _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
      return INFINITY;

    } /* end switch */

    /* accept or reject */
    V = _unur_call_urng(urng);

    /* squeeze rejection */
    if (V <= iv->sq)
      	return X;

    /* evaluate hat at X:
       get uniform random number between 0 and hat(X) */
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      V *= iv->fx * exp(iv->dTfx*(X - iv->x)); break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      V *= 1./(Thx*Thx); break;
    case TDR_VAR_T_POW:
      /** TODO **/
    default:
      return INFINITY;
    } /* end switch */

    /* evaluate PDF at X */
    fx = PDF(X);

    /* main rejection */
    if (V <= fx)
      return X;

    /* being above squeeze is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_ps_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }

} /* end of _unur_tdr_ps_sample() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_ps_sample_check( struct unur_gen *gen )
     /*----------------------------------------------------------------------*/
     /* sample from generator and verify results (proportional squeeze)      */
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
  struct unur_tdr_interval *iv;
  double U, V;                 /* uniform random number                      */
  double X;                    /* generated point                            */
  double fx, sqx, hx;          /* values of density, squeeze, and hat at X   */
  int squeeze_rejection = FALSE; /* indicates squeeze rejection              */
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
    X = _unur_tdr_ps_eval_invcdfhat( gen, U, &hx, &fx, &sqx, &iv );

    /* accept or reject */
    V = _unur_call_urng(urng);

    /* squeeze rejection */
    if (V <= iv->sq)
      squeeze_rejection = TRUE;

    /* get uniform random number between 0 and hat(X) */
    V *= hx;

    /* check result */
    if (_unur_FP_less(X, DISTR.BD_LEFT) || _unur_FP_greater(X, DISTR.BD_RIGHT) ) {
      _unur_warning(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"generated point out of domain");
      error = 1;
    }
    if (_unur_FP_greater(fx, hx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF > hat. Not T-concave!");
      error = 1;
    }
    if (_unur_FP_less(fx, sqx)) {
      _unur_warning(gen->genid,UNUR_ERR_GEN_CONDITION,"PDF < squeeze. Not T-concave!");
      error = 1;
    }

#ifdef UNUR_ENABLE_LOGGING
    /* write info into LOG file (in case error) */
    if (error && (gen->debug & TDR_DEBUG_SAMPLE)) 
      _unur_tdr_ps_debug_sample( gen, iv, X, fx, hx, sqx ); 
#endif

    /* squeeze rejection */
    if (squeeze_rejection)
      return X;

    /* main rejection */
    if (V <= fx)
      return X;

    /* being above squeeze is bad. improve the situation! */
    if (GEN->n_ivs < GEN->max_ivs) {
      if ( (_unur_tdr_ps_improve_hat( gen, iv, X, fx) != UNUR_SUCCESS)
	   && (gen->variant & TDR_VARFLAG_PEDANTIC) )
	return UNUR_INFINITY;
    }

    /* else reject and try again */

    /* use the auxilliary generator the next time
       (it can be the same as the main generator) */
    urng = gen->urng_aux;

  }
} /* end of _unur_tdr_ps_sample_check() */

/*---------------------------------------------------------------------------*/

double
_unur_tdr_ps_eval_invcdfhat( const struct unur_gen *gen, double U,
			     double *hx, double *fx, double *sqx,
			     struct unur_tdr_interval **ivl )
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
     /*                                                                      */
     /* return:                                                              */
     /*   inverse hat CDF.                                                   */
     /*                                                                      */
     /*   values of hat, density, squeeze (computation is suppressed if      */
     /*   corresponding pointer is NULL).                                    */
     /*                                                                      */
     /*   ivl can be NULL. then no pointer is stored.                        */
     /*----------------------------------------------------------------------*/
{
  struct unur_tdr_interval *iv;             /* pointer to intervals          */
  double X;                                 /* inverse of hat CDF at U       */  
  double Thx;                               /* transformed squeeze           */
  double t;                                 /* aux variable                  */

  /** -1- Compute inverse of hat CDF at U **/ 

  /* look up in guide table and search for segment */
  iv =  GEN->guide[(int) (U * GEN->guide_size)];
  U *= GEN->Atotal;
  while (iv->Acum < U) {
    iv = iv->next;
  }

  /* transform u such that u in (-A_hatl, A_hatr) */
  /* (reuse of uniform random number)             */
  U -= iv->Acum - iv->Ahatr;

  /* inverse of CDF (random variate) */
  switch (gen->variant & TDR_VARMASK_T) {

  case TDR_VAR_T_LOG:
    if (_unur_iszero(iv->dTfx))
      X = iv->x + U / iv->fx;
    else {
      t = iv->dTfx * U / iv->fx;
      if (fabs(t) > 1.e-6)
	/* x = iv->x + log(t + 1.) / iv->dTfx; is cheaper but numerical unstable */
	X = iv->x + log(t + 1.) * U / (iv->fx * t);
      else if (fabs(t) > 1.e-8)
	/* use Taylor series */
	X = iv->x + U / iv->fx * (1 - t/2. + t*t/3.);
      else
	X = iv->x + U / iv->fx * (1 - t/2.);
    }
    break;

  case TDR_VAR_T_SQRT:
    if (_unur_iszero(iv->dTfx))
      X = iv->x + U / iv->fx;
    else {
      /* it would be less expensive to use:
	 X = iv->x + iv->Tfx/iv->dTfx * (1. - 1./(1. + iv->dTfx * iv->Tfx * U) )
	 however, this is unstable for small iv->dTfx */
      X = iv->x + (iv->Tfx*iv->Tfx*U) / (1.-iv->Tfx*iv->dTfx*U);  
      /* It cannot happen, that the denominator becomes 0 ! */
    }
    break;
    
  case TDR_VAR_T_POW:
    /** TODO **/

  default:
    _unur_error(gen->genid,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
    return INFINITY;

  } /* end switch */

  /** -2- Evaluate hat at X **/
  if (hx != NULL) 
    switch (gen->variant & TDR_VARMASK_T) {
    case TDR_VAR_T_LOG:
      *hx = iv->fx * exp(iv->dTfx*(X - iv->x));
      break;
    case TDR_VAR_T_SQRT:
      Thx = iv->Tfx + iv->dTfx * (X - iv->x);     /* transformed hat at X */ 
      *hx = 1./(Thx*Thx);
      break;
    case TDR_VAR_T_POW:
      /** TODO **/
    default:
      *hx = INFINITY;
    } /* end switch */

  /** -3- Evaluate density at X **/
  if (fx != NULL) {
    *fx = PDF(X);
  }

  /** -4- Evaluate squeeze at X **/
  if (sqx != NULL && hx != NULL) {
    *sqx = *hx * iv->sq;
  }

  /* store interval pointer */
  if (ivl) *ivl = iv;

  return X;
} /* end of _unur_tdr_ps_eval_invcdfhat() */


/*---------------------------------------------------------------------------*/

int
_unur_tdr_ps_improve_hat( struct unur_gen *gen, struct unur_tdr_interval *iv,
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
  result = _unur_tdr_ps_interval_split(gen, iv, x, fx);
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
} /* end of _unur_tdr_ps_improve_hat() */

/*****************************************************************************/
